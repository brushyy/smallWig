/*
 * merge block address files for random_acess compression in parallel
 */		
 

	
int fmergeadd(int nproc,char * argv[]) {

	char FileName[200];//name of processed Wig file
	
	//clock_t start = clock();
	//printf("=================\nfmergecount.c start!\n");

	strcpy(FileName, argv[1]);
	int N = 5;
	char Type[5][40]={"DiffSeqMatlabChrmAdd","LenDiffSeqMatlabChrmAdd","BlockAux","DiffSeqMatlabBlockAdd","LenDiffSeqMatlabBlockAdd"};
	char c[10000];
	char str[1000];
	FILE **fp = (FILE**)malloc(nproc * sizeof(FILE*));
	FILE *fpout;
	int i,j;
	int *tmpint = (int *)malloc((nproc+1) * sizeof(int));
	int addr; //address
//	int exceed; //exceed for fpExceed
	int tmp = 0; //first 0 to be added to beginning of address files
	
	
	for (j=0; j<N; j++){
		//output count file
		strcpy(c,FileName);
		strcat(c,Type[j]);
		fpout = fopen(c,"w+");
		if (fpout == NULL){
			fprintf(stderr, "Can't open file %s!\n",c);
			exit(1);
		}
		if (j<3)
			fprintf(fpout,"0 ");//first addr
		else if (j>=3)//BlockAdd, also write first addr as tmp=0
			fwrite(&tmp, sizeof(int), 1, fpout);
		//input prll addr files 
		tmpint[0] = 0; //previous files size
		for (i=0;i<nproc;i++){
			strcpy(c,FileName);
			strcat(c,Type[j]);
			sprintf(str,"%02d",i);
			strcat(c,str);
			fp[i] = fopen(c,"r");
			if (fp[i] == NULL){
				fprintf(stderr, "Can't open file %s!\n",c);
				exit(1);
			}
			if (j<3){
				fscanf(fp[i], "%d ",&addr); //skip first addr=0
				while(!feof(fp[i])){
					fscanf(fp[i], "%d ",&addr);
					addr += tmpint[i];
					fprintf(fpout, "%d ",addr);
				}
				tmpint[i+1] = addr; 
			}//when j=2, or BlockAux, tmpint[i] stores cummulated number of blocks up to the first chrm corresponding to the i-th processor  
			else{//BlockAdd
				fseek (fp[i] , 0 , SEEK_END);
				long int size = ftell (fp[i]);
				rewind (fp[i]);
				fread(&addr, sizeof(int), 1, fp[i]); //skip first addr=0
				int printcount;
				for (printcount=1; printcount<size/sizeof(int); printcount++){
					fread(&addr, sizeof(int), 1, fp[i]);
					addr += tmpint[i];
					fwrite(&addr, sizeof(int), 1, fpout);
				}
				tmpint[i+1] = addr; 
			}
			fclose(fp[i]);
			remove(c);
		}
		fclose(fpout);
	}

	
	free(tmpint);
	free(fp);
	
	//clock_t end = clock();
	//printf("fmergeadd takes %2.3f sec.\n", (double)(end - start)/CLOCKS_PER_SEC);
		
	
	return 0;
	
}
