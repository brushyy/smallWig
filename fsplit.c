/*split file for decompress
 * options:
 * 0--regular
 * 1--lpaq
 * 2--random access
 */
 


	
int fsplit(int *opt, char * argv[]) {
	
	char FileName[400];//name of processed Wig file
	
	clock_t start = clock();
	SPAM(("=================\nfsplit.cpp start!\n"));
	
	strcpy(FileName, argv[1]);
	
	
	int N=20;
	
	char Type[20][40];
	char c[20][10000];
	FILE *fp[20];
	FILE *fpout;
	long int sz[20];
	int i,j,tmp;
	char *buffer;
	
	//output file
	i=0;
	strcpy(c[i],FileName);
	fpout = fopen(c[0], "r");
	if (fpout  == NULL) {
		fprintf(stderr, "Can't open file %s!\n",c[0]);
		exit(1);
	}
	
	//compression option
	fscanf(fpout,"%d ",opt);
	SPAM(("opt=%d\n",*opt));
	switch (*opt){
		case 2://encoded by blocks, 
			fscanf(fpout,"%d ",&BlockSize);
			//int FloatSize;
			//fscanf(fpout,"%d ",&FloatSize);
			//SetFloatSize = pow10(FloatSize);
			break;
	}
	
	for (i=0; i<20; i++){
		for (j=0; j<40; j++){
			Type[i][j] = list[*opt][i][j];
		}
	}
	
	i=0;
	while (i<20 && Type[i][0]){
		strcpy(c[i],FileName);
		strcat(c[i],Type[i]);
		i++;
	}
	N = i;
	
	
	//file sizes in first line of output file
	for (i=1; i<N; i++){
		fp[i] = fopen(c[i], "w"); 
		if (fp[i] == NULL) {
			fprintf(stderr, "Can't open file %s!\n",c[i]);
			exit(1);
		}
		//file size
		fscanf(fpout,"%ld ",&sz[i]);
	}
	
	
	//copy every file
	for (i=1; i<N; i++){
		buffer = (char*)malloc(sz[i]*sizeof(char));
		fread(buffer,sizeof(char),sz[i],fpout);
		fwrite(buffer,sizeof(char),sz[i],fp[i]);
		fclose(fp[i]);
		free(buffer);
	}
	fclose(fpout);
	
	
	clock_t end = clock();
	SPAM(("fsplit takes %2.3f sec.\n", (double)(end - start)/CLOCKS_PER_SEC));
		
	
	return 1;
	
}


