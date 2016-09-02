/*
 * merge sub files
 */		
 
	
int fmergesub(int nproc,char *c) {
	
	
	//clock_t start = clock();
	//printf("=================\nfmergesub.c start! %s\n",c);
	

	char str[2],tmpchar[1000];
	FILE **fp = (FILE**)malloc(nproc * sizeof(FILE*));
	FILE *fpout;
	int i,j;
	char *buffer;
	long int sz;
	long int maxsz = 100000000; //max memory allocation size
	long int tmpsz; //size for each allocation
	
	for (i=0;i<nproc;i++){
		strcpy(tmpchar,c);
		sprintf(str,"%02d",i);
		strcat(tmpchar,str);
		if (i==0){//output count file
			fpout = fopen(tmpchar,"a");
			if (fpout == NULL){
				fprintf(stderr, "Can't open file %s!\n",tmpchar);
				exit(1);
			}
		}
		else{ //input prll count files 
			fp[i] = fopen(tmpchar,"r");
			if (fp[i] == NULL){
				fprintf(stderr, "Can't open file %s!\n",tmpchar);
				exit(1);
			}
			sz = fsize(fp[i]);
			
			for (j=0; j<(sz+maxsz-1)/maxsz; j++){ //ceil(sz/maxsz) = floor((sz+maxsz-1)/maxsz)
				if (j == ((sz+maxsz-1)/maxsz) - 1)
					tmpsz = sz%maxsz;
				else
					tmpsz = maxsz;
				if (j==0){
					buffer = (char*)malloc(tmpsz*sizeof(char));
					if (buffer==NULL){
						printf("memory allocation not successful!\n");
						exit(1);
					}
				}
				fread(buffer,sizeof(char),tmpsz,fp[i]);
				fwrite(buffer,sizeof(char),tmpsz,fpout);
				if (j == ((sz+maxsz-1)/maxsz) - 1) 
					free(buffer);
			} 
			fclose(fp[i]);
			remove(tmpchar);
		}
	}
	fclose(fpout);
	free(fp);
	
	//rename output file 
	strcpy(tmpchar,c);
	sprintf(str,"00");
	strcat(tmpchar,str);		
	rename(tmpchar, c);
	
	//clock_t end = clock();
	//printf("fmergecount takes %2.3f sec.\n", (double)(end - start)/CLOCKS_PER_SEC);
		
	return 0;
	
}
