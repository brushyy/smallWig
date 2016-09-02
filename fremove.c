/*
 * remove files after decompression
 * options: 
 * 0--regular
 * 1--lpaq
 * 2--random access
 * 3--random access with parallel, but in this case opt=2 is written in compressed file
 */	

int fremove(int opt, char * argv[], int nproc) {
	
	char FileName[1000];//name of processed Wig file
	char tmpchar[1000], str[2];
	
	clock_t start = clock();
	SPAM(("=================\nfremove start! opt=%d\n",opt));
	
	strcpy(FileName, argv[1]);
	
	
	
	
	char Type[20][40];
	char c[20][10000];
	int i,j;
	
	
	for (i=0; i<20 && list[opt][i][0]; i++){
		for (j=0; j<40; j++){
			Type[i][j] = list[opt][i][j];
		}
	}
	int N=i;
	//extra files to delete
	char dlist[2][40]={"DiffSeqMatlabDecode","LenDiffSeqMatlabDecode"};
	for (i=0; i<2; i++){
		for (j=0; j<40; j++){
			Type[i+N][j] = dlist[i][j];
		}
		//printf("type %d is %s\n",i+N,Type[i+N]);
	}

	
	i=1;
	while (i<20 && Type[i][0]){
		strcpy(c[i],FileName);
		strcat(c[i],Type[i]);
		if (opt==3 && i>=N){
			for (j=0;j<nproc;j++){
				strcpy(tmpchar,c[i]);
				sprintf(str,"%02d",j);
				strcat(tmpchar,str);
				remove(tmpchar);
			}
		}
		else		
			remove(c[i]);
		//printf("remove file %s\n",Type[i]);
		i++;
	}
	
	i=0;
	if (opt==3){ //mergesub of the output wig file
		strcpy(c[i],argv[2]);
		fmergesub(nproc,c[i]);
	}
	
	
	clock_t end = clock();
	SPAM(("fremove takes %2.3f sec.\n", (double)(end - start)/CLOCKS_PER_SEC));
		
	
	return 1;
	
}
 
