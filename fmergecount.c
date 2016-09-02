/*
 * merge count files for encoder in parallel
 */		
 
int multimerge(
    FILE **fp,      				 // count files
    int number_of_arrays,            // number of files
    FILE *fpout               		 // file for output alphabet and count
){
    int i = 0;       // output cursor
    int j = 0;       // index for minimum search
    double min;         // minimum in this iteration
    int minposition; // position of the minimum
    int tmpcount;	 // count for alphabet=min
    double *arrays = (double *)malloc(sizeof(double)*number_of_arrays);	//alphabet
    int *count = (int *)malloc(sizeof(int)*number_of_arrays);	//count

    // cursor for the arrays
    int * cursor = (int *)malloc(number_of_arrays*sizeof(int));
	
	// initialize arrays to first items in fp[]
	for (j=0; j<number_of_arrays; j++){
		//fscanf(fp[j],"%d\t%d\n", &arrays[j],&count[j]);
		fread(&arrays[j],sizeof(double),1,fp[j]);
		fread(&count[j],sizeof(int),1,fp[j]);
	}

    while(1){
        min = 1<<30;
        minposition = -1; // invalid position
        tmpcount = 0;

        // Go through the current positions and get the minimum
        for(j = 0; j < number_of_arrays; ++j){
            if (arrays[j] < min){ // the element is smaller
				min = arrays[j];  // save the minimum ...
                minposition = j;  // ... and its position
            }
        }

        // if there is no minimum, then the position will be invalid
        if(minposition == -1)
            break;
		
		//find all poisitions that are equal to min
		for (j=0; j < number_of_arrays; j++){
			if (arrays[j]-min<1e-10 &&
				arrays[j]-min>-1e-10){	// the element arrays[j] is min
				tmpcount += count[j];// count for alphabet=min
				if (fread(&arrays[j],sizeof(double),1,fp[j])==1){
					//fscanf(fp[j],"%d\t%d\n", &arrays[j],&count[j]);// move to next elt
					fread(&count[j],sizeof(int),1,fp[j]);
				}
				else
					arrays[j] = (1<<30) + 1; // too large a value
		   }
		}
		
        //fprintf(fpout,"%d\t%d\n", min,tmpcount);//output to file   
        fwrite(&min,sizeof(double),1,fpout);
		fwrite(&tmpcount,sizeof(int),1,fpout);         
        i++;
    }
    free(cursor);
    free(arrays);
    free(count);

    return 0;
}
	
int fmergecount(int nproc,char * argv[]) {

	char FileName[200];//name of processed Wig file
	
	clock_t start = clock();
	SPAM(("=================\nfmergecount.c start!\n"));
	
	strcpy(FileName, argv[1]);
	
	char Type[2][40]={"Diff","LenDiff"};
	char **c;
	char str[1000];
	FILE **fp = (FILE**)malloc(nproc * sizeof(FILE*));
	FILE *fpout;
	int i,j;
	
	
	c = (char **)malloc(nproc * sizeof(char *));
	for(i = 0; i < nproc; i++){
		c[i] = (char *)malloc(10000 * sizeof(char));
	}
	
	for (j=0; j<2; j++){
		//output count file
		i=0; //c[0] will be replaced with sub file later
		strcpy(c[i],FileName);
		strcat(c[i],Type[j]);
		fpout = fopen(c[i],"w+");
		if (fpout == NULL){
			fprintf(stderr, "Can't open file %s!\n",c[i]);
			exit(1);
		}
		//input prll count files 
		for (i=0;i<nproc;i++){
			strcpy(c[i],FileName);
			strcat(c[i],Type[j]);
			sprintf(str,"%02d",i);
			strcat(c[i],str);
			fp[i] = fopen(c[i],"r");
			if (fp[i] == NULL){
				fprintf(stderr, "Can't open file %s!\n",c[i]);
				exit(1);
			}
			fgets(str,1000,fp[i]);
		}
		fputs(str,fpout);//first line of annotations
		multimerge(fp, nproc, fpout);
		
		for (i=0; i<nproc; i++){
			fclose(fp[i]);
			remove(c[i]);
		}
		fclose(fpout);
	}

	
	clock_t end = clock();
	SPAM(("fmergecount takes %2.3f sec.\n", (double)(end - start)/CLOCKS_PER_SEC));
	
	for(i = 0; i < nproc; i++){
		free(c[i]); 
	}
	free(c);	
	free(fp);
	
	return 0;
	
}
