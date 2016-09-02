/* 
 * get matlabseq from seq
 * get sequences for matlab arithmetic coding use
 * DiffMatlab is the count vector w. all counts >=0
 * DiffSeqMatlab is consecutive values from 0 to size DiffMatlab
 * 
 * use this after getblock.c (and fmergecount.c if parallel)
 */
	
int getmatlabseq(char * argv[],int nproc,int iproc){	
	
	int Nd = 18;
	int SetMaxDiff = exp2(Nd); //set max value difference 
	double DiffTable[SetMaxDiff]; //alphabet
	memset(DiffTable,0,SetMaxDiff*sizeof(double));
	int DiffSize;
	char tmpFileName[1000], str[1000];
	int i, tmpint;
	double Diff;
	int DiffMatlab;
	char counttype[2][10]={"Diff","LenDiff"}, seqtype[2][20]={"DiffSeq", "LenDiffSeq"};
	char seqmatlabtype[2][20]={"DiffSeqMatlab","LenDiffSeqMatlab"};
	int type;
	int indicator;
	
	for (type=0; type<2; type++){
		
		//read from count files
		strcpy(tmpFileName,argv[1]);
		strcat(tmpFileName,counttype[type]);
		FILE *fpDiff = fopen(tmpFileName,"r");
		if (fpDiff == NULL) {
			fprintf(stderr,"Can't open output Count file!\n");
			exit(1);
		}
		//skip annotation
		fgets(str,sizeof(str),fpDiff);
		
		i = 0;
		while (fread(&DiffTable[i],sizeof(double),1,fpDiff)==1 && i < SetMaxDiff){
			//fscanf(fpDiff,"%d\t%*u\n", &DiffTable[i]);
			fread(&tmpint,sizeof(int),1,fpDiff);
			i++;
		}
		DiffSize = i;
		if (DiffSize > SetMaxDiff){
			fprintf(stderr,"Too large alphabet!\n");
			exit(1);
		}
		fclose(fpDiff);
		  
		// diffseq
		strcpy(tmpFileName,argv[1]);
		strcat(tmpFileName,seqtype[type]);
		if (nproc>0){
			tmpint = sprintf(str,"%02d",iproc);
			strcat(tmpFileName,str);
		}
		FILE *fpDiffSeq = fopen(tmpFileName, "r");
		if (fpDiffSeq == NULL) {
			fprintf(stderr, "Can't open file %s!\n",tmpFileName);
			exit(1);
		}
		
		//matlabseq files
		strcpy(tmpFileName,argv[1]);
		strcat(tmpFileName,seqmatlabtype[type]);
		if (nproc>0){
			tmpint = sprintf(str,"%02d",iproc);
			strcat(tmpFileName,str);
		}
		FILE *fpDiffSeqMatlab = fopen(tmpFileName,"w");
		if (fpDiffSeqMatlab == NULL) {
			fprintf(stderr,"Can't open output Count file!\n");
			exit(1);
		}
		
		//write to matlab seq files
		while (fread(&Diff,sizeof(double),1,fpDiffSeq)==1){
			//fscanf(fpDiffSeq,"%d ",&Diff);
			DiffMatlab = binfind(DiffSize,DiffTable,Diff,&indicator); //binary find
			if (!indicator)
				fprintf(fpDiffSeqMatlab,"%d ",DiffMatlab);
			else{
				fprintf(stderr,"Alphabet=%lg not found while getting matlab sequence!\n",Diff);
				exit(1);
			}
		}
		
		fclose(fpDiffSeq);
		fclose(fpDiffSeqMatlab);
    }
    return 0;
}	
	
	
