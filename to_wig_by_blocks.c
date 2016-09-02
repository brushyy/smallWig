/*
 * Change a matlab seq_in_blocks back to wig format
 * Also support parallel processing
*/


int to_wig_by_blocks(char * argv[], int nproc, int iproc){
	
	clock_t start = clock();
	SPAM(("================\nto_wig_by_blocks.c start!\n"));
	
	char FileName[400];//name of processed Wig file
	strcpy(FileName, argv[1]);
	
	
	//for parallel
	int First; //first chrm to read
	//int Last;
	//int pflag; //indicate whether should start read
	if (nproc<=0){
		First = 0;
		//Last = 32;//assume no more than 32 chrm
		//pflag=1;
	}
	else{
		First = (24/nproc)*iproc; 
		/*if (iproc < nproc-1)
			Last = (24/nproc)*(iproc+1); 
		else
			Last = 32;*/
		//pflag=0; 
	}
	
	int Nd = 18;
	int SetMaxDiff = exp2(Nd); //set max alphabet size
	double Map[SetMaxDiff]; //maps each consecutive integer to the actual alphabet
	memset(Map,0,SetMaxDiff*sizeof(double)); 
	double LenMap[SetMaxDiff]; 
	memset(LenMap,0,SetMaxDiff*sizeof(double));  
	char ChrmNames[8*32]; //every chrm has a name less than 8 characters, assume no more than 32 chrm 
	memset(ChrmNames,0,8*32*sizeof(char));
	int ChrmCount=0;
	
	//int type;
	//char Type[40];
	
	double PreValue=-2;
	double Value=0;
	double Diff;//difference of values of contigs
	//int NextDiff=1;
	double PreLocation=-1;
	double Location=-1;
	double Contig = 1; //length of current contig
	double PreContig = 1; //length of previous contig
	double LenDiff; //length difference of two consecutive contigs
	char str[200];
	char tmpFileName[400];
	int i,j;
	int tmpint,tmp;
	int pcount=0;
	long long lcount=0; //count number of lines in wig file
	
	
	//for blocks use
	int SetBlockSize = exp2(BlockSize); //number of diff that are encoded together
	int IndexAux[32]; //(the last location of a chrm)>>SearchBlock, assume no more than 32 chrm
	memset(IndexAux,0,32*sizeof(unsigned int));
	int bcount=0; //counter inside each block
	//int BlockCount=0; //counter for number of blocks
	
	
	//read from count files
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"Diff");
	FILE *fpDiff = fopen(tmpFileName,"r");
	if (fpDiff == NULL) {
		fprintf(stderr,"Can't open output Count file!\n");
		exit(1);
	}
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"LenDiff");
	FILE *fpLenDiff = fopen(tmpFileName,"r");
	if (fpLenDiff == NULL) {
		fprintf(stderr,"Can't open output LenDiff file!\n");
		exit(1);
	}
	//skip annotation
	fgets(str,sizeof(str),fpDiff);
	fgets(str,sizeof(str),fpLenDiff);
	//change consecutive alphabet to the actual alphabet
	i=0;
	while (fread(&Map[i],sizeof(double),1,fpDiff)==1){
		//fscanf(fpDiff,"%d\t%*u\n", &Map[i]);
		fread(&tmpint,sizeof(int),1,fpDiff);
		i ++;
	}
	i=0;
	while (fread(&LenMap[i],sizeof(double),1,fpLenDiff)==1){
		//fscanf(fpLenDiff,"%d\t%*u\n", &LenMap[i]);
		fread(&tmpint,sizeof(int),1,fpLenDiff);
		i ++;
	}
			
		
	fclose(fpDiff);
	fclose(fpLenDiff);
	/*printf("map done!\n");
	for (j=0; j<80; j++)
		printf("j=%d, %d %d\n",j, Map[j],LenMap[j]);*/
	
	
	///////////////matlab seq to seq	
	//open matlab seq
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"DiffSeqMatlabDecode");
	if (nproc>0){
		tmpint = sprintf(str,"%02d",iproc);
		strcat(tmpFileName,str);
	}
	FILE *fpDiffSeqMatlab = fopen(tmpFileName, "r");
	if (fpDiffSeqMatlab  == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}
	//open matlab len seq
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"LenDiffSeqMatlabDecode");
	if (nproc>0){
		tmpint = sprintf(str,"%02d",iproc);
		strcat(tmpFileName,str);
	}
	FILE *fpLenDiffSeqMatlab = fopen(tmpFileName, "r");
	if (fpLenDiffSeqMatlab  == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}
	//output wig
	strcpy(tmpFileName,argv[2]);
	if (nproc>0){
		tmpint = sprintf(str,"%02d",iproc);
		strcat(tmpFileName,str);
	}
	FILE *fp;
	fp = fopen(tmpFileName, "w");
	if (fp == NULL) {
		fprintf(stderr, "Can't open output file %s!\n",tmpFileName);
		exit(1);
	}
	
	//chrm names
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"ChrmName");
	FILE *fpName = fopen(tmpFileName, "r");
	if (fpName == NULL) {
		fprintf(stderr, "Can't open input file %s!\n",tmpFileName);
		exit(1);
	}
    while(!feof(fpName)){
        fscanf(fpName,"%8s",&ChrmNames[8*ChrmCount]);
        ChrmCount++;
    }
	
	
	int StartBlockIndex=0;
	//int EndBlockIndex=0;
	
	//fpBlockAux
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"BlockAux");
	FILE *fpBlockAux = fopen(tmpFileName, "r");
	if (fpBlockAux == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}
	for(i=0; i<First+1 && !feof(fpBlockAux); i++)
		fscanf(fpBlockAux,"%d ",&tmp);
	StartBlockIndex = tmp ;
	/*for (i; i<Last+1 && !feof(fpBlockAux); i++)
		fscanf(fpBlockAux,"%d ",&tmp);
	EndBlockIndex = tmp + 1;*/
	fclose(fpBlockAux);
	
	
	j=0;
	i = First; //Chrm count
	pcount = 1024; //count to print annotation
	bcount = 0;
	Diff=1;
	LenDiff=1;
	while(!feof(fpDiffSeqMatlab) && !feof(fpLenDiffSeqMatlab)){
			
		/*if (Location==3839085){ //(j<100){ 
			printf("PreLocation=%lu Location=%lu Contig=%ld LenDiff=%ld Value=%d Diff=%d\n", PreLocation,Location, Contig, LenDiff,Value,Diff);
		    printf("bcount=%d\n",bcount);
		}*/
		
		if (bcount==SetBlockSize)
			bcount = 0;
		bcount++;
		
		if (pcount==1024){
			fprintf(fp,"variableStep chrom=%s span=1\n", &ChrmNames[i*8]);
			pcount = 0;
			lcount++;
		}
		fscanf(fpDiffSeqMatlab,"%d ",&tmpint);
		Diff = Map[tmpint];
		fscanf(fpLenDiffSeqMatlab,"%d ",&tmpint);
		LenDiff = LenMap[tmpint];
	 
		
						
		if (Diff == -(1<<30)){//end of a chrm
			i++;
			pcount = 1024;
			PreValue = -2;
			PreContig = 1;
			PreLocation = -1;
			//printf("chrm %d finished\n",i);
			//remove padded 0s
			tmpint=(SetBlockSize-bcount)%SetBlockSize; //number of padded 0s
			for (j=0; j<tmpint; j++){
				fscanf(fpDiffSeqMatlab,"%d ",&tmp);
				fscanf(fpLenDiffSeqMatlab,"%d ",&tmp);
			}
			bcount=0;
			continue;
		}
				
		//write to wig	
		Value = PreValue + Diff;
		Contig = PreContig + LenDiff;
		
		
		if (Contig < 0){
			printf("Contig=%lg out of range! Check input file please!\n", Contig);
			printf("Location=%lg PreLocation=%lg Value=%lg PreValue=%lg ChrmName=%s\n",Location,PreLocation,Value,PreValue,&ChrmNames[i*8]);
			printf("Line # in wig is %lld\n",lcount);
			exit(1);
		}
			
		if ((Value>1e-8 || Value<-1e-8) && PreLocation!=-1){
			for (Location=PreLocation; Location < PreLocation + Contig; Location++ ){
				//fprintf(fp,"%ld\t%lg\n", Location,(double)(Value)/SetFloatSize);
				fprintf(fp,"%ld\t%lg\n", (long int)Location,Value);
				lcount ++;
				pcount ++;
			}
		}
		else
			Location = PreLocation + Contig;
		
		//reset
		PreValue = Value;
		PreContig = Contig;
		PreLocation = Location;
		//Diff = NextDiff;
		
	}
	clock_t end = clock();

    FILE *fpResult = fopen("result","a");
	if (fpResult == NULL) {
		fprintf(stderr,"Can't open output Count file!\n");
		exit(1);
	}
    SPAMR((fpResult,"//%s to_wig\n",FileName));
	SPAMR((fpResult,"Time for processing file in decode %2.3f sec.\n", (double)(end-start)/CLOCKS_PER_SEC));
	
	if (!feof(fpDiffSeqMatlab) || !feof(fpLenDiffSeqMatlab)){
		printf("two matlab seq files not same sizes!\n");
		exit(1);
	}
	
	fclose(fp);
	fclose(fpDiffSeqMatlab);
	fclose(fpLenDiffSeqMatlab);
	fclose(fpName);
	fclose(fpResult);

	return 1;
}
