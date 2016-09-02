/*
 * Change a matlab seq back to wig format
 * argv stores filenames
 * for subsequence query, Qchrm,Qstart,Qend are the query chrm name, start location, end location 
*/

int to_wig_sub(char * argv[], char *Qchrm, int Qstart, int Qend){
	
	if (Qstart<0 || Qend<Qstart){
		printf("must have Qend > Qstart >=0!\n");
		exit(1);
	}
		
	clock_t start = clock();
	SPAM(("================\nto_wig_sub.c start!\n"));
	
	char FileName[400];//name of processed Wig file
	strcpy(FileName, argv[1]);
	
	
	int Nd = 18;
	int SetMaxDiff = exp2(Nd); //set max alphabet size
	unsigned int DiffCount; //frequncy of each difference value of contigs
	unsigned int LenDiffCount;
	int Map[SetMaxDiff]; //maps each consecutive integer to the actual alphabet
	memset(Map,0,SetMaxDiff*sizeof(int)); 
	int LenMap[SetMaxDiff]; 
	memset(LenMap,0,SetMaxDiff*sizeof(int));  
	char ChrmNames[8*32]; //every chrm has a name less than 8 characters, assume no more than 32 chrm 
	memset(ChrmNames,0,8*32*sizeof(char));
	int ChrmCount=0;
	
	//int type;
	//char Type[40];
	
	int PreValue=-2;
	int Value=0;
	int Diff;//difference of values of contigs
	//int NextDiff=1;
	long unsigned int PreLocation=0;
	long unsigned int Location=0;
	long int Contig = 1; //length of current contig
	long int PreContig = 1; //length of previous contig
	long int LenDiff; //length difference of two consecutive contigs
	char str[200];
	char tmpFileName[400];
	int i,j;
	int tmpint;
	int pcount=0;
	long long lcount=0; //count number of lines in wig file
	
	//for subsequence query
	int Qnum=-1; // query chrm number
	
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
	int tmpSub = 0;
	int tmpLenSub = 0;
	/*for (i=0; i<SetMaxDiff; i++){
		fscanf(fpDiff,"%d\t%u\n", &Diff,&DiffCount);
		if (DiffCount==0)
			tmpSub++;
		else
			Map[i-tmpSub]=Diff;
		fscanf(fpLenDiff,"%ld\t%u\n", &LenDiff,&LenDiffCount);
		if (LenDiffCount==0)
			tmpLenSub++;
		else
			LenMap[i-tmpLenSub]=LenDiff;
	}*/
	i=-1;
	while (!feof(fpDiff)){
		fscanf(fpDiff,"%d\t%u\n", &Diff,&DiffCount);
		tmpSub += mydiff(Diff) - i -1;
		Map[mydiff(Diff)-tmpSub]=Diff;
		i = mydiff(Diff);
	}
	i=-1;
	while (!feof(fpLenDiff)){
		fscanf(fpLenDiff,"%ld\t%u\n", &LenDiff,&LenDiffCount);
		tmpLenSub += mydiff(LenDiff) - i -1;
		LenMap[mydiff(LenDiff)-tmpLenSub]=LenDiff;
		i = mydiff(LenDiff);
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
	FILE *fpDiffSeqMatlab = fopen(tmpFileName, "r");
	if (fpDiffSeqMatlab  == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}
	//open matlab len seq
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"LenDiffSeqMatlabDecode");
	FILE *fpLenDiffSeqMatlab = fopen(tmpFileName, "r");
	if (fpLenDiffSeqMatlab  == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}
	//output wig
	strcpy(tmpFileName,argv[2]);
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
        if (memcmp(&ChrmNames[8*ChrmCount],Qchrm,8)==0)
			Qnum = ChrmCount;			
        ChrmCount++;
    }
    if (Qnum==-1){
		printf("No query chrom with name %s found!\n",Qchrm);
		exit(1);
	}
	
	
	// exceed
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"Exceed");
	FILE *fpExceed = fopen(tmpFileName, "r");
	if (fpExceed == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}
	// lenexceed
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"LenExceed");
	FILE *fpLenExceed = fopen(tmpFileName, "r");
	if (fpLenExceed == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}
	
	j=0;
	i = 0; //Chrm count
	pcount = 1024; //count to print annotation
	Diff=1;
	LenDiff=1;
	
	while(!feof(fpDiffSeqMatlab) && !feof(fpLenDiffSeqMatlab)){
		/*if (Location==16458){ //(j<100){ 
			printf("Location=%lu Contig=%ld LenDiff=%ld Value=%d Diff=%d\n", Location, Contig, LenDiff,Value,Diff);
		    j++;
		}*/
		
		if (pcount==1024 && Qnum==i){
			fprintf(fp,"variableStep chrom=%s span=1\n", &ChrmNames[i*8]);
			pcount = 0;
		}
		fscanf(fpDiffSeqMatlab,"%d ",&tmpint);
		Diff = Map[tmpint];
		fscanf(fpLenDiffSeqMatlab,"%d ",&tmpint);
		LenDiff = LenMap[tmpint];
	
						
		if (mydiff(Diff) == SetMaxDiff-1){//end of a chrm
			i++;
			pcount = 1024;
			PreValue = -2;
			PreContig = 1;
			PreLocation = 0;
			//printf("chrm %d finished\n",i);
			if (Qnum >= i)
				continue;
			else  //Qnum<i
				break;
		}
		if (i==Qnum){
			//exceed max value
			if (mydiff(Diff) == SetMaxDiff-2)
				fscanf(fpExceed, "%d ", &Diff);
			if (mydiff(LenDiff) == SetMaxDiff-2){
				fscanf(fpLenExceed, "%ld ", &LenDiff);
				//printf("Exceed LenDiff=%ld\n",LenDiff);
			}
			
			//write to wig	
			Value = PreValue + Diff;
			Contig = PreContig + LenDiff;
			
			
			if (Contig < 0 || (Contig>SetMaxDiff && Value)){
				printf("Contig=%lu out of range! Check input file please!\n", Contig);
				printf("Location=%lu PreLocation=%lu Value=%d PreValue=%d ChrmName=%s\n",Location,PreLocation,Value,PreValue,&ChrmNames[i*8]);
				printf("Line # in wig is %lld\n",lcount);
				exit(1);
			}
				
			if (Value>0){
				for (Location=PreLocation; Location < PreLocation + Contig; Location++ ){
					if (Location>= Qstart && Location <= Qend){
						fprintf(fp,"%lu\t%d\n", Location,Value);
						lcount ++;
						pcount ++;
					}
				}
			}
			Location = PreLocation + Contig;
			if (Location > Qend)
				break;
			//reset
			PreValue = Value;
			PreContig = Contig;
			PreLocation = Location;
			//Diff = NextDiff;
		}
		
	}
	clock_t end = clock();
	
	
    FILE *fpResult = fopen("result","a");
	if (fpResult == NULL) {
		fprintf(stderr,"Can't open output Count file!\n");
		exit(1);
	}
    SPAMR((fpResult,"//%s to_wig\n",FileName));
	SPAMR((fpResult,"Time for processing file in decode %2.3f sec.\n", (double)(end-start)/CLOCKS_PER_SEC));
	
/*	if (!feof(fpDiffSeqMatlab) || !feof(fpLenDiffSeqMatlab)){
		printf("two matlab seq files not same sizes!\n");
		exit(1);
	}
*/	
	fclose(fp);
	fclose(fpDiffSeqMatlab);
	fclose(fpLenDiffSeqMatlab);
	fclose(fpExceed);
	fclose(fpLenExceed);
	fclose(fpName);
	fclose(fpResult);
	
	return 0;
}

