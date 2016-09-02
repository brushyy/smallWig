//
//  getseqblock.c
//
//  Created by Zhi Ying on 12/5/13.
//  Copyright (c) 2013 Zhi Ying. All rights reserved.
//
//	Get a sequence of difference values, and difference run lengths, together with their counts
//  To use them in matlab, the diff sequences are converted to values within 1~size(DiffCount)
//  Each pair (Diff,LenDiff) corresponds to the difference of a run
//	First pair is assumed to be (Value=-1,Contig=1), and first location is assumed to be 1.
//
//  Also to support random access, encode blocks of 2^BlockSize diff together
//
//  fpStart stores the starting location (=PreLocation) of each block
//  fpPre stores the prevalue, and precontig at the beginning of a block; 
//  fpSmry stores min, max, mean=sum(values), std=sum(values^2), coverage=(#non-NaN basis), totalloc=(# basis)
//  fpSmryAll stores for every chromosome: min, max, mean=sum(values)/coverage, std=sum(values^2)/coverage, coverage=(#non-NaN basis), totalloc=(# basis)
//  fpBlockAux stores the accumulated number of code blocks so far before each chrm (first item=0, last item=total # blocks, #items = #chrm+1)
//  fpChrmAdd stores the address of beginning of each chrm, to be found in arithematic coding, also uses info from fpBlockAux (1st item=0, last item=total file size, #items = #chrm+1)
//  fpBlockAddress, address of coded data for each block, to be found in arithematic coding (1st item=0, last item=total file size, #items = #blocks+1)
//
//  Also include the function that when a totoal of n processors is used, only process 24/n 
//  if nproc<=0, then no parallel


//record and update summary info
void summary(double *min, double *max, double *mean, double *std, double *coverage, double *totalloc, double *smry, int i, FILE *fpSmry){
	fwrite(min, sizeof(double),1,fpSmry);
	fwrite(max, sizeof(double),1,fpSmry);
	fwrite(mean, sizeof(double),1,fpSmry);
	fwrite(std, sizeof(double),1,fpSmry);
	fwrite(coverage, sizeof(double),1,fpSmry);
	fwrite(totalloc, sizeof(double),1,fpSmry);
	if (*min<smry[6*i])
		smry[6*i] = *min;
	if (*max>smry[6*i+1])
		smry[6*i+1] = *max;
	smry[6*i+2] = (smry[6*i+2]*smry[6*i+4]+ *mean)/(smry[6*i+4]+ *coverage);
	smry[6*i+3] = (smry[6*i+3]*smry[6*i+4]+ *std)/(smry[6*i+4]+ *coverage);
	smry[6*i+4] += *coverage;
	smry[6*i+5] += *totalloc; 
	//reset
	*min=1<<30;
	*max=-(1<<30);
	*mean=0;
	*std=0;
	*coverage=0;	
	*totalloc=0;
}


//record an item
void record(
double Value, double PreValue, 
double PreLocation, double Location, 
int SetMaxDiff, int SetBlockSize, 
double *Diff, double *NextDiff, 
double *Contig, double *PreContig, 
double *LenDiff,
int *bcount, int *BlockCount,
unsigned int *DiffCount, unsigned int *LenDiffCount,
double *DiffTable, double *LenDiffTable, int *DiffSize, int *LenDiffSize,
FILE *fpLenDiffSeq, FILE *fpDiffSeq,
FILE *fpStart, FILE *fpPre, FILE *fpSmry,
long unsigned int *BlockStart, long unsigned int *NextBlockStart, 
double *min, double *max, double *mean, double *std, double *coverage, double *totalloc,
double *smry, int ChrmCount){
	

	*Diff = *NextDiff; 
	*NextDiff = Value-PreValue;	
	*LenDiff = (double)(*Contig)-*PreContig;
	
	(*bcount)++;
	if (*bcount==1){//beginning of a block
		(*BlockCount)++;
	}
	if (*bcount==SetBlockSize){//end of a block
		*bcount=0;
		*BlockStart=*NextBlockStart;
		*NextBlockStart=PreLocation+1;

		//store statistics for current block
		summary(min,max,mean,std,coverage,totalloc,smry, ChrmCount, fpSmry);
		//store prevalue, change-value location, contig for next block use
		double tmp; //change-value locatino
		if (Value==0)
			tmp = PreLocation+1;
		else
			tmp = Location;
		fwrite(&PreValue, sizeof(double),1,fpPre);
		fwrite(Contig, sizeof(double),1,fpPre);
		fwrite(&tmp, sizeof(double),1,fpStart);
		//if (*BlockCount < 10) printf("Start Location %d\n",tmp);
	}

	//find Diff and LenDiff in the frequency talble
	int i, j, indicator;
	
	i = binfind(*DiffSize,DiffTable,*Diff,&indicator);
	if (indicator != 0){//not found then insert
		(*DiffSize) ++;
		if (*DiffSize > SetMaxDiff){
			fprintf(stderr,"Too large alphabet size of diff!\n");
			exit(1);
		}
		for (j=*DiffSize-2; j>=i; j--){
			DiffTable[j+1] = DiffTable[j];
			DiffCount[j+1] = DiffCount[j];
		}
		DiffTable[i] = *Diff;
		DiffCount[i] = 1;
	}
	else
		DiffCount[i] ++;
	/*printf("Diff=%lg\tTable\tCount\n",*Diff);
	for (i=0; i<*DiffSize; i++)
		printf("%lg\t%d\n",DiffTable[i],DiffCount[i]);*/
		
	//lendiff
	i = binfind(*LenDiffSize,LenDiffTable,*LenDiff,&indicator);
	if (indicator != 0){//not found then insert
		(*LenDiffSize) ++;
		if (*LenDiffSize > SetMaxDiff){
			fprintf(stderr,"Too large alphabet size of lendiff!\n");
			exit(1);
		}
		for (j=*LenDiffSize-2; j>=i; j--){
			LenDiffTable[j+1] = LenDiffTable[j];
			LenDiffCount[j+1] = LenDiffCount[j];
		}
		LenDiffTable[i] = *LenDiff;
		LenDiffCount[i] = 1;
	}
	else
		LenDiffCount[i] ++;
	
	//record
	//fprintf(fpLenDiffSeq,"%ld ", *LenDiff);
	//fprintf(fpDiffSeq,"%d ",*Diff);
	fwrite(LenDiff,sizeof(double),1,fpLenDiffSeq);
	fwrite(Diff,sizeof(double),1,fpDiffSeq);
	
	
	*PreContig = *Contig; //reset
	*Contig = 1;

}

//////////////////// main////////////////////////////////////////////////////
int getseqblock(char * argv[],int nproc,int iproc)
{
	clock_t start = clock();
	//printf("let's begin\n");
	//printf("sizeof(long unsinged int)=%ld,sizeof(int)=%ld\n",sizeof(long unsigned int),sizeof(int));
	char FileName[400];//name of output swig file

	strcpy(FileName, argv[1]);
	
	if (BlockSize >=8 && BlockSize <=32){
		SPAM(("BlockSize=%d\n",BlockSize));
	}
	else{
		printf("Wrong range of block size or search size! Should be from 8 to 32!\n");
		exit(1);
	}
	
	//for parallel
	int First; //first chrm to read
	int Last;
	int pflag; //indicate whether should start read
	if (nproc<=0){
		First = 0;
		Last = 32;//assume no more than 32 chrm
		pflag=1;
	}
	else{
		First = (24/nproc)*iproc; 
		if (iproc < nproc-1)
			Last = (24/nproc)*(iproc+1); 
		else
			Last = 32;
		pflag=0; 
	}
	
	
	//int N = 16;
	int Nd = 18;
	//int SetMaxLength = exp2(Nd-1); //set max contig length
	int SetMaxDiff = exp2(Nd); //set max value difference 
	int SetBlockSize = exp2(BlockSize); //number of diff that are encoded together
	//int SearchInd[(int)(1)<<IndexSize]; //auxilary table indicating how many times a code block is used for search, assume average contig length is 2^4, has the same number of elements as OffSet
	//memset(SearchInd, 0, ((int)(1)<<IndexSize)*sizeof(int));
	char ChrmNames[8*32]; //every chrm has a name less than 8 characters, assume no more than 32 chrm 
	memset(ChrmNames,0,8*32*sizeof(char));
	char *tmpChrmName;
	unsigned int DiffCount[SetMaxDiff]; //frequncy of each difference value of contigs
	memset(DiffCount,0,SetMaxDiff*sizeof(unsigned int));
	unsigned int LenDiffCount[SetMaxDiff]; //frequncy of each difference length of contigs
	memset(LenDiffCount,0,SetMaxDiff*sizeof(unsigned int));
	double DiffTable[SetMaxDiff]; //difference alphabet table
	memset(DiffTable,0,SetMaxDiff*sizeof(double));
	double LenDiffTable[SetMaxDiff]; //len diff alphabet table
	memset(LenDiffTable,0,SetMaxDiff*sizeof(double));
	int DiffSize=2, LenDiffSize=2; //Diff and LenDiff alphabet sizes, Talbe[0]=end_marker, Table[1]=exceed_marker
	DiffTable[0] = -(1<<30); //end marker
	LenDiffTable[0] = -(1<<30);
	DiffTable[1] = -(1<<29);//exceed marker
	LenDiffTable[1] = -(1<<29);

	
	double PreValue=-1;
	double Value=0;
	double Diff=1;//difference of values of contigs
	//int PreDiff=1;
	double NextDiff=1;
	double PreLocation=-1;
	double Location=-1;
	long unsigned int MaxLocation =10000;//the largest location. Not necessarily the last one, because the chromozones are not in order 
	double Contig = 1; //length of current contig
	double PreContig = 1; //length of previous contig
	double LenDiff=1; //length difference of two consecutive contigs
	//long int PreLenDiff=1; //length diff of previous two contigs
	int bcount=0; //counter inside each block
	int BlockCount=0; //counter for number of blocks
	long unsigned int BlockStart; //starting chrm location of a encode block
	long unsigned int NextBlockStart; // previous BlockStart
	int SearchInd; //number of search blocks in each code block
	//long unsigned int ChrmStart[100]; //assume at most 100 chromosomes, ChrmStart[i] stores the starting location of chrm i
	int ChrmCount=-1; //temp chrm count
 	//int ChrmNum=0; //number of chrms
	char tmp;
	char str[1000];  //annotation line
	char Prestr[1000]; //previous annotation line
    char tmpFileName[400];	
    fpos_t position;
	int i,j=0;
	int k=0;
	int tmpcount=0;
	int tmpint=0;
	int tmpflag=0;
	double tmpdouble=0;
	int warnflag=0; //if warning of float size has been printed
	
	//statistics of each block
	double min=1<<30, max=-(1<<30), mean=0; 
	double coverage=0, totalloc=0; //min and max value, sum of non-NaN basis, number of covered base pairs, number of total base pairs
	double std=0; //sum of squre of non-NaN basis
	double smry[32*6];//stats of every chrm
	for (i=0; i<32*6; i++)
		smry[i] = 0; //all set of 0
	for (i=0; i<32; i++){
		smry[6*i] = 1<<30; //min is a large number
		smry[6*i+1]= -(1<<30); //max is set to be small
	}

	
	
	//printf("let's begin\n");
	//printf("test Count[0~4] %lu %lu %lu %lu %lu", Count[0],Count[1],Count[2],Count[3],Count[4]);
	//printf("sizeof Location= %lu\n", sizeof(Location));

	//input
	strcpy(tmpFileName,FileName);
	FILE *fp;
	fp = fopen(tmpFileName, "r");
	if (fp == NULL) {
		fprintf(stderr, "Can't open input file %s!\n",tmpFileName);
		exit(1);
	}
	//output diffseq
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"DiffSeq");
	if (nproc>0){
		tmpint = sprintf(str,"%02d",iproc);
		strcat(tmpFileName,str);
	}
	FILE *fpDiffSeq = fopen(tmpFileName, "w+");
	if (fpDiffSeq == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}
	//output lendiffseq
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"LenDiffSeq");
	if (nproc>0){
		tmpint = sprintf(str,"%02d",iproc);
		strcat(tmpFileName,str);
	}
	FILE *fpLenDiffSeq = fopen(tmpFileName, "w+");
	if (fpLenDiffSeq == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}

	
	//output start, for preloc, 
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"Start");
	if (nproc>0){
		tmpint = sprintf(str,"%02d",iproc);
		strcat(tmpFileName,str);
	}
	FILE *fpStart = fopen(tmpFileName, "wb");
	if (fpStart  == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}
	//output prevalue, (Pre)contig
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"Pre");
	if (nproc>0){
		tmpint = sprintf(str,"%02d",iproc);
		strcat(tmpFileName,str);
	}
	FILE *fpPre = fopen(tmpFileName, "wb");
	if (fpPre == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}
	//output summary info: min, max, mean, std, coverage, totalloc
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"Smry");
	if (nproc>0){
		tmpint = sprintf(str,"%02d",iproc);
		strcat(tmpFileName,str);
	}
	FILE *fpSmry = fopen(tmpFileName, "wb");
	if (fpSmry  == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}
	//output blockaux
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"BlockAux");
	if (nproc>0){
		tmpint = sprintf(str,"%02d",iproc);
		strcat(tmpFileName,str);
	}
	FILE *fpBlockAux = fopen(tmpFileName, "w");
	if (fpBlockAux  == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}
	//put a dummy 0 at the beginning
	fprintf(fpBlockAux,"0 ");
	
	
	//output result
	FILE *fpResult=fopen("result","a");
	if (fpResult == NULL) {
		fprintf(stderr,"Can't open output result file!\n");
		exit(1);
	}
	SPAMR((fpResult,"\n//combine/getseqblock.c\n//%s\n",FileName));
	
	
	
	
	SPAM(("===============================\ngetseqblock.c\n"));
	
	//start processing input file
	strcpy(Prestr, "empty");
	while (!feof(fp)){		
		tmp = fgetc(fp);
		if  (tmp >= '0' && tmp <= '9'){
			if (!pflag)
				fgets(str,sizeof(str),fp);
			else{
				ungetc(tmp,fp);
				fgetpos (fp, &position);
				fscanf(fp,"%lf %lf\n", &Location,&Value);
				//Value = (round)(tmpdouble*SetFloatSize); //in case Value is stored in a floating point manner
				/*if (Value - tmpdouble*SetFloatSize > (double)1/SetfloatSize){
					if (warnflag == 0){
						fprintf(stderr,"WARNING: THE FILE CONTAINS VALUE=%g WITH PRECISION HIGHER THAN FLOATSIZE=%d! INCREASE FLOATSIZE! diff=%g\n",tmpdouble,(int)(log10(SetFloatSize)),Value - tmpdouble*SetFloatSize);
						warnflag=1;
					}
				}*/
				
				if (Location > MaxLocation)
					MaxLocation = Location;
				if (Location != PreLocation+1){//for contig before 0s
					Value = 0;
				}
				if (Value != PreValue){
					//stats
					if (PreValue != 0 && PreLocation > 0){
						if (PreValue < min)
							min = PreValue;
						if (PreValue > max)
							max = PreValue;
						mean += PreValue*Contig;
						std += (double)PreValue*PreValue*Contig;
						coverage += Contig;
					}
					totalloc += Contig;  
					/*if (k<5){
						printf("line %d: Value %d, PreLocation %d, Location %d, Contig %d min1 %d, max1 %d, mean1 %d, std1 %ld, coverage1 %d, totoalloc1 %d\n",
						        k,       Value, PreLocation, Location, Contig, min,max,mean,std,coverage,totalloc);
					}
					k++;*/
					//record
					record(Value, PreValue, PreLocation, Location, SetMaxDiff, SetBlockSize, 
						&Diff, &NextDiff, &Contig, &PreContig, &LenDiff, 
						&bcount, &BlockCount, DiffCount, LenDiffCount, DiffTable, LenDiffTable, &DiffSize, &LenDiffSize, fpLenDiffSeq, fpDiffSeq, 
						fpStart, fpPre, fpSmry, &BlockStart, &NextBlockStart, 
						&min, &max, &mean, &std, &coverage, &totalloc, smry, ChrmCount);
				}
				else{//Diff == 0
						Contig++;
				}
				
				if (Location != PreLocation+1){
					//record skipped 0s, push back the line
					fsetpos (fp, &position);
					PreValue = 0;
					if (Location<PreLocation){
						printf("Location %lg < PreLocation %lg!\n",Location, PreLocation);
						exit(1);
					}
					Contig = Location - PreLocation-1;
					PreLocation = Location-1;		
				}				
				else{
					PreValue = Value;
					PreLocation = Location;			
				}		
			}
		}
		else{ //get a line of annotations
			ungetc(tmp,fp);
			fgets(str,sizeof(str),fp); 
			if (strcmp(str, Prestr)){ //different chrm
				//old chrm
				if (ChrmCount>=0 && pflag){ 
					//insert a virtual data point and then record
					PreValue = Value;
					Value++;
					PreLocation = Location;
					Location++;
					record(Value, PreValue, PreLocation, Location, SetMaxDiff, SetBlockSize, 
						&Diff, &NextDiff, &Contig, &PreContig, &LenDiff, 
						&bcount, &BlockCount, DiffCount, LenDiffCount, DiffTable, LenDiffTable, &DiffSize, &LenDiffSize, fpLenDiffSeq, fpDiffSeq, 
						fpStart, fpPre, fpSmry, &BlockStart, &NextBlockStart, 
						&min, &max, &mean, &std, &coverage, &totalloc, smry, ChrmCount);
					//end marker of a chrom
					//fprintf(fpLenDiffSeq,"%d ", -(1<<30));
					//fprintf(fpDiffSeq,"%d ",-(1<<30));
					tmpdouble = -(1<<30);
					fwrite(&tmpdouble,sizeof(double),1,fpLenDiffSeq);
					fwrite(&tmpdouble,sizeof(double),1,fpDiffSeq);
	
					bcount++;
					if (bcount==1){//beginning of a block
						BlockCount++;
					}
					//pad 0s to fill a code block
					tmpint=(SetBlockSize-bcount)%SetBlockSize; //number of padded 0s
					for (j=0; j<tmpint; j++){
						//fprintf(fpDiffSeq,"1 ");
						//fprintf(fpLenDiffSeq,"0 ");
						tmpdouble=0;
						fwrite(&tmpdouble,sizeof(double),1,fpLenDiffSeq);
						tmpdouble = 1;
						fwrite(&tmpdouble,sizeof(double),1,fpDiffSeq);
					}
					DiffCount[mydiff(1)]+=tmpint;
					LenDiffCount[mydiff(0)]+=tmpint;
					if(tmpint){//search index of the last block
						BlockStart=NextBlockStart;
						NextBlockStart=Location;
					}
					//number of blocks so far
					fprintf(fpBlockAux,"%d ",BlockCount);
					//stats of the last block of the chrm
					summary(&min,&max,&mean,&std,&coverage,&totalloc,smry, ChrmCount, fpSmry);
				}	
				//new chrm
				//printf("chrm %d, %s", ChrmCount, str);
				strcpy(Prestr, str);
				ChrmCount++;
				if (ChrmCount>=32){
					printf("Too many chromosomes!\n");
					exit(1);
				}
				if (ChrmCount == Last)
					break;
				if (ChrmCount == First)
					pflag = 1;	
				if (pflag){
					tmpChrmName = strstr(str, "chrom=");
					if (tmpChrmName == NULL){
						printf("Wrong annotation line: %s, must contain chromosome name!\n", str);
						exit(1);
					}
					j=6;
					while (tmpChrmName[j]!=' ' && j<14){
						ChrmNames[8*ChrmCount+j-6]=tmpChrmName[j];
						j++;
					}
					
					//first (PreValue,PreContig,PreLocation) for each chrm
					
					tmpdouble=-2;//prevalue
					fwrite(&tmpdouble, sizeof(double),1,fpPre);
					tmpdouble=1;//precontig
					fwrite(&tmpdouble, sizeof(double),1,fpPre);
					tmpdouble=-1;//prelocation
					fwrite(&tmpdouble, sizeof(double),1,fpStart);
					//printf("fpstart first elt %d\n",tmpstart);
					
					
					//starting index
					//reset parameters
					PreValue=-1;
					Value=0;
					Diff=1;
					NextDiff=1;
					PreLocation=-1;
					Location=-1;
					Contig = 1; 
					PreContig = 1; 
					LenDiff=1;
					bcount=0;
					BlockStart=0;
					NextBlockStart=0; 
				}
			}
		}
	}
	//ending part of last chrm, only process it if nonparelle, or if it is the last processor
	if (Location>=0 && (nproc<=0 || (nproc>0 && iproc==nproc-1))){ 
		//printf("ending of chrm %d\n",ChrmCount);
		//insert a virtual data point and then record
		PreValue = Value;
		Value++;
		PreLocation = Location;
		Location++;
		record(Value, PreValue, PreLocation, Location, SetMaxDiff, SetBlockSize, 
			&Diff, &NextDiff, &Contig, &PreContig, &LenDiff, 
			&bcount, &BlockCount, DiffCount, LenDiffCount, DiffTable, LenDiffTable, &DiffSize, &LenDiffSize, fpLenDiffSeq, fpDiffSeq, 
			fpStart, fpPre, fpSmry, &BlockStart, &NextBlockStart, 
			&min, &max, &mean, &std, &coverage, &totalloc, smry, ChrmCount);
		//end marker of a chrom, Table[0]=marker=-(1<<30);
		//fprintf(fpLenDiffSeq,"%d ", -(1<<30));
		//fprintf(fpDiffSeq,"%d ",-(1<<30));
		tmpdouble = -(1<<30);
		fwrite(&tmpdouble,sizeof(double),1,fpLenDiffSeq);
		fwrite(&tmpdouble,sizeof(double),1,fpDiffSeq);
		bcount++;
		if (bcount==1){//beginning of a block
			BlockCount++;
		}
		//pad 0s to fill a code block
		tmpint=(SetBlockSize-bcount)%SetBlockSize; //number of padded 0s
		for (j=0; j<tmpint; j++){
			//fprintf(fpDiffSeq,"1 ");
			//fprintf(fpLenDiffSeq,"0 ");
			tmpdouble=0;
			fwrite(&tmpdouble,sizeof(double),1,fpLenDiffSeq);
			tmpdouble = 1;
			fwrite(&tmpdouble,sizeof(double),1,fpDiffSeq);
		}
		DiffCount[mydiff(1)]+=tmpint;
		LenDiffCount[mydiff(0)]+=tmpint;
		if(tmpint){//search index of the last block
			BlockStart=NextBlockStart;
			NextBlockStart=Location;
		}
		//number of blocks so far
		fprintf(fpBlockAux,"%d ",BlockCount);
		//stats of the last block of the chrm
		summary(&min,&max,&mean,&std,&coverage,&totalloc,smry, ChrmCount, fpSmry);
	}
	ChrmCount++;
	if (Last>ChrmCount)
		Last = ChrmCount;
	
	//open output files 1
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"Diff");
	if (nproc>0){
		tmpint = sprintf(str,"%02d",iproc);
		strcat(tmpFileName,str);
	}
	FILE *fpDiff = fopen(tmpFileName,"w");
	if (fpDiff == NULL) {
		fprintf(stderr,"Can't open output Count file!\n");
		exit(1);
	}
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"LenDiff");
	if (nproc>0){
		tmpint = sprintf(str,"%02d",iproc);
		strcat(tmpFileName,str);
	}
	FILE *fpLenDiff = fopen(tmpFileName,"w");
	if (fpLenDiff == NULL) {
		fprintf(stderr,"Can't open output Count file!\n");
		exit(1);
	}
	//write to output files 1
	if (nproc<0 || iproc==0){//if non parallel or if first processor
		DiffCount[1] = 1;//exceed marker
		LenDiffCount[1] = 1; 
		DiffCount[0] = 1;//end marker of a chrm
		LenDiffCount[0] = 1; 
	}
	else{
		DiffCount[1] = 0;//exceed marker
		LenDiffCount[1] = 0; 
		DiffCount[0] = 0;//end marker of a chrm
		LenDiffCount[0] = 0; 
	}

	
	fprintf(fpDiff,"#Difference\tCount\n");
	fprintf(fpLenDiff,"#Difference of Length\tCount\n");
	for (i=0; i<DiffSize; i++){
		//fprintf(fpDiff,"%d\t%u\n", DiffTable[i],DiffCount[i]);
		fwrite(&DiffTable[i],sizeof(double),1,fpDiff);
		fwrite(&DiffCount[i],sizeof(int),1,fpDiff);
	}
	for (i=0; i<LenDiffSize; i++){
		//fprintf(fpLenDiff,"%d\t%u\n", LenDiffTable[i],LenDiffCount[i]);
		fwrite(&LenDiffTable[i],sizeof(double),1,fpLenDiff);
		fwrite(&LenDiffCount[i],sizeof(int),1,fpLenDiff);
	}
	
	//output file chrmname
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"ChrmName");
	if (nproc>0){
		tmpint = sprintf(str,"%02d",iproc);
		strcat(tmpFileName,str);
	}
	FILE *fpName = fopen(tmpFileName,"w");
	if (fpName == NULL) {
		fprintf(stderr,"Can't open output file %s!\n",tmpFileName);
		exit(1);
	}
	//for (i=0; i<ChrmCount; i++){
	for (i=First; i<Last; i++){
		fprintf(fpName,"%s ",&ChrmNames[8*i]);
	}

	fclose(fpName);
	
	
	//output file SmryAll
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"SmryAll");
	if (nproc>0){
		tmpint = sprintf(str,"%02d",iproc);
		strcat(tmpFileName,str);
	}
	FILE *fpSmryAll = fopen(tmpFileName,"wb");
	if (fpSmryAll == NULL) {
		fprintf(stderr,"Can't open output file %s!\n",tmpFileName);
		exit(1);
	}
	fwrite(&smry[First*6],sizeof(double), (Last-First)*6,fpSmryAll);
	fclose(fpSmryAll);	
	

	
	

	
	clock_t end = clock();
	SPAMR((fpResult,"Time for processing file in seconds\n%2.3f\nProcess index iproc=%d, nproc=%d\n", (double)(end-start)/CLOCKS_PER_SEC,iproc,nproc));
	//printf("Time for processing file %2.3f sec.\n", (double)(end-start)/CLOCKS_PER_SEC);
	//fprintf(fpResult,"NoExceedValue=%d, NoExceededLength=%d\n",ExceededValue,ExceededLength);
	//fprintf(fpResult,"BlockSize=%d, SearchSize=%d\n",BlockSize, SearchSize);
	
	
	
	
	fclose(fp);
	fclose(fpDiffSeq);
	fclose(fpLenDiffSeq);
	fclose(fpDiff);
	fclose(fpLenDiff);
	fclose(fpResult);
	fclose(fpBlockAux);
	fclose(fpStart);
	fclose(fpSmry);
	fclose(fpPre);
	
	return 0;


}
