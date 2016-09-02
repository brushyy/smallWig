//
//  getseq.c regular
//
//  Created by Zhiying on 12/5/13.
//  Copyright (c) 2013 Zhi Ying. All rights reserved.
//
//	Get a sequence of difference values, and difference run lengths, together with their counts
//  To use them in matlab, the diff sequences are converted to values within 1~size(DiffCount)
//  Each pair (Diff,LenDiff) corresponds to the difference of a run
//	First pair is assumed to be (Value=-1,Contig=1), and first location is assumed to be 1.
//





//record an item
void record0(
int Value, int PreValue, 
long unsigned int PreLocation,
int SetMaxDiff, 
int *Diff, int *NextDiff, 
int *ExceededValue, int *ExceededLength, 
long unsigned int *Contig, long unsigned int *PreContig, 
long int *LenDiff, 
unsigned int *DiffCount, unsigned int *LenDiffCount,
FILE *fpLenDiffSeq, FILE *fpDiffSeq, FILE *fpExceed, FILE *fpLenExceed){
	

	*Diff = *NextDiff; 
	*NextDiff = Value-PreValue;	
	*LenDiff = *Contig-*PreContig;
	
	//if exceed Max Value or Length
	if (mydiff(*Diff) >= SetMaxDiff-2){	
		ExceededValue ++;
		fprintf(fpExceed, "%d ", *Diff);
		*Diff = getdiff(SetMaxDiff-2);
		//printf("Diff is %d\n", *Diff);
	}
	if  (mydiff(*LenDiff) >= SetMaxDiff-2){
		//printf("LenDiff=%ld\n",*LenDiff);
		ExceededLength ++;
		fprintf(fpLenExceed, "%ld ", *LenDiff);
		*LenDiff = getdiff(SetMaxDiff-2);
		//printf("Max LenDiff is %ld\n", *LenDiff);
	}
	
	//record
	fprintf(fpLenDiffSeq,"%ld ", *LenDiff);
	fprintf(fpDiffSeq,"%d ",*Diff);
	
	DiffCount[mydiff(*Diff)]++;
	LenDiffCount[mydiff(*LenDiff)]++;
	*PreContig = *Contig; //reset
	*Contig = 1;
		
}

//////////////////// main////////////////////////////////////////////////////
int getseq(char * argv[])
{
	clock_t start = clock();
	SPAM(("=================\ngetseq.c start!\n"));
	//printf("sizeof(long unsinged int)=%ld,sizeof(int)=%ld\n",sizeof(long unsigned int),sizeof(int));

	char FileName[400];//name of processed Wig file
	

	strcpy(FileName, argv[1]);
		
	
	//int N = 16;
	int Nd = 18;
	//int SetMaxLength = exp2(Nd-1); //set max contig length
	int SetMaxDiff = exp2(Nd); //set max value difference 
	char ChrmNames[8*32]; //every chrm has a name less than 8 characters, assume no more than 32 chrm 
	memset(ChrmNames,0,8*32*sizeof(char));
	char *tmpChrmName;
	unsigned int DiffCount[SetMaxDiff]; //frequncy of each difference value of contigs
	memset(DiffCount,0,SetMaxDiff*sizeof(unsigned int));
	unsigned int LenDiffCount[SetMaxDiff]; //frequncy of each difference length of contigs
	memset(LenDiffCount,0,SetMaxDiff*sizeof(unsigned int));
	int DiffSub[SetMaxDiff]; //subtract this vector from the Diff Seq
	memset(DiffSub,0,SetMaxDiff*sizeof(int));
	int LenDiffSub[SetMaxDiff]; //subtract this vector from the LenDiff Seq
	memset(LenDiffSub,0,SetMaxDiff*sizeof(int));
  

	int PreValue=-1;
	int Value=0;
	int Diff=1;//difference of values of contigs
	//int PreDiff=1;
	int NextDiff=1;
	long unsigned int PreLocation=0;
	long unsigned int Location=0;
	long unsigned int MaxLocation =10000;//the largest location. Not necessarily the last one, because the chromozones are not in order 
	int ExceededValue = 0; //number of data points that exceed SetMaxDiff 
	int ExceededLength = 0; //number of contig length diff that exceed SetMaxDiff
	long unsigned int Contig = 1; //length of current contig
	long unsigned int PreContig = 1; //length of previous contig
	long int LenDiff=1; //length difference of two consecutive contigs
	//long int PreLenDiff=1; //length diff of previous two contigs
	int ChrmCount=-1; //temp chrm count
	char tmp;
	char str[1000];  //annotation line
	char Prestr[1000]; //previous annotation line
    char tmpFileName[400];	
    fpos_t position;
    int i=0;
	int j=0;
	//int k=0;
	//int tmpcount=0;
	//int tmpint=0;
	//int tmpflag=0;
	double tmpdouble=0;
	
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
	FILE *fpDiffSeq = fopen(tmpFileName, "w+");
	if (fpDiffSeq == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}
	//output lendiffseq
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"LenDiffSeq");
	FILE *fpLenDiffSeq = fopen(tmpFileName, "w+");
	if (fpLenDiffSeq == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}
	//output exceed
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"Exceed");
	FILE *fpExceed = fopen(tmpFileName, "w+");
	if (fpExceed == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}
	//output lenexceed
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"LenExceed");
	FILE *fpLenExceed = fopen(tmpFileName, "w+");
	if (fpLenExceed == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}
	
	//output result
	FILE *fpResult=fopen("result","a");
	if (fpResult == NULL) {
		fprintf(stderr,"Can't open output result file!\n");
		exit(1);
	}
	SPAMR((fpResult,"\n//getseq.c\n//%s\n",FileName));
	
	
	
	
	SPAM(("\n===============================\nstart processing...\n"));

	//start processing input file
	strcpy(Prestr, "empty");
	while (!feof(fp)){		
		tmp = fgetc(fp);
		if  (tmp >= '0' && tmp <= '9'){
			ungetc(tmp,fp);
			fgetpos (fp, &position);
			fscanf(fp,"%lu %lf\n", &Location,&tmpdouble);
			Value = (int)(tmpdouble); //in case Value is stored in a floating point manner
			//if (tmpflag==1)
			//	fprintf(fpResult,"%lu %d\n", Location,Value);
			if (Location > MaxLocation)
				MaxLocation = Location;
			if (Location != PreLocation+1){//for contig before 0s
				Value = 0;
			}
			if (Value != PreValue){
				record0(Value, PreValue, PreLocation, SetMaxDiff, &Diff, &NextDiff, &ExceededValue, &ExceededLength, 
				&Contig, &PreContig, &LenDiff, DiffCount, LenDiffCount, fpLenDiffSeq, fpDiffSeq, fpExceed, fpLenExceed);
			}
			else{//Diff == 0
					Contig++;
			}
			
			if (Location != PreLocation+1){
				//record skipped 0s, push back the line
				fsetpos (fp, &position);
			    PreValue = 0;
			    if (Location<PreLocation){
					printf("Location %lu < PreLocation %lu!",Location, PreLocation);
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
		else{ //get a line of annotations
			ungetc(tmp,fp);
			fgets(str,sizeof(str),fp); 
			if (strcmp(str, Prestr)){ //different chrm
				//old chrm
				if (ChrmCount>=0){ 
					//insert a virtual data point and then record
					PreValue = Value;
					Value++;
					PreLocation = Location;
					Location++;
					record0(Value, PreValue, PreLocation, SetMaxDiff, &Diff, &NextDiff, &ExceededValue, &ExceededLength, 
					&Contig, &PreContig, &LenDiff, DiffCount, LenDiffCount, fpLenDiffSeq, fpDiffSeq, fpExceed, fpLenExceed);
					//end marker of a chrom
					fprintf(fpLenDiffSeq,"%ld ", getdiff(SetMaxDiff-1));
					fprintf(fpDiffSeq,"%ld ", getdiff(SetMaxDiff-1));
				}	
				//new chrm
				//printf("chrm %d, %s", ChrmCount, str);
				ChrmCount++;
				if (ChrmCount>=32){
					printf("Too many chromosomes!\n");
					exit(1);
				}	
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
				strcpy(Prestr, str);
				
				//get first location, value
				/*fgetpos(fp,&position);
				fscanf(fp,"%lu %lf\n", &Location,&tmpdouble);
				Value = (int)(tmpdouble); 
				fsetpos(fp,&position);*/
				//starting index
				//IndexAux[ChrmCount*2] = (((BlockStart)-1)>>BlockSize)+1;
				//reset parameters
				PreValue=-1;
				Value=0;
				Diff=1;
				NextDiff=1;
				PreLocation=0;
				Location=0;
				Contig = 1; 
				PreContig = 1; 
				LenDiff=1;
				
			}
		}
	}
	//ending part of last chrm
	if (Location>=0){ 
		//printf("ending of chrm %d\n",ChrmCount);
		//insert a virtual data point and then record
		PreValue = Value;
		Value++;
		PreLocation = Location;
		Location++;
		record0(Value, PreValue, PreLocation, SetMaxDiff, &Diff, &NextDiff, &ExceededValue, &ExceededLength, 
		&Contig, &PreContig, &LenDiff, DiffCount, LenDiffCount, fpLenDiffSeq, fpDiffSeq, fpExceed, fpLenExceed);											
		//end marker of a chrom
		fprintf(fpLenDiffSeq,"%ld ", getdiff(SetMaxDiff-1));
		fprintf(fpDiffSeq,"%ld ",getdiff(SetMaxDiff-1));
	}
	ChrmCount++;
	
	
	//open output files 1
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"Diff");
	FILE *fpDiff = fopen(tmpFileName,"w");
	if (fpDiff == NULL) {
		fprintf(stderr,"Can't open output Count file!\n");
		exit(1);
	}
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"LenDiff");
	FILE *fpLenDiff = fopen(tmpFileName,"w");
	if (fpLenDiff == NULL) {
		fprintf(stderr,"Can't open output Count file!\n");
		exit(1);
	}
	//write to output files 1
	DiffCount[SetMaxDiff-2] = 1;//exceed marker
	LenDiffCount[SetMaxDiff-2] = 1; 
	DiffCount[SetMaxDiff-1] = 1;//end marker of a chrm
	LenDiffCount[SetMaxDiff-1] = 1; 
	fprintf(fpDiff,"#Difference\tCount\n");
	fprintf(fpLenDiff,"#Difference of Length\tCount\n");
	for (i=0; i<SetMaxDiff; i++){
		if (DiffCount[i])
			fprintf(fpDiff,"%ld\t%u\n", getdiff(i),DiffCount[i]);
		if (LenDiffCount[i])
			fprintf(fpLenDiff,"%ld\t%u\n", getdiff(i),LenDiffCount[i]);
	}
	
	
	// output file ChrmName
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"ChrmName");
	FILE *fpName = fopen(tmpFileName,"w");
	if (fpName == NULL) {
		fprintf(stderr,"Can't open output file!\n");
		exit(1);
	}
    for (i=0; i<ChrmCount; i++){
		fprintf(fpName,"%s ",&ChrmNames[8*i]);
	}
//	fprintf(fpName,"\n");
	
	

	//get sequences for matlab arithmetic coding use
	//DiffMatlab is the count vector w. all counts >=0
	//DiffSeqMatlab is consecutive values from 1 to size DiffMatlab
	int tmpSub = 0;
	int tmpLenSub = 0;
	for (i=0; i<SetMaxDiff; i++){
	if (DiffCount[i]==0)
		tmpSub++;
	else
		DiffSub[i]=tmpSub;
	if (LenDiffCount[i]==0)
		tmpLenSub++;
	else
		LenDiffSub[i]=tmpLenSub;
	}

	//matlab seq files
	rewind(fpDiffSeq);
	rewind(fpLenDiffSeq);	
	
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"DiffSeqMatlab");
	FILE *fpDiffSeqMatlab = fopen(tmpFileName,"w");
	if (fpDiffSeqMatlab == NULL) {
		fprintf(stderr,"Can't open output Count file!\n");
		exit(1);
	}
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"LenDiffSeqMatlab");
	FILE *fpLenDiffSeqMatlab = fopen(tmpFileName,"w");
	if (fpLenDiffSeqMatlab == NULL) {
		fprintf(stderr,"Can't open output Count file!\n");
		exit(1);
	}
	//write to matlab seq files
	while (!feof(fpDiffSeq)){
		fscanf(fpDiffSeq,"%d ",&Diff);
		fscanf(fpLenDiffSeq,"%ld ", &LenDiff);
		fprintf(fpDiffSeqMatlab,"%d ",mydiff(Diff)-DiffSub[mydiff(Diff)]);
		fprintf(fpLenDiffSeqMatlab,"%d ",mydiff(LenDiff)-LenDiffSub[mydiff(LenDiff)]);
	}
    
	
	
	clock_t end = clock();
	SPAMR((fpResult,"Time for processing file in seconds\n%2.3f\n", (double)(end-start)/CLOCKS_PER_SEC));
	SPAM(("Time for processing file %2.3f sec.\n", (double)(end-start)/CLOCKS_PER_SEC));
	SPAMR((fpResult,"NoExceedValue=%d, NoExceededLength=%d\n",ExceededValue,ExceededLength));

	
	fclose(fp);
	fclose(fpDiffSeq);
	fclose(fpLenDiffSeq);
	fclose(fpDiffSeqMatlab);
	fclose(fpLenDiffSeqMatlab);
	fclose(fpDiff);
	fclose(fpLenDiff);
    fclose(fpName);
	fclose(fpResult);
	fclose(fpExceed);
	fclose(fpLenExceed);
	
	return 1;


}
