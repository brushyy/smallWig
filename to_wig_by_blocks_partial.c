/*
 * Change a matlab seq_in_blocks back to wig format
 * Bstart is the code block index for Qstart (counting all blocks before Qchrm also), Bend is the last block index
*/


int to_wig_by_blocks_partial(char * argv[],char *Qchrm,int Qstart,int Qend, int *Bstart, int *Bend){
	
	clock_t start = clock();
	SPAM(("================\nto_wig_by_blocks_partial.c start!\n"));
	
	char FileName[400];//name of processed Wig file
	strcpy(FileName, argv[1]);
	
	
	int Nd = 18;
	int SetMaxDiff = exp2(Nd); //set max alphabet size
	//unsigned int DiffCount; //frequncy of each difference value of contigs
	//unsigned int LenDiffCount;
	double Map[SetMaxDiff]; //maps each consecutive integer to the actual alphabet
	memset(Map,0,SetMaxDiff*sizeof(double)); 
	double LenMap[SetMaxDiff]; 
	memset(LenMap,0,SetMaxDiff*sizeof(double));  
	//char ChrmNames[8*32]; //every chrm has a name less than 8 characters, assume no more than 32 chrm 
	//memset(ChrmNames,0,8*32*sizeof(char));
	//int ChrmCount=0;
	
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
	//int SetBlockSize = exp2(BlockSize); //number of diff that are encoded together
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
		
	
	//fpstart                 
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"Start");
	FILE *fpStart = fopen(tmpFileName, "r");
	if (fpStart == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}
	//fpPre                
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"Pre");
	FILE *fpPre = fopen(tmpFileName, "r");
	if (fpPre == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}
	//fpSmry               
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"Smry");
	FILE *fpSmry = fopen(tmpFileName, "r");
	if (fpSmry == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}
	
	SPAM(("start and end block index %d %d\n",*Bstart, *Bend));
	//get the starting info
	fseek(fpStart, sizeof(double)*(*Bstart), SEEK_SET);
	fread(&PreLocation,sizeof(double),1,fpStart);
	double Location1=0, Location2=0; //first location (PreLocation) of *Bstart+1 and *Bend
	fread(&Location1,sizeof(double),1,fpStart);
	fseek(fpStart, sizeof(double)*(*Bend), SEEK_SET);
	fread(&Location2,sizeof(double),1,fpStart);
	
	//Pre info
	fseek(fpPre, (sizeof(double)*2)*(*Bstart), SEEK_SET);
	fread(&PreValue,sizeof(double),1,fpPre);
	fread(&PreContig,sizeof(double),1,fpPre);
	fclose(fpStart);
	fclose(fpPre);
	SPAM(("PreValue=%lg, PreLocation=%ld, PreContig=%ld, Location1=%ld, Location2=%ld\n",
		PreValue, (long int)PreLocation, (long int)PreContig, (long int)Location1, (long int)Location2));
	
	//Smry info
	double min1, max1, mean1, coverage1, totalloc1; //result from each block
	double std1;
	double min=1<<30, max=-(1<<30),  coverage=0, totalloc=0; //overall result
	double mean=0, std=0;
	
	if (*Bend-*Bstart>=2){ //at least 3 blocks
		fseek(fpSmry, (6*sizeof(double))*(*Bstart + 1), SEEK_SET);
		for (i=*Bstart+1; i<*Bend; i++){
			fread(&min1, sizeof(double), 1, fpSmry);
			fread(&max1, sizeof(double), 1, fpSmry);
			fread(&mean1, sizeof(double), 1, fpSmry);
			fread(&std1, sizeof(double), 1, fpSmry);
			fread(&coverage1, sizeof(double), 1, fpSmry);
			fread(&totalloc1, sizeof(double), 1, fpSmry);
		
			if (min1 < min)
				min = min1;
			if (max1 > max)
				max = max1;
			mean = (mean*coverage + (double)mean1)/(coverage + coverage1);
			std = (std*coverage + (double)std1)/(coverage + coverage1);
			coverage += coverage1;
			totalloc += totalloc1; 
			//printf("middle blocks min1 %d, max1 %d, mean1 %d, std1 %ld, coverage1 %d, totoalloc1 %d\n",min1,max1,mean1,std1,coverage1,totalloc1);
		}
	}
	fclose(fpSmry);
	
	
	j=0;
	i = 0; //Chrm count
	pcount = 1024; //count to print annotation
	bcount = 0;
	Diff=1;
	LenDiff=1;
	int tmpcount=0; //tmp debug output lines count
	
	while(!feof(fpDiffSeqMatlab) && !feof(fpLenDiffSeqMatlab)){
		bcount++;
		/*if (j<100){ 
			printf("Location=%d Contig=%ld LenDiff=%ld Value=%d Diff=%d\n", Location, Contig, LenDiff,Value,Diff);
		    j++;
		}*/
		
		if (pcount==1024){
			fprintf(fp,"variableStep chrom=%s span=1\n", Qchrm);
			pcount = 0;
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
			continue;
		}
				
		//write to wig	
		Value = PreValue + Diff;
		Contig = PreContig + LenDiff;
		
		
		if (Contig < 0){
			printf("Contig=%ld out of range! Check input file please!\n", (long int)Contig);
			printf("Location=%ld PreLocation=%ld Value=%lg PreValue=%lg ChrmName=%s\n",(long int)Location,(long int)PreLocation,Value,PreValue,Qchrm);
			printf("Line # in wig is %lld\n",lcount);
			exit(1);
		}
		
		//first or last block overlapping with query
		if ((Location1 <= Location2 && PreLocation < Location1 && PreLocation+Contig-1 >= (double)Qstart) 
			|| (PreLocation+Contig-1 >= Location2 && PreLocation <= (double)Qend && PreLocation+Contig-1 >= (double)Qstart)){
			
			/*if (tmpcount<10){
				printf("PreLocation=%d, LastLocation=%d, Value=%d, Contig=%d\n",PreLocation, PreLocation+Contig-1,Value,Contig);
				printf("min %d max %d mean %.2lf std %.2lf coverage %d totalloc %d\n",min,max,mean,std,coverage,totalloc);
				tmpcount++;
			}*/
			
			double tmpdouble=Contig; //temp contig
			if (PreLocation+Contig-1 > (double)Qend)
				tmpdouble = tmpdouble - (PreLocation+Contig-1 - Qend);
			if (PreLocation < (double)Qstart)
				tmpdouble = tmpdouble - (Qstart-PreLocation);
			
			if ((Value>1e-8 || Value<-1e-8) && PreLocation!=-1){
				if (Value < min)
					min = Value;
				if (Value > max)
					max = Value;
				mean = (mean*coverage + (double)Value*tmpdouble)/(coverage + tmpdouble);
				std = (std*coverage + (double)Value*Value*tmpdouble)/(coverage + tmpdouble);
				coverage += tmpdouble;
				//printf("Value %d, coverage %d, contig %d, std %lf\n",Value,coverage,tmp,std);
			}
			totalloc += tmpdouble;  
			
		}
		
			
		if ((Value>1e-8 || Value<-1e-8) && PreLocation!=-1){
			for (Location=PreLocation; Location < PreLocation + Contig; Location++ ){
				if (Location >= (double)Qstart && Location <= (double)Qend){
					//fprintf(fp,"%ld\t%lg\n", Location,(double)(Value)/SetFloatSize);
					fprintf(fp,"%ld\t%lg\n", (long int)Location,Value);
					lcount ++;
					pcount ++;
				}
			}
		}
		else{
			Location = PreLocation + Contig;
		}
		//reset
		PreValue = Value;
		PreContig = Contig;
		PreLocation = Location;
		//Diff = NextDiff;
		
	}
	clock_t end = clock();
	
	/*printf("minimum\t\t%g\n"
	       "maximum\t\t%g\n"
	       "average\t\t%lf\n"
	       "std dev\t\t%lf\n"
	       "coverage\t%lf\n",
	       (double)(min)/SetFloatSize,
	       (double)(max)/SetFloatSize,
	       (double)(mean)/SetFloatSize,
	       sqrt((std-mean*mean)*coverage/(coverage-1))/SetFloatSize,
	       double(coverage)/totalloc);*/
	printf("minimum\t\t%g\n"
	       "maximum\t\t%g\n"
	       "average\t\t%lf\n"
	       "std dev\t\t%lf\n"
	       "coverage\t%lf\n",
	       min, max, mean, sqrt((std-mean*mean)*coverage/(coverage-1)),
	       double(coverage)/totalloc);
	
	
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
	fclose(fpResult);

	return 1;
}
