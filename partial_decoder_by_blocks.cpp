//partially decode data by blocks of size 2^BlockSize using arithmetic code
//returns BStart, BEnd, which are starting and ending block indices in BlockAdd[]

void partial_decoder_by_blocks(char * argv[], char *Qchrm, int QStart, int QEnd, int *BStart, int *BEnd) {

	char FileName[400];//name of processed Wig file
	char Type[40]; // type of Diff or LenDiff
	char T[40];
	//int QStart, QEnd; //start and end locations of query
	
	strcpy(FileName, argv[1]);
	
	SPAM(("=============\npartial_decoder_by_blocks.cpp start!\n"));
	
	
	
	char tmpFileName[400];
	//parameters that can be changed for testing
	//the RANGE_SIZE_IN_BITS depends on your alphabet size
	int alphabet_size;//256; //Please notice the size of the alphabet. 
	//int data_size;// = 39432012;
	int RANGE_SIZE_IN_BITS=30;
	//ranges_in_bits<32, it cannot be more than sizeof(int), which is 4 bytes=32 bits, otherwise makeRanges does not work
	//range_in_bits is the number of bits used to represent the probabilities
	// L=64-range_in_bits is the number of bits used to store the current range of low and high
	//Need for any probability p of a symbol, p>= 1/2^(L-3). So L-3>=range_in_bits, or range_in_bits<=30
	int index; //indexd data 
	//int cm_size;
	int *cm;
	//int FileSize;
	int SetBlockSize = (int)(1)<<BlockSize;
	int BlockNum=0;
	//int BlockAdd[BlockNum];
	//int BMASK = ((int)(1)<<BlockSize) - 1; //mask used to do modulo 2^BlockSize
	char ChrmNames[32*8];
	//char *tmpName;
	int ChrmChosen; //chrm index that is chosen
	int ChrmEnd; //ending search block index of the chrm (start=0)
	//int StartLocation; //start chrm locatioin of first block overlapping the range
	int StartAdd[2]; //start address in bytes of first encoded block overlapping the range
	int EndAdd[2]; //end address in bytes of first encoded block overlapping the range
	
	bool isOK;
	long int i=0;
	int j=0;
	int tmpint,tmp;
	int tmpint2;
	int BlockCount; 
	int pcount=0; //tmp print count
	//FILE *fpDiffSeqMatlab;
	//long int underflow; //count number of times low and high are too close
	


	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"ChrmName");
	FILE *fpName = fopen(tmpFileName, "r");
	if (fpName  == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}

	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"BlockAux");
	FILE *fpBlockAux = fopen(tmpFileName, "r");
	if (fpBlockAux  == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}

	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"Start");
	FILE *fpStart = fopen(tmpFileName, "rb");
	if (fpStart  == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}

	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"DiffSeqMatlab");
	strcat(tmpFileName,"BlockAdd");
	FILE *fpBlockAdd0 = fopen(tmpFileName, "rb");
	if (fpBlockAdd0  == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}
	strcpy(tmpFileName,FileName);
	strcat(tmpFileName,"LenDiffSeqMatlab");
	strcat(tmpFileName,"BlockAdd");
	FILE *fpBlockAdd1 = fopen(tmpFileName, "rb");
	if (fpBlockAdd1  == NULL) {
		fprintf(stderr, "Can't open file %s!\n",tmpFileName);
		exit(1);
	}

	//find chrm
	int ChrmCount = 0;
	ChrmChosen = -1;
	while(!feof(fpName)){
        fscanf(fpName,"%s",ChrmNames);
        if (strcmp(ChrmNames,Qchrm)==0){
			ChrmChosen = ChrmCount;		
			break;
		}	
        ChrmCount++;
        //SPAM(("ChrmNames %s Qchrm %s\n",ChrmNames,Qchrm));
    }
    if (ChrmChosen==-1){
		printf("No query chrom with name %s found!\n",Qchrm);
		exit(1);
	}
	fclose(fpName);

	//get search blocks indices of the chrm
	//find block indices for chrmchosen
	int BlockAux[2]; 
	for (i=0; i<ChrmChosen+1; i++)
		fscanf(fpBlockAux,"%d ", &BlockAux[0]);
	fscanf(fpBlockAux,"%d ", &BlockAux[1]);
	fclose(fpBlockAux);
	//get starting locations for chrmchosen
	double *Start;
	Start = (double*)malloc((BlockAux[1]-BlockAux[0])*sizeof(double));
	fseek(fpStart, sizeof(double)*BlockAux[0], SEEK_SET);
	fread(Start,sizeof(double),BlockAux[1]-BlockAux[0],fpStart);
	fclose(fpStart);
	//printf("BlockAux[0,1]=%d %d, start[0,1]=%d %d\n",BlockAux[0],BlockAux[1],Start[0],Start[1]);
	//check range
	if (QStart<0){
		QStart = 0;
		printf("Query Start too small, changed  to 0!\n");
	}
	if (QEnd<QStart){
		QEnd = QStart+1;
		printf("must have Qend>Qstart, changed Qend to %d\n",QEnd);
	}
	//binary search	
	*BStart = binsearch(BlockAux[1]-BlockAux[0],Start,(double)QStart);
	*BEnd = binsearch(BlockAux[1]-BlockAux[0],Start,(double)QEnd);

	if (*BStart < 0){
		QStart = Start[0];
		printf("Query start too large, changed  to %d!\n",QStart);
	}
	if (*BEnd < 0){
		QEnd = Start[BlockAux[1]-BlockAux[0]-1];
		printf("Query end too large, changed  to %d!\n",QEnd);
	}
	free(Start);
	*BStart +=  BlockAux[0];
	*BEnd +=  BlockAux[0];
	BlockNum = *BEnd -  *BStart + 1;
	//find Start(End)Add 
	fseek(fpBlockAdd0, sizeof(int)*(*BStart), SEEK_SET);
	fread(&StartAdd[0],sizeof(int), 1, fpBlockAdd0);
	fseek(fpBlockAdd0, sizeof(int)*(*BEnd + 1), SEEK_SET);
	fread(&EndAdd[0],sizeof(int), 1, fpBlockAdd0);
	fseek(fpBlockAdd1, sizeof(int)*(*BStart), SEEK_SET);
	fread(&StartAdd[1],sizeof(int), 1, fpBlockAdd1);
	fseek(fpBlockAdd1, sizeof(int)*(*BEnd + 1), SEEK_SET);
	fread(&EndAdd[1],sizeof(int), 1, fpBlockAdd1);
	fclose(fpBlockAdd0);
	fclose(fpBlockAdd1);
	
	
	SPAM(("BStart %d, *BEnd %d, blocknum %d\n", *BStart,*BEnd,BlockNum));
	SPAM(("type 0 startadd %d, endadd %d\n", StartAdd[0], EndAdd[0]));
	SPAM(("type 1 startadd %d, endadd %d\n", StartAdd[1], EndAdd[1]));




	//////////////////////////////start decoding////////////////////////


	//for (int type=1; type<2; type++){
	for (int type=1; type>=0; type--){
		if (type==0){
			strcpy(Type,"DiffSeqMatlab");
			strcpy(T,"Diff");
			//alphabet_size = 6134;
		}
		else{
			strcpy(Type,"LenDiffSeqMatlab");
			strcpy(T,"LenDiff");
			//alphabet_size = 67955;
		}	
		
		
		
		clock_t start_decoding = clock();
		//read g_buffer from file
		strcpy(tmpFileName,FileName);
		strcat(tmpFileName,Type);
		strcat(tmpFileName,"Code");
		FILE *fpCode = fopen(tmpFileName, "rb");
		if (fpCode  == NULL) {
			fprintf(stderr, "Can't open file %s!\n",tmpFileName);
			exit(1);
		}	
				
		/////////load g_buffer from Code file
		
		fseek(fpCode, StartAdd[type], SEEK_SET);
		g_buffer = (unsigned char*)malloc(EndAdd[type]-StartAdd[type]);
		fread(g_buffer,1,EndAdd[type]-StartAdd[type],fpCode);
/*		for (j=0; j<20; j++)
			printf("%hhX ",g_buffer[j]);
		printf("\nabove is g_buffer from file\n");
		printf("g_buffer size is %d\n", EndAdd[type]-StartAdd[type]);
*/		
		
		
		//load count table file
		strcpy(tmpFileName,FileName);
		strcat(tmpFileName,T);
		FILE *fpDiff = fopen(tmpFileName, "r");
		if (fpDiff  == NULL) {
			fprintf(stderr, "Can't open file %s!\n",tmpFileName);
			exit(1);
		}

		absize(&alphabet_size, fpDiff);
		cm = (int*)malloc((alphabet_size+2)*sizeof(int));
		makeRanges(cm, alphabet_size, RANGE_SIZE_IN_BITS,fpDiff);
		fclose(fpDiff);

		
		//write decode to a file
		strcpy(tmpFileName,FileName);
		strcat(tmpFileName,Type);
		strcat(tmpFileName,"Decode");
		FILE *fpDecode = fopen(tmpFileName, "w+");
		if (fpDecode  == NULL) {
			fprintf(stderr, "Can't open file %s!\n",tmpFileName);
			exit(1);
		}

		
		
		isOK = true;	
		current_byte = 0;
		//BlockCount =  *BStart;
		BlockCount = 0;
		//i =  *BStart<<BlockSize;
		i = 0;
		
		//printf("==========\nBegin to decode.\n");
		
		//printf("data[4096]=%d\n",data[4097]);
		while (BlockCount<BlockNum && isOK){
			//if (BlockCount<2)
			//	printf("current_byte %d\n",current_byte);
			TMPREAD_COUNT = 0;
			TMPREAD = 0x00;
			RangeMapper* rm_decode = new RangeMapper(RANGE_SIZE_IN_BITS);
			rm_decode->init();
			while (true) {
				int midpoint = rm_decode->getMidPoint();
				//next is binary search algorithm that does not need having lookup array
				index = findInterval(cm, alphabet_size + 2, midpoint);
/*				if (pcount<20){
					printf("%d ",index);
					pcount++;
				}				
				if (i<=10){
					printf("i=%d, midpoint=%d,index=%d\n",i, midpoint,index);
					printf("i=%d, LOW is %llX, HIGH is %llX, midpoint is %d\n",i,rm_decode->getLOW(),rm_decode->getHIGH(),rm_decode->getMidPoint());
					printf("current_byte=%d, tmpread_count=%hhX, TMPWREAD=%hhX \n", current_byte,TMPREAD_COUNT,TMPREAD);					
				}
*/				if (index == alphabet_size){//end of data marker
					//finish reading leftover bits
					break; 
				}
/*				if (index != data[i]) {
					printf("Data mismatch at i=%d, decoded element=%d, data element=%d\n", i-1, index, data[i-1]);
					isOK = false;
					break;
				}
*/				i++;
				fprintf(fpDecode,"%d ",index);
				rm_decode->decodeRange(cm[index], cm[index+1]);		
			}
			BlockCount++;
			//printf("i %d\n",i);
			if (i != (BlockCount<<BlockSize)){
				printf("Decoded block not a multiple of SetBlockSize %d! i=%ld, BlockCount=%d\n",SetBlockSize,i,BlockCount);
				isOK = false;
			}		
			delete rm_decode;
		}
		
//		if (i != data_size) isOK = false;
	
	
		clock_t end_decoding = clock();
		//end decoding
		
		
		///////////////////////////compare//////////////////////////////
		//load data from file for correctness check
/*		strcpy(tmpFileName,FileName);
		strcat(tmpFileName,Type);
		fpDiffSeqMatlab = fopen(tmpFileName, "r");
		if (fpDiffSeqMatlab  == NULL) {
			fprintf(stderr, "Can't open file %s!\n",tmpFileName);
			exit(1);
		}
		for (i=0; i<( *BStart<<BlockSize) && !feof(fpDiffSeqMatlab); i++){//skip blocks before start block
			fscanf(fpDiffSeqMatlab,"%d ", &tmpint2);
		}
		rewind(fpDecode);
		while (!feof(fpDiffSeqMatlab) && !feof(fpDecode)){
			fscanf(fpDiffSeqMatlab,"%d ", &tmpint2);
			fscanf(fpDecode,"%d ",&tmpint);
			if (tmpint != tmpint2){
				printf("Decode not matched at symbol %ld! data is %d, decoded data is %d\n",i,tmpint2, tmpint);
				isOK = false;
				break;
			}
			i++;
		}
		if (i != (*BEnd+1)<<BlockSize){//different lengths of files
			printf("Wrong seq sizes %ld!",(*BEnd- *BStart+1)<<BlockSize);
			isOK = false;
		}
		data_size = (*BEnd- *BStart+1)<<BlockSize;
		fclose(fpDiffSeqMatlab);
*/		
		
		
		if (isOK){
			//printf("Round trip is OK\n");
			//write to summary result file
			FILE *fpResult=fopen("result","a");
			if (fpResult == NULL) {
				fprintf(stderr,"Can't open output result file!\n");
				exit(1);
			}
			SPAMR((fpResult,"\n\n//random_access/partial_decoder_by_blocks\n//File is %s%s\n", FileName, Type));
			SPAMR((fpResult,"Time for decoding %2.6f sec.\n", (double)(end_decoding - start_decoding)/CLOCKS_PER_SEC));
			//fprintf(fpResult,"Diff Seq length is %d, cm_size is %d.\n",data_size, cm_size);
			SPAMR((fpResult,"Query length is %d.\n",QEnd-QStart));
			fclose(fpResult);
			
			//printf("\n\n//random_access/partial_decoder_by_blocks\n//File is %s%s\n", FileName, Type);
			//printf("Time for decoding %2.6f sec.\n", (double)(end_decoding - start_decoding)/CLOCKS_PER_SEC);
			//printf("Diff Seq length is %d, cm_size is %d.\n",data_size, cm_size);
		}
		/*else{
			printf("\n\nData mismatch, decoded size=%ld, data size=%d\n",i-( *BStart<<BlockSize),data_size);
			printf("Diff Seq length is %d, cm_size is %d.\n",data_size, cm_size);
		}*/
		free(g_buffer);
		free(cm);
		fclose(fpCode);
		fclose(fpDecode);
	}
	
	//return( *BStart);
	SPAM(("BStart=%d BEnd=%d  \n",*BStart,*BEnd));
}
