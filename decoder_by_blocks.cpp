//decode data by blocks of size 2^BlockSize using arithmetic code
//do not need BlockAdd table to determine beginning of a new block  
//also supports parallel processing


int decoder_by_blocks(char * argv[],int nproc, int iproc) {

	char FileName[400];//name of processed Wig file
	char Type[40];
	char T[40];

	strcpy(FileName, argv[1]);

	//for parallel
	char str[2];
	int First; //first chrm to read
	int Last;
	//int pflag; //indicate whether should start read
	if (nproc<=0){
		First = 0;
		Last = 32;//assume no more than 32 chrm
		//pflag=1;
	}
	else{
		First = (24/nproc)*iproc; 
		if (iproc < nproc-1)
			Last = (24/nproc)*(iproc+1); 
		else
			Last = 32;
		//pflag=0; 
	}
	long int StartAdd=0, EndAdd=0;


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
	//int *data; //data read from Code file for comparison
	int index; //indexd data 
	//int cm_size;
	int *cm;
	int buffer_size;
	int SetBlockSize = (int)(1)<<BlockSize;
	//int BlockNum=400000;
	//int BlockAdd[BlockNum];
	//int BMASK = ((int)(1)<<BlockSize) - 1; //mask used to do modulo 2^BlockSize
	bool isOK;
	int i=0;
	int j=0;
	int BlockCount; 
	int pcount=0; //tmp print count
	//FILE *fpDiffSeqMatlab;
	int tmpint, tmpint2;
	//long int underflow; //count number of times low and high are too close
	
	SPAM(("=============\ndecoder_by_blocks.cpp start!\n"));

	//for (int type=0; type<1; type++){
	for (int type=0; type<2; type++){
		if (type==1){
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
		
		//read chrmadd from file
		strcpy(tmpFileName,FileName);
		strcat(tmpFileName,Type);
		strcat(tmpFileName,"ChrmAdd");
		FILE *fpChrmAdd = fopen(tmpFileName, "r");
		if (fpChrmAdd == NULL) {
			fprintf(stderr, "Can't open file %s!\n",tmpFileName);
			exit(1);
		}
		for (i=0; i<First+1 && !feof(fpChrmAdd); i++){
			fscanf(fpChrmAdd,"%ld ", &StartAdd);				
		}
		for (; i<Last+1 && !feof(fpChrmAdd); i++){
			fscanf(fpChrmAdd,"%ld ", &EndAdd);
		}
		fclose(fpChrmAdd);
		
		//read g_buffer from file
		strcpy(tmpFileName,FileName);
		strcat(tmpFileName,Type);
		strcat(tmpFileName,"Code");
		FILE *fpCode = fopen(tmpFileName, "rb");
		if (fpCode  == NULL) {
			fprintf(stderr, "Can't open file %s!\n",tmpFileName);
			exit(1);
		}		

		fseek(fpCode, StartAdd, SEEK_SET);
		g_buffer = (unsigned char*)malloc(EndAdd-StartAdd);
		fread(g_buffer,1,EndAdd-StartAdd,fpCode);
		buffer_size = EndAdd-StartAdd;
			
		
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
		
		//write index to a file
		strcpy(tmpFileName,FileName);
		strcat(tmpFileName,Type);
		strcat(tmpFileName,"Decode");
		if (nproc>0){
			tmpint = sprintf(str,"%02d",iproc);
			strcat(tmpFileName,str);
		}
		FILE *fpDecode = fopen(tmpFileName, "w+");
		if (fpDecode  == NULL) {
			fprintf(stderr, "Can't open file %s!\n",tmpFileName);
			exit(1);
		}

		
		pcount=0;
		isOK = true;	
		current_byte = 0;
		BlockCount = 0;
		i=0;
		
		//printf("==========\nBegin to decode. BlockNum=%d\n",BlockNum);
		//printf("data[4096]=%d\n",data[4097]);
		//while (BlockCount<BlockNum && isOK){
		while (isOK && (current_byte<buffer_size)){
			/*if (current_byte==buffer_size){
				printf("current_byte exceeds buffer_size %d, BlockCount=%d\n",buffer_size, BlockCount);
				exit(1);
			}*/
			//current_byte is automatically the next block address
			//current_byte = BlockAdd[BlockCount];
			/*if (BlockCount>=BlockNum-2 || BlockCount<=4)
				printf("BlockCount %d, BlockAdd %d\n",BlockCount,BlockAdd[BlockCount]);*/
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
				if (i==4096 || i==4095){
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
				printf("Decoded block not a multiple of SetBlockSize %d! i=%d, BlockCount=%d\n",SetBlockSize,i,BlockCount);
				isOK = false;
			}		
			delete rm_decode;
		}
		
//		if (i != data_size) isOK = false;
	
	
		clock_t end_decoding = clock();
		//end decoding
		
		
		/*///////////////////////////compare//////////////////////////////
		//load data from file for correctness check
		strcpy(tmpFileName,FileName);
		strcat(tmpFileName,Type);
		fpDiffSeqMatlab = fopen(tmpFileName, "r");
		if (fpDiffSeqMatlab  == NULL) {
			fprintf(stderr, "Can't open file %s!\n",tmpFileName);
			exit(1);
		}
		i=0;//counter inside a block
		rewind(fpDecode);
		while (!feof(fpDiffSeqMatlab) && !feof(fpDecode)){
			fscanf(fpDiffSeqMatlab,"%d ", &tmpint2);
			fscanf(fpDecode,"%d ",&tmpint);
			if (tmpint != tmpint2){
				printf("Decode not matched at symbol %d!",i);
				isOK = false;
			}
			i++;
		}
		if (!feof(fpDiffSeqMatlab) || !feof(fpDecode)){//different lengths of files
			printf("Wrond seq sizes %d!",i);
			isOK = false;
		}
		data_size = i;
		fclose(fpDiffSeqMatlab);
		if (isOK)
			printf("Round trip is OK\n");
		*/
		
		
		//if (isOK){
			//printf("Round trip is OK\n");
			//write to summary result file
			FILE *fpResult=fopen("result","a");
			if (fpResult == NULL) {
				fprintf(stderr,"Can't open output result file!\n");
				exit(1);
			}
			SPAMR((fpResult,"\n\n//random_access/decoder_by_blocks\n//File is %s%s\n", FileName, Type));
			SPAMR((fpResult,"Time for decoding %2.3f sec.\n", (double)(end_decoding - start_decoding)/CLOCKS_PER_SEC));
			//fprintf(fpResult,"Diff Seq length is %d, cm_size is %d.\n",data_size, cm_size);
			fclose(fpResult);
			
			/*printf("\n\n//random_access/decoder_by_blocks\n//File is %s%s\n", FileName, Type);
			printf("Time for decoding %2.3f sec.\n", (double)(end_decoding - start_decoding)/CLOCKS_PER_SEC);
			printf("Diff Seq length is %d, cm_size is %d.\n",data_size, cm_size);*/
		//}
		/*else{
			printf("\n\nData mismatch, decoded size=%d, data size=%d\n",i,data_size);
			printf("Diff Seq length is %d, cm_size is %d.\n",data_size, cm_size);
		}*/
		free(g_buffer);
		free(cm);
		
		fclose(fpCode);
		fclose(fpDecode);
	}
	return 0;
	
}
