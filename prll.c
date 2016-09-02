/*  run some function f using parallel processes,
*   each process run on a sperate processor
*/


int prll( int *nproc, char **argv, int (*f)(char**,int,int) )
{
	clock_t start=clock();
	pid_t mypid;
	
	/* Determine the actual number of processors */
	int NUM_PROCS = sysconf(_SC_NPROCESSORS_CONF);
  	//printf("System has %i processor(s).\n", NUM_PROCS);

    /* Check for user specified parameters */
	if (*nproc < 1){
		printf("WARNING: Must utilize at least 1 cpu. Spinning "
		" all %i cpu(s) instead...\n", NUM_PROCS);
		*nproc = NUM_PROCS;
	}
	else if (*nproc > NUM_PROCS){
		printf("WARNING: %i cpu(s), are not "
		"available on this system, spinning all %i cpu(s) "
		"instead...\n", *nproc, NUM_PROCS);
		*nproc = NUM_PROCS;
	}
	else{	
		SPAM(("Maxing computation on %i cpu(s)...\n", *nproc));
	}

	/* Kick off the actual work of spawning threads and computing */
	int created_thread = 0;

	/* We need a thread for each cpu we have... */
	while ( created_thread < *nproc - 1 )
	{
		mypid = fork();
		
		if (mypid < 0){
			printf("Failed to fork process 2\n");
			exit(1);
		}
		
		else if (mypid == 0) /* Child process */
 		{
			SPAM(("\tCreating Child Thread: #%i\n", created_thread));
			break;
		}

		else /* Only parent executes this */
		{ 
			/* Continue looping until we spawned enough threads! */
			created_thread++;
		} 
	}

	/* NOTE: All threads execute code from here down! */

	cpu_set_t mask;

	/* CPU_ZERO initializes all the bits in the mask to zero. */ 
	CPU_ZERO( &mask ); 	

	/* CPU_SET sets only the bit corresponding to cpu. */
	CPU_SET( created_thread, &mask );  
        
	/* sched_setaffinity returns 0 in success */
    if( sched_setaffinity( 0, sizeof(mask), &mask ) == -1 ){
		printf("WARNING: Could not set CPU Affinity, continuing...\n");
	}
	/* sched_setaffinity sets the CPU affinity mask of the process denoted by pid. 
	If pid is zero, then the current process is used.*/


	/* Now we have a single thread bound to each cpu on the system */
	int NPROC = *nproc; //actual number of processors used
	(*f)(argv,NPROC,created_thread);
	
	/*cpu_set_t mycpuid;	
	sched_getaffinity(0, sizeof(mycpuid), &mycpuid);*/
	
	if (mypid==0)//child process exits
		exit(0);
	else{//parent wait till all are finished
		int status;
		pid_t pid;
		int i;
		for (i=0;i<*nproc-1;i++){
			pid = wait(&status);
			//printf("\tParent detects process %d was done\n", pid);
		}
		clock_t end=clock();
		SPAM(("parallel time is  %2.3f sec.\n", (double)(end-start)/CLOCKS_PER_SEC));
		return 0;
	}
	
}





