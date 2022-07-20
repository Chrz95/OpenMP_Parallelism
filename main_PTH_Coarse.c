#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#define MAXCHAR 100
#include <pthread.h>
#include <assert.h>


// Global Variables

double time_E = 0 ;
double time_F = 0;
long long int NumOfCells = 0 ;
pthread_mutex_t ReadMutex;
int pair = 0 ;
long int Traceback_Steps = 0 ;
int Q_min ;
int Q_max;
int D_Sz_all;
char* output_filename ;
char* path ;
int match ;
int mismatch ;
int gap ;
FILE * fp ;
FILE * fp1 ;
int NumOfThreads ;
int Q_min ;
int Q_max  ;
int D_Sz_all  ;

double gettime(void)
{
	struct timeval ttime;

	gettimeofday(&ttime, NULL);
	return ttime.tv_sec+ttime.tv_usec * 0.000001;
}

int totalMax (int A,int B,int C,int D)
{
	int max = 0 ;

	if (A > B)
		max = A ;
	else
		max = B ;

	if (C > max)
		max = C ;

	if (D > max)
		max = D ;

	return max ;
}

char *strrev(char *str)
{
      char *p1, *p2;

      if (! str || ! *str)
            return str;
      for (p1 = str, p2 = str + strlen(str) - 1; p2 > p1; ++p1, --p2)
      {
            *p1 ^= *p2;
            *p2 ^= *p1;
            *p1 ^= *p2;
      }
      return str;
}

 long int Smith_Waterman (char* FirstString,char* SecondString,int* ScoringMatrix)
 {
	 	int length1 = Q_max; // Q
	 	int length2 = D_Sz_all; // D
	 	int StartPointI[length2 + 1];
		int StartPointJ[length1 + 1];
	 	int MaxValue = 0 , cnt = 0 ;
		int index = 0;
		int index_hor = 0 ;
		int index_vert = 0 ;
		int index_diag = 0 ;
		double timeE_0 = gettime();
		double timeE_1 = 0 ;
		double timeF_0 , timeF_1 ;

		// printf ("The Scoring Matrix has a size of %d bytes \n",(length2 + 1) * (length1 + 1) * sizeof(int)) ;

		// Scoring Matrix Calculation

		for (int i = 0 ; i <  length2 + 1 ; i++) // Initialize first column to 0
		{
			index = i * (length1 + 1) ;
			ScoringMatrix[index] = 0 ;
		}

		for (int j = 0 ; j <  length1 + 1 ; j++) // Initialize first row to 0
		{
			index = j ;
			ScoringMatrix[index] = 0 ;
		}

		index = 0 ;

		for (int i = 1 ; i < length2 + 1  ; i++) // Calculating the marix
		{
				for (int j = 1 ; j < length1 + 1; j++)
				{
							index = j + i * (length1 + 1) ;
							index_hor = (j-1) + i * (length1 + 1) ;
							index_vert =  j + (i-1) * (length1 + 1) ;
							index_diag = (j-1) + (i-1) * (length1 + 1) ;

							if (SecondString[i-1] == FirstString[j - 1]) ScoringMatrix[index] = totalMax(0,ScoringMatrix[index_vert] + gap,ScoringMatrix[index_hor] + gap,ScoringMatrix[index_diag] + match);
							else ScoringMatrix[index] = totalMax(0,ScoringMatrix[index_vert] + gap,ScoringMatrix[index_hor] + gap,ScoringMatrix[index_diag] + mismatch);

							if ((ScoringMatrix[index] > MaxValue) && (ScoringMatrix[index] !=0))  // Each time a new max value is found 
							{
								cnt = 0;

								for (int k = 0 ; k < length2 + 1 ; k++) // Re Initialize the Start Point matrix
								{
										StartPointI[k] = 0 ;
										StartPointJ[k] = 0 ;
								}

								MaxValue = ScoringMatrix[index] ;
								StartPointI[cnt] = i ;
								StartPointJ[cnt] = j ;
							}
							else if (ScoringMatrix[index] == MaxValue) // If the same max value is found, place its coordinates in the StartPoint matricrs
							{
									cnt++ ;
									if (cnt < length2 + 1)
									{
											StartPointI[cnt] = i  ;
											StartPointJ[cnt] = j  ;
									}
							}
					}
			}

	timeE_1 = gettime();
	time_E += timeE_1 - timeE_0 ;

 	// Traceback

	timeF_0 = gettime();

	int I , J ;
	int trace_cnt = 0 ;
	int Start = 0, Stop = 0 ;
	int traceback_length = 0 ;

	traceback_length = length2/100 + 2 ; // Εmperical equation 

	char* Q_traceback = (char *) malloc(traceback_length*length1) ;
	assert(Q_traceback!=NULL);
	char* D_traceback = (char *) malloc(traceback_length*length1) ;
	assert(D_traceback!=NULL);

	fprintf (fp1,"\nPair %d :\n",pair);
	//printf ("\nPair %d :\n",pair);

	for (int i=0 ; i <= cnt ; i++)
	{
		trace_cnt = 0 ;

		I = StartPointI[i]  ;
		J = StartPointJ[i] ;

		index = J + I * (length1 + 1) ;
		index_hor = (J-1) + I * (length1 + 1) ;
		index_vert =  J + (I-1) * (length1 + 1) ;
		index_diag = (J-1) + (I-1) * (length1 + 1) ;

		for (int i=0 ; i < traceback_length*length1 ; i++) // Initialize traceback matrices
		{
			Q_traceback[i] = '\0';
			D_traceback[i] = '\0';
		}

		while ((ScoringMatrix[index_vert] != 0) || (ScoringMatrix[index_diag] != 0) || (ScoringMatrix[index_hor] != 0)) // All surrounding elements must be zero to stop
		{
		//	printf ("(I,J) = (%d,%d) \n",I,J);
			index = J + I * (length1 + 1) ;
			index_hor = (J-1) + I * (length1 + 1) ;
			index_vert =  J + (I-1) * (length1 + 1) ;
			index_diag = (J-1) + (I-1) * (length1 + 1) ;
			//printf ("HELP \n ") ;

			if ((ScoringMatrix[index_diag] >= ScoringMatrix[index_hor]) && (ScoringMatrix[index_diag] >= ScoringMatrix[index_vert]))
			{
				Q_traceback[trace_cnt] = FirstString[J - 1] ;
				D_traceback[trace_cnt] = SecondString[I -  1] ;
				trace_cnt++;
				Traceback_Steps ++ ;
				I-- ;
				J-- ;
			}
			else if (ScoringMatrix[index_hor] >= ScoringMatrix[index_vert]) // Horizontal
			{
				Q_traceback[trace_cnt] = FirstString[J - 1] ;
				D_traceback[trace_cnt] = '-' ;
				trace_cnt++;
				Traceback_Steps++ ;
				J-- ;
			}
			else // if (ScoringMatrix[I -1][J] >= ScoringMatrix[I][J -1]) // Vertical
			{
				Q_traceback[trace_cnt] = '-' ;
				D_traceback[trace_cnt] =  SecondString[I - 1] ;
				trace_cnt++;
				Traceback_Steps++ ;
				I-- ;
			}
		}

			Stop = StartPointI[i] - 1 ;
			Start = I ;
			fprintf (fp1,"\nMatch %d [Score: %d,Start: %d,Stop: %d]\n",i+1,MaxValue,Start,Stop);
			fprintf (fp1,"D: %s \n",strrev(D_traceback));
		    fprintf (fp1,"Q: %s \n",strrev(Q_traceback));
	}

		timeF_1 = gettime();
		time_F += timeF_1 - timeF_0 ;

		//free (Q_traceback);
		//free (D_traceback);

	return Traceback_Steps;
 }

 void * Algorithm_Per_Pair (void * Args) // Thread function
 {
 	///////////////// Read a pair of strings /////////////////////////////////////
 	pthread_mutex_lock (&ReadMutex);

	int* ScoringMatrix = (int* ) Args ;

 	char * Q_string = (char *) malloc(Q_max + 2); // Initialize Q_string
	assert(Q_string!=NULL);
 	char * D_string = (char *) malloc(D_Sz_all + 2); // Initialize D_string
	assert(D_string!=NULL);
    char* ch_str = (char *) malloc (2) ;
	assert(ch_str!=NULL);
 	char ch = fgetc(fp) ;

 		pair ++ ;
 		while((ch != EOF) && (ch != 'D'))
 		{
 			ch = fgetc(fp) ;
 			if ((ch != ':') && (ch != ' ') && (ch != '\n') && (ch != '\t')  && (ch  != EOF) && (ch  != 'Q') && (ch  != 'D'))
 			{
 				ch_str[0] = ch ;
 				strcat(Q_string,ch_str);
 			}
 		}

 		fseek (fp,-1,SEEK_CUR); // Move the pointer left of D

 		while((ch != EOF) && (ch != 'Q'))
 		{
 			ch = fgetc(fp) ;
 			if ((ch != ':') && (ch != ' ') && (ch != '\n') && (ch != '\t')  && (ch  != EOF) && (ch  != 'D') && (ch  != 'Q'))
 			{
 				ch_str[0] = ch ;
 				strcat(D_string,ch_str);
 			}
 		}

 		fseek (fp,-1,SEEK_CUR);

 		NumOfCells += (strlen(Q_string))*(strlen(D_string));
 		pthread_mutex_unlock (&ReadMutex);


 	// Execute the Smith-Waterman algorithm for each pair of strings
 	Traceback_Steps += Smith_Waterman (Q_string,D_string,ScoringMatrix);

 		free (Q_string);
  	free (D_string);
 }

 void Read_Dataset ()
 {
    char* str= (char *) malloc(MAXCHAR);
	assert(str!=NULL);
    int count = 0 ;
    char* delim = (char *) malloc(2) ;
	assert(delim!=NULL);
    char* token = (char *) malloc(MAXCHAR);
	assert(token!=NULL);
	int pairs ;

    strcpy(delim,":");

	fp = fopen(path,"r"); // Input file

    if (fp == NULL)
    {
        printf("Could not open file %s",path);
        return ;
    }

	fp1 = fopen(output_filename,"w"); // Output file

	if (fp1 == NULL)
	{
			printf("Could not open file %s",output_filename);
			return ;
	}

    while ((fgets(str, MAXCHAR, fp) != NULL) && (count <=3)) // For starters , read only the variable lines
    {
	   	// get the first token
   		token = strtok(str,delim);

  		// walk through other tokens
	   	if (token != NULL)  // Extract the number from each line
	   	{
	    	token = strtok(NULL,delim);
	    	if (count == 0) pairs = atoi (token) ;
			if (count == 1) Q_min = atoi (token) ;
			if (count == 2) Q_max = atoi (token) ;
			if (count == 3) D_Sz_all = atoi (token) ;
	   	}
	   	count++;
	}

	printf ("A) Total Number of Pairs : %d \n",pairs);
	free(delim);
	free(str);

	int NeededThreads = NumOfThreads;
	if (pairs < NumOfThreads) NeededThreads = pairs ;

	pthread_t threads[NeededThreads] ;
	int* ScoringMatrix  ; 

	for (int i=0 ; i< pairs ; i++) // Create thread for each pair to calculate it
	{
			ScoringMatrix = (int *) malloc((D_Sz_all + 1) * (Q_max + 1) * sizeof(int));
			assert (ScoringMatrix != NULL) ; 
			int rc = pthread_create(&threads[i % NeededThreads], NULL, &Algorithm_Per_Pair, (void*) ScoringMatrix) ;

			if (rc)
		   	{
						printf("Problem \n");
						exit(-1);
			}

			if (pthread_join(threads[i % NeededThreads],NULL) )
			{
						printf("error joining thread.");
						abort();
			}

			free (ScoringMatrix);
	}

	printf ("B) Total Number of Cells that take value : %lld \n",NumOfCells);
	printf ("C) Total Number of Traceback Steps : %ld \n",Traceback_Steps);

	fclose(fp1);
	fclose(fp);
	return;
}

int main(int argc, char *argv[])
{
	double time0 = gettime();
	double time1 = 0 ;
	double elapsed_time = 0 ;
	match = atoi(argv[3]);
	mismatch = atoi(argv[4]);
	gap = atoi(argv[5]);
	NumOfThreads = atoi(argv[6]);
	output_filename = (char *) malloc(MAXCHAR);
	assert(output_filename!=NULL);
	path = (char *) malloc(MAXCHAR);
	assert(path!=NULL);
	strcpy(output_filename,"Report_");
	output_filename = strcat(output_filename,argv[1]);
	output_filename = strcat(output_filename,"_PTH_");
	output_filename = strcat(output_filename,argv[6]);
	output_filename = strcat(output_filename,".txt");
	path = argv[2];

	// Read the Dataset
	Read_Dataset();

	// Print the performance indicators
	time1  = gettime();
	elapsed_time = time1 - time0 ;
	printf("D) Program Execution Time : %lf sec\n",elapsed_time);
	printf("E) Scoring Matrix Calculation Time : %lf sec\n",time_E/NumOfThreads);
	printf("F) Traceback Time : %lf sec\n",time_F/NumOfThreads);
	printf("G) CUPS (Program Execution Time) : %lf CUPS\n",NumOfCells/elapsed_time);
	printf("H) CUPS (Cell Calculation Time) : %lf CUPS\n",NumOfCells/time_E);
	pthread_exit(NULL);
}
