#include <stdio.h>
#include <math.h>
#include <omp.h>

#define N 729
#define reps 100


double a[N][N], b[N][N], c[N];
int jmax[N];
int ckk = 4;  // chunk size for schedule. 4 is best for loop2 as found before
int num_threads;  //number of threads
int num_iter_thread[24];  // number of iterations for each thread. 24 cause that is the maximum in this exercise
int num_chunks;  // chunk size for each thread
int iter_thread[24][N][3];  // individual iter count for each thread & bool values
omp_lock_t locks[24][N];  // omp lock variable


void init1(void);
void init2(void);
void loop1(int start_p, int end_p);
void loop2(int start_p, int end_p);
void valid1(void);
void valid2(void);


int main(int argc, char *argv[]) {

  double start1,start2,end1,end2;
  double l1_t_def, l2_t_def;
  int r;
  int i, j, k;
  int thread_iter = 0;
  int remaining = 0;
  int start_location = 0;
  int origin_iter = 0;
  int chunk = 0;

  // assign the points to threads
  #pragma omp parallel default(none) shared(num_threads)
  {
    #pragma omp master
    {
      num_threads = omp_get_num_threads();
    }
  }

  // loop over the max number of threads
  for (i=0; i<24; i++)
  {
    // loop over the total number of iterations
    for (j=0; j<N; j++)
    {
      omp_init_lock(&locks[i][j]);  // initialise the locks
      // get over the start, end, and bool val
      for (k=0; k<3; k++)
      {
        // set everything to 0
        iter_thread[i][j][k] = 0;
      }
    }
  }


    thread_iter = N / num_threads;  // get how many iterations are per threads equally
    remaining = N % num_threads;  // get the remaining value of iterations if N does not divide nicely by num_threads

    for (i=0; i<num_threads; i++)
    {
        num_iter_thread[i] = thread_iter;  // initialise everything to the first batch of iterations
        if (i < remaining)
        {
          num_iter_thread[i] += 1;  // add one to the first threads, so the remaining is distributed nicely
        }
        origin_iter = num_iter_thread[i];
        printf("ori = %d\n", origin_iter);

        for (j=0; origin_iter>0; j++)
        {
          chunk = origin_iter / num_threads + 1;

          iter_thread[i][j][0] = start_location;  // start location for thread i after j steps
          iter_thread[i][j][1] = start_location + chunk - 1; // end position. -1 so it will actually go there
          iter_thread[i][j][2] = 1;
          start_location += chunk;  // updates the overall start location. it will extend over the threads
          origin_iter -= chunk; // original will go to 0 for each thread and stop
          printf("chunk=%d, i=%d, j=%d origin_iter=%d\n", chunk, i, j, origin_iter);
          printf("1=%d, 2=%d\n", iter_thread[i][j][0], iter_thread[i][j][1]);
        }
    }
    num_chunks = j;  // initial number of different sizes of the chunk
    printf("num_chunks = %d\n", num_chunks);


    // first loop
    init1();
    start1 = omp_get_wtime();

    for (int r=0; r<reps; r++){
      loop1_main();
    }

    end1  = omp_get_wtime();
    valid1();
    printf("Total time for %d reps of loop 1 = %f\n",reps, (float)(end1-start1));


    // second loop
    init2();
    start2 = omp_get_wtime();

    for (r=0; r<reps; r++)
    {
      loop2_main();
    }

    end2  = omp_get_wtime();
    valid2();
    printf("Total time for %d reps of loop 2 = %f\n",reps, (float)(end2-start2));

// final bracket for main
}

void init1(void){
  int i,j;

  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      a[i][j] = 0.0;
      b[i][j] = 3.142*(i+j);
    }
  }
}

void init2(void){
  int i,j, expr;
  for (i=0; i<N; i++){
    expr =  i%( 3*(i/30) + 1);
    if ( expr == 0) {
      jmax[i] = N;
    }
    else {
      jmax[i] = 1;
    }
    c[i] = 0.0;
  }

  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      b[i][j] = (double) (i*j+1) / (double) (N*N);
    }
  }
}

// the initial loops
void loop1(int start_p, int end_p)
{
  int i,j;

  for (i=start_p; i<end_p; i++){
    //printf("thread = %d i = %d \n", omp_get_thread_num(), i);
    for (j=N-1; j>i; j--){
      a[i][j] += cos(b[i][j]);
    }
  }
}

void loop2(int start_p, int end_p)
  {
  int i,j,k;
  double rN2;

  rN2 = 1.0 / (double) (N*N);

  for (i=start_p; i<end_p; i++){
    for (j=0; j < jmax[i]; j++){
      for (k=0; k<j; k++){
	       c[i] += (k+1) * log (b[i][j]) * rN2;
      }
    }
  }
}


//valid part
void valid1(void) {
  int i,j;
  double suma;

  suma= 0.0;
  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      suma += a[i][j];
    }
  }
  //printf("Loop 1 check: Sum of a is %lf\n", suma);
}

void valid2(void) {
  int i;
  double sumc;

  sumc= 0.0;
  for (i=0; i<N; i++){
    sumc += c[i];
  }
  //printf("Loop 2 check: Sum of c is %f\n", sumc);
}
