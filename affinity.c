#include <stdio.h>
#include <math.h>
#include <omp.h>

#define N 729
#define reps 100
#define MAX_THREAD 24

double a[N][N], b[N][N], c[N];
int jmax[N];
int num_threads;  //number of threads
int iter_start_loc[MAX_THREAD];
int iter_stop_loc[MAX_THREAD];
int start_loc[MAX_THREAD], stop_loc[MAX_THREAD];

void init1(void);
void init2(void);
void loop1(int start_p, int end_p);
void loop2(int start_p, int end_p);
void create_chunk(int thread);
void loop1_main(void);
void loop2_main(void);
void valid1(void);
void valid2(void);


int main(int argc, char *argv[]) {

  double start1,start2,end1,end2;
  int r;
  int i, j, k;
  int thread_iter = 0;
  int remaining = 0;

  // assign the points to threads
  #pragma omp parallel default(none) shared(num_threads)
  {
    #pragma omp master
    {
      num_threads = omp_get_num_threads();
    }
  }

  thread_iter = N / num_threads;  // get how many iterations are per threads equally
  remaining = N % num_threads;  // get the remaining value of iterations if N does not divide nicely by num_threads
  start_loc[0] = 0;
  stop_loc[0] = thread_iter;
  if (remaining != 0)
    {
      stop_loc[0] += 1;
    }
  for (i=1; i<num_threads; i++)
    {
        start_loc[i] = stop_loc[i-1];
        stop_loc[i] = start_loc[i] + thread_iter;
        if (i < remaining)
        {
          stop_loc[i] += 1;  // add one to the first threads, so the remaining is distributed nicely
        }
    }

    // first loop
    init1();
    start1 = omp_get_wtime();

    for (int r=0; r<reps; r++){
      for (int i=0; i<num_threads; i++)
      {
        // make the initial positions
        iter_start_loc[i] = start_loc[i];
        iter_stop_loc[i] = start_loc[i] + ceil((stop_loc[i] - start_loc[i]) / (double)num_threads);
      }
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
      for (int i=0; i<num_threads; i++)
      {
        // make the initial starting positions
        iter_start_loc[i] = start_loc[i];
        iter_stop_loc[i] = start_loc[i] + ceil((stop_loc[i] - start_loc[i]) / (double)num_threads);
      }
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

void create_chunk(int thread)
{
  int chunk;
  chunk = ceil((stop_loc[thread] - iter_stop_loc[thread]) / (double)omp_get_num_threads());
  iter_start_loc[thread] = iter_stop_loc[thread];
  iter_stop_loc[thread] = iter_stop_loc[thread] + chunk;
}

// the initial loops
void loop1(int start_p, int end_p)
{
  int i,j;
  for (i=start_p; i<end_p; i++)
  {
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

// the affinity loop functions

void loop1_main(void)
{
  int i, j;
  int thread_sum = 0;


  #pragma omp parallel private(i,j) shared(thread_sum)
  {
    int thread_num = omp_get_thread_num();
    int max_num_threads = omp_get_num_threads();

    int min_thread=0;
    int min_iter=0;
    int ub;
    int lb;

    lb = iter_start_loc[thread_num];
    ub = iter_stop_loc[thread_num];

    while (thread_sum < max_num_threads)
    {
      loop1(lb, ub);
      #pragma omp critical
      {
        if (iter_start_loc[thread_num] != stop_loc[thread_num])
        {
          create_chunk(thread_num);
          lb = iter_start_loc[thread_num];
          ub = iter_stop_loc[thread_num];
        }
        else
        {
          min_iter = 0;
          for (i=0; i<max_num_threads; i++)
          {
              if (min_iter < stop_loc[i] - iter_stop_loc[i])
              {
                min_thread = i;
                min_iter = stop_loc[i] - iter_stop_loc[i];
              }
          }
          if (min_iter != 0)
          {
            create_chunk(min_thread);
            lb = iter_start_loc[min_thread];
            ub = iter_stop_loc[min_thread];
          }
          else
          {
            thread_sum = max_num_threads;
          }
        }
      }
    }
  }
}

void loop2_main(void)
{
  int i, j;
  int thread_sum = 0;


  #pragma omp parallel private(i,j) shared(thread_sum)
  {
    int thread_num = omp_get_thread_num();
    int max_num_threads = omp_get_num_threads();

    int min_thread=0;
    int min_iter=0;
    int ub;
    int lb;

    lb = iter_start_loc[thread_num];
    ub = iter_stop_loc[thread_num];

    while (thread_sum < max_num_threads)
    {
      loop2(lb, ub);
      #pragma omp critical
      {
        if (iter_start_loc[thread_num] != stop_loc[thread_num])
        {
          create_chunk(thread_num);
          lb = iter_start_loc[thread_num];
          ub = iter_stop_loc[thread_num];
        }
        else
        {
          min_iter = 0;
          for (i=0; i<max_num_threads; i++)
          {
              if (min_iter < stop_loc[i] - iter_stop_loc[i])
              {
                min_thread = i;
                min_iter = stop_loc[i] - iter_stop_loc[i];
              }
          }
          if (min_iter != 0)
          {
            create_chunk(min_thread);
            lb = iter_start_loc[min_thread];
            ub = iter_stop_loc[min_thread];
          }
          else
          {
            thread_sum = max_num_threads;
          }
        }
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
  printf("Loop 1 check: Sum of a is %lf\n", suma);
}

void valid2(void) {
  int i;
  double sumc;

  sumc= 0.0;
  for (i=0; i<N; i++){
    sumc += c[i];
  }
  printf("Loop 2 check: Sum of c is %f\n", sumc);
}
