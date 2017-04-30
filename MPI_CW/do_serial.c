#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "pgmio.h"

double average(double *times, int r, int iter);
double mySum(void *myArray, int size);
double myAverage(void *myArray, int size);
float maxValue(float *myArray, int size);
double avg_time(double *myArray, int size);
float** make_2d_dyn(int rows, int cols);


//#define M 256
//#define N 192

#define MAXITER 5000
//#define DELTA_FREQ 100
//#define MAX_DELTA 0.05
//#define AVG_FREQ 200


int main(int argc, char** argv)
{
  int i,j,iter=0;
  int rank=0, size=0;
  float local_sum, local_avg;
  char *ptr;
  double times[MAXITER];

  // set defaults
  int M=256, N=192, DELTA_FREQ=100, AVG_FREQ=200;
  float MAX_DELTA=0.05;

  //read the inputs
  for (int arg = 0; arg < argc; arg++)
  {
    if (strcmp(argv[arg], "M") == 0)
    {
      M = (int)strtol(argv[arg+1], &ptr, 10);
    } else if (strcmp(argv[arg], "N") == 0)
    {
      N = (int)strtol(argv[arg+1], &ptr, 10);
    } else if (strcmp(argv[arg], "Df") == 0)
    {
      DELTA_FREQ = (int)strtol(argv[arg+1], &ptr, 10);
    } else if (strcmp(argv[arg], "Af") == 0)
    {
      AVG_FREQ = (int)strtol(argv[arg+1], &ptr, 10);
    } else if (strcmp(argv[arg], "MD") == 0)
    {
      MAX_DELTA = (double)strtol(argv[arg+1], &ptr, 10);
    }
  }

  float max_delta = MAX_DELTA + 1;
  fflush(stdout);
  float **masterbuf = make_2d_dyn(M, N);

  fflush(stdout);
  float **edge = make_2d_dyn(M+2, N+2);
  fflush(stdout);
  float **old = make_2d_dyn(M+2, N+2);
  fflush(stdout);
  float **new = make_2d_dyn(M+2, N+2);
  fflush(stdout);
  float **delta = make_2d_dyn(M, N);
  fflush(stdout);

  clock_t start_time = clock();
  char filename[16], filename_end[22];
  sprintf(filename, "edge%dx%d.pgm", M, N);
  sprintf(filename_end, "edge%dx%d_%.3f.pgm", M, N, MAX_DELTA);

  printf("Reading %s\n", filename);
  pgmread(filename, *masterbuf, M, N);
  printf("Finished reading\n");

  for (int i=0; i<M+2; i++)
  {
    for (int j=0; j<N+2; j++)
    {
      edge[i][j] = 255.0;
    }
  }
  for (i=1; i<M+1; i++)
  {
    for (j=1; j<N+1; j++)
    {
      edge[i][j] = masterbuf[i-1][j-1];
    }
  }

  for (int i=0; i<M+2; i++)
  {
    for (int j = 0; j < N+2; j++)
    {
      old[i][j] = edge[i][j];
    }
  }



  while (iter < MAXITER && MAX_DELTA < max_delta)
  {
    //printf("iter=%d\n",iter);
    //communicate(&old, prev,next, MP, NP);
    times[iter] = clock();
    for (int i=1; i<M+1; i++)
    {
      for (int j=1; j<N+1; j++)
      {
        new[i][j] =  (old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - edge[i][j]) * 0.25;
      }
    }


    // Calculate the maximum DELTA here
    if (iter % DELTA_FREQ == 0)
    {
      for (int i=0; i<M; i++)
      {
        for (int j=0; j<N; j++)
        {
          delta[i][j] = fabsf(old[i+1][j+1] - new[i+1][j+1]);
        }
      }
      max_delta = maxValue(*delta, M*N);
      printf("max_delta=%.10f iter=%d rank=%d size=%d limit_delta=%f\n", max_delta, iter, rank, size, MAX_DELTA);
    }

    // Calculate the average here
    if (iter % AVG_FREQ == 0)
    {
      local_avg = myAverage(*new, (M+2) * (N+2));
      printf("local_avg=%.10f iter=%d size=%d limit_delta=%f\n", local_avg, iter, size, MAX_DELTA);
    }

    for (int i=1; i<M+1; i++)
    {
      for (int j=1; j<N+1; j++)
      {
        old[i][j] = new[i][j];
      }
    }

    iter++;
    times[iter-1] -= clock();
  }

  for (int i=0; i<M; i++)
  {
    for (int j=0; j<N; j++)
    {
      masterbuf[i][j] = old[i+1][j+1];
    }
  }

  printf("Writing\n");
  pgmwrite(filename_end, *masterbuf, M, N);
  printf("Finished writing\n");

  clock_t end_time = clock();

  printf("iter=%d overall_time=%f iter_time=%f total_loop=%f max_delta=%f\n", iter, (double)(end_time-start_time)/(double)CLOCKS_PER_SEC, avg_time(times, iter), avg_time(times, iter) * iter, max_delta);

  free(masterbuf);
  free(edge);
  free(old);
  free(new);
  free(delta);
}
float** make_2d_dyn(int rows, int cols)
{
  float **myArray = (float **) malloc(rows * sizeof(float *));
	myArray[0] = (float *) malloc(rows * cols * sizeof(float));
	for (int i = 1; i<rows; i++)
		myArray[i] = myArray[0] + i * cols;
  return myArray;
}

double mySum(void *myArray, int size)
{
  double sum = 0;
  float *x = (float *)myArray;
  for (int i = 0; i<size; i++)
  {
    sum += x[i];
  }
  return sum;
}

double myAverage(void *myArray, int size)
{
  return mySum(myArray, size)/ (double)size;
}

float maxValue(float *myArray, int size)
{
  int i;
  float max = 0;
  for (i = 0; i<size; i++)
  {
    if (max < myArray[i]) max = myArray[i];
  }
  return max;
}

double avg_time(double *myArray, int size)
{
  double sum = 0;
  double *x = (double *)myArray;
  for (int i=0; i<size; i++)
  {
    sum += x[i];
  }
  sum = sum / (double)(size * CLOCKS_PER_SEC);
  return sum;
}
