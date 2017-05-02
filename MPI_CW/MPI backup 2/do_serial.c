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

#define MAXITER 100000

int main(int argc, char** argv)
{
  /*
    The main part of the serial code. It executes the conversion of the pictures
    with just edges defined to the mode detailed picture. It uses just a single
    route (process) to achieve this goal.
  */
  int i,j,iter=0;
  int rank=0, size=0;
  float local_sum, local_avg;
  char *ptr;
  double times[MAXITER];

  // set defaults, in case inputs are not given.
  int M=256, N=192, DELTA_FREQ=100, AVG_FREQ=200;
  float MAX_DELTA=0.05;
  int val_array[10] = {1, 2, 5, 10, 25, 50 , 75, 100, 150, 200};
  float m_delta[9] = {1.0, 0.75, 0.5, 0.25, 0.1, 0.075, 0.05, 0.025, 0.01};
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
      DELTA_FREQ = val_array[(int)strtol(argv[arg+1], &ptr, 10)-1];
	printf("!!!%d\n",(int)strtol(argv[arg+1], &ptr, 10));
    } else if (strcmp(argv[arg], "Af") == 0)
    {
      AVG_FREQ = val_array[(int)strtol(argv[arg+1], &ptr, 10)-1];
    } else if (strcmp(argv[arg], "MD") == 0)
    {
      MAX_DELTA = m_delta[(int)strtol(argv[arg+1], &ptr, 10)-1];
	printf("!!%f\n", MAX_DELTA);
    }
  }

  //create the arrays which will be used for image processing
  float max_delta = MAX_DELTA + 1;
  float **masterbuf = make_2d_dyn(M, N);
  float **edge = make_2d_dyn(M+2, N+2);
  float **old = make_2d_dyn(M+2, N+2);
  float **new = make_2d_dyn(M+2, N+2);
  float **delta = make_2d_dyn(M, N);
  fflush(stdout);

  //start measuring the time & create the files to be read & written to
  clock_t start_time = clock();
  char filename[16], filename_end[22];
  sprintf(filename, "edge%dx%d.pgm", M, N);
  sprintf(filename_end, "edge%dx%d_%.3f.pgm", M, N, MAX_DELTA);

  printf("Reading %s\n", filename);
  pgmread(filename, *masterbuf, M, N);
  printf("Finished reading\n");

  //print statement used for post-processing
  printf("init size=-1 MP=-1 NP=-1 M=%d N=%d max_delta=%f delta_freq=%d avg_freq=%d\n", M, N, MAX_DELTA, DELTA_FREQ, AVG_FREQ);

  //fill the initial arrays with the edge values & the actual values from the
  // picture
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

  // start the main loop. Iterate over either the maximum number MAXITER
  // or until the maximum change (max_delta) is smaller than the threashold
  while (iter < MAXITER && MAX_DELTA < max_delta)
  {
    times[iter] = clock();
    // calculate the new values based on the old ones
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
/*
    // Calculate the average here
    if (iter % AVG_FREQ == 0)
    {
      local_avg = myAverage(*new, (M+2) * (N+2));
      printf("local_avg=%.10f iter=%d size=%d limit_delta=%f\n", local_avg, iter, size, MAX_DELTA);
    }
*/
    // reassign the new values to the old so the iteration can restart
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

  //write the final values to the masterbuf
  for (int i=0; i<M; i++)
  {
    for (int j=0; j<N; j++)
    {
      masterbuf[i][j] = old[i+1][j+1];
    }
  }

  //write the latest values generated to the file
  printf("Writing\n");
  pgmwrite(filename_end, *masterbuf, M, N);
  printf("Finished writing\n");

  clock_t end_time = clock();

  // print values of the time for performance
  printf("iter=%d overall_time=%f iter_time=%f total_loop=%f max_delta=%f\n", iter, (double)(end_time-start_time)/(double)CLOCKS_PER_SEC, avg_time(times, iter), avg_time(times, iter) * iter, max_delta);

  // free the memory
  free(masterbuf);
  free(edge);
  free(old);
  free(new);
  free(delta);
}


// ---------------------------
// Additional functions
// ---------------------------

float** make_2d_dyn(int rows, int cols)
{
  /*
    Function that makes the dynamic allocation of memory for a 2D array of
    M rows & N columns.
  */
  float **myArray = (float **) malloc(rows * sizeof(float *));
	myArray[0] = (float *) malloc(rows * cols * sizeof(float));
	for (int i = 1; i<rows; i++)
		myArray[i] = myArray[0] + i * cols;
  return myArray;
}

double mySum(void *myArray, int size)
{
  /*
    Function that calculates the sum of an array of a given size.
  */
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
  // Function that calculates the average of an array of a given size
  return mySum(myArray, size)/ (double)size;
}

float maxValue(float *myArray, int size)
{
  // Function that calculates the maximum value from an array.
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
  /*
  Function that calculates the average time for an iteration. Different from
  myAverage & mySum because a different type of array was given.
  */

  double sum = 0;
  double *x = (double *)myArray;
  for (int i=0; i<size; i++)
  {
    sum += x[i];
  }
  sum = sum / (double)size;
  return sum;
}
