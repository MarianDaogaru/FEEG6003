#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#include "pgmio.h"
double average(double *times, int r, int iter);


/*
void pgmread(char *filename, void *vx, int nx, int ny);
void pgmwrite(char *filename, void *vx, int nx, int ny);
//void communicate(void *oldbuf, int prev, int next, int MP, int NP);
void datread(char *filename, void *vx, int nx, int ny);
*/
#define M 256
#define N 192
#define P 6

#define MAXITER 50000
#define DELTA_FREQ 100
#define MAX_DELTA 0.05
#define AVG_FREQ 200

float masterbuf[M][N];
float buf[M][N];
float edge[M+2][N+2];
float old[M+2][N+2];
float new[M+2][N+2];

int main(void)
{

  int i,j,iter;

  char filename[16], filename_end[22];
  sprintf(filename, "edge%dx%d.pgm", M, N);
  sprintf(filename_end, "edge%dx%d_%.3f.pgm", M, N, MAX_DELTA);

  printf("Reading\n");
  pgmread(filename, masterbuf, M, N);
  printf("Finished reading\n");

  for (int i=0; i<M; i++)
  {
    for (int j=0; j<N; j++)
    {
      edge[i][j] = 255;
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



  while (iter < MAXITER)
  {

    //communicate(&old, prev,next, MP, NP);

    for (int i=1; i<M+1; i++)
    {
      for (int j=1; j<N+1; j++)
      {
        new[i][j] =  (old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - edge[i][j]) * 0.25;
      }
    }

    for (int i=1; i<M+1; i++)
    {
      for (int j=1; j<N+1; j++)
      {
        old[i][j] = new[i][j];
      }
    }
    for (int i=0; i<M; i++)
    {
      for (int j=0; j<N; j++)
      {
        buf[i][j] = old[i+1][j+1];
      }
    }
  }

  for (int i=0; i<M; i++)
  {
    for (int j=0; j<N; j++)
    {
      masterbuf[i][j] = buf[i][j];
    }
  }
  printf("Writing\n");
  pgmwrite(filename_end, masterbuf, M, N);
  printf("Finished writing\n");

}

double average(double *times, int r, int iter)
{
  double sum = 0;
  //double *timeZ = (double *)times;

  for (int i = 0; i<iter; i++)
  {
    //printf("%f\n", times[i]);
    sum -= times[i];
  }

  return sum / (double)iter;
}
