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
#define P 4

#define MP M/P
#define NP N

#define ITER 1000

float masterbuf[M][N];
float buf[MP][NP];
float edge[MP+2][NP+2];
float old[MP+2][NP+2];
float new[MP+2][NP+2];
double times[P][ITER];

int main(void)
{

  int i,j,k;
  int tag, size, rank, next, prev;
  int start_time, stop_time;
  MPI_Comm comm;
  MPI_Status status;
  MPI_Request *requests;
  MPI_Status *statuses;

  char *filename;
  filename = "edge256x192.pgm";

  comm = MPI_COMM_WORLD;
  tag = 1;

  MPI_Init(NULL,NULL);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  next = rank + 1;
  prev = rank - 1;

  if (size != P)
  {
    if (rank ==0)
    {
      printf("ERROR, size != P\n");
    }
    MPI_Finalize();
    exit(-1);
  }

  if (next >= size)
    {
      next = MPI_PROC_NULL;
    }

  if (prev < 0)
    {
      prev = MPI_PROC_NULL;
    }

    if (rank == 0)
    {
      printf("Reading\n");
      pgmread(filename, masterbuf, M, N);
      printf("Finished reading\n");

    }

  MPI_Scatter(masterbuf, MP * NP, MPI_FLOAT, &buf, MP * NP, MPI_FLOAT, 0, comm);


  for (int i=0; i<MP+2; i++)
  {
    for (int j=0; j<NP+2; j++)
    {
      edge[i][j] = 255;
    }
  }
  for (i=1; i<MP+1; i++)
  {
    for (j=1; j<NP+1; j++)
    {
      edge[i][j] = buf[i-1][j-1];
    }
  }

  for (int i=0; i<MP+2; i++)
  {
      printf("i=%d, rank=%d\n", i, rank);
    for (int j = 0; j < NP+2; j++)
    {
      old[i][j] = edge[i][j];
    }
  }
    printf("finished rank %d \n", rank);



  for (int k=0; k<ITER; k++)
  {
    times[rank][k] = MPI_Wtime();
    printf("in for loop %d\n", rank);
    MPI_Sendrecv(&old[MP][1], NP, MPI_FLOAT, next, 1, &old[0][1],  NP, MPI_FLOAT, prev, 1, MPI_COMM_WORLD, &status);
    printf("past first sendrecv %d\n", rank);
    MPI_Sendrecv(&old[1][1], NP, MPI_FLOAT, prev, 2, &old[MP+1][1], NP, MPI_FLOAT, next, 2, MPI_COMM_WORLD, &status);
    printf("past second sendrecv %d\n", rank);

    //communicate(&old, prev,next, MP, NP);

    for (int i=1; i<MP+1; i++)
    {
      for (int j=1; j<NP+1; j++)
      {
        new[i][j] =  (old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - edge[i][j]) * 0.25;
      }
    }

    for (int i=1; i<MP+1; i++)
    {
      for (int j=1; j<NP+1; j++)
      {
        old[i][j] = new[i][j];
      }
    }

    for (int i=0; i<MP; i++)
    {
      for (int j=0; j<NP; j++)
      {
        buf[i][j] = old[i+1][j+1];
      }
    }

      times[rank][k] -= MPI_Wtime();
  }
  printf("finished %d n=%d, p=%d\n", rank, next, prev);

  MPI_Gather(buf, MP * NP, MPI_FLOAT, masterbuf, MP * NP, MPI_FLOAT, 0, comm);
  printf("finished AFTER gather %d n=%d, p=%d\n", rank, next, prev);
  if (rank == 0)
    {
      filename = "edge256x192_1.pgm";
      printf("Writing\n");
      pgmwrite(filename, masterbuf, M, N);
      printf("Finished writing\n");
    }

  MPI_Finalize();

  for (int i=0; i<size; i++)
  {
    printf("Tread %d %f r=%d\n", i, average(times[i], i, ITER), rank);
  }

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
