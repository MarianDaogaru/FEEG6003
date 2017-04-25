#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

void pgmread(char *filename, void *vx, int nx, int ny);
void pgmwrite(char *filename, void *vx, int nx, int ny);


#define M 192
#define N 128
#define P 4

#define MP M/P
#define NP N

float masterbuf[M][N];
float buf[MP][NP];
float edge[MP+2][NP+2];
float old[MP+2][NP+2];
float new[MP+2][NP+2];

int main(void)
{

  int iter = 50;
  int i,j,k;
  int tag, size, rank, next, prev;
  MPI_Comm comm;
  MPI_Status status;
  MPI_Request request;

  char *filename;
  filename = "edge192x128.pgm";

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

  MPI_Scatter(masterbuf, MP * NP, MPI_FLOAT, buf, MP * NP, MPI_FLOAT, 0, comm);

  for (int i=0; i<MP+2; i++)
  {
    for (int j=0; j<NP+2; j++)
    {
      edge[i][j] = 255.0;
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



  for (int k=0; k<iter; k++)
  {
    printf("in for loop %d\n", rank);
    MPI_Sendrecv(&old[(rank+1)*MP][1], NP, MPI_FLOAT, next, 1, &old[rank*MP][1],  NP, MPI_FLOAT, prev, 1, MPI_COMM_WORLD, &status);
    printf("past first sendrecv %d\n", rank);
    MPI_Sendrecv(&old[1+rank*MP][1], NP, MPI_FLOAT, prev, 2, &old[(rank+1)*MP+1][1], NP, MPI_FLOAT, next, 2, MPI_COMM_WORLD, &status);
    printf("past second sendrecv %d\n", rank);
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
  }
  MPI_Gather(buf, MP * NP, MPI_FLOAT, masterbuf, MP * NP, MPI_FLOAT, 0, comm);
  if (rank == 0)
    {
      filename = "edge192x128_1.pgm";
      printf("Writing\n");
      pgmwrite(filename, buf, M, N);
      printf("Finished writing\n");
    }

  MPI_Finalize();
}
