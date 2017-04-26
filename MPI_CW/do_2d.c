#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#include "pgmio.h"
void choose_MN(void *myArray);
double mySum(void *myArray, int size);
double myAverage(void *myArray, int size);
float maxValue(float myArray[], int size);

/*
//void communicate(void *oldbuf, int prev, int next, int MP, int NP);
*/
#define M 256
#define N 192
#define P 24

#define MAXITER 10000
#define DELTA_FREQ 10
#define MAX_DELTA 0.05
#define AVG_FREQ 200

float masterbuf[M][N];

double times[P][MAXITER];
float max_delta_thread[P];
float max_delta = MAX_DELTA + 1;

int main(void)
{
  int MN[2];
  int MP, NP;
  int i,j;
  int iter=0;
  int size, rank, next, prev;
  float local_sum, global_sum[P], local_avg, global_avg;
  MPI_Comm comm;
  MPI_Status status;
  MPI_Request *requests;
  MPI_Status *statuses;

  char *filename;
  filename = "edge256x192.pgm";

  fflush(stdout);
  choose_MN(MN);
  MP = MN[0];
  NP = MN[1];
  printf("MP = %d, NP=%d\n", MP, NP);
  float buf[MP][NP];
  float edge[MP+2][NP+2];
  float old[MP+2][NP+2];
  float new[MP+2][NP+2];
  float delta[MP*NP];

  comm = MPI_COMM_WORLD;
  fflush(stdout);
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
      //printf("i=%d, rank=%d\n", i, rank);
    for (int j = 0; j < NP+2; j++)
    {
      old[i][j] = edge[i][j];
    }
  }

  printf("finished rank %d, before while, iter = %d, max_detal = %f \n", rank, iter, max_delta);



  while (iter < MAXITER && MAX_DELTA < max_delta)
  {
    times[rank][iter] = -MPI_Wtime();
    //printf("in for loop %d\n", rank);
    MPI_Sendrecv(&old[MP][1], NP, MPI_FLOAT, next, 1, &old[0][1],  NP, MPI_FLOAT, prev, 1, MPI_COMM_WORLD, &status);
    //printf("past first sendrecv %d\n", rank);
    MPI_Sendrecv(&old[1][1], NP, MPI_FLOAT, prev, 2, &old[MP+1][1], NP, MPI_FLOAT, next, 2, MPI_COMM_WORLD, &status);
    //printf("past second sendrecv %d\n", rank);

    //communicate(&old, prev,next, MP, NP);

    for (int i=1; i<MP+1; i++)
    {
      for (int j=1; j<NP+1; j++)
      {
        new[i][j] =  (old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - edge[i][j]) * 0.25;
      }
    }

    // AICI CALC DELTA
    if (iter % DELTA_FREQ == 0)
    {
      for (int i=0; i<MP; i++)
      {
        for (int j=0; j<NP; j++)
        {
          delta[i*NP + j] = fabs(old[i+1][j+1] - new[i+1][j+1]);
        }
      }
      MPI_Allreduce(delta, max_delta_thread, P, MPI_FLOAT, MPI_MAX, comm);

      for (int i=0; i<P; i++) printf("mD = %f, i=%d\n", max_delta_thread[i], i);
      //printf("maxValue = %f \n", maxValue(max_delta_thread, P));
      max_delta = maxValue(max_delta_thread, P);
      printf("max_delta = %f\n", max_delta);
    }

    // AVG AICI
    if (iter % AVG_FREQ == 0)
    {
      printf("in avg, iter = %d", iter);
      local_avg = myAverage(buf, MP * NP);
      local_sum = mySum(buf, MP * NP);
      printf("local_sum = %f, local_avg = %f, rank=%d, iter=%d\n", local_sum, local_avg, rank, iter);
      MPI_Reduce(&local_sum, &global_sum, 1, MPI_FLOAT, MPI_SUM, 0, comm);
      printf("dupa reduce rank = %d \n", rank);
      if (rank == 0)
      {
          printf("in if\n");
        for (int i = 0; i<P; i++)
        {
          printf("gs = %f\n", global_sum[i]);
        }
        global_avg = mySum(global_sum, P) / ((double) M * N);
        printf("dupa avb \n");
        printf("At iter = %d, avg = %f, sum = %f, rank = %d\n", iter, global_avg, mySum(global_sum, P), rank);
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

      times[rank][iter] += MPI_Wtime();
      iter++;
  }
  printf("finished %d n=%d, p=%d\n", rank, next, prev);
  if (MAX_DELTA > max_delta)
  {
    printf("Finished earlier, iter=%d, max_delta=%f \n", iter, max_delta);
  }
  MPI_Gather(buf, MP * NP, MPI_FLOAT, masterbuf, MP * NP, MPI_FLOAT, 0, comm);
  printf("finished AFTER gather %d n=%d, p=%d\n", rank, next, prev);
  if (rank == 0)
    {
      printf("global avg = %f\n", global_avg);
      filename = "edge256x192_1.pgm";
      printf("Writing\n");
      pgmwrite(filename, masterbuf, M, N);
      printf("Finished writing\n");
    }
  MPI_Finalize();


  for (int i=0; i<size; i++)
  {
    if (i == rank) printf("Tread %d %f r=%d\n", i, myAverage(times[i], iter), rank);
  }
}


// ---------------------------
// Additional functions
// ---------------------------


void choose_MN(void *myArray)
{
  int M_i, N_i;
  int MP, NP;
  int *MN = (int *)myArray;
  int done_choosing_MN = 0, cont = 1;

  if (fmod(P, sqrt(P)) == 0)
  {
    // is square of another number
    MN[0] = M / (int)sqrt(P);
    MN[1] = N / (int)sqrt(P);
    printf("MP=%d, NP=%d\n", MN[0], MN[1]);
  }
  else
  {
    // start mapping the 2 closest intergers that multiply give P
    while (done_choosing_MN == 0)
    {
      M_i = P / sqrt(P) + cont;
      printf("in while Mi=%d, Ni=%d, cont=%d\n", M_i, N_i, cont);
      while(P % M_i != 0)
      {
        cont++;
        M_i = P / sqrt(P) + cont;
        printf("Mi=%d, Ni=%d, cont=%d\n", M_i, N_i, cont);
      }
      printf("after 2nd while\n");
      N_i = P / M_i;
      cont++; //if they divide, but not divide in the next if, this must be increased so you actually get the next number

      if (N_i == 1)
      {
        // N_i will go down as M_i is increased, so we might reach the end here
        if (M % M_i == 0)
        {
          // see if we actually divide this part, otherwise go to the next axis and divide that
          MN[0] = M / M_i;
          MN[1] = N;
        }
        else
        {
          MN[0] = M;
          MN[1] = N / M_i;
        }
        done_choosing_MN = 1; // we are actually done choosing now
      }
      else if(fmod(M, M_i) == 0 && fmod(N, N_i) == 0)
      { //if the numbers are not 0 and they both divide these axis, we will chose these ones, if not, we will swap the division,
        MN[0] = M / M_i;
        MN[1] = N / N_i;
        printf("MP=%d, NP=%d 22 \n", MN[0], MN[1]);
        done_choosing_MN = 1;
      }
      else if (fmod(M, N_i) == 0 && fmod(N, M_i) == 0)
      {
        MN[0] = M / N_i;
        MN[1] = N / M_i;
        printf("MP=%d, NP=%d 23 \n", MN[0], MN[1]);
        done_choosing_MN = 1;
      }
    }
  }
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

float maxValue(float myArray[], int size)
{
  int i;
  float max = 0;
  for (i = 0; i<size; i++)
  {
    if (max < myArray[i]) max = myArray[i];
  }
  return max;
}
