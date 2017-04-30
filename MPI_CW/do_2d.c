#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#include "pgmio.h"
#include "comms.h"
void choose_MN(void *myArray, int size);
void my_Scatter(void *buf, int m, int n, float masterbuf[m][n], int rank, int MP, int NP);
void my_Gather(void *masterbuf, int MP, int NP, float buf[MP][NP], int rank, int size);
void my_Gather_process(float *lmbuf, int MP, int NP, float buf[MP][NP], int rank);
double mySum(void *myArray, int size);
double myAverage(void *myArray, int size);
float maxValue(float myArray[], int size);
double avg_time(double *myArray, int size);

/*
//void communicate(void *oldbuf, int prev, int next, int MP, int NP);
*/
#define M 256
#define N 192
#define P 1

#define MAXITER 50000
#define DELTA_FREQ 100
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
  int size, rank, left, right, up, down, MP_fact;
  float local_sum, global_sum, local_avg, global_avg;
  double start_time, make_MP_time, choose_neighbours, make_buff, reconstruct_time, barrier_time;
  MPI_Comm comm;
  MPI_Status status;

  char filename[16], filename_end[22];
  sprintf(filename, "edge%dx%d.pgm", M, N);
  sprintf(filename_end, "edge%dx%d_%.3f.pgm", M, N, MAX_DELTA);
  //filename = "edge256x192.pgm";
  //filename_end = "edge256x192_end.pgm";
  //printf("Reading\n");
  pgmread(filename, masterbuf, M, N);
  //printf("Finished reading\n");

  comm = MPI_COMM_WORLD;

  MPI_Init(NULL,NULL);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  start_time = -MPI_Wtime();
  //printf("before choose\n");
  make_MP_time = -MPI_Wtime();
  choose_MN(MN, size);
  MP = MN[0];
  NP = MN[1];
  make_MP_time += MPI_Wtime();
  //printf("choose after\n");
  //printf("MP = %d, NP=%d\n", MP, NP);
  float buf[MP][NP];
  float edge[MP+2][NP+2];
  float old[MP+2][NP+2];
  float new[MP+2][NP+2];
  float delta[MP*NP];

  choose_neighbours = -MPI_Wtime();
  MP_fact = M / MP;
  right = rank + 1;
  left = rank - 1;
  up = rank - MP_fact;
  down = rank + MP_fact;
  if (rank % MP_fact== 0)
  {
    left = MPI_PROC_NULL;
  }
  if ((rank + 1) % MP_fact == 0)
  {
    right = MPI_PROC_NULL;
  }
  if (down >= size)
    {
      down = MPI_PROC_NULL;
    }
  if (up < 0)
    {
      up = MPI_PROC_NULL;
    }
  choose_neighbours += MPI_Wtime();
  //printf("!!!rank=%d, left=%d, right=%d, up=%d, down=%d, MP_fact=%d\n", rank, left, right, up, down, MP_fact);

  if (rank==0) printf("init size=%d MP=%d NP=%d M=%d N=%d max_delta=%f delta_freq=%d avg_freq=%d\n", size, MP, NP, M, N, MAX_DELTA, DELTA_FREQ, AVG_FREQ);

  //MPI_Scatter(masterbuf, MP * NP, MPI_FLOAT, &buf, MP * NP, MPI_FLOAT, 0, comm);
  make_buff = -MPI_Wtime();
  my_Scatter(buf, M, N, masterbuf, rank, MP, NP);

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
      //printf("rank=%d, i=%d, j=%d, buf=%f, buf2=%f\n",rank, i, j, ((float *)buf)[(i-1)*NP+j-1], buf[i-1][j-1]);
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
  //printf("finished rank %d, before while, iter = %d, max_detal = %f \n", rank, iter, max_delta);
  make_buff += MPI_Wtime();

  while (iter < MAXITER && MAX_DELTA < max_delta)
  {
    times[rank][iter] = -MPI_Wtime();

    //printf("in for loop %d, rank=%d\n", iter, rank);
    /*MPI_Sendrecv(&old[MP][1], NP, MPI_FLOAT, right, 1, &old[0][1],  NP, MPI_FLOAT, left, 1, MPI_COMM_WORLD, &status);
    printf("past first sendrecv %d\n", rank);
    MPI_Sendrecv(&old[1][1], NP, MPI_FLOAT, left, 2, &old[MP+1][1], NP, MPI_FLOAT, right, 2, MPI_COMM_WORLD, &status);
    printf("past second sendrecv %d\n", rank);
    */
    communicate_lr(old, left, right, MP, NP);
    //printf("before comms UD, rank=%d\n", rank);
    communicate_ud(old, up, down, MP, NP);

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
      MPI_Allreduce(delta, max_delta_thread, size, MPI_FLOAT, MPI_MAX, comm);

      //for (int i=0; i<size; i++) printf("mD = %f, i=%d\n", max_delta_thread[i], i);
      //printf("maxValue = %f \n", maxValue(max_delta_thread, P));

      max_delta = maxValue(max_delta_thread, size);
      if (rank==0) printf("max_delta=%.10f iter=%d rank=%d size=%d limit_delta=%f\n", max_delta, iter, rank, size, MAX_DELTA);
    }

    // AVG AICI
    if (iter % AVG_FREQ == 0)
    {
      //printf("in avg, iter = %d\n", iter);
      local_avg = myAverage(buf, MP * NP);
      local_sum = mySum(buf, MP * NP);
      //printf("local_sum = %f, local_avg = %f, rank=%d, iter=%d\n", local_sum, local_avg, rank, iter);
      MPI_Reduce(&local_sum, &global_sum, 1, MPI_FLOAT, MPI_SUM, 0, comm);
      //printf("dupa reduce rank = %d \n", rank);
      if (rank == 0)
      {
        //printf("in if\n");
        /*for (int i = 0; i<size; i++)
        {
          printf("gs = %f, rank=%d\n", global_sum,i);
        }*/
        global_avg = global_sum / (double)(M * N);
        //printf("dupa avb \n");
        //printf("At iter = %d, avg = %f, sum = %f, rank = %d\n", iter, global_avg, global_sum, rank); //mySum(global_sum, P)
        printf("local_avg=%.10f iter=%d size=%d limit_delta=%f\n", global_avg, iter, size, MAX_DELTA);
      }
    }


    for (int i=1; i<MP+1; i++)
    {
      for (int j=1; j<NP+1; j++)
      {
      //printf("rank=%d, i=%d, j=%d, old=%f, new=%f\n",rank, i, j, old[i][j], new[i][j]);
      //if (fabs(new[i][j]) > 255) printf("!!!!! rank=%d, i=%d, j=%d, new=%f\n", rank, i, j, new[i][j]);
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
  //printf("finished %d n=%d, p=%d\n", rank, right, left);
  /*if (MAX_DELTA > max_delta)
  {
    printf("Finished earlier, iter=%d, max_delta=%f \n", iter, max_delta);
  }*/
  //printf("SIZEEE = %d, buf[1][1]=%f\n", (int)(sizeof(masterbuf)/sizeof(float)),buf[24][24]);
  reconstruct_time = -MPI_Wtime();
  my_Gather(masterbuf, MP, NP, buf, rank, size);
  //printf("afte my gather, rank=%d\n", rank);
  reconstruct_time += MPI_Wtime();
  barrier_time = -MPI_Wtime();
  MPI_Barrier(comm);
  barrier_time += MPI_Wtime();
  //printf("afte barrier, rank=%d\n", rank);

  local_sum = mySum(buf, MP * NP);
  MPI_Reduce(&local_sum, &global_sum, 1, MPI_FLOAT, MPI_SUM, 0, comm);

  printf("avg_time=%.10f rank=%d overall_time=%.10f total_loop=%.10f make_MP_time=%.10f choose_neighbours=%.10f make_buff=%.10f reconstruct_time=%.10f barrier_time=%.10f\n",
  avg_time(times[rank], iter), rank, start_time+MPI_Wtime(), avg_time(times[rank], iter) * iter, make_MP_time, choose_neighbours, make_buff, reconstruct_time, barrier_time);



  //printf("finished AFTER gather %d n=%d, p=%d\n", rank, right, left);

  if (rank == 0)
    {
      global_avg = global_sum / (double)(M * N);
      printf("global_avg=%.10f max_delta=%.10f iter=%d\n", global_avg, max_delta, iter);
      //printf("Writing\n");
      pgmwrite(filename_end, masterbuf, M, N);
      //printf("Finished writing\n");
    }
    /*
  for (int i=0; i<size; i++)
  {
    if (i == rank) printf("avg_time=%f rank=%d overall_time=%f total_loop=%f make_MP_time=%f choose_neighbours=%f make_buff=%f reconstruct_time=%f barrier_time=%f\n",
    avg_time(times[i], iter), rank, start_time+MPI_Wtime(), avg_time(times[i], iter) * iter, make_MP_time, choose_neighbours, make_buff, reconstruct_time, barrier_time);
  }
  */
  MPI_Barrier(comm);
  MPI_Finalize();
}


// ---------------------------
// Additional functions
// ---------------------------


void choose_MN(void *myArray, int size)
{
  /*
    Function that calculates the the optimum MP & NP based on the number of
    processes acting. The basic concept is to find the closest 2 integers
    that multiplied equal to size. It then return an array to be allocated
    to MP & NP.
  */
  int M_i, N_i;
  int MP, NP;
  int *MN = (int *)myArray;
  int done_choosing_MN = 0, cont = 1;

  if (fmod(size, sqrt(size)) == 0)
  {
    // is square of another number
    MN[0] = M / (int)sqrt(size);
    MN[1] = N / (int)sqrt(size);
    //printf("MP=%d, NP=%d\n", MN[0], MN[1]);
  }
  else
  {
    // start mapping the 2 closest intergers that multiply give P
    while (done_choosing_MN == 0)
    {
      M_i = size / sqrt(size) + cont;
      //printf("in while Mi=%d, Ni=%d, cont=%d\n", M_i, N_i, cont);
      while(size % M_i != 0)
      {
        cont++;
        M_i = size / sqrt(size) + cont;
        //printf("Mi=%d, Ni=%d, cont=%d\n", M_i, N_i, cont);
      }
      //printf("after 2nd while\n");
      N_i = size / M_i;
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
        //printf("MP=%d, NP=%d 22 \n", MN[0], MN[1]);
        done_choosing_MN = 1;
      }
      else if (fmod(M, N_i) == 0 && fmod(N, M_i) == 0)
      {
        MN[0] = M / N_i;
        MN[1] = N / M_i;
        //printf("MP=%d, NP=%d 23 \n", MN[0], MN[1]);
        done_choosing_MN = 1;
      }
    }
  }
}

void my_Scatter(void *buf, int m, int n, float masterbuf[m][n], int rank, int MP, int NP)
{
  /*
    Function to scatter the data between processes from the initial buffer based
    on the current proces rank.
  */
  int i, j;
  int x, y, chunks;
  float *lbuf = (float *)buf;

  chunks = M / MP;
  x = rank % chunks;
  y = rank / chunks;

  for (i=0; i<MP; i++)
  {
    for (j=0; j<NP; j++)
    {
      lbuf[i*NP+j] = masterbuf[x*MP+i][y*NP+j];
    }
  }
}


void my_Gather(void *masterbuf, int MP, int NP, float buf[MP][NP], int rank, int size)
{ //I'm making my own Gather, with blackjack & hookers

  int i, j, k;
  int x, y, chunks;
  float *lmbuf = (float *)masterbuf;
  chunks = M / MP;
  float localbuf[MP][NP];
  if (rank != 0)
  {
    //printf("my_gather, rank=%d, buf[2][2]=%f\n", rank, buf[24][24]);
    communicate_chunk(buf, rank, 0, MP, NP);

    //MPI_Send(&buf, MP*NP, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
    //printf("my_gather 2, rank=%d\n", rank);
  }
  else
  {
    my_Gather_process(lmbuf, MP, NP, buf, rank);
    for (k=1; k<size; k++)
    {
      //printf("my_gather, rank=%d, i=%d, doining first process, MP=%d, NP=%d\n", rank, k, MP, NP);
      communicate_chunk(localbuf, 0, k, MP, NP);
      //printf("after RECVS i=%d, localbuf[2[2]]=%f, MP=%d, NP=%d\n", k, localbuf[24][24], MP, NP);
      my_Gather_process(lmbuf, MP, NP, localbuf, k);
    }
  }
}

void my_Gather_process(float *lmbuf, int MP, int NP, float buf[MP][NP], int rank)
{ //I'm making my own Gather, with blackjack & hookers
  int i, j;
  int x, y, chunks;
  // float *lmbuf = (float *)masterbuf;
  fflush(stdout);
  //printf("in gather_procs rank=%d, MP=%d, NP=%d, M=%d, N=%d\n", rank, MP, NP, M, N);
  chunks = M / MP;

  //printf("in gather_procs rank=%d, chunks=%d\n", rank, chunks);
  x = rank % chunks;
  y = rank / chunks;
  fflush(stdout);
  //printf("after chunks \n");
  for (i=0; i<MP; i++)
  {
    for (j=0; j<NP; j++)
    {
      //printf("rank=%d, i=%d, j=%d, x=%d, y=%d, lx=%d, x=%d, y=%d, lbuf=%d, size = %f, STUFF1=%f  ",rank, i, j, x, y,(x*MP+i)*N+y*NP+j, (x*MP+i), (x*MP+i)*N, y*NP+j, buf[i][j], lmbuf[(x*MP+i)*N+y*NP+j]);
      lmbuf[(x*MP+i)*N+y*NP+j] = buf[i][j];
      //printf("STUFF 2=%f \n", lmbuf[(x*MP+i)*N+y*NP+j]);
      fflush(stdout);
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

double avg_time(double *myArray, int size)
{
  double sum = 0;
  double *x = (double *)myArray;
  for (int i=0; i<size; i++)
  {
    sum += x[i];
  }
  sum = sum / (double)size;
  return sum;
}
