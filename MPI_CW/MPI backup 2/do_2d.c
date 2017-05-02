#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#include "pgmio.h"
#include "comms.h"
#include "do_2d.h"

#define MAXITER 100000
double times[MAXITER];

int main(int argc, char** argv)
{
  int MN[2];
  int MP, NP;
  int i,j;
  int iter=0;
  int size, rank, left, right, up, down, MP_fact;
  float local_sum, global_sum, local_avg, global_avg;
  double start_time, make_MP_time, choose_neighbours, make_buff, reconstruct_time, barrier_time;
  char *ptr;

  MPI_Comm comm;
  MPI_Status status;

  // set defaults, in case inputs are not given.
  int M=768, N=768, DELTA_FREQ=100, AVG_FREQ=200;
  float MAX_DELTA=0.05;
  int val_array[10] = {1, 2, 5, 10, 25, 50, 75, 100, 150, 200};
  float m_delta[9] = {1., 0.75, 0.5, 0.25, 0.1, 0.075, 0.05, 0.025, 0.01};
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
    } else if (strcmp(argv[arg], "Af") == 0)
    {
      AVG_FREQ = val_array[(int)strtol(argv[arg+1], &ptr, 10)-1];
    } else if (strcmp(argv[arg], "MD") == 0)
    {
      MAX_DELTA = m_delta[(int)strtol(argv[arg+1], &ptr, 10)-1];
    }
  }

  //create the arrays which will be used for image processing
  float max_delta = MAX_DELTA + 1;
  float masterbuf[M][N];
  fflush(stdout);

  //generate the name of the files to be read & written to
  char filename[16], filename_end[22];
  sprintf(filename, "edge%dx%d.pgm", M, N);
  sprintf(filename_end, "edge%dx%d_%.3f.pgm", M, N, MAX_DELTA);

  //read the initial data
  printf("Reading\n");
  pgmread(filename, masterbuf, M, N);
  printf("Finished reading\n");


  // MPI STARTS FROM HERE
  comm = MPI_COMM_WORLD;

  MPI_Init(NULL,NULL);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  start_time = -MPI_Wtime(); //timing the entire MPI section

  // choosing the optimum split for the big array
  make_MP_time = -MPI_Wtime();
  choose_MN(MN, size, M, N);
  MP = MN[0];
  NP = MN[1];
  make_MP_time += MPI_Wtime();

  // create the smaller arrays
  float max_delta_thread[size];
  float buf[MP][NP];
  float edge[MP+2][NP+2];
  float old[MP+2][NP+2];
  float new[MP+2][NP+2];
  float delta[MP*NP];

  // map the neighbours based on current rank, left righ up & down
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

  // initial print detailing all the current variables. To be used in post processing
  if (rank==0) printf("init size=%d MP=%d NP=%d M=%d N=%d max_delta=%f delta_freq=%d avg_freq=%d\n", size, MP, NP, M, N, MAX_DELTA, DELTA_FREQ, AVG_FREQ);

  // Scatter data and assign proper values to corresponding arrays
  make_buff = -MPI_Wtime();
  my_Scatter(buf, M, N, masterbuf, rank, MP, NP); // scatter the value to appropriate buf according to the rank

  for (int i=0; i<MP+2; i++)
  {
    for (int j=0; j<NP+2; j++)
    {
      edge[i][j] = 255; // the edge matrix is given 255, including the halo
      old[i][j] = 255; // same for the old matrix
    }
  }
  for (i=1; i<MP+1; i++)
  {
    for (j=1; j<NP+1; j++)
    {
      edge[i][j] = buf[i-1][j-1]; // edge matrix will record the just the intial edge picture
      old[i][j] = buf[i-1][j-1]; // old initially looks like the edge, as it is the initial iteration
    }
  }
  make_buff += MPI_Wtime();

  // start the main loop
  // loop finishes when either the maximum number of iterations has been achieved
  // or the maximum change in a pixel past which it becomes redundent is computed.
  // if no max_delta, we just go to max iterations.
  while (iter < MAXITER && MAX_DELTA < max_delta)
  {
    times[iter] = -MPI_Wtime(); // record the time for each loop

    // send the halos to the corresponding neighbours
    communicate_lr(old, left, right, MP, NP);
    communicate_ud(old, up, down, MP, NP);

    // compute the new pixel value here
    for (int i=1; i<MP+1; i++)
    {
      for (int j=1; j<NP+1; j++)
      {
        new[i][j] =  (old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - edge[i][j]) * 0.25;
      }
    }

    // compute the value of the maximum delta change after a certain number of iterations
    if (iter % DELTA_FREQ == 0)
    {
      for (int i=0; i<MP; i++)
      {
        for (int j=0; j<NP; j++)
        {
          delta[i*NP+j] = fabsf(old[i+1][j+1] - new[i+1][j+1]); // for the current rank
        }
      }
      // after all ranks calculated their change, from each array, the maximum will be send to another array that hold just the maximum value
      MPI_Allreduce(delta, max_delta_thread, size, MPI_FLOAT, MPI_MAX, comm);
      max_delta = maxValue(max_delta_thread, size); //get the new maximum for the entire picture
      if (rank==0) printf("max_delta=%.10f iter=%d rank=%d size=%d limit_delta=%f\n", max_delta, iter, rank, size, MAX_DELTA); // for post processing
    }
/*
    // Calculate the average pixel value
    if (iter % AVG_FREQ == 0)
    {
      // get the local average and the local sum
      local_avg = myAverage(buf, MP * NP);
      local_sum = mySum(buf, MP * NP);
      MPI_Reduce(&local_sum, &global_sum, 1, MPI_FLOAT, MPI_SUM, 0, comm); // send the local sum to a global sum stored in MASTER
      if (rank == 0)
      {
        // the MASTER now calculates the actual global average prints it for post processing.
        global_avg = global_sum / (double)(M * N);
        printf("local_avg=%.10f iter=%d size=%d limit_delta=%f\n", global_avg, iter, size, MAX_DELTA);
      }
    }
*/
    // the new becomes the old for the next iteration & the cycle begins again
    for (int i=1; i<MP+1; i++)
    {
      for (int j=1; j<NP+1; j++)
      {
        old[i][j] = new[i][j];
        buf[i-1][j-1] = old[i][j];
      }
    }

    times[iter] += MPI_Wtime();
    iter++;
  }

  // star reconstructing the MASTERBUF mastrix & time everything
  reconstruct_time = -MPI_Wtime();
  my_Gather(masterbuf, MP, NP, buf, rank, size, M, N);
  reconstruct_time += MPI_Wtime();
  barrier_time = -MPI_Wtime();
  MPI_Barrier(comm); // wait for all processes to send their buf
  barrier_time += MPI_Wtime();

  // calculate the current value of the global sum & average
  local_sum = mySum(buf, MP * NP);
  MPI_Reduce(&local_sum, &global_sum, 1, MPI_FLOAT, MPI_SUM, 0, comm);

  // print all the times for each process for post processing
  printf("avg_time=%.10f rank=%d overall_time=%.10f total_loop=%.10f make_MP_time=%.10f choose_neighbours=%.10f make_buff=%.10f reconstruct_time=%.10f barrier_time=%.10f\n",
  avg_time(times, iter), rank, start_time+MPI_Wtime(), avg_time(times, iter) * iter, make_MP_time, choose_neighbours, make_buff, reconstruct_time, barrier_time);

  // if rank is MASTER, then print the global avg, max_delta achieved & the iteration at which the loop stopped
  // also, write the final picture back to memory
  if (rank == 0)
  {
    global_avg = global_sum / (double)(M * N);
    printf("global_avg=%.10f max_delta=%.10f iter=%d\n", global_avg, max_delta, iter);
    printf("Writing\n");
    pgmwrite(filename_end, masterbuf, M, N);
    printf("Finished writing\n");
  }

  // finalise the MPI
  MPI_Barrier(comm);
  MPI_Finalize();
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

void choose_MN(void *myArray, int size, int M, int N)
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
  }
  else
  {
    // start mapping the 2 closest intergers that multiply give P
    while (done_choosing_MN == 0)
    {
      M_i = size / sqrt(size) + cont;
      while(size % M_i != 0)
      {
        cont++;
        M_i = size / sqrt(size) + cont;
      }
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
        done_choosing_MN = 1;
      }
      else if (fmod(M, N_i) == 0 && fmod(N, M_i) == 0)
      {
        MN[0] = M / N_i;
        MN[1] = N / M_i;
        done_choosing_MN = 1;
      }
    }
  }
}

void my_Scatter(void *buf, int M, int N, float masterbuf[M][N], int rank, int MP, int NP)
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


void my_Gather(void *masterbuf, int MP, int NP, float buf[MP][NP], int rank, int size, int M, int N)
{ //I'm making my own Gather, with blackjack & hookers
  /*
  Function that does the MPI_GATHER but for the 2D case.
  It takes the the buff of each thread going in and sending to MASTER then wait
  for completion. MASTER will write its own buff to masterbuf, then start
  to receive stuff from other processes in order & writing them accordingly.
  */
  int i, j, k;
  int x, y, chunks;
  float *lmbuf = (float *)masterbuf;
  chunks = M / MP;
  float localbuf[MP][NP];

  if (rank != 0)
  { //if not the MASTER, send the buff
    communicate_chunk(buf, rank, 0, MP, NP);
  }
  else
  {
    my_Gather_process(lmbuf, MP, NP, buf, rank, M, N); // write the initial buff
    for (k=1; k<size; k++)
    { // start to write the other buffs to the main one
      communicate_chunk(localbuf, 0, k, MP, NP);
      my_Gather_process(lmbuf, MP, NP, localbuf, k, M, N);
    }
  }
}

void my_Gather_process(float *lmbuf, int MP, int NP, float buf[MP][NP], int rank, int M, int N)
{ //I'm making my own Gather, with blackjack & hookers
  /*
  Second part of gathering. This function actually writes the data to the master buff.
  Takes the buff & the masterbuf and writes the buff into the master.
  */
  int i, j;
  int x, y, chunks;
  fflush(stdout);
  chunks = M / MP;

  x = rank % chunks;
  y = rank / chunks;
  fflush(stdout);

  for (i=0; i<MP; i++)
  {
    for (j=0; j<NP; j++)
    {
      lmbuf[(x*MP+i)*N+y*NP+j] = buf[i][j];
      fflush(stdout);
    }
  }
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
