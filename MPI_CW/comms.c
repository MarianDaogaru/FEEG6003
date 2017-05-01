#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>


void communicate_lr(void *oldbuf, int left, int right, int MP, int NP)
{
  /*
    Function that makes the communication to the left and right neighbours.
    It takes the buf, the left & right neightbours. Then sends the appropriate
    values to them. In this case, the 2nd column to left & the second to last column to right.
    Then it receives from the the halo values.
  */
  float *buf = (float *)oldbuf;
  int prev = left;
  int next = right;
  MPI_Status status_p;
  MPI_Request request_p;
  MPI_Status status_n;
  MPI_Request request_n;

  MPI_Issend(&buf[(NP+2)*MP+1], NP, MPI_FLOAT, next, 0, MPI_COMM_WORLD, &request_p);
  MPI_Issend(&buf[(NP+2)+1], NP, MPI_FLOAT, prev, 0, MPI_COMM_WORLD, &request_n);
  MPI_Recv(&buf[1], NP, MPI_FLOAT, prev, 0, MPI_COMM_WORLD, &status_p);
  MPI_Recv(&buf[(NP+2)*(MP+1)+1], NP, MPI_FLOAT, next, 0, MPI_COMM_WORLD, &status_n);

  MPI_Wait(&request_p, &status_p);
  MPI_Wait(&request_n, &status_n);
}

void communicate_ud(void *oldbuf, int up, int down, int MP, int NP)
{
  /*
    Function that sends & receives the halos from the neighbours positioned
    above and below the current process in the topology.
    This function uses the Type_vector MPI method to actually map the values
    needed to be sent.
  */
  float *buf = (float *) oldbuf;
  int prev = up;
  int next = down;
  MPI_Status status_p;
  MPI_Request request_p;
  MPI_Status status_n;
  MPI_Request request_n;

  MPI_Datatype new_array;

  MPI_Type_vector(MP, 1, NP+2, MPI_FLOAT, &new_array);
  MPI_Type_commit(&new_array);

  MPI_Issend(&buf[(NP+2)+1], 1, new_array, up, 0, MPI_COMM_WORLD, &request_p);
  MPI_Issend(&buf[(NP+2)+NP], 1, new_array, down, 0, MPI_COMM_WORLD, &request_n);
  MPI_Recv(&buf[NP+2], 1, new_array, up, 0, MPI_COMM_WORLD, &status_p);
  MPI_Recv(&buf[(NP+2)+NP+1], 1, new_array, down, 0, MPI_COMM_WORLD, &status_n);

  MPI_Wait(&request_p, &status_p);
  MPI_Wait(&request_n, &status_n);
}

void communicate_chunk(void *buf, int rank, int recv_rank, int MP, int NP)
{
  /*
    Function that send the entire buf to rank 0. MASTER then in turns receives
    the entire buf.
  */
  MPI_Status status;
  MPI_Request request;
  float *buff =(float *)buf;
  if (rank!=0)
  {
    MPI_Send(&buff[0], MP*NP, MPI_FLOAT, recv_rank, 0, MPI_COMM_WORLD);
  }
  else
  {
    MPI_Recv(&buff[0], MP*NP, MPI_FLOAT, recv_rank, 0, MPI_COMM_WORLD, &status);
  }
}
