#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>


void communicate(void *oldbuf, int left, int right, int MP, int NP)
{
  int count = MP * NP;
  float *buf = (float *)oldbuf;
  int prev = left;
  int next = right;
  MPI_Status status_p;
  MPI_Request request_p;
  MPI_Status status_n;
  MPI_Request request_n;
  //MPI_Comm comm;

  printf("in comms, before send %d %d\n", prev, next);
  MPI_Ibsend(&buf[(NP+2)*MP+1], NP, MPI_FLOAT, prev, 0, MPI_COMM_WORLD, &request_p);
  printf("in comms, after first send %d %d\n", prev, next);
  MPI_Ibsend(&buf[(NP+2)+1], NP, MPI_FLOAT, next, 0, MPI_COMM_WORLD, &request_n);
  printf("in comms, after second send %d %d\n", prev, next);
  MPI_Recv(&buf[1], NP, MPI_FLOAT, prev, 0, MPI_COMM_WORLD, &status_p);
  printf("in comms, after first rec %d %d\n", prev, next);
  MPI_Recv(&buf[(NP+2)*(MP+1)+1], NP, MPI_FLOAT, next, 0, MPI_COMM_WORLD, &status_n);
  printf("in comms, after second recv%d %d\n", prev, next);

  MPI_Wait(&request_p, &status_p);
  printf("in comms, after first wait %d %d\n", prev, next);
  MPI_Wait(&request_n, &status_n);
  printf("in comms, after second wait %d %d\n", prev, next);
}
