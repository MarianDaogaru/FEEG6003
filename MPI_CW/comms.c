#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>


void communicate_lr(void *oldbuf, int left, int right, int MP, int NP)
{
  float *buf = (float *)oldbuf;
  int prev = left;
  int next = right;
  MPI_Status status_p;
  MPI_Request request_p;
  MPI_Status status_n;
  MPI_Request request_n;
  //MPI_Comm comm;

  printf("in comms, before send %d %d\n", prev, next);
  MPI_Ibsend(&buf[(NP+2)*MP+1], NP, MPI_FLOAT, next, 0, MPI_COMM_WORLD, &request_p);
  printf("in comms, after first send %d %d\n", prev, next);
  MPI_Ibsend(&buf[(NP+2)+1], NP, MPI_FLOAT, prev, 0, MPI_COMM_WORLD, &request_n);
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

void communicate_ud(void *oldbuf, int up, int down, int MP, int NP)
{
  float *buf = (float *) oldbuf;
  int prev = up;
  int next = down;
  MPI_Status status_p;
  MPI_Request request_p;
  MPI_Status status_n;
  MPI_Request request_n;

  MPI_Datatype array_up_send, array_up_recv;
  MPI_Datatype array_down_send, array_down_recv;
  MPI_Datatype test;

  /*MPI_Type_vector(MP, 1, NP+2, MPI_FLOAT, &array_up_send);
  MPI_Type_commit(&array_up_send);
  MPI_Type_vector(MP+2, 1, NP+2, MPI_FLOAT, &array_down_send);
  MPI_Type_commit(&array_down_send);
  MPI_Type_vector(MP+2, 1, NP+2, MPI_FLOAT, &array_up_send);
  MPI_Type_commit(&array_up_send);
  MPI_Type_vector(MP+2, 1, NP+2, MPI_FLOAT, &array_down_recv);
  MPI_Type_commit(&array_down_recv);
  MPI_Type_vector(MP+2, 1, NP+2, MPI_FLOAT, &array_up_recv);
  MPI_Type_commit(&array_up_recv);
  MPI_Type_vector(MP+2, 1, NP+2, MPI_FLOAT, &array_down_recv);
  MPI_Type_commit(&array_down_recv); */
  MPI_Type_vector(MP, 1, NP+2, MPI_FLOAT, &test);
  MPI_Type_commit(&test);

  printf("AAAin comms, before send %d %d\n", prev, next);
  MPI_Ibsend(&buf[(NP+2)+1], 1, test, up, 0, MPI_COMM_WORLD, &request_p);
  printf("AAAin comms, after first send %d %d\n", prev, next);
  MPI_Ibsend(&buf[(NP+2)+NP], 1, test, down, 0, MPI_COMM_WORLD, &request_n);
  printf("AAAin comms, after second send %d %d\n", prev, next);
  MPI_Recv(&buf[NP+2], 1, test, up, 0, MPI_COMM_WORLD, &status_p);
  printf("AAAin comms, after first rec %d %d\n", prev, next);
  MPI_Recv(&buf[(NP+2)+NP+1], 1, test, down, 0, MPI_COMM_WORLD, &status_n);
  printf("AAAin comms, after second recv%d %d\n", prev, next);

  MPI_Wait(&request_p, &status_p);
  printf("AAAin comms, after first wait %d %d\n", prev, next);
  MPI_Wait(&request_n, &status_n);
  printf("AAAin comms, after second wait %d %d\n", prev, next);
}

void communicate_chunk(void *buf, int rank, int recv_rank, int MP, int NP)
{
  MPI_Status status;
  MPI_Request request;
  float *buff =(float *)buf;
  if (rank!=0)
  {
    printf("in comms_chunk, SEND rank=%d, buff=%f\n", rank, buff[24*NP + 24]);
    MPI_Send(&buff[0], MP*NP, MPI_FLOAT, recv_rank, 0, MPI_COMM_WORLD);
  }
  else
  {
    printf("in comms_chunk, RECV rank=%d\n", rank);
    MPI_Recv(&buff[0], MP*NP, MPI_FLOAT, recv_rank, 0, MPI_COMM_WORLD, &status);
    printf("in comms_chunk, WAIT rank=%d, buff=%f\n", rank,buff[24*NP + 24]);
    //MPI_Wait(&request, &status);
    printf("dupa wait\n");
  }
  printf("afara din if \n");
}
