#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int main()
{
  int rank, size, tag;
  int start,stop,i;
  int l, r, ad, pass, sum;
  MPI_Comm comm;
  MPI_Status status;
  MPI_Request request;

  comm = MPI_COMM_WORLD;
  tag = 1;
  MPI_Init(NULL,NULL);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  l = rank - 1;
  r = rank + 1;

  if (r == size)
  {
    r = 0;
  }
  if (rank == 0)
  {
    l = size - 1;
  }
  sum = 0;

  pass = (rank + 1) * (rank + 1);
  for (int i =0; i<size; i++)
  {
    MPI_Issend(&pass, 1, MPI_INT, r, tag, comm, &request);
    MPI_Recv(&ad, 1, MPI_INT, l, tag, comm, &status);
    MPI_Wait(&request, &status);

    sum += ad;
    pass = ad;
  }

  printf("rank = %d, sum =%d \n", rank, sum);

  MPI_Finalize();
}
