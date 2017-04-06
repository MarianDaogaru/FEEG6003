#include <stdio.h>
#include <math.h>


#define N 729
#define reps 100
#include <omp.h>

double a[N][N], b[N][N], c[N];
int jmax[N];


void init1(void);
void init2(void);
void loop1(void);
void loop2(void);
void loop1_static(void);
void loop2_static(void);
void valid1(void);
void valid2(void);


int main(int argc, char *argv[]) {

  double start1,start2,end1,end2;
  double l1_t_def, l2_t_def;
  int r;

  init1();

  #pragma opm parallel default(none) private(start1, start2, end1, end2) shared(a, b)
  {
  start1 = omp_get_wtime();

  for (r=0; r<reps; r++){
    loop1();
  }

  end1  = omp_get_wtime();

  valid1();

  printf("Total time for %d reps of loop 1 = %f\n",reps, (float)(end1-start1));
  l1_t_def = end1 - start1;
// STATIC
init1();
  start1 = omp_get_wtime();

  for (r=0; r<reps; r++){
    loop1_static();
  }

  end1  = omp_get_wtime();

  valid1();

  printf("Total time for %d reps of loop 1 with STATIC = %f\n",reps, (float)(end1-start1));
  printf("Time dif for loop1 for STATIC = %f\n", (float)(l1_t_def - (end1 - start1)));


// section for loop2

init2();

start2 = omp_get_wtime();

for (r=0; r<reps; r++){
  loop2();
}

end2  = omp_get_wtime();

valid2();

printf("Total time for %d reps of loop 2 = %f\n",reps, (float)(end2-start2));
l2_t_def = end2 - start2;
//STATIC
init2();

start2 = omp_get_wtime();

for (r=0; r<reps; r++){
  loop2_static();
}

end2  = omp_get_wtime();

valid2();

printf("Total time for %d reps of loop 2 STATIC = %f\n",reps, (float)(end2-start2));
printf("Time dif for loop2 for STATIC = %f\n", (float)(l2_t_def - (end2 - start2)));
}
}

void init1(void){
  int i,j;

  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      a[i][j] = 0.0;
      b[i][j] = 3.142*(i+j);
    }
  }

}

void init2(void){
  int i,j, expr;

  for (i=0; i<N; i++){
    expr =  i%( 3*(i/30) + 1);
    if ( expr == 0) {
      jmax[i] = N;
    }
    else {
      jmax[i] = 1;
    }
    c[i] = 0.0;
  }

  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      b[i][j] = (double) (i*j+1) / (double) (N*N);
    }
  }

}

void loop1(void) {
  int i,j;

  for (i=0; i<N; i++){
    for (j=N-1; j>i; j--){
      a[i][j] += cos(b[i][j]);
    }
  }

}



void loop2(void) {
  int i,j,k;
  double rN2;

  rN2 = 1.0 / (double) (N*N);

  for (i=0; i<N; i++){
    for (j=0; j < jmax[i]; j++){
      for (k=0; k<j; k++){
	c[i] += (k+1) * log (b[i][j]) * rN2;
      }
    }
  }

}

void loop1_static(void) {
  /* the static parallelisation of the first loop.
  it uses STATIC as the parallel for kind. */
  int i,j;

#pragma omp for schedule(static) private(i,j)
  for (i=0; i<N; i++){
    for (j=N-1; j>i; j--){
      a[i][j] += cos(b[i][j]);
    }
  }
}



void loop2_static(void) {
  int i,j,k;
  double rN2;

  rN2 = 1.0 / (double) (N*N);
#pragma omp fir schedule(static) private(i, j, k)
  for (i=0; i<N; i++){
    for (j=0; j < jmax[i]; j++){
      for (k=0; k<j; k++){
	       c[i] += (k+1) * log (b[i][j]) * rN2;
      }
    }
  }

}

void valid1(void) {
  int i,j;
  double suma;

  suma= 0.0;
  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      suma += a[i][j];
    }
  }
  printf("Loop 1 check: Sum of a is %lf\n", suma);

}


void valid2(void) {
  int i;
  double sumc;

  sumc= 0.0;
  for (i=0; i<N; i++){
    sumc += c[i];
  }
  printf("Loop 2 check: Sum of c is %f\n", sumc);
}
