#include <stdio.h>
#include <math.h>


#define N 729
#define reps 100
#include <omp.h>


double a[N][N], b[N][N], c[N];
int jmax[N];
int ckk;

void init1(void);
void init2(void);
void loop1(void);
void loop2(void);
void loop1_static(void);
void loop2_static(void);
void loop1_auto(void);
void loop2_auto(void);
void loop1_static_n(int n);
void loop2_static_n(int n);
void loop1_dynamic(int n);
void loop2_dynamic(int n);
void loop1_guided(int n);
void loop2_guided(int n);

void valid1(void);
void valid2(void);


int main(int argc, char *argv[]) {

  double start1,start2,end1,end2;
  double l1_t_def, l2_t_def;
  int r;


  ckk = 4;

  init1();

  #pragma omp parallel default(none) shared(ckk, a, b, c, start1, start2, end1, end2, l1_t_def, l2_t_def)
  {
    #pragma omp single
    {
      start1 = omp_get_wtime();

      for (int r=0; r<reps; r++){
        loop1();
      }

      end1  = omp_get_wtime();

      valid1();

      printf("Total time for %d reps of loop 1 = %f\n",reps, (float)(end1-start1));
      l1_t_def = end1 - start1;

    // STATIC
    init1();
      start1 = omp_get_wtime();
    }

    #pragma omp barrier
    {
    for (int r=0; r<reps; r++){
      loop1_static();
    }
    }
  #pragma omp single
  {
    end1  = omp_get_wtime();
    valid1();

  printf("Total time for %d reps of loop 1 with STATIC  in thread %d = %f\n",reps, omp_get_thread_num(), (float)(end1-start1));
  printf("Time dif for loop1 for STATIC = %f\n", (float)(l1_t_def - (end1 - start1)));

  // AUTO
  init1();
    start1 = omp_get_wtime();
  }
    #pragma omp barrier
    {
    for (int r=0; r<reps; r++){
      loop1_auto();
    }
    }
    #pragma omp single
    {
    end1  = omp_get_wtime();

    valid1();

    printf("Total time for %d reps of loop 1 with AUTO = %f\n",reps, (float)(end1-start1));
    printf("Time dif for loop1 for AUTO = %f\n", (float)(l1_t_def - (end1 - start1)));


// STATIC N

    init1();
      start1 = omp_get_wtime();
    }
  #pragma omp barrier
  {
    for (int r=0; r<reps; r++){
      loop1_static_n(ckk);
    }

  }
  #pragma omp single
  {
  end1  = omp_get_wtime();

  valid1();

  printf("Total time for %d reps of loop 1 with STATIC_%d = %f\n",reps, ckk, (float)(end1-start1));
  printf("Time dif for loop1 for STATIC_%d = %f\n", ckk, (float)(l1_t_def - (end1 - start1)));


// DYNAMIC

    init1();
      start1 = omp_get_wtime();
    }
  #pragma omp barrier
  {
    for (int r=0; r<reps; r++){
      loop1_dynamic(ckk);
    }

  }
  #pragma omp single
  {
  end1  = omp_get_wtime();

  valid1();

  printf("Total time for %d reps of loop 1 with DYNAMIC_%d = %f\n",reps, ckk, (float)(end1-start1));
  printf("Time dif for loop1 for DYNAMIC_%d = %f\n", ckk, (float)(l1_t_def - (end1 - start1)));



// GUIDED

init1();
  start1 = omp_get_wtime();
}
#pragma omp barrier
{
for (int r=0; r<reps; r++)
{
  loop1_guided(ckk);
}

}
#pragma omp single
{
end1  = omp_get_wtime();

valid1();

printf("Total time for %d reps of loop 1 with GUIDED_%d = %f\n",reps, ckk, (float)(end1-start1));
printf("Time dif for loop1 for GUIDED_%d = %f\n", ckk, (float)(l1_t_def - (end1 - start1)));

// section for loop2

init2();

start2 = omp_get_wtime();

for (int r=0; r<reps; r++){
  loop2();
}

end2  = omp_get_wtime();

valid2();

printf("Total time for %d reps of loop 2 = %f\n",reps, (float)(end2-start2));
l2_t_def = end2 - start2;


//STATIC
  init2();


start2 = omp_get_wtime();
}
#pragma omp barrier
{
for (int r=0; r<reps; r++){
  loop2_static();
}
}
#pragma omp single
{
end2  = omp_get_wtime();

valid2();

printf("Total time for %d reps of loop 2 STATIC = %f\n",reps, (float)(end2-start2));
printf("Time dif for loop2 for STATIC = %f\n", (float)(l2_t_def - (end2 - start2)));

//AUTO
init2();

start2 = omp_get_wtime();
}
#pragma omp barrier
{
for (int r=0; r<reps; r++){
  loop2_auto();
}
}
#pragma omp single
{
end2  = omp_get_wtime();

valid2();

printf("Total time for %d reps of loop 2 AUTO = %f\n",reps, (float)(end2-start2));
printf("Time dif for loop2 for AUTO = %f\n", (float)(l2_t_def - (end2 - start2)));


// STATIC N

    init2();
      start2 = omp_get_wtime();
    }
  #pragma omp barrier
  {
    for (int r=0; r<reps; r++)
    {
      loop2_static_n(ckk);
    }

  }
  #pragma omp single
  {
  end2  = omp_get_wtime();

  valid2();

  printf("Total time for %d reps of loop 2 with STATIC_%d = %f\n",reps, ckk, (float)(end2 - start2));
  printf("Time dif for loop2 for STATIC_%d = %f\n", ckk, (float)(l2_t_def - (end2 - start2)));


// DYNAMIC

init2();
  start2 = omp_get_wtime();
}
#pragma omp barrier
{
for (int r=0; r<reps; r++)
{
  loop2_dynamic(ckk);
}

}
#pragma omp single
{
end2  = omp_get_wtime();

valid2();

printf("Total time for %d reps of loop 2 with DYNAMIC_%d = %f\n",reps, ckk, (float)(end2 - start2));
printf("Time dif for loop2 for DYNAMIC_%d = %f\n", ckk, (float)(l2_t_def - (end2 - start2)));


// GUIDED
init2();
  start2 = omp_get_wtime();
}
#pragma omp barrier
{
for (int r=0; r<reps; r++)
{
  loop2_guided(ckk);
}

}
#pragma omp single
{
end2  = omp_get_wtime();

valid2();

printf("Total time for %d reps of loop 2 with GUIDED_%d = %f\n",reps, ckk, (float)(end2 - start2));
printf("Time dif for loop2 for GUIDED_%d = %f\n", ckk, (float)(l2_t_def - (end2 - start2)));
}
}
}

// init section
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
#pragma omp for schedule(static)
  for (i=0; i<N; i++){
    for (j=N-1; j>i; j--){
      a[i][j] += cos(b[i][j]);
    }
  }
}



void loop2_static(void) {
  /* loop2 using STATIC */
  int i,j,k;
  double rN2;

  rN2 = 1.0 / (double) (N*N);
#pragma omp for schedule(static)
  for (i=0; i<N; i++){
    for (j=0; j < jmax[i]; j++){
      for (k=0; k<j; k++){
	       c[i] += (k+1) * log (b[i][j]) * rN2;
      }
    }
  }

}


//AUTO
void loop1_auto(void) {
  int i,j;

#pragma omp for schedule(auto) private(i,j)
  for (i=0; i<N; i++){
    for (j=N-1; j>i; j--){
      a[i][j] += cos(b[i][j]);
    }
  }
}

void loop2_auto(void) {
  int i,j,k;
  double rN2;

  rN2 = 1.0 / (double) (N*N);

#pragma omp for schedule(static) private(i,j,k)
  for (i=0; i<N; i++){
    for (j=0; j < jmax[i]; j++){
      for (k=0; k<j; k++){
	c[i] += (k+1) * log (b[i][j]) * rN2;
      }
    }
  }
}


//STATIC N
void loop1_static_n(int n) {
  /* the static parallelisation of the first loop.
  it uses STATIC as the parallel for kind. */
  int i,j;
#pragma omp for schedule(static, n)
  for (i=0; i<N; i++){
    for (j=N-1; j>i; j--){
      a[i][j] += cos(b[i][j]);
    }
  }
}

void loop2_static_n(int n) {
  /* loop2 using STATIC */
  int i,j,k;
  double rN2;

  rN2 = 1.0 / (double) (N*N);
#pragma omp for schedule(static, n)
  for (i=0; i<N; i++){
    for (j=0; j < jmax[i]; j++){
      for (k=0; k<j; k++){
	       c[i] += (k+1) * log (b[i][j]) * rN2;
      }
    }
  }
}

//DYNAMICS

void loop1_dynamic(int n) {
  /* the dynamic parallelisation of the first loop.
  it uses DYNAMIC as the parallel for kind. */
  int i,j;
#pragma omp for schedule(dynamic, n)
  for (i=0; i<N; i++){
    for (j=N-1; j>i; j--){
      a[i][j] += cos(b[i][j]);
    }
  }
}

void loop2_dynamic(int n) {
  /* loop2 using DYNAMIC */
  int i,j,k;
  double rN2;

  rN2 = 1.0 / (double) (N*N);
#pragma omp for schedule(dynamic, n)
  for (i=0; i<N; i++){
    for (j=0; j < jmax[i]; j++){
      for (k=0; k<j; k++){
	       c[i] += (k+1) * log (b[i][j]) * rN2;
      }
    }
  }
}

//GUIDED
void loop1_guided(int n) {
  /* the guided parallelisation of the first loop.
  it uses GUIDED as the parallel for kind. */
  int i,j;
#pragma omp for schedule(guided, n)
  for (i=0; i<N; i++){
    for (j=N-1; j>i; j--){
      a[i][j] += cos(b[i][j]);
    }
  }
}

void loop2_guided(int n) {
  /* loop2 using GUIDED */
  int i,j,k;
  double rN2;

  rN2 = 1.0 / (double) (N*N);

#pragma omp for schedule(guided, n)
  for (i=0; i<N; i++){
    for (j=0; j < jmax[i]; j++){
      for (k=0; k<j; k++){
	       c[i] += (k+1) * log (b[i][j]) * rN2;
      }
    }
  }
}


//valid part
void valid1(void) {
  int i,j;
  double suma;

  suma= 0.0;
  for (i=0; i<N; i++){
    for (j=0; j<N; j++){
      suma += a[i][j];
    }
  }
}


void valid2(void) {
  int i;
  double sumc;

  sumc= 0.0;
  for (i=0; i<N; i++){
    sumc += c[i];
  }
}
