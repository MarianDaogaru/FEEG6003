float** make_2d_dyn(int rows, int cols);
void choose_MN(void *myArray, int size, int M, int N);
void my_Scatter(void *buf, int M, int N, float masterbuf[M][N], int rank, int MP, int NP);
void my_Gather(void *masterbuf, int MP, int NP, float buf[MP][NP], int rank, int size, int M);
void my_Gather_process(float *lmbuf, int MP, int NP, float buf[MP][NP], int rank, int M);
double mySum(void *myArray, int size);
double myAverage(void *myArray, int size);
float maxValue(float *myArray, int size);
double avg_time(double *myArray, int size);
