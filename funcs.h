#include <stdlib.h>
#include <stdio.h>
#include <math.h>

static double gilbert(int i, int j, int i1, int j1, int m)
{
    return 1. / (1 + i1 + (i * m) + j1 + (j * m));
}

static double single(int i, int j, int i1, int j1, int m)
{
    (void) m;
    return ((i == j)&&(i1 == j1)) ? 1. : 0.;
}

static double difference(int i, int j, int i1, int j1, int m)
{
    return fabs((i1 + (i * m)) - (j1 + (j * m)));
}

static int file_read(double &a, FILE *fin)
{
    if (!fscanf(fin, "%lf", &a))
    {
        printf("File error\n");
        fclose(fin);
        return -1;
    }
    return 1;
}

static int print(double &a, FILE *fin)
{
    fin = NULL;
    printf("%7.8f ", a);
    return 1;
}

class maximum
{
public:
    int    max_i;
    int    max_j;
    double min;
    double max;
    int    error;
    bool   isHere;
    maximum() : max_i(0), max_j(0), min(0), max(0), error(0), isHere(true) {}
};

class args_t
{
public:
    int thread;
    int p;
    int m;
    int l;
    int s;
    double *A;
    double *B;
    double **U;
    double **E;
    maximum *max;
    int *index;
    int *index_nb;
    double pr_time;
    args_t() : thread(0), p(0), m(0), l(0), s(0), A(NULL), B(NULL), U(NULL), E(NULL), pr_time(0) {}
};

int      block_filling        (int argc, char* argv[], double *A, int m, int l, int s, int &inputMode);
void     block_print          (double *A, int n, int m, int l, int s);
double **memory               (double **A, int m, int l, int s);
void     calc_param           (char* argv[], int &n, int &m, int &l, int &s, int &p, bool &isWrongPrm);
void    *block_discrepancy    (void *prt_args);
void     clean_memory         (double **A, int l);
void    *blockJordanInversion (void *ptr_arg);
void     multiplication       (double *A, double *B, double *C, int m, int n, int r);
int      restore_A            (char *argv[], double **&A, int m, int l, int s, int &inputMode);
void     block_multiplication (double **A, double **B, double *U, int m, int l, int s);
double  *filling(int m, int s, int l, double *A, double (*p)(int i, int j, int i1, int j1, int m));
int      Jordan_Inversion     (double *A, double *B, int n, int *index);
void     subtract             (double *A, double *B, int n, int m);
double  *add                  (double *A, double *B, int n, int m);
int      RowPass              (double *A, int m, int l, int s, FILE *fin, int (p)(double &a, FILE *fin));
void     part_add             (double *A, double *B, int n, int s, int m);
void     RowTransposition     (double *A, int num, int m, int k, int s, int l, int t, int max_i, int p);
void     ColumnTransposition  (double *A, int num, int m, int k, int s, int l, int t, int max_j, int p);
void     DivideMainRow        (double *A, double **E, double **U, int w, int m, int s, int l, int k, int t, args_t *ptr);
void     SubtractMainRow      (double *B, double *A, double **U, int w, int m, int s, int l, int k, int t, args_t *ptr);
void     CopyRow              (double *A, double *B, int m, int s, int l, int copy_index_A, int copy_index_B);

extern pthread_barrier_t barr;
