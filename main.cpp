#include <iostream>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "funcs.h"
#include "get_time.h"

using namespace std;
pthread_barrier_t barr;
int count;

int main(int argc, char* argv[])
{
    pthread_t tid;
    int n = 0; int m = 0;
    int l = 0; int s = 0;
    int p = 0;
    bool isWrongPrm = false;
    int inputMode   = -1;

    calc_param(argv, n, m, l, s, p, isWrongPrm);

    if (isWrongPrm == true)
    {
        printf("Parmeters of command line are not correct, finish the program\n");
        return 0;
    }

    pthread_barrier_init(&barr, NULL, p);

/******************************************************************/

    double *A = new double [n * n];
    double *B = new double [n * n];

    double **E     = new double* [p];
    double **U     = new double* [p];
    int **index_nb = new int*    [p];

    for (int i = 0; i < p; i++)
    {
        U[i]        = new double [m * m];
        E[i]        = new double [m * m];
        index_nb[i] = new int [m];
    }

    block_filling(argc, argv, A, m, l, s, inputMode);
    B = filling(m, s, l, B, single);

/******************************************************************/

    args_t *args = new args_t [p];
    if(!args)
    {
        printf("Not enought memory\n");
        return -1;
    }

    maximum *max = new maximum [p];
    int *index   = new int [l];

    for (int i = 0; i < l; i++)
        index[i] = i;

    for (int i = 0; i < p; i++)
    {

        args[i].p = p;
        args[i].m = m;
        args[i].l = l;
        args[i].s = s;
        args[i].A = A;
        args[i].B = B;
        args[i].U = U;
        args[i].E = E;

        args[i].thread   = i;
        args[i].max      = max;
        args[i].index    = index;
        args[i].index_nb = index_nb[i];
    }

    double time = get_full_time();

    for (int k = 1; k < p; k++)
    {
        if (pthread_create(&tid, 0, &blockJordanInversion, args + k))
        {
            printf ("Error %d\n", k);
            return 100;
        }
    }
    blockJordanInversion((void*) args);

    time = get_full_time() - time;
    printf("Time is %f\n", time);


    block_print(B, n, m, l, s);
    block_filling(argc, argv, A, m, l, s, inputMode);

/******************************************************************/

    for (int k = 0; k < p; k++)
        printf("Thread %d time is %f\n", k, args[k].pr_time);

    for (int i = 0; i < p; i++)
    {
        args[i].A = A;
        args[i].max[i].max = -1.;

    }

    for (int k = 1; k < p; k++)
    {
        if (pthread_create(&tid, 0, &block_discrepancy, args + k))
        {
            printf ("Error with thread %d\n", k);
            return 200;
        }
    }
    block_discrepancy((void*) args);

/******************************************************************/

    for (int i = 0; i < p; i++)
    {
        delete [] U[i];
        delete [] E[i];
        delete [] index_nb[i];
    }
    /*delete [] U;
    delete [] E;*/
    delete [] index_nb;
    delete [] index;
    delete [] max;
    delete [] args;

    return 0;
}
