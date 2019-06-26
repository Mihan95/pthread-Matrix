#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "funcs.h"
#include "get_time.h"

#define MAX_OUTPUT_SIZE 10

using namespace std;

void clean_memory(double **A, int l)
{
    #define A(i,j) A[(i) * (l + 1) + (j)]

    for (int i = 0; i < l; i++)
    {
        delete [] A(i,l);
        delete [] A(l,i);
        for (int j = 0; j < l; j++)
            delete [] A(i,j);
    }
    delete [] A(l,l);
    #undef A

    delete [] A;
}


double *filling(int m, int s, int l, double *A, double (*p)(int i, int j, int i1, int j1, int m))
{
    int i  = 0;
    int j  = 0;
    int i1 = 0;
    int j1 = 0;

    int row    = m * m * l + s * m;
    int col_m  = m * m;
    int col_sm = s * m;

    for (i = 0; i < l; i++)
    {
        for (j = 0; j < l; j++){
            for (i1 = 0; i1 < m; i1++)
                for (j1 = 0; j1 < m; j1++)
                    A[i * row + j * col_m + (i1 * m + j1)] = (*p)(i, j, i1, j1, m);
        }


            for (i1 = 0; i1 < m; i1++)
                for (j1 = 0; j1 < s; j1++)
                    A[i * row + l * col_m + (i1 * s + j1)] = (*p)(i, j, i1, j1, m);
    }

    for (j = 0; j < l; j++){
        for (i1 = 0; i1 < s; i1++)
            for (j1 = 0; j1 < m; j1++)
                A[l * row + j * col_sm + (i1 * m + j1)] = (*p)(l, j, i1, j1, m);
    }

        for ( i1 = 0; i1 < s; i1++)
            for ( j1 = 0; j1 < s; j1++)
                A[l * row + l * col_sm + (i1 * s + j1)] = (*p)(l, l, i1, j1, m);

    return A;
}

int block_filling(int argc, char *argv[], double *A, int m, int l, int s, int &inputMode)  ///Заполнение матрицы
{
    if (argc == 5)
    {
        FILE *fin = fopen(argv[4],"r");
        if (!fin)
        {
            cout << "File error" << endl;

            exit(0);
        }

        if (RowPass(A, m, l, s, fin, file_read) == -1)
            exit(0);

        fclose(fin);
    }
    if (argc == 4)
    {
      if (inputMode == -1)
      {
        cout << "Select way to assignment the matrix: \n"
        "1. Gilbert matrix\n"
        "2. Single matrix\n"
        "3. Matrix of difference\n";
        /*cin >> inputMode;
        cout << endl;*/
        inputMode = 3;
    }
        switch (inputMode)
        {
        case (1) :
            A = filling(m, s, l, A, gilbert);
            break;
        case(2) :
            A = filling(m, s, l, A, single);
            break;
        case (3) :
            A = filling(m, s, l, A, difference);
            break;
        }
    }
    return 1;
   #undef A
}

void block_print(double *A, int n, int m, int l, int s)
{
    if (n < MAX_OUTPUT_SIZE)
    {
        RowPass(A, m, l, s, NULL, print);
    }
    else
    {
        if ((n >= MAX_OUTPUT_SIZE) && (m >= MAX_OUTPUT_SIZE))
        {
            for (int k = 0; k < MAX_OUTPUT_SIZE; k++)
            {
                for (int r = 0; r < MAX_OUTPUT_SIZE; r++)
                    printf("%7.8f ", A[k * m + r]);

                cout << endl;
            }
        }
        else
        {
            if ((n >= MAX_OUTPUT_SIZE) && (m < MAX_OUTPUT_SIZE))
            {
                int p = MAX_OUTPUT_SIZE / m;       ///Сколько m влезает в MAX_OUTPUT_SIZE
                int q = MAX_OUTPUT_SIZE % m;       ///Сколько еще надо взять из соседних блоков

                RowPass(A, m, p, q, NULL, print);
            }
        }
    }
}

void FindMaxRow(double *A, double *B, double *U, double *E, int m, int w, int s, int l, int i, args_t *ptr, int col)
{
    #define E(i,j) E[(i) * m + (j)]

    int row    = m * m * l + s * m;
    int col_m  = m * m;
    int col_sm = s * m;

    for (int j = 0; j < l; j++)
    {
        for (int k = 0; k < l; k++)
        {
            multiplication(A + (i * row + k * col), B + (k * row + j * col_m), U, w, m, m);
            add(E, U, w, m);
        }
        multiplication(A + (i * row + l * col), B + (l * row + j * col_sm), U, w, s, m);
        add(E, U, w, m);
    }

        for (int k = 0; k < l; k++)
        {
            multiplication(A + (i * row + k * col), B + (k * row + l * col_m), U, w, m, s);
            part_add(E, U, w, s, m);
        }
        multiplication(A + (i * row + l * col), B + (l * row + l * col_sm), U, w, s, s);
        part_add(E, U, w, s, m);

    for (int j = 0; j < m; j++)
    {
        E(j,0) = fabs(E(j,0));
        for(int k = 1; k < m; k++)
             E(j,0) += fabs(E(j,k));

        if (E(j,0) > ptr->max[ptr->thread].max)
            ptr->max[ptr->thread].max = E(j,0);
    }

    for(int q = 0; q < m * m; q++)
    {
        U[q] = 0;
        E[q] = 0;
    }

    #undef E
}

void *block_discrepancy(void *ptr_arg)
{
    args_t *ptr = (args_t*) ptr_arg;

    double  *A = ptr->A;
    double  *B = ptr->B;
    double  *U = ptr->U[ptr->thread];
    double  *E = ptr->E[ptr->thread];

    int m = ptr->m;
    int l = ptr->l;
    int s = ptr->s;

    double norm = -1.;

    for (int i = 0; i < m * m; i++)
    {
        U[i] = 0.;
        E[i] = 0.;
    }

    for (int i = ptr->thread; i < l; i += ptr->p)
        FindMaxRow(A, B, U, E, m, m, s, l, i, ptr, m*m);

    if (ptr->thread == 0)
        FindMaxRow(A, B, U, E, m, s, s, l, l, ptr, m*s);

    pthread_barrier_wait(&barr);

    if (ptr->thread == 0)
    {
        for (int i = 0; i < ptr->p; i++)
        {
            if (ptr->max[i].max > norm)
                norm = ptr->max[i].max;
        }
        norm = fabs(norm - 1);
        printf("Discrepancy is %e\n\n", norm);
    }
    pthread_barrier_wait(&barr);

    return NULL;
}

void pFindMax(double *A, double **U, double **E, int m, int s, int l, maximum *max, int thread,
              int p, int *index_nb, int k)
{
    int row    = m * m * l + s * m;
    int col_m  = m * m;

    double sum  = 0.;
    double norm = 0.;
    bool isAllNormsZero = true;

    ///Ищем главный элемент
    for (int i = k+thread; i < l; i += p)
    {
        for (int j = k; j < l; j++)
        {
            for (int i1 = 0; i1 < m; i1++)
            {
                for (int j1 = 0; j1 < m; j1++)
                {
                    U[thread][i1 * m + j1] = A[i * row + j * col_m + (i1 * m +j1)];
                    E[thread][i1 * m + j1] = (i1 == j1) ? 1. : 0.;
                }
            }

            if (Jordan_Inversion(U[thread], E[thread], m, index_nb) == -1)
                continue;

            ///Считаем норму
            norm = 0.;
            for (int j2 = 0; j2 < m; j2++)
                norm += fabs(E[thread][j2]);

            for (int i2 = 0; i2 < m; i2++)
            {
                sum = 0.;
                for (int j2 = 0; j2 < m; j2++)
                    sum += fabs(E[thread][i2*m+j2]);

                norm = (sum > norm) ? sum : norm;
            }

            if ((isAllNormsZero == true) && (norm > 0))
            {
                max[thread].min   = norm;
                max[thread].max_i = i;
                max[thread].max_j = j;
                isAllNormsZero = false;
            }
            else if (isAllNormsZero == false)
                 {
                    if ((norm < max[thread].min) && (norm > 0)){
                        max[thread].min   = norm;
                        max[thread].max_i = i;
                        max[thread].max_j = j;
                    }
                 }
        }
    }
}

int blockFindMax(double *A, double *B, double **U, double **E, int m, int l, int s,/*шаг алгоритма*/ int k, maximum *max,
                /*номер потока*/int thread, int p, int *index, int *index_nb)
{
    int row    = m * m * l + s * m;
    int col_m  = m * m;

    ///U - временная матрица, E - единичная
    int j = 0; int t = (l + 1) % p;
    max[thread].min   = 0.;
    max[thread].max_i = k + thread;
    max[thread].max_j = k;

    if (k + thread < l)
    {
        max[thread].isHere = true;
        pFindMax(A, U, E, m, s, l, max, thread, p, index_nb, k);
    }
    else max[thread].isHere = false;



    pthread_barrier_wait(&barr);

    int max_i  = max[0].max_i;
    int max_j  = max[0].max_j;
    double min = max[0].min;

    for (int i = 1; i < p; i++)
    {
        if (max[i].min < min)
        {
            if (max[i].isHere == true)
            {
                min   = max[i].min;
                max_i = max[i].max_i;
                max_j = max[i].max_j;
            }
            else continue;

        }
        else if (max[i].min == min)
             {
                if (max[i].max_i * m + max[i].max_j > max_i * m + max_j)
                {
                    if (max[i].isHere == true)
                    {
                        max_i = max[i].max_i;
                        max_j = max[i].max_j;
                    }
                    else continue;
                }
             }
    }

    pthread_barrier_wait(&barr);

    RowTransposition(A, thread, m, k, s, l, t, max_i, p);
    RowTransposition(B, thread, m, k, s, l, t, max_i, p);

    pthread_barrier_wait(&barr);

    ColumnTransposition(A, thread, m, k, s, l, t, max_j, p);

    pthread_barrier_wait(&barr);

    for (int i1 = 0; i1 < m; i1++)
    {
        for (int j1 = 0; j1 < m; j1++)
        {
            U[thread][i1 * m + j1] = A[k * row + k * col_m + (i1 * m + j1)];
            E[thread][i1 * m + j1] = (i1 == j1) ? 1. : 0.;
        }
    }

    if (thread == 0)
    {
        j            = index[k];
        index[k]     = index[max_j];
        index[max_j] = j;
    }

    if (Jordan_Inversion(U[thread], E[thread], m, index_nb) == -1){
        if (thread == 0)
            cout << "Jordan method is not possible" << endl;

        max[thread].error = -1;
    }

    pthread_barrier_wait(&barr);

    return 1;
}


void *blockJordanInversion(void *ptr_arg)
{
    args_t *ptr  = (args_t*) ptr_arg;
    ptr->pr_time = get_time();

    int m = ptr->m;
    int l = ptr->l;
    int s = ptr->s;
    int t = (l + 1) % ptr->p;    ///Номер потока для остатков

    int row    = m * m * l + s * m;
    int col_m  = m * m;
    int col_sm = s * m;

    double *A  = ptr->A;
    double *B  = ptr->B;
    double **E = ptr->E;
    double **U = ptr->U;

    for (int k = 0; k < l; k++)
    {


        blockFindMax(A, B, U, E, m, l, s, k, ptr->max, ptr->thread, ptr->p, ptr->index, ptr->index_nb);

        if (ptr->max[ptr->thread].error == -1)
            exit(0);

        for (int i1 = 0; i1 < m; i1++)
        {
            for (int j1 = 0; j1 < m; j1++)
            {
                E[ptr->thread][i1 * m + j1] = (i1 == j1) ? 1 : 0;
                U[ptr->thread][i1 * m + j1] = A[k * row + k * col_m + (i1 * m + j1)];
            }
        }

        Jordan_Inversion(U[ptr->thread], E[ptr->thread], m, ptr->index_nb);

        DivideMainRow(A, E, U, k + 1, m, s, l, k, t, ptr);
        DivideMainRow(B, E, U, 0    , m, s, l, k, t, ptr);

        pthread_barrier_wait(&barr);

        SubtractMainRow(A, A, U, k + 1, m, s, l, k, t, ptr);
        SubtractMainRow(B, A, U, 0    , m, s, l, k, t, ptr);

        pthread_barrier_wait(&barr);
    }

    //Работаем с последним(нижним правым) блоком
    for (int i1 = 0; i1 < s; i1++)
    {
        for (int j1 = 0; j1 < s; j1++)
        {
            E[ptr->thread][i1 * s + j1] = (i1 == j1) ? 1 : 0;
            U[ptr->thread][i1 * s + j1] = A[l * row + l * col_sm + (i1 * s + j1)];
        }
    }

    if (Jordan_Inversion(U[ptr->thread], E[ptr->thread], s, ptr->index_nb) == -1){
        if (ptr->thread == 0)
            cout << "Jordan method is not possible" << endl;

        pthread_barrier_wait(&barr);
        exit(0);
    }

    pthread_barrier_wait(&barr);

    for (int j = ptr->thread; j < l; j += ptr->p)
    {
       multiplication(E[ptr->thread], B + (l * row + j * col_sm), U[ptr->thread], s, s, m);

       for (int i1 = 0; i1 < s; i1++)
           for (int j1 = 0; j1 < m; j1++)
               B[l * row + j * col_sm + (i1 * m + j1)] = U[ptr->thread][i1 * m + j1];

    }
    /*block_print(B, 5, m, l, s);
    printf("\n");*/
    if (ptr->thread == t)
    {
        multiplication(E[ptr->thread], B + (l * row + l * col_sm), U[ptr->thread], s, s, s);

        for (int i1 = 0; i1 < s; i1++)
           for (int j1 = 0; j1 < s; j1++)
               B[l * row + l * col_sm + (i1 * s + j1)] = U[ptr->thread][i1 * s + j1];

    }

    pthread_barrier_wait(&barr);

    for (int i = ptr->thread; i < l; i += ptr->p)
    {
        for (int j = 0; j < l; j++)
        {
            multiplication(A + (i * row + l * col_m), B + (l * row + j * col_sm), U[ptr->thread], m, s, m);
            subtract(B + (i * row + j * col_m), U[ptr->thread], m, m);
        }
        multiplication(A + (i * row + l * col_m), B + (l * row + l * col_sm), U[ptr->thread], m, s, s);
        subtract(B + (i * row + l * col_m), U[ptr->thread], m, s);
    }

    pthread_barrier_wait(&barr);

    ///Переставляем строки в самом верном порядке
    if (ptr->thread == 0)
    {
        for (int i = 0; i < l; i++)
            CopyRow(A, B, m, s, l, ptr->index[i], i);

        for (int i = 0; i < l; i++)
            CopyRow(B, A, m, s, l, i, i);
    }

   ptr->pr_time = get_time() - ptr->pr_time;
   pthread_barrier_wait(&barr);

   return NULL;
}

