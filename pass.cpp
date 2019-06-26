#include "funcs.h"
#include <iostream>

using namespace std;

int RowPass(double *A, int m, int l, int s, FILE *fin, int (p)(double &a, FILE *fin))
{
    int row    = m * m * l + s * m;
    int col_m  = m * m;
    int col_sm = s * m;

    for (int i = 0; i < l; i++)
    {
        for (int k = 0; k < m; k++)
        {
            for (int j = 0; j < l; j++)
            {
                for (int r = 0; r < m; r++)
                {
                    if (p(A[i * row + j * col_m + (k * m + r)], fin) == -1)
                        return -1;
                }
            }


            ///Остаточные столбцы
            for (int r = 0; r < s; r++)
            {
                if (p(A[i * row + l * col_m + (k * s + r)], fin) == -1)
                    return -1;
            }

            if (fin == NULL)
                printf("\n");
        }
    }
    ///Остаточные строки
    for (int k = 0; k < s; k++)
    {
        for (int j = 0; j < l; j++)
        {
            for (int r = 0; r < m; r++)
            {
                if (p(A[l * row + j * col_sm + (k * m + r)], fin) == -1)
                    return -1;
            }
        }

        ///Остаточные столбцы
        for (int r = 0; r < s; r++)
        {
            if (p(A[l * row + l * col_sm + (k * s + r)], fin) == -1)
                return -1;
        }

        if (fin == NULL){
            printf("\n");
        }
    }
    return 1;
}

void RowTransposition(double *A, /*номер потока*/ int thread, int m,/*шаг алгоритма*/ int k, int s, int l, int t, int max_i,/*число потоков*/ int p)
{
    int row    = m * m * l + s * m;
    int col_m  = m * m;

    double tmp = 0.;
    for (int j = thread; j < l; j += p)
    {
        for (int i1 = 0; i1 < m; i1++)
        {
            for (int j1 = 0; j1 < m; j1++)
            {
                tmp = A[max_i * row + j * col_m + (i1 * m + j1)];
                A[max_i * row + j * col_m + (i1 * m + j1)] = A[k * row + j * col_m + (i1 * m +j1)];
                A[k * row + j * col_m +(i1 * m +j1)] = tmp;
            }
        }
    }

    if (thread == t)
    {
        for (int i1 = 0; i1 < m; i1++)
        {
            for (int j1 = 0; j1 < s; j1++)
            {
                tmp = A[max_i * row + l * col_m + (i1 * s + j1)];
                A[max_i * row + l * col_m + (i1 * s + j1)] = A[k * row + l * col_m + (i1 * s + j1)];
                A[k * row + l * col_m + (i1 * s + j1)] = tmp;
            }
        }
    }
}

void ColumnTransposition(double *A, int thread, int m, int k, int s, int l, int t, int max_j, int p)
{
    int row    = m * m * l + s * m;
    int col_m  = m * m;
    int col_sm = s * m;

    double tmp = 0.;
    for (int i = thread; i < l; i += p)
    {
        for (int i1 = 0; i1 < m; i1++)
        {
            for (int j1 = 0; j1 < m; j1++)
            {
                tmp = A[i * row + max_j * col_m + (i1 * m +j1)];
                A[i * row + max_j * col_m + (i1 * m +j1)] = A[i * row + k * col_m +(i1 * m +j1)];
                A[i * row + k * col_m +(i1 * m +j1)] = tmp;
            }
        }
    }
    if (thread == t)
    {
        for (int i1 = 0; i1 < s; i1++)
        {
            for (int j1 = 0; j1 < m; j1++)
            {
                tmp = A[l * row + max_j * col_sm + (i1 * m + j1)];
                A[l * row + max_j * col_sm + (i1 * m + j1)] = A[l * row + k * col_sm + (i1 * m + j1)];
                A[l * row + k * col_sm + (i1 * m + j1)] = tmp;
           }
       }
    }
}

void DivideMainRow(double *A, double **E, double **U, int w, int m, int s, int l, int k, int t, args_t *ptr)
{
    int row    = m * m * l + s * m;
    int col_m  = m * m;

    for (int j = w + ptr->thread; j < l; j += ptr->p)
    {
        multiplication(E[ptr->thread], A + (k * row + j * col_m), U[ptr->thread], m, m, m);
        for (int i1 = 0; i1 < m; i1++)
            for (int j1 = 0; j1 < m; j1++)
                A[k * row + j * col_m + (i1 * m + j1)] = U[ptr->thread][i1 * m + j1];

    }
    if (ptr->thread == t)
    {
        multiplication(E[ptr->thread], A + (k * row + l * col_m), U[ptr->thread], m, m, s);

        for (int i1 = 0; i1 < m; i1++)
            for (int j1 = 0; j1 < s; j1++)
                A[k * row + l * col_m + (i1 * s + j1)] = U[ptr->thread][i1 * s + j1];
    }
}

void SubtractBloks(double *B, double *A, double **U, int w, int x, int m, int s, int l, int k, int i, args_t *ptr, int col)
{
    int row    = m * m * l + s * m;
    int col_m = m*m;

    for (int j = w; j < l; j++)
    {
        multiplication(A + (i * row + k * col), B + (k * row + j * col_m), U[ptr->thread], x, m, m);
        subtract(B + (i * row + j * col), U[ptr->thread], x, m);
    }
    multiplication(A + (i * row + k * col), B + (k * row + l * col_m), U[ptr->thread], x, m, s);
    subtract(B + (i * row + l * col), U[ptr->thread], x, s);
}

void SubtractMainRow(double *B, double *A, double **U, int w, int m, int s, int l, int k, int t, args_t *ptr)
{
    for (int i = ptr->thread; i < k; i += ptr->p)
    {
        SubtractBloks(B, A, U, w, m, m, s, l, k, i, ptr, m*m);
    }

    for (int i = k+1 + ptr->thread; i < l; i += ptr->p){
        SubtractBloks(B, A, U, w, m, m, s, l, k, i, ptr, m*m);
    }

    if (ptr->thread == t){
        SubtractBloks(B, A, U, w, s, m, s, l, k, l, ptr, s*m);
    }
}

void CopyRow(double *A, double *B, int m, int s, int l, int copy_index_A, int copy_index_B)
{
    int row    = m * m * l + s * m;
    int col_m  = m * m;

    for (int j = 0; j < l; j++)
        for (int i1 = 0; i1 < m; i1++)
            for (int j1 = 0; j1 < m; j1++)
                A[copy_index_A * row + j * col_m + (i1 * m + j1)] = B[copy_index_B * row + j * col_m + (i1 *m + j1)];

    for (int i1 = 0; i1 < m; i1++)
        for (int j1 = 0; j1 < s; j1++)
            A[copy_index_A * row + l * col_m + (i1 * s + j1)] = B[copy_index_B * row + l * col_m + (i1 * s + j1)];
}
