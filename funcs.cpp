#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#define EPS 2e-15

using namespace std;

void calc_param(char* argv[], int &n, int &m, int &l, int &s, int &p, bool &isWrongParam)
{
    n = atoi(argv[1]);  ///Размер матрицы
    m = atoi(argv[2]);  ///Размер блоков
    p = atoi(argv[3]);  ///Число потоков

    if ((n == 0) || (m == 0) || (p == 0))
    {
        isWrongParam = true;
        return;
    }

    if (n < m)
    {
        cout << "Wrong parametrs" << endl;
        exit(0);
    }
    ///Рассчитаем количество блоков и выделим память
    l = n / m;      ///Количество блоков размера m
    s = n % m;      ///Размер крайнего нижнего блока

    if (s == 0)
    {
        s = m;
        l--;
    }
    if (n == m)
    {
        s = 0;
        l++;
    }
}

int find_max(double *A, double *B, int n, int k, int &max_i, int &max_j)
{
    ///A и B размера n
    ///k - шаг алгоритма

    #define A(i,j) A[(i) * n + (j)]
    #define B(i,j) B[(i) * n + (j)]
    double max = fabs(A(k,k)); max_i = k; max_j = k;

    ///Находим главный элемент
    for (int i = k; i < n; i++)
    {
        for (int j = k; j < n; j++)
        {
            if (fabs(A(i,j)) > max)
            {
                  max = fabs(A(i,j));
                max_i = i;
                max_j = j;
            }
        }
    }

    if (max < EPS)
    {
        return -1;
    }

    double temp = 0.;

    ///Переставляем строки A
    for (int j = 0; j < n; j++)
    {
        temp        = A(max_i, j);
        A(max_i, j) = A(k, j);
        A(k, j)     = temp;
    }

    ///Переставляем строки B
    for (int j = 0; j < n; j++)
    {
        temp        = B(max_i, j);
        B(max_i, j) = B(k, j);
        B(k, j)     = temp;
    }

    ///Переставляем столбцы A
    for (int i = 0; i < n; i++)
    {
        temp        = A(i, max_j);
        A(i, max_j) = A(i, k);
        A(i, k)     = temp;
    }

    #undef A
    #undef B
    return 1;
}

int Jordan_Inversion(double *A, double *B, int n, int *index)
{
    #define A(i,j) A[(i) * n + (j)]
    #define B(i,j) B[(i) * n + (j)]

    int i = 0; int j = 0; int k = 0;

    for (int i = 0; i < n; ++i)
        index[i] = i;

    double temp = 0.;
    int max_j   = 0;
    int max_i   = 0;
    int max     = 0;

    for (k = 0; k < n; k++)
    {
        max = 0;
        max = find_max(A, B, n, k, max_i, max_j);
        if (max == -1)
          return -1;

        j            = index[k];
        index[k]     = index[max_j];
        index[max_j] = j;

        temp = A(k,k);

        ///Делим "первую"(в данном шаге алгоритма) строку на первый элемент
        for (j = k + 1; j < n; j++)
            A(k, j) = A(k, j) / temp;

        A(k,k) = 1;

        ///Делим "первую"(в данном шаге алгоритма) строку на первый элемент
        for (j = 0; j < n; j++)
            B(k, j) = B(k, j) / temp;

        ///Вычитаем "первую" строку из нижних и верхних
        for (i = 0; i < k; i++)
            for (j = k + 1; j < n; j++)
                A(i,j) = A(i,j) - A(i,k) * A(k,j);

        for (i = k + 1; i < n; i++)
            for (j = k + 1; j < n; j++)
                A(i,j) = A(i,j) - A(i,k) * A(k,j);

        ///Вычитаем "первую" строку из нижних и верхних
        for (i = 0; i < k; i++)
            for (j = 0; j < n; j++)
                B(i,j) = B(i,j) - A(i,k) * B(k,j);

        for (i = k+1; i < n; i++)
            for ( j = 0; j < n; j++)
                B(i,j) = B(i,j) - A(i,k)*B(k,j);

    }

    ///Переставляем строки B обратно
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            A(index[i],j) = B(i,j);

    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            B(i,j) = A(i,j);

    #undef A
    #undef B
    return 1;
}

void multiplication(double *A, double *B, double *C, int n, int m, int l){
    #define A(i,j) A[(i) * m + (j)]
    #define B(i,j) B[(i) * l + (j)]
    #define C(i,j) C[(i) * l + (j)]

    int i, j, k;

    double s0, s1, s2, s3, s4, s5, s6, s7, s8;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < l; j++)
            C(i,j) = 0;

    int on = n - n % 3, ol = l - l % 3;

    for (i = 0; i < on; i+=3)
    {
        for (j = 0; j < ol; j+=3)
        {
        s0 = s1 = s2 = s3 = s4 = s5 = s6 = s7 = s8 = 0.;
            for (k = 0; k < m; k++)
            {
                s0 += A(i,k)   * B(k,j);
                s1 += A(i+1,k) * B(k,j);
                s2 += A(i+2,k) * B(k,j);
                s3 += A(i,k)   * B(k,j+1);
                s4 += A(i+1,k) * B(k,j+1);
                s5 += A(i+2,k) * B(k,j+1);
                s6 += A(i,k)   * B(k,j+2);
                s7 += A(i+1,k) * B(k,j+2);
                s8 += A(i+2,k) * B(k,j+2);
            }
            C(i,j)      += s0;
            C(i+1,j)    += s1;
            C(i+2,j)    += s2;
            C(i,j+1)    += s3;
            C(i+1,j+1)  += s4;
            C(i+2, j+1) += s5;
            C(i,j+2)    += s6;
            C(i+1,j+2)  += s7;
            C(i+2,j+2)  += s8;
        }
    }

    for (i = 0; i < n; i++)
    {
        for (j = ol; j < l; j++)
        {
            C(i,j) = 0;
            for (k = 0; k < m; k++)
                C(i,j) += A(i,k) * B(k,j);
        }
    }

    for (i = on; i < n; i++)
    {
        for (j = 0; j < l; j++)
        {
            C(i,j) = 0;
            for (k = 0; k < m; k++)
                C(i,j) += A(i,k) * B(k,j);
        }
    }

    #undef A
    #undef B
    #undef C
}


void subtract(double *A, double *B, int n, int m)
{
    #define A(i,j) A[(i) * m + (j)]
    #define B(i,j) B[(i) * m + (j)]

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            A(i,j) -= B(i,j);

    #undef A
    #undef B
}

void add(double *A, double *B, int n, int m)
{
    #define A(i,j) A[(i) * m + (j)]
    #define B(i,j) B[(i) * m + (j)]

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            A(i,j) += B(i,j);

    #undef A
    #undef B
}

void part_add(double *A, double *B, int n, int s, int m)
{
    m += 0;
    #define A(i,j) A[(i) * m + (j)]
    #define B(i,j) B[(i) * s + (j)]

    for (int i = 0; i < n; i++)
        for (int j = 0; j < s; j++)
            A(i,j) += B(i,j);

    #undef A
    #undef B
}


