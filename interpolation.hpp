#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <assert.h>

class Polynomials
{
protected:
    std::vector <double> x;
    std::vector <double> f;
    std::vector <double> xx;
    unsigned int n; // n is the size of an array with input points,
    unsigned int N; // N is the size of an array with points to interpolate, please, don't confuse them
public:
    Polynomials (std::vector <double> &points, std::vector <double> &func, std::vector <double> &new_pts);
    std::vector <double> ff;
    virtual ~Polynomials ()
    {
    }
};

Polynomials::Polynomials (std::vector <double> &points, std::vector <double> &func, std::vector <double> &new_pts)
{
    if (points.size() != func.size())
    {
        std::cout<<"Input vectors must be the same size!"<<std::endl;
        exit (-1);
    }
    if (new_pts.size() == 0)
    {
        std::cout<<"Add new points"<<std::endl;
        exit (-1);
    }
    n = points.size ();
    N = new_pts.size ();
    x = points;
    f = func;
    xx = new_pts;
    ff = std::vector <double> (xx.size (), 0.0);
}

class Lagrange_polynomial: public Polynomials
{
private:
    double basis_fun;
public:
    Lagrange_polynomial (std::vector <double> &points, std::vector <double> &func, std::vector <double> &new_pts) : 
            Polynomials (points, func, new_pts)
            {
                basis_fun = 1;
            }
    void Interpolate ();     
};

void Lagrange_polynomial::Interpolate ()
{
    for (unsigned int k = 0; k < N; k++)
    {
        for (unsigned int i = 0; i < n; i++)
        {
            for (unsigned j = 0; j < n; j++)
                if (i != j)
                    basis_fun *= (xx [k] - x [j]) / (x [i] - x [j]);
            ff [k] += f [i] * basis_fun;
            basis_fun = 1;
        }
    }
}

class Newton_polynomial: public Polynomials
{
private:
    std::vector <double> div_dif;
    double tmp;
    unsigned int index (unsigned int i, unsigned int j, unsigned int size);
public:
    Newton_polynomial (std::vector <double> &points, std::vector <double> &func, std::vector <double> &new_pts) :  
    Polynomials (points, func, new_pts)
    {
        tmp = 1;
        div_dif = std::vector <double> ((n - 1) * n / 2, 0.0); //divided differnces are stored in an uppper-triangular matrix
    }                                                          //since 0-order divided differnces are stored in x vector, the matrix is (n-1)*(n-1)
    void Interpolate ();
};

unsigned int Newton_polynomial::index (unsigned int i, unsigned int j, unsigned int size)
{
    if (size <= 0)
    {
        std::cout << "Incorrect size" << std::endl;
        exit (-1);
    }
    if (i >= size || j >= size)
    {
        std::cout << "Incorrect indices" << std::endl;
        exit (-1);
    }
    return size * (size - 1) / 2 - (size - i) * (size - i - 1) / 2 + j - i + 1;
}

void Newton_polynomial::Interpolate ()
{
    for (unsigned int i = 0; i < n - 1; i ++)
        for (unsigned int j = 0; j < n - i - 1; j++)
        {
            if (!i)
                div_dif [index (i, j, n - 1)] = (f [j + 1] - f [j]) / (x [j + 1] - x [j]);
            else
                div_dif [index (i, j + i, n - 1)] = (div_dif [index (i - 1, j + i, n - 1)] - 
                                                 div_dif [index (i - 1, j + i - 1, n - 1)]) / (x [j + i + 1] - x [j]); 
        }
    for (unsigned int i = 0; i < N; i++)
    {
        ff [i] += f[0];
        for (unsigned int j = 0; j < n - 1; j++)
        {
            for (unsigned int k = 0; k <= j; k++)
                tmp *= xx[i] - x[k];
            ff [i] += div_dif [index (j, j, n - 1)] * tmp;
            tmp = 1;
        }
    }
}

class Spline: public Polynomials
{
private:
    std::vector <double> A; //A is a tridiagonal symmetric matrix
    std::vector <double> F; 
    std::vector <double> coef;
    std::vector <double> p;
    std::vector <double> q;
    std::vector <double> b;
    std::vector <double> d;
    unsigned int k_u;
    void Sweep ();
public:
    Spline (std::vector <double> &points, std::vector <double> &func, std::vector <double> &new_pts) : 
    Polynomials (points, func, new_pts)
    {
        if (n < 4)
        {
            std::cout << "Not enough points" << std::endl;
            exit (-1);
        }
        k_u = 0;
        A = std::vector <double> (2 * n - 5, 0.0); //first n - 2 elements form the main diagonal
        F = std::vector <double> (n - 2, 0.0);
        coef = std::vector <double> (n, 0.0);
        for (unsigned int i = 0; i < n - 2; i++)
        {
            A [i] = (x [i + 2] - x [i]) / 3;
            if (i < n - 3)
                A [i + n - 2] = (x [i + 1] - x [i]) / 6;
            F [i] = (f [i + 2] - f [i + 1]) / (x [i + 2] - x [i + 1]) - (f [i + 1] - f [i]) / (x [i + 1] - x [i]); 
        }
        p = std::vector <double> (n - 3, 0.0);
        q = std::vector <double> (n - 3, 0.0);
        b = std::vector <double> (n - 1, 0.0);
        d = std::vector <double> (n - 1, 0.0);
    }
    void Interpolate ();
};

void Spline::Sweep ()
{
    p [0] = A [n - 2] /(-A [0]);
    q [0] = F [0] / A [0]; 
    for (unsigned int i = 1; i < n - 3; i++)
    {
        p [i] = A [n - 2 + i] / ((-A [i]) - A [n - 2 + i - 1] * p [i - 1]);
        q [i] = (A [n - 2 + i - 1] * q [i - 1] - F [i]) / ((-A [i]) - A [n - 2 + i - 1] * p [i - 1]);
    }
    coef [n - 2] = (F [n - 3] - A [2 * n - 6] * q [n - 4]) / (A [2 * n - 6] * p [n - 4] + A [n - 3]);
    for (int i = n - 3; i > 0; i--)
        coef [i] = p [i - 1] * coef [i + 1] + q [i - 1];
}

void Spline::Interpolate ()
{
    Sweep ();
    for (unsigned int i = 0; i < N; i++)
    {
        if (xx [i] < x [0] || xx [i] > x [n - 1])
        {
            std::cout << "The new point is out of segment" << std::endl;
            exit (-1);
        }
        k_u = std::upper_bound (x.begin (), x.end (), xx [i]) - x.begin ();
        if (k_u == x.size ())
            k_u--;
        b [n - 2] = (f [n - 1] - f [n - 2]) / (x [n - 1] - x [n - 2]) - 2 / 3 * (x [n - 1] - x [n - 2]) * coef [n - 2];
        d [n - 2] = (-coef [n - 2]) / 3 / (x [n - 1] - x [n - 2]);
        for (unsigned int j = 0; j < n - 2; j++)
        {
            b [j] = (f [j + 1] - f [j]) / (x [j + 1] - x [j]) - (x [j + 1] - x [j]) / 3 * (coef [j + 1] + 2 * coef [j]);
            d [j] = (coef [j + 1] - coef [j]) / 3 / (x [j + 1] - x [j]);
        }
        ff [i] = f [k_u - 1] + b [k_u - 1] * (xx [i] - x [k_u - 1]) + coef [k_u - 1] * pow((xx [i] - x[k_u - 1]), 2) + d[k_u - 1] * pow ((xx [i] - x [k_u - 1]), 3); 
    }
}

#endif
