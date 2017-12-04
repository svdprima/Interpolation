#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include "interpolation.hpp"

/*
double func (double x)
{
    return cos (x);
    //return (x != 0) ? sin (x) / x : 1;
    //return (1 + 0.05 * cos (0.5 * x)) * cos (x);
}
*/

std::vector <double> Range (double from, double bin, double to)
{
    if (bin >= to - from)
    {
        std::cout << "Cannot create linear array" << std::endl;
        exit (-1);
    }
    std::vector <double> result;
    while (from <= to)
    {
        result.push_back (from);
        from += bin;
    }
    return result;
}

int main ()
{
    const unsigned int N = 6;
    std::vector <double> x = std::vector <double> (N);
    std::vector <double> f = std::vector <double> (N);
    for (unsigned int i = 0; i < N; i++)
    {
        x [i] = i;
        //f [i] = func (x [i]);
        f [i] = pow ((-1), i);
    }
    for (unsigned int i = 0; i < N; i++)
    {
        std::cout << x[i] << " ";
    }
    std::cout<<std::endl;
    for (unsigned int i = 0; i < N; i++)
    {
        std::cout << f[i] << " ";
    }
    std::cout<<std::endl;
    std::vector <double> xx = Range (x[0], 0.1, x [N - 1]); 
#ifdef LAGRANGE
    Lagrange_polynomial L (x, f, xx);
    L.Interpolate ();
    {
        std::ofstream output;
        output.open ("Lagrange.txt");
        for (unsigned int i = 0; i < xx.size (); i++)
            output << xx [i] << " " << L.ff [i] << std::endl;
        output.close ();
    }
#endif
#ifdef NEWTON
    Newton_polynomial Np (x, f, xx);
    Np.Interpolate ();
    {
        std::ofstream output;
        output.open ("Newton.txt");
        for (unsigned int i = 0; i < xx.size (); i++)
            output << xx [i] << " " << Np.ff [i] << std::endl;
        output.close ();
    }
#endif
#ifdef SPLINE
    Spline S (x, f, xx);
    S.Interpolate ();
    {
        std::ofstream output;
        output.open ("Spline.txt");
        for (unsigned int i = 0; i < xx.size (); i++)
            output << xx [i] << " " << S.ff [i] << std::endl;
        output.close ();
    }
#endif
    return 0;
}
