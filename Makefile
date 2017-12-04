all: Lagrange Newton Spline
Lagrange: interpolation.hpp main.cpp
	g++ main.cpp -Wall -Wextra -DLAGRANGE -o Lagrange
Newton: interpolation.hpp main.cpp
	g++ main.cpp -Wall -Wextra -DNEWTON -o Newton
Spline:interpolation.hpp main.cpp
	g++ main.cpp -Wall -Wextra -DSPLINE -g -o Spline
