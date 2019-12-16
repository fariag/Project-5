In the programs folder you can find all the codes for Project 5 In the results folder you can find the results gathered for the project

Here's a list of how to compile and execute each program:

5c)
g++ -O3 5_c.cpp lib.cpp -o [Program name] -larmadillo -std=c++11

./5_c.exe [Min alpha] [Max alpha] [Omega/w]

Output:

Alpha | Acceptance rate | Mean Energy | Mean Variance | Average distance

Example:

g++ -O3 5_c.cpp lib.cpp -o 5_c.exe -larmadillo -std=c++11

./5_c.exe 0.75 1 1.0

Output:

0.80 0.505182 3.78686 0.274853 1.78306

5d)
g++ -O3 5_d.cpp lib.cpp -o [Program name] -larmadillo -std=c++11

./5_d.exe [Min alpha] [Max alpha] [Min beta] [Max beta] [w/Omega]

Output:

Beta | Acceptance rate | Mean Energy | Mean Variance | Average distance

Example:

g++ -O3 5_d.cpp lib.cpp -o 5_d.exe -larmadillo -std=c++11

./5_d.exe 0.75 1.0 0.2 0.3 1

Output:

0.30 0.53162 2.00321 0.00328482 2.57145

5e)
g++ -O3 5_e.cpp lib.cpp -o [Program name] -larmadillo -std=c++11

./5_e.exe [alpha] [beta] [Min omega/w] [Max omega/w]

Output:

Omega | Mean Energy | Mean Variance | Mean Kinetic energy | Mean Potential energy

Example:

g++ -O3 5_e.cpp lib.cpp -o 5_e.exe -larmadillo -std=c++11

./5_e.exe 0.99 0.29 0.75 1

Output:

1.00 3.73022 0.000263698 1.35555 2.37467
