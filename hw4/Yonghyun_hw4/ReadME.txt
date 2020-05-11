Problem 1
How to compile:
gcc -c mat_vec.c -Wall -pedantic
gcc -c eval.c -Wall -pedantic
gcc -o testeval testeval.c eval.o mat_vec.o -Wall -pedantic -lm -llapack

Problem 2
How to compile:
gcc -Wall -pedantic -o RmathGamma RmathGamma.c -lRmath
How to run:
./RmathGamma

Problem 3
How to compile:
R CMD SHLIB skewkt_Call.c
How to run a program:
Rscript run_skewkt_Call.R