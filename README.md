# aks-primality-test

A C implementation of Agrawal–Kayal–Saxena (AKS) primality test algorithm. 

Requires GMP library to handle big numbers.

A Makefile for linux-gnu is included. The code should work in multiple additional platforms and compilers.

## Author

Ting Liu (tingianliu@gmail.com)

## Building

To use OpenMP support:
```
make
make test
```

To disable OpenMP, and build single-threaded:
```
PARALLEL=false make clean all
make test
```

To experiment with different thread counts, build with OpenMP, and set the following environment 
variable when running:
```
OMP_NUM_THREADS=4 ./aks selected.txt
```

When `OMP_NUM_THREADS` is not set, it defaults to the number returned by the command:
```
grep -c ^processor /proc/cpuinfo
```

## References

Manindra Agrawal, Neeraj Kayal, Nitin Saxena, [PRIMES is in P](https://www.cse.iitk.ac.in/users/manindra/algebra/primality_v6.pdf), Annals of Mathematics 160(2): 781-793, 2004. The original version of the paper is [here](https://www.cse.iitk.ac.in/users/manindra/algebra/primality_original.pdf).

[OpenMP](https://www.openmp.org/) :: Resources :: [Tutorials](https://www.openmp.org/resources/tutorials-articles/)