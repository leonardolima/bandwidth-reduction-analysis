# Peak Memory Optimization

This project is an implementation of an adapted version of the [bandwidth minimization problem](docs/references/cuthill1969.pdf), namely the peak memory optimization problem. You can find a more detailed explanation [here](docs/slides/presentation-10-07-2019.pdf).

## Getting Started

These are the basic steps to take if you need to run the project on your local machine.

### Prerequisites

[Eigen](http://eigen.tuxfamily.org), a C++ template library for linear algebra, is used for specifying and manipulating matrices.

### Running

At this point, only the bandwidth minimization solution found [here](docs/references/cuthill1969.pdf) is already implemented. In the [bandwidth_minimization.cpp](src/bandwidth_minimization.cpp) file, you can find a function named *execute_algorithm()*, from which it is possible to change the matrix that you are currently working on.

In the [src](src/) directory there is a **Makefile**. You can run 

```
make
```

And

```
./main
```

You may need to change the permissions of **main** in order to execute it. When finished, just type

```
make clean
```

## License

This project is licensed under the GPL-3.0 license - see [LICENSE.md](LICENSE.md) for details.

