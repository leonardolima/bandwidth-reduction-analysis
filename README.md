# Peak Memory Minimization

This project aims to adapt, implement and analyze a variety of solving methods for the [bandwidth minimization problem](docs/references/cuthill1969.pdf) in order to minimize peak memory consumption of programs. You can find a more detailed explanation [here](docs/slides/presentation-10-07-2019.pdf).

## Getting Started

These are the basic steps to take if you need to run the project on your local machine.

### Prerequisites

[Eigen](http://eigen.tuxfamily.org), a C++ template library for linear algebra, is used for specifying and manipulating matrices.

### Running

Currently, just the motivation of the work is fully working, more specifically the impact of bandwidth reduction when numerically solving the heat equation using the Crank-Nicolson method. In the [main.cpp](src/bandwidth_analysis/main.cpp) file, you can change the parameters accordingly.

Inside the [src/bandwidth_analysis](src/bandwidth_analysis) and [src/](src/) directories there are **Makefiles**. You can run 

```
$ make
```

And

```
$ ./main
```

You may need to change the permissions of **main** in order to execute it. When finished, just type

```
$ make clean
```

## License

This project is licensed under the GPL-3.0 license - see [LICENSE](LICENSE) for details.

