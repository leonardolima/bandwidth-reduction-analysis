CXX = g++
CXXFLAGS = -I /usr/local/include/eigen3 -std=c++17 -Wall -pedantic-errors -g

SRCS = main.cpp cuthill_mckee.cpp gauss_jordan.cpp tridiagonal.cpp crank_nicolson.cpp lu_decomposition.cpp ../io.cpp ../basic.cpp compare_solvers.cpp
OBJS = ${SRCS:.cpp=.o}
HEADERS = cuthill_mckee.h gauss_jordan.h tridiagonal.h crank_nicolson.h lu_decomposition.h ../io.h ../basic.h compare_solvers.h

MAIN = main

all: ${MAIN}

${MAIN}: ${OBJS}
	${CXX} ${CXXFLAGS} ${OBJS} -o ${MAIN}

.cpp.o:
	${CXX} ${CXXFLAGS} -c $< -o $@

clean:
	${RM} ${PROGS} ${OBJS} *.o *~.
