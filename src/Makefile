CXX = g++
CXXFLAGS = -I /usr/local/include/eigen3 -std=c++0x -Wall -pedantic-errors -g

SRCS = main.cpp bandwidth_minimization.cpp gauss_jordan.cpp tridiagonal.cpp crank_nicolson.cpp
OBJS = ${SRCS:.cpp=.o}
HEADERS = bandwidth_minimization.h gauss_jordan.h tridiagonal.h crank_nicolson.h

MAIN = main

all: ${MAIN}

${MAIN}: ${OBJS}
	${CXX} ${CXXFLAGS} ${OBJS} -o ${MAIN}

.cpp.o:
	${CXX} ${CXXFLAGS} -c $< -o $@

clean:
	${RM} ${PROGS} ${OBJS} *.o *~.
