CXX = g++
CXXFLAGS = -I /usr/local/include/eigen3 -std=c++17 -Wall -pedantic-errors -g

SRCS = main.cpp peak_mem.cpp bandwidth_analysis/cuthill_mckee.cpp adapted_cuthill_mckee.cpp basic.cpp io.cpp level.cpp topological.cpp
OBJS = ${SRCS:.cpp=.o}
HEADERS = peak_mem.h bandwidth_analysis/cuthill_mckee.h adapted_cuthill_mckee.h basic.h io.h level.h topological.h

MAIN = main

all: ${MAIN}

${MAIN}: ${OBJS}
	${CXX} ${CXXFLAGS} ${OBJS} -o ${MAIN}

.cpp.o:
	${CXX} ${CXXFLAGS} -c $< -o $@

clean:
	${RM} ${PROGS} ${OBJS} *.o *~.
