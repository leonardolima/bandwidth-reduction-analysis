CXX = g++
CXXFLAGS = -I /usr/local/include/eigen3 -std=c++17 -Wall -pedantic-errors -g

SRCS = main.cpp peak_mem.cpp adapted_cuthill_mckee.cpp basic.cpp bandwidth_analysis/cuthill_mckee.cpp io.cpp depth.cpp
OBJS = ${SRCS:.cpp=.o}
HEADERS = peak_mem.h adapted_cuthill_mckee.h basic.h bandwidth_analysis/cuthill_mckee.h io.h depth.h

MAIN = main

all: ${MAIN}

${MAIN}: ${OBJS}
	${CXX} ${CXXFLAGS} ${OBJS} -o ${MAIN}

.cpp.o:
	${CXX} ${CXXFLAGS} -c $< -o $@

clean:
	${RM} ${PROGS} ${OBJS} *.o *~.
