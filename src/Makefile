CXX = g++
CXXFLAGS = -I /usr/local/include/eigen3 -std=c++17 -Wall -pedantic-errors -g

SRCS = main.cpp peak_mem.cpp cuthill_mckee.cpp basic.cpp
OBJS = ${SRCS:.cpp=.o}
HEADERS = peak_mem.h cuthill_mckee.h basic.h

MAIN = main

all: ${MAIN}

${MAIN}: ${OBJS}
	${CXX} ${CXXFLAGS} ${OBJS} -o ${MAIN}

.cpp.o:
	${CXX} ${CXXFLAGS} -c $< -o $@

clean:
	${RM} ${PROGS} ${OBJS} *.o *~.
