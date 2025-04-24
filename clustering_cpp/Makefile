# Use GUROBI_HOME from environment if set, otherwise allow override on the command line
GUROBI_HOME ?= /opt/gurobi1100/linux64

sources = main.cpp cg_reduction.cpp cg_stabilization.cpp cg_util.cpp cg_pricing.cpp CGStruct.cpp Kmeans2D.cpp cg_pricing_multi.cpp cg_branching.cpp branch_and_price.cpp JobQueue.cpp ThreadPool.cpp
objects = $(sources:.cpp=.o)

all: bb

bb: $(objects)
	g++ -L $(GUROBI_HOME)/lib/ -o $@ $^ -pthread -lgurobi_c++ -lgurobi110 -Wl,-rpath,$(GUROBI_HOME)/lib/

%.o: %.cpp
		g++ -g -Wall -O2 -o $@ -c $< -I $(GUROBI_HOME)/include/


.PHONY: clean
		
clean:
		rm -f bb *.o
