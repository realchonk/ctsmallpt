
CXX=g++

CEXPR_DEPTH=8290000000

CXXFLAGS=-O3 -Wall -Wextra -std=c++20 -fopenmp

ifeq ($(CXX),g++)
CXXFLAGS += -fconstexpr-ops-limit=$(CEXPR_DEPTH)
else ifeq ($(CXX),clang++)
CXXFLAGS += -fconstexpr-depth=$(CEXPR_DEPTH) -fconstexpr-steps=$(CEXPR_DEPTH)
endif


smallpt: smallpt.cpp
	$(CXX) -o $@ $< $(CXXFLAGS)

clean:
	rm -f smallpt

run: smallpt
	./$< 1
	sxiv image.ppm

.PHONY: all clean run
