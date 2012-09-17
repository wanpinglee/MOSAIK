# ==========================
# MOSAIK Makefile
# (c) 2012 Wan-Ping Lee
# ==========================

FLAGS = -Wall -O3


all: xc_pe.cpp xc_se.cpp
	@$(CXX) $(FLAGS) xc_pe.cpp -o xc_pe
	@$(CXX) $(FLAGS) xc_se.cpp -o xc_se

.PHONY: all

clean:
	@echo "Cleaning up attachXC."
	@rm xc_pe
	@rm xc_se
