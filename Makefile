# ==========================
# MOSAIK Makefile
# (c) 2009 Michael Stromberg
# ==========================

# define our object and binary directories
export OBJ_DIR = ../obj
export BIN_DIR = ../bin

# define our common source directories
export ASSEMBLY_DIR = CommonSource/AssemblyFormats
export CONFIG_DIR = CommonSource/Config
export DATA_STRUCT_DIR = CommonSource/DataStructures
export EXTERNAL_DIR = CommonSource/ExternalReadFormats
export MOSAIKREAD_DIR = CommonSource/MosaikReadFormat
export PAIRWISE_DIR = CommonSource/PairwiseAlignment
export UTILITIES_DIR = CommonSource/Utilities

# define some default flags
export CFLAGS ?= -Wall -O3 -Wno-char-subscripts
export CXXFLAGS ?= $(CFLAGS)
export LDFLAGS ?= -Wl,-s
export CXX ?= g++

# define our platform
export BLD_PLATFORM ?= linux
include includes/$(BLD_PLATFORM).inc

# define our source subdirectories
SUBDIRS = CommonSource MosaikBuild MosaikAligner MosaikSort MosaikMerge MosaikAssembler
UTIL_SUBDIRS = CommonSource MosaikJump MosaikCoverage MosaikText MosaikDupSnoop

all:
	@echo "Building MOSAIK for the following platform: $$BLD_PLATFORM"
	@echo "========================================================="
	@./UpdateConfigH
	@for dir in $(SUBDIRS); do \
		echo "- Building in $$dir"; \
		$(MAKE) --no-print-directory -C $$dir; \
		echo ""; \
	done

utils:
	@echo "Building MOSAIK for the following platform: $$BLD_PLATFORM"
	@echo "========================================================="
	@./UpdateConfigH
	@for dir in $(UTIL_SUBDIRS); do \
		echo "- Building in $$dir"; \
		$(MAKE) --no-print-directory -C $$dir; \
		echo ""; \
	done

.PHONY: all

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean
