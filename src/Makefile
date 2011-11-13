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
FLAGS = -Wall -Wno-char-subscripts -ansi -O3
CFLAGS = -std=C99
CXXFLAGS =
#CXXFLAGS = -ansi -pedantic -Wextra -Weffc++
CFLAGS += $(FLAGS)
CXXFLAGS += $(FLAGS)
export CFLAGS
export CXXFLAGS
export LDFLAGS = -Wl
export CXX ?= g++

# define our platform
export BLD_PLATFORM ?= linux
include includes/$(BLD_PLATFORM).inc

# define our source subdirectories
SUBDIRS = fann-2.1.0 CommonSource MosaikBuild MosaikAligner MosaikJump MosaikText
#UTIL_SUBDIRS = MosaikJump

all:
	@test -d $(OBJ_DIR) || mkdir $(OBJ_DIR)
	@test -d $(BIN_DIR) || mkdir $(BIN_DIR)
	@echo "Building MOSAIK for the following platform: $$BLD_PLATFORM"
	@echo "========================================================="
#	@./UpdateConfigH
	@for dir in $(SUBDIRS); do \
		echo "- Building in $$dir"; \
		$(MAKE) --no-print-directory -C $$dir; \
		echo ""; \
	done

.PHONY: all

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean
