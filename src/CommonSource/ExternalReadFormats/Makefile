# =====================================
# MOSAIK External Read Formats Makefile
# (c) 2008 Michael Stromberg
# =====================================

# -------------------
# define our includes
# -------------------
TD = ../../
INCLUDES = -I$(TD)$(DATA_STRUCT_DIR) -I$(TD)$(CONFIG_DIR) -I$(TD)$(UTILITIES_DIR) -I$(TD)$(MOSAIKREAD_DIR) -I.

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= Fasta.cpp Fastq.cpp SRF.cpp BamWriter.cpp
OBJECTS= $(SOURCES:.cpp=.o)
BUILT_OBJECTS= $(patsubst %,$(TD)$(OBJ_DIR)/%,$(OBJECTS))

all: $(BUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(CXXFLAGS) $(PLATFORM_FLAGS) $(INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(TD)$(OBJ_DIR)/* $(TD)$(BIN_DIR)/*

.PHONY: clean
