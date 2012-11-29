# ==================================
# MOSAIK Pairwise Alignment Makefile
# (c) 2008 Michael Stromberg
# ==================================

# -------------------
# define our includes
# -------------------
TD = ../../
INCLUDES = -I$(TD)$(DATA_STRUCT_DIR) -I$(TD)$(CONFIG_DIR) -I$(TD)$(UTILITIES_DIR) -I.

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= SmithWatermanGotoh.cpp BandedSmithWaterman.cpp ssw_cpp.cpp
OBJECTS= $(SOURCES:.cpp=.o)
BUILT_OBJECTS= $(patsubst %,$(TD)$(OBJ_DIR)/%,$(OBJECTS))

CSOURCES = ssw.c
COBJECTS= $(CSOURCES:.c=.o)
CBUILT_OBJECTS= $(patsubst %,$(TD)$(OBJ_DIR)/%,$(COBJECTS))

all: $(BUILT_OBJECTS) $(CBUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(CXXFLAGS) $(PLATFORM_FLAGS) $(INCLUDES)

$(CBUILT_OBJECTS): $(CSOURCES)
	@echo "  * compiling" $(*F).c
	@$(CC) -c -o $@ $(*F).c -O3 $(PLATFORM_FLAGS) $(INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(TD)$(OBJ_DIR)/* $(TD)$(BIN_DIR)/*

.PHONY: clean
