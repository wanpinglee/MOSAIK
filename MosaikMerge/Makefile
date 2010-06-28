# ============================
# MOSAIK Mosaik Merge Makefile
# (c) 2008 Michael Stromberg
# ============================

# -------------------
# define our includes
# -------------------
TD = ../
INCLUDES = -I$(TD)$(CONFIG_DIR) -I$(TD)$(DATA_STRUCT_DIR) -I$(TD)$(UTILITIES_DIR) -I$(TD)$(MOSAIKREAD_DIR) -I.

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= MergeMain.cpp MosaikMerge.cpp
OBJECTS= $(SOURCES:.cpp=.o)
_EXT_OBJECTS= AlignmentReader.o AlignmentQuality.o AlignmentWriter.o Benchmark.o ConsoleUtilities.o fastlz.o FastLZIO.o FileUtilities.o MemoryUtilities.o MosaikString.o Options.o NaiveAlignmentSet.o SequenceUtilities.o SHA1.o TimeSupport.o
EXT_OBJECTS=$(patsubst %,$(TD)$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(TD)$(OBJ_DIR)/%,$(OBJECTS))
PROGRAM= MosaikMerge
LIBS= $(PLATFORM_LIBS) -lpthread

all: $(PROGRAM)

.PHONY: all

$(PROGRAM): $(BUILT_OBJECTS) $(EXT_OBJECTS)
	@echo "  * linking $(PROGRAM)"
	@$(CXX) $(LDFLAGS) $(CXXFLAGS) $(PLATFORM_FLAGS) -o $(TD)$(BIN_DIR)/$@ $^ $(LIBS)

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(CXXFLAGS) $(PLATFORM_FLAGS) $(INCLUDES)

$(EXT_OBJECTS):
	@$(MAKE) --no-print-directory -C $(TD)$(MOSAIKREAD_DIR)
	@$(MAKE) --no-print-directory -C $(TD)$(DATA_STRUCT_DIR)
	@$(MAKE) --no-print-directory -C $(TD)$(UTILITIES_DIR)

clean:
	@echo "Cleaning up."
	@rm -f $(TD)$(OBJ_DIR)/* $(TD)$(BIN_DIR)/*

.PHONY: clean
