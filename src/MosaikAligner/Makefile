# =============================
# MOSAIK MosaikAligner Makefile
# (c) 2008 Michael Stromberg
# =============================

# -------------------
# define our includes
# -------------------
TD = ../
INCLUDES = -I$(TD)$(CONFIG_DIR) -I$(TD)$(DATA_STRUCT_DIR) -I$(TD)$(MOSAIKREAD_DIR) -I$(TD)$(PAIRWISE_DIR) -I$(TD)$(UTILITIES_DIR) -I$(TD)$(EXTERNAL_DIR) -I/home/wanping/tools/fann-2.1.0/src/include/

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= AlignerMain.cpp AlignmentThread.cpp MosaikAligner.cpp
OBJECTS= $(SOURCES:.cpp=.o)
_EXT_OBJECTS= AbstractDnaHash.o AlignmentQuality.o AlignmentReader.o AlignmentWriter.o ArchiveMerge.o ArchiveSort.o BamWriter.o BandedSmithWaterman.o Benchmark.o BestNSecondBestSelection.o CigarTager.o ColorspaceUtilities.o ConsoleUtilities.o disorder.o DnaHash.o fastlz.o FastLZIO.o FileUtilities.o HashRegionTree.o JumpDnaHash.o MdTager.o MosaikString.o MemoryUtilities.o MultiDnaHash.o NaiveAlignmentSet.o Options.o PairwiseUtilities.o ReadReader.o ReferenceSequenceReader.o RegexUtilities.o SequenceUtilities.o SHA1.o SmithWatermanGotoh.o SortThread.o StatisticsMaps.o TimeSupport.o UbiqDnaHash.o ZaTager.o floatfann.o
EXT_OBJECTS=$(patsubst %,$(TD)$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(TD)$(OBJ_DIR)/%,$(OBJECTS))
PROGRAM= MosaikAligner
LIBS=$(PLATFORM_LIBS) -lz -lpthread

all: $(PROGRAM)

.PHONY: all

$(PROGRAM): $(BUILT_OBJECTS) $(EXT_OBJECTS)
	@echo "  * linking $(PROGRAM)"
	@$(CXX) $(LDFLAGS) $(CXXFLAGS) $(PLATFORM_FLAGS) -o $(TD)$(BIN_DIR)/$@ $^ $(LIBS)

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(CXXFLAGS) $(PLATFORM_FLAGS) $(INCLUDES)

$(EXT_OBJECTS):
	@$(MAKE) --no-print-directory -C $(TD)$(DATA_STRUCT_DIR)
	@$(MAKE) --no-print-directory -C $(TD)$(MOSAIKREAD_DIR)
	@$(MAKE) --no-print-directory -C $(TD)$(PAIRWISE_DIR)
	@$(MAKE) --no-print-directory -C $(TD)$(UTILITIES_DIR)

clean:
	@echo "Cleaning up."
	@rm -f $(TD)$(OBJ_DIR)/* $(TD)$(BIN_DIR)/*

.PHONY: clean
