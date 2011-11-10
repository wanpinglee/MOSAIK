# ==========================
# MOSAIK Utilities Makefile
# (c) 2008 Michael Stromberg
# ==========================

# -------------------
# define our includes
# -------------------
TD = ../../
INCLUDES = -I$(TD)$(CONFIG_DIR) -I$(TD)$(DATA_STRUCT_DIR) -I$(TD)$(MOSAIKREAD_DIR) -I$(TD)$(EXTERNAL_DIR) -I$(TD)/fann-2.1.0/include/ -I.

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= ArchiveMerge.cpp \
		ArchiveSort.cpp \
		Benchmark.cpp \
		BestNSecondBestSelection.cpp \
		CigarTager.cpp \
		ColorspaceUtilities.cpp \
		ConsoleUtilities.cpp \
		Entropy.cpp \
		FastLZIO.cpp \
		FileUtilities.cpp \
		MdTager.cpp \
		MemoryUtilities.cpp \
		PairwiseUtilities.cpp \
		QualityNeuralNetwork.cpp \
		RegexUtilities.cpp \
		SequenceUtilities.cpp \
		SHA1.cpp \
		SortNMergeUtilities.cpp \
		SortThread.cpp \
		StatisticsMaps.cpp \
		TimeSupport.cpp \
		ZaTager.cpp
OBJECTS= $(SOURCES:.cpp=.o)
BUILT_OBJECTS= $(patsubst %,$(TD)$(OBJ_DIR)/%,$(OBJECTS))

CSOURCES = fastlz.c md5.c
COBJECTS= $(CSOURCES:.c=.o)
CBUILT_OBJECTS= $(patsubst %,$(TD)$(OBJ_DIR)/%,$(COBJECTS))

all: $(BUILT_OBJECTS) $(CBUILT_OBJECTS)

.PHONY: all

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(CXXFLAGS) $(PLATFORM_FLAGS) $(INCLUDES)

$(CBUILT_OBJECTS): $(CSOURCES)
	@echo "  * compiling" $(*F).c
	@$(CC) -c -o $@ $(*F).c -O3 -w -DSQLITE_OMIT_LOAD_EXTENSION $(PLATFORM_FLAGS) $(INCLUDES)

clean:
	@echo "Cleaning up."
	@rm -f $(TD)$(OBJ_DIR)/* $(TD)$(BIN_DIR)/*

.PHONY: clean
