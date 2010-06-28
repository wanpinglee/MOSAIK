# ================================
# MOSAIK Mosaik Assembler Makefile
# (c) 2008 Michael Stromberg
# ================================

# -------------------
# define our includes
# -------------------
TD = ../
INCLUDES = -I$(TD)$(ASSEMBLY_DIR) -I$(TD)$(CONFIG_DIR) -I$(TD)$(DATA_STRUCT_DIR) -I$(TD)$(MOSAIKREAD_DIR) -I$(TD)$(UTILITIES_DIR)

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES= AssemblerMain.cpp MosaikAssembler.cpp
OBJECTS= $(SOURCES:.cpp=.o)
_EXT_OBJECTS= AbstractAssemblyFormat.o Ace.o AlignmentReader.o Benchmark.o ConsoleUtilities.o fastlz.o FastLZIO.o FileUtilities.o GigaBayesFormat.o MemoryUtilities.o MosaikString.o Options.o ReferenceSequenceReader.o SequenceUtilities.o SHA1.o
EXT_OBJECTS=$(patsubst %,$(TD)$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(TD)$(OBJ_DIR)/%,$(OBJECTS))
PROGRAM= MosaikAssembler
LIBS=$(PLATFORM_LIBS) -lpthread -lz

all: $(PROGRAM)

.PHONY: all

$(PROGRAM): $(BUILT_OBJECTS) $(EXT_OBJECTS)
	@echo "  * linking $(PROGRAM)"
	@$(CXX) $(LDFLAGS) $(CXXFLAGS) $(PLATFORM_FLAGS) -o $(TD)$(BIN_DIR)/$@ $^ $(LIBS)

$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX) -c -o $@ $(*F).cpp $(CXXFLAGS) $(PLATFORM_FLAGS) $(INCLUDES)

$(EXT_OBJECTS):
	@$(MAKE) --no-print-directory -C $(TD)$(ASSEMBLY_DIR)
	@$(MAKE) --no-print-directory -C $(TD)$(DATA_STRUCT_DIR)
	@$(MAKE) --no-print-directory -C $(TD)$(MOSAIKREAD_DIR)
	@$(MAKE) --no-print-directory -C $(TD)$(UTILITIES_DIR)

clean:
	@echo "Cleaning up."
	@rm -f $(TD)$(OBJ_DIR)/* $(TD)$(BIN_DIR)/*

.PHONY: clean
