
TD = ../

CSOURCES= floatfann.c
COBJECTS= $(CSOURCES:.c=.o)
CBUILT_OBJECTS= $(patsubst %,$(TD)$(OBJ_DIR)/%,$(COBJECTS))

INCLUDES= include/

all: $(CBUILT_OBJECTS)
.PHONY: all

$(CBUILT_OBJECTS): $(CSOURCES)
	@echo "  * compiling" $(*F).c
	@$(CC) -c -o $@ $(*F).c -O3 $(PLATFORM_FLAGS) -I$(INCLUDES) 

clean:
	@echo "Cleaning up."
	@rm -f $(TD)$(OBJ_DIR)/* $(TD)$(BIN_DIR)/*

.PHONY: clean
