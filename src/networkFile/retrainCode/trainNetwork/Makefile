# ==========================
# MOSAIK Makefile
# (c) 2012 Wan-Ping Lee
# ==========================

OBJ_DIR = ./obj

FLAGS =  -Wall -O3
CFLAGS = -O3
FANN=../../../fann-2.1.0

SOURCES = sam_parser_float.cpp parameter_parser_float.cpp mq_train_float.cpp
PROGRAM = train_mq

all: $(FANN)/floatfann.c $(SOURCES)
	@test -d $(OBJ_DIR) || mkdir $(OBJ_DIR)
	@echo "  * compiling ......"
	@$(CC) $(CFLAGS) -c -o $(OBJ_DIR)/floatfann.o $(FANN)/floatfann.c -I$(FANN)/include
	@$(CXX) -c $(FLAGS) -o $(OBJ_DIR)/sam_parser_float.o sam_parser_float.cpp
	@$(CXX) -c $(FLAGS) -o $(OBJ_DIR)/parameter_parser_float.o parameter_parser_float.cpp
	@$(CXX) -c $(FLAGS) -o $(OBJ_DIR)/mq_train_float.o mq_train_float.cpp -I$(FANN)/include
	@echo "  * linking ......"
	@$(CXX) $(FLAGS) $(OBJ_DIR)/*.o -o $(PROGRAM)

.PHONY: all

clean:
	@echo "Cleaning up trainNetwork."
	@rm -rf $(OBJ_DIR)
	@rm $(PROGRAM)
