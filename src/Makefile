BIN=../bin
OBJ=../obj
HTS_ROOT=$(HOME)/src/htslib
CFLAGS=-g -D_FILE_OFFSET_BITS=64
#CFLAGS=-O3
CC=gcc

LIB= fehh.o \
     n_choose_2.o

.SUFFIXES:

OBJS=$(addprefix $(OBJ)/, $(LIB))


PROG=ehh

LIST=$(addprefix $(BIN)/, $(PROG))

all: $(LIST)
	@mkdir -p $(OBJ)
	@mkdir -p $(BIN)

$(OBJ)/%.o: %.c 
	$(CC) -c $(CFLAGS) -o $@ $< \
		-I$(HTS_ROOT)


$(BIN)/%: %.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ \
		-I$(HTS_ROOT) \
		-DSAMTOOLS=1 \
		$(HTS_ROOT)/libhts.a \
		-ll -ldl -lz -lm  -pthread
