BIN=bin
OBJ=obj

all: 
	@mkdir -p $(OBJ)
	@mkdir -p $(BIN)
	cd src; make

clean:
	rm -rf $(BIN)/*
	rm -rf $(OBJ)/*
