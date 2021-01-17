CC=mpicc

DEBUG=-g -Wall -pedantic -std=c99 -lm


OPTIM=-O3



CFLAGS= $(DEBUG) $(OPTIM)


EXE = multiscale_parallele
SRC = elements_code.c

.SUFFIXES: .c .o

OBJS =  $(SRC:.c=.o)

all: $(EXE)

.c.o :
	$(CC)  -c $< $(CFLAGS)

$(EXE): $(OBJS)
	$(CC) -o $(EXE) $(OBJS) $(CFLAGS) $(LDFLAGS)

clean:
	rm -f $(OBJS) $(EXE)

mrproper: clean
	 rm -f $(EXE)
