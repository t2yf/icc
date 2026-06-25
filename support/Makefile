# OTPIM ?= 0
# include "likwid.h"

TARGET = prova
CC = gcc 
CFLAGS = -Wall -Werror -DLIKWID_PERFMON -I$(LIKWID_INCLUDE) #primeiro -DOPTIM=$(OPTIM) entre -werroe e likwid -O3 -march=native -maxv -fopt-info-vec 
LFLAGS = -L$(LIKWID_LIB) -llikwid -lm 

SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@ 

clean: 
	rm -f *.o $(TARGET)