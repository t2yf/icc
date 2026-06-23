V ?= 0 

TARGET = prova
GCC = gcc 
CFLAGS = -DV=$(V) -Wall -Werror -O3 -mavx2 -march=native -DLIKWID_PERFMON -I$(LIKWID_INCLUDE) 
LFLAGS = -L$(LIKWID_LIB) -llikwid -lm 

SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

%.o:%.c 
	$(CC) $(CFLAGS) -c $< -o $@ 

clean:
	rm -f *.o $(TARGET)