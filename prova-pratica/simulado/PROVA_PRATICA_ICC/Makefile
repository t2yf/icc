V ?= 0

TARGET = prova
CC = gcc
CFLAGS = -Wall -Werror -DLIKWID_PERFMON -I$(LIKWID_INCLUDE)
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