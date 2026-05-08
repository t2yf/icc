       PROG   = matmult

           CC = gcc -std=c11 -g
         OBJS = matriz.o

       CFLAGS =  -I${LIKWID_HOME}/include -DLIKWID_PERFMON -O3 -mavx2 -march=native
       LFLAGS = -L${LIKWID_HOME}/lib -llikwid -lm 

# Lista de arquivos para distribuição
DISTFILES = *.c *.h *.txt Makefile Resultados perfctr 
DISTDIR = `basename ${PWD}`

.PHONY: all debug clean limpa purge faxina dist

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

all: $(PROG)

debug:         CFLAGS += -DDEBUG
debug:         $(PROG)

$(PROG):  $(PROG).o

$(PROG): $(OBJS) 
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

clean limpa:
	@echo "Limpando ...."
	@rm -f *~ *.bak *.tmp

purge faxina:   clean
	@echo "Faxina ...."
	@rm -f  $(PROG) *.o core a.out
	@rm -f *.png marker.out *.log

dist: purge
	@echo "Gerando arquivo de distribuição ($(DISTDIR).tar) ..."
	@ln -s . $(DISTDIR)
	@tar -chvf $(DISTDIR).tar $(addprefix ./$(DISTDIR)/, $(DISTFILES))
	@rm -f $(DISTDIR)

