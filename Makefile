#Makefile

CC = g++
CFLAGS = -g -Wall
CFLAGS += -I ./seqan
LIBS =

all:action read_input

action:
	@echo Compiling...

read_input:Main.cpp build_chains.cpp read_fasta.cpp table_entry.cpp RNA_seq.cpp
	${CC} ${CFLAGS} -o $@ $^ ${LIBS}

clean:
	@echo Cleaning...
	rm -f read_input
