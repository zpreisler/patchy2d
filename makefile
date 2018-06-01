CC=gcc
CFLAGS=-O2 -Wall -Wextra
CLIBS=-I $(HOME)/include -L $(HOME)/lib -lutils -lSDL2 -lSDL2main -lGL -lGLEW 
EXEC=sdl
INSTALL=install -m 111
BINDIR=$(HOME)/bin/

all: main.c
	$(CC) $(CFLAGS) main.c -o $(EXEC) $(CLIBS)

install:
	$(INSTALL) $(EXEC) $(BINDIR)
