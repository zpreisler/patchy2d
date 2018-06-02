CC=gcc
CFLAGS=-O2 -Wall -Wextra
CLIBS=-I $(HOME)/include -L $(HOME)/lib -lutils -lGL -lGLEW 
SDL_LIBS=-lSDL2 -lSDL2main
EXEC=sdl
INSTALL=install -m 111
BINDIR=$(HOME)/bin/

MAIN=main.c
SOURCES=program_gl.c mySDL.c
OBJS=$(SOURCES:.c=.o)
DEPS=$(SOURCES:.c=.h)

TARGET=main

all: $(TARGET)

main: $(MAIN) $(OBJS)
	$(CC) $(CFLAGS) $(MAIN) $(OBJS) -o $(EXEC) $(CLIBS) $(SDL_LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@ $(CLIBS) $(SDL_LIBS)

install:
	$(INSTALL) $(EXEC) $(BINDIR)
