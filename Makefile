CFLAGS=-O3
DEBUGFLAGS=-Wall -ggdb3
INCLUDES=-I/opt/X11/include
LDFLAGS=-L/opt/X11/lib -lX11 -lm

galsim: graphics.o galsim.c
	gcc $(DEBUGFLAGS) -o galsim galsim.c graphics.o $(LDFLAGS)

galsim_final: graphics.o galsim.c
	gcc $(CFLAGS) -o galsim galsim.c graphics.o $(LDFLAGS)

graphics.o: graphics/graphics.c graphics/graphics.h
	gcc $(CFLAGS) $(INCLUDES) -c graphics.c

clean:
	rm -f ./galsim graphics/*.o
