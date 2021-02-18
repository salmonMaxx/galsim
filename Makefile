# CFLAGS=-O3 	# improvement, speedup of ~3
CFLAGS2=-Ofast # pretty good improvement!
CFLAGS=-O3 -march=native -ffast-math # A little faster.
DEBUGFLAGS=-Wall -ggdb3
INCLUDES=-I/opt/X11/include
LDFLAGS=-L/opt/X11/lib -lX11 -lm

galsim: graphics.o galsim.c
	gcc $(CFLAGS) -o galsim galsim.c graphics.o $(LDFLAGS)

debug: graphics.o galsim.c
	gcc $(DEBUGFLAGS) -o galsim galsim.c graphics.o $(LDFLAGS)

galsim_2: graphics.o galsim.c
	gcc $(CFLAGS2) -o galsim_2 galsim.c graphics.o $(LDFLAGS)

graphics.o: graphics/graphics.c graphics/graphics.h
	gcc $(CFLAGS) $(INCLUDES) -c graphics/graphics.c

clean:
	rm -f galsim galsim_2 *.o graphics/*.o *.gal
