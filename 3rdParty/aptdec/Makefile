CC = g++
BIN = /usr/bin
INCLUDES = -I.
CFLAGS = -g -O0 -ggdb $(INCLUDES)
OBJS = main.o image.o dsp.o filter.o reg.o fcolor.o

aptdec:	$(OBJS)
	$(CC) -o $@ $(OBJS) -lm -lsndfile -lpng

main.o:	  main.c palette.h offsets.h messages.h
dsp.o:	  dsp.c filtercoeff.h filter.h
filter.o: filter.c filter.h
image.o:  image.c satcal.h offsets.h messages.h
fcolor.o: fcolor.c offsets.h

clean:
	rm -f *.o aptdec

install:
	install -m 755 aptdec $(BIN)

uninstall:
	rm $(BIN)/aptdec
