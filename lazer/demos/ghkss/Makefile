.PHONY: all clean

LAZER_DIR = ../..

LIBS = $(LAZER_DIR)/liblazer.a -lmpfr -lgmp -lm ../../third_party/hexl-development/build/hexl/lib/libhexl.a -lstdc++ #XXX
CFLAGS = -Wall -Wextra -Wshadow -Wundef -O3 -g

all: ghkss

ghkss: params.h ghkss.c
	$(CC) $(CFLAGS) -I$(LAZER_DIR) -o $@ ghkss.c $(LIBS)

example1: params.h example1.c
	$(CC) $(CFLAGS) -I$(LAZER_DIR) -o $@ example1.c $(LIBS)

example2: params.h example2.c
	$(CC) $(CFLAGS) -I$(LAZER_DIR) -o $@ example2.c $(LIBS)

params.h: params.py
	cd ../../scripts && sage lin-codegen.sage ../demos/ghkss/params.py > ../demos/ghkss/params.h

clean:
	rm -f ghkss