.PHONY: all clean

all: _dualcode_params_cffi.o

_dualcode_params_cffi.o: dualcode_params.h ../_lazer_cffi.o
	python3 ../params_cffi_build.py dualcode_params.h ../..

dualcode_params.h: dualcode_params.py
	cd ../../scripts && sage lin-codegen.sage ../python/dualcode/dualcode_params.py > ../python/dualcode/dualcode_params.h

../_lazer_cffi.o:
	cd .. && make

clean:
	rm -rf __pycache__
	rm -f *.c *.o *_cffi.* dualcode_params.h
	find . -maxdepth 1 -name "*.so" ! -name "liblazer.so" -delete
