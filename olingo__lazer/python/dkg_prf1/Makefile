.PHONY: all clean

all: _ctxr_params_cffi.o

_ctxr_params_cffi.o: ctxr_params.h ../_lazer_cffi.o
	python3 ../params_cffi_build.py ctxr_params.h ../..

ctxr_params.h: ctxr_params.py
	cd ../../scripts && sage lin-codegen.sage ../python/dkg_prf1/ctxr_params.py > ../python/dkg_prf1/ctxr_params.h

../_lazer_cffi.o:
	cd .. && make

clean:
	rm -rf __pycache__
	rm -f *.c *.o *_cffi.* ctxr_params.h
	find . -maxdepth 1 -name "*.so" ! -name "liblazer.so" -delete
