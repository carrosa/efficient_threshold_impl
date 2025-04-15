# 1. Create a python file xxx_params.py
# Do this lazer/python dir
# 2. Generate a header file xxx_params.h (run the following from lazer/scripts)
sage lin-codegen.sage ../python/xxx_params.py > ../python/xxx_params-h
# Generate C code + create a python module (Run this from lazer/python)
python3 params_cffi_build.py xxx_params.h
# Run python code
python3 xxx.py
