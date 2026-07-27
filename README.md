This is a Fortran code to generate length vs time trajectory for a single actin filament undergoing polymerization, severing and hydrolysis. Three microscopic monomer states are there, namely ADP-pi (undecorated, denoted by +1 value in the code), ADP (undecorated, denoted by 0 value in the code), and ADP-cof (decorated, denoted by -1 value in the code).
Details are annotated in the code, along with the parameter values.
To successfully run the code, you need gfortran or ifort as the compiler, installed on your machine.
Command for compilation & code running: gfortran -O3 code_file_name.f90; ./a.out
