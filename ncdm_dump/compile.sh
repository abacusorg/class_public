cd ..
make
cd ncdm_dump
g++ -o dump_ncdm -I../include -L../build -L../build/lib ../build/growTable.o ../build/dei_rkck.o ../build/sparse.o ../build/evolver_rkck.o ../build/evolver_ndf15.o ../build/arrays.o ../build/parser.o ../build/quadrature.o ../build/hyperspherical.o ../build/common.o ../build/input.o ../build/background.o ../build/thermodynamics.o ../build/perturbations.o ../build/primordial.o ../build/nonlinear.o ../build/transfer.o ../build/spectra.o ../build/lensing.o -fopenmp -g -fPIC -O4 -ffast-math -I../hyrec ../build/hyrectools.o ../build/helium.o ../build/hydrogen.o ../build/history.o ../build/output.o -DHYREC dump_ncdm_perturb.cpp -I/usr/include/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_cpp -lhdf5
