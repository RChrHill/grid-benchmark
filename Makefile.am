AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I .buildutils/m4
 
bin_PROGRAMS =                \
  Benchmark_Grid              \
  Benchmark_IO

Benchmark_Grid_SOURCES              = Benchmark_Grid.cpp \
                                      Instantiation/Sp4Fund/Wilson.cpp \
                                      Instantiation/Sp4Fund/CayleyFermion5D.cc \
                                      Instantiation/Sp4Fund/WilsonFermion5DInstantiationSpWilsonImplD.cc \
                                      Instantiation/Sp4Fund/WilsonKernelsInstantiationSpWilsonImplD.cc \
                                      Instantiation/Sp4TwoIndexAntiSymmetric/Wilson.cpp \
                                      Instantiation/Sp4TwoIndexAntiSymmetric/CayleyFermion5DInstantiationSpWilsonTwoIndexAntiSymmetricImplD.cc \
                                      Instantiation/Sp4TwoIndexAntiSymmetric/WilsonFermion5DInstantiationSpWilsonTwoIndexAntiSymmetricImplD.cc \
                                      Instantiation/Sp4TwoIndexAntiSymmetric/WilsonKernelsInstantiationSpWilsonTwoIndexAntiSymmetricImplD.cc

Benchmark_IO_SOURCES                = Benchmark_IO.cpp
