#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/implementation/WilsonKernelsImplementation.h>
#include <Grid/qcd/action/fermion/implementation/WilsonKernelsHandImplementation.h>
#include <Grid/qcd/action/fermion/implementation/WilsonKernelsAsmImplementation.h>
#include "Sp4TwoIndexAntiSymmetric.hpp"

NAMESPACE_BEGIN(Grid);

template class WilsonKernels<Sp4TwoIndexAntiSymmetricWilsonImplD>;

NAMESPACE_END(Grid);
