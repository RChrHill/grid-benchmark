#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/implementation/WilsonKernelsImplementation.h>
#include <Grid/qcd/action/fermion/implementation/WilsonKernelsHandImplementation.h>
#include <Grid/qcd/action/fermion/implementation/WilsonKernelsAsmImplementation.h>
#include "Sp4Fund.hpp"

NAMESPACE_BEGIN(Grid);

template class WilsonKernels<Sp4FundWilsonImplF>;

NAMESPACE_END(Grid);
