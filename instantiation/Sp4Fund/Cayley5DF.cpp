#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/implementation/CayleyFermion5DImplementation.h>
#include <Grid/qcd/action/fermion/implementation/CayleyFermion5Dcache.h>
#include "Sp4Fund.hpp"

NAMESPACE_BEGIN(Grid);

template class CayleyFermion5D<Sp4FundWilsonImplF>;

NAMESPACE_END(Grid);
