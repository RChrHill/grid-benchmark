#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/implementation/WilsonFermion5DImplementation.h>
#include "Sp4Fund.hpp"

NAMESPACE_BEGIN(Grid);

template class WilsonFermion5D<Sp4FundWilsonImplF>;

NAMESPACE_END(Grid);
