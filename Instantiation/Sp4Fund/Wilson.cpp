#include <Grid/Grid.h>

#include <Grid/qcd/action/fermion/implementation/WilsonFermionImplementation.h>
#include "../KernelImplementation.hpp"
#include "Implementation.hpp"

NAMESPACE_BEGIN(Grid);

template class WilsonFermion<Sp4FundWilsonImplF>;
template class WilsonFermion<Sp4FundWilsonImplD>;

NAMESPACE_END(Grid);