#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/implementation/WilsonFermion5DImplementation.h>
#include "Sp4TwoIndexAntiSymmetric.hpp"

NAMESPACE_BEGIN(Grid);

template class WilsonFermion5D<Sp4TwoIndexAntiSymmetricWilsonImplD>;

NAMESPACE_END(Grid);
