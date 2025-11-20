#pragma once

#include <Grid/Grid.h>
#include <Grid/qcd/action/fermion/WilsonImpl.h>


NAMESPACE_BEGIN(Grid);

typedef TwoIndexRep<4, TwoIndexSymmetry::AntiSymmetric, Grid::GroupName::Sp> TwoIndexAntiSymmetricRepresentationSp4;

typedef WilsonImpl<Grid::vComplexF, 
                   TwoIndexAntiSymmetricRepresentationSp4, 
                   Grid::CoeffReal> Sp4TwoIndexAntiSymmetricWilsonImplF;
typedef WilsonImpl<Grid::vComplexD, 
                   TwoIndexAntiSymmetricRepresentationSp4, 
                   Grid::CoeffReal> Sp4TwoIndexAntiSymmetricWilsonImplD;

typedef DomainWallFermion<Sp4TwoIndexAntiSymmetricWilsonImplF> DomainWallFermionSp4TwoIndexAntiSymmetricF;
typedef DomainWallFermion<Sp4TwoIndexAntiSymmetricWilsonImplD> DomainWallFermionSp4TwoIndexAntiSymmetricD;

NAMESPACE_END(Grid);
