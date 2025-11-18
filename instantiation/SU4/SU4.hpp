#pragma once

#include <Grid/Grid.h>


NAMESPACE_BEGIN(Grid);

typedef Grid::FundamentalRep<4,Grid::GroupName::SU> FundamentalRepresentationSU4;

typedef WilsonImpl<Grid::vComplexF, 
                   FundamentalRepresentationSU4,
                   Grid::CoeffReal> SU4FundWilsonImplF;
typedef WilsonImpl<Grid::vComplexD, 
                   FundamentalRepresentationSU4,
                   Grid::CoeffReal> SU4FundWilsonImplD;

typedef DomainWallFermion<SU4FundWilsonImplF> DomainWallFermionSU4F;
typedef DomainWallFermion<SU4FundWilsonImplD> DomainWallFermionSU4D;

NAMESPACE_END(Grid);
