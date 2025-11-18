#pragma once

#include <Grid/Grid.h>


NAMESPACE_BEGIN(Grid);

typedef Grid::FundamentalRep<4,Grid::GroupName::Sp> FundamentalRepresentationSp4;

typedef WilsonImpl<Grid::vComplexF, 
                   FundamentalRepresentationSp4, 
                   Grid::CoeffReal> Sp4FundWilsonImplF;
typedef WilsonImpl<Grid::vComplexD, 
                   FundamentalRepresentationSp4, 
                   Grid::CoeffReal> Sp4FundWilsonImplD;

typedef DomainWallFermion<Sp4FundWilsonImplD> DomainWallFermionSp4F;
typedef DomainWallFermion<Sp4FundWilsonImplD> DomainWallFermionSp4D;

NAMESPACE_END(Grid);
