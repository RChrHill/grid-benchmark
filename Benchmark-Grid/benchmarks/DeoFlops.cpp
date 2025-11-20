#include "DeoFlops.hpp"


template<>
void update_representation<DWFSp4TwoIndexASInfo::Action>(
  typename DWFSp4TwoIndexASInfo::RepresentedGaugeField &Uas,
  const typename DWFSp4TwoIndexASInfo::FundamentalGaugeField &Uin)
{
  using namespace Grid;

  constexpr int Nc = DWFSp4TwoIndexASInfo::Nc;
  std::cout << GridLogDebug << "Updating TwoIndex representation\n";
  // Uas is in the TwoIndex antisymmetric representation
  // Uin is in the fundamental representation
  // get the U in TwoIndexRep
  // (U)_{(ij)(lk)} = tr [ adj(e^(ij)) U e^(lk) transpose(U) ]
  conformable(Uas, Uin);
  Uas = Zero();

  GeneralLinkField<vComplexD, Nc> tmp(Uin.Grid());

  const int Dimension = GaugeGroupTwoIndex<Nc, AntiSymmetric, GroupName::Sp>::Dimension;
  std::vector<typename GaugeGroup<Nc, GroupName::Sp>::MatrixD> eij(Dimension);

  for (int a = 0; a < Dimension; a++)
    GaugeGroupTwoIndex<Nc, AntiSymmetric, GroupName::Sp>::base(a, eij[a]);

  for (int mu = 0; mu < Nd; mu++) {
    auto Uin_mu = peekLorentz(Uin, mu);
    auto U_mu = peekLorentz(Uas, mu);
    for (int a = 0; a < Dimension; a++) {
      tmp = transpose(Uin_mu) * adj(eij[a]) * Uin_mu;
      for (int b = 0; b < Dimension; b++) {
        pokeColour(U_mu, trace(tmp * eij[b]), a, b);
      }
    }
    pokeLorentz(Uas, U_mu, mu);
  }
}
