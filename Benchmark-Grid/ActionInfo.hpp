#pragma once

#include <Grid/Grid.h>


/* Grid's template types for actions and representations do not expose template parameters and attributes that we require.
   Here we forcibly expose them via the template system. */

/* Type aliases because Grid's built-in ones are hard-coded to the configure-time Nc */
template<typename Simd, int dimension>
using GeneralGaugeField = Grid::Lattice<Grid::iVector<Grid::iScalar<Grid::iMatrix<Simd, dimension>>, Grid::Nd>>;

template<typename Simd, int dimension>
using GeneralLinkField = Grid::Lattice<Grid::iScalar<Grid::iScalar<Grid::iMatrix<Simd, dimension>>>>;


// ## getRepresentation #########################################
template<typename Impl>
struct getRepresentation
{
  static_assert(!std::is_same_v<Impl,Impl>, "getRepresentation undefined for given type");
};

template<typename S, typename Representation, typename Options>
struct getRepresentation<Grid::WilsonImpl<S, Representation, Options>>
{
  typedef Representation type;
};

template<typename S, typename Representation>
struct getRepresentation<Grid::StaggeredImpl<S, Representation>>
{
  typedef Representation type;
};


// ## getActionRepresentation ###################################
template<typename Action>
struct getActionRepresentation
{
  typedef typename getRepresentation<typename Action::Impl_t>::type type;
};


// ## groupNameStr ##############################################
template<typename groupName_>
struct groupNameStr {
  static_assert(!std::is_same_v<groupName_, groupName_>, "groupNameStr undefined for given type");
};

template<>
struct groupNameStr<Grid::GroupName::SU> {
  static constexpr const char* value = "SU";
};

template<>
struct groupNameStr<Grid::GroupName::Sp> {
  static constexpr const char* value = "Sp";
};


// ## getGroupInfo ##############################################
template<typename Representation>
struct getGroupInfo
{
  static_assert(!std::is_same_v<Representation, Representation>, "getGroupInfo undefined for given type");
};

template<int Nc_, typename groupName_>
struct getGroupInfo<Grid::FundamentalRep<Nc_, groupName_>>
{
  static constexpr int Nc = Nc_;
  typedef groupName_ label;
  typedef Grid::GaugeGroup<Nc_, groupName_> group;
  static constexpr const char* getRepresentationName()
  {
    return "Fundamental";
  }
  static constexpr const char* groupName = groupNameStr<groupName_>::value;
};

template<int Nc_, Grid::TwoIndexSymmetry S_, typename groupName_>
struct getGroupInfo<Grid::TwoIndexRep<Nc_, S_, groupName_>>
{
  static constexpr int Nc = Nc_;
  typedef groupName_ label;
  typedef Grid::GaugeGroupTwoIndex<Nc_, S_, groupName_> group;
  static constexpr const char* getRepresentationName()
  {
    if constexpr (S_ == Grid::Symmetric)
    {
	    return "Two-Index Symmetric";
    }
    else if constexpr (S_ == Grid::AntiSymmetric)
    {
	    return "Two-Index Antisymmetric";
    }
    else
    {
      static_assert(!std::is_same_v<groupName_, groupName_>, "Unrecognised TwoIndexSymmetry");
    }
  }
  static constexpr const char* groupName = groupNameStr<groupName_>::value;
};


// ## ActionFactory #############################################
#define DEFINE_ACTION_FACTORY_TYPES(ActionType)\
typedef ActionType Action;\
typedef typename getActionRepresentation<Action>::type Representation;\
typedef getGroupInfo<Representation> GroupInfo;\
static constexpr int Nc = GroupInfo::Nc;\
static constexpr int dim_rep = Action::Dimension;\
typedef GeneralGaugeField<typename Action::Simd, Action::Dimension> RepresentedGaugeField;\
typedef GeneralGaugeField<typename Action::Simd, GroupInfo::Nc> FundamentalGaugeField;

template<typename ActionType>
struct ActionFactory
{
  DEFINE_ACTION_FACTORY_TYPES(ActionType)
  
  static auto create(RepresentedGaugeField& Umu, Grid::GridCartesian* Grid4D, Grid::GridRedBlackCartesian* RbGrid4D, Grid::GridCartesian* Grid5D, Grid::GridRedBlackCartesian* RbGrid5D)
  {
    using namespace Grid;
    static_assert(!std::is_same_v<Action,Action>, "create is not defined for provided Action");
  }
};

template<typename Impl>
struct ActionFactory<Grid::DomainWallFermion<Impl>>
{
  DEFINE_ACTION_FACTORY_TYPES(Grid::DomainWallFermion<Impl>)
  
  static auto create(RepresentedGaugeField& Umu, Grid::GridCartesian* Grid4D, Grid::GridRedBlackCartesian* RbGrid4D, Grid::GridCartesian* Grid5D, Grid::GridRedBlackCartesian* RbGrid5D)
  {
    using namespace Grid;
    RealD mass = 0.1;
    RealD M5 = 1.8;
    return Action(Umu, *Grid5D, *RbGrid5D, *Grid4D, *RbGrid4D, mass, M5);
  }
  static double fps()
  {
    using namespace Grid;
    // Nc=3 gives
    // 1344= 3*(2*8+6)*2*8 + 8*3*2*2 + 3*4*2*8
    // 1344 = Nc* (6+(Nc-1)*8)*2*Nd + Nd*Nc*2*2  + Nd*Nc*Ns*2
    //	double flops=(1344.0*volume)/2;
    #if 0
      double fps = Nc* (6+(Nc-1)*8)*Ns*Nd + Nd*Nc*Ns  + Nd*Nc*Ns*2;
    #else
      return  dim_rep * (6 + (dim_rep - 1) * 8) * Ns * Nd + 2 * Nd * dim_rep * Ns + 2 * Nd * dim_rep * Ns * 2;
    #endif
  }
  static std::string name() { return "DWF"; }
};

template<typename Impl>
struct ActionFactory<Grid::ImprovedStaggeredFermion<Impl>>
{
  DEFINE_ACTION_FACTORY_TYPES(Grid::ImprovedStaggeredFermion<Impl>)

  static auto create(RepresentedGaugeField& Umu, Grid::GridCartesian* Grid4D, Grid::GridRedBlackCartesian* RbGrid4D, Grid::GridCartesian* Grid5D, Grid::GridRedBlackCartesian* RbGrid5D)
  {
    using namespace Grid;
    RealD mass = 0.1;
    RealD c1 = 9.0 / 8.0;
    RealD c2 = -1.0 / 24.0;
    RealD u0 = 1.0;
    typename Action::ImplParams params;
    return Action(Umu, Umu, *Grid4D, *RbGrid4D, mass, c1, c2, u0, params);
  }
  static double fps()
  {
    using namespace Grid;
    if constexpr (Nc != 3)
    {
      static_assert(!std::is_same_v<Action,Action>, "Nc!=3 is not supported for ImprovedStaggered");
    }
    return 1146.0;
  }
  static std::string name() { return "ImprovedStaggered"; }
};

#undef DEFINE_ACTION_FACTORY_TYPES


// ## actionPrec ################################################
template<typename Action>
static std::string actionPrec()
{
  using namespace Grid;
  typedef typename Action::Simd::Real Real_t;
  if constexpr (std::is_same_v<Real_t, float>)
  {
    return "SINGLE";
  }
  else if constexpr (std::is_same_v<Real_t, double>)
  {
    return "DOUBLE";
  }
  else
  {
    static_assert(!std::is_same_v<Action,Action>, "Unknown precision for provided action");
  }
}
