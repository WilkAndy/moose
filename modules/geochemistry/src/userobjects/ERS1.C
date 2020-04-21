// TODO: when activity is fixed, remove from algebraic system, but never swap out of basis, then
// calculate like a mineral
// TODO: option for enforcing charge neutrality at every step
// TODO: presumably, fixing pH means that m_H is replaced with 10^(-pH)/a_H
//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ERS1.h"
#include "GeochemistryConstants.h"
#include "GeochemistryIonicStrength.h"
#include "GeochemistryActivity.h"

registerMooseObject("GeochemistryApp", ERS1);

defineLegacyParams(ERS1);

InputParameters
ERS1::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addRequiredParam<UserObjectName>(
      "model_definition",
      "The name of the GeochemicalModelDefinition user object.  Only equilibrium reactions are "
      "solved by ERS1, so the model_definition must not contain any kinetic "
      "species");
  params.addParam<std::vector<std::string>>(
      "swap_out_of_basis",
      "Species that should be removed from the model_definition's basis and be replaced with the "
      "swap_into_basis species.  There must be the same number of species in swap_out_of_basis and "
      "swap_into_basis.  If this list contains more than one species, the swapping is performed "
      "one-by-one, starting with the first pair (swap_out_of_basis[0] and swap_into_basis[0]), "
      "then the next pair, etc");
  params.addParam<std::vector<std::string>>(
      "swap_into_basis",
      "Species that should be removed from the model_definition's "
      "equilibrium species list and added to the basis");
  MultiMooseEnum constraint_meaning("moles_bulk_water kg_solvent_water moles_bulk_species "
                                    "free_molality free_moles_mineral_species fugacity activity");
  params.addRequiredParam<MultiMooseEnum>(
      "constraint_meaning",
      constraint_meaning,
      "Meanings of the numerical values given in constrain_value");
  params.addRequiredParam<std::vector<std::string>>(
      "constraint_species",
      "Names of the species that have their values fixed to constraint_value with meaning "
      "constraint_meaning.  All basis species (after swap_into_basis and swap_out_of_basis) must "
      "be provided with exactly one constraint");
  params.addRequiredParam<std::vector<Real>>(
      "constraint_value", "Numerical value of the containts on constraint_species");
  params.addParam<unsigned int>(
      "precision",
      4,
      "Precision for printing values.  Also, if the absolute value of a stoichiometric coefficient "
      "is less than 10^(-precision) then it is set to zero");
  params.addParam<std::string>("output_species",
                               "",
                               "Only output results for this species.  If not provided, results "
                               "for all species will be outputted");
  params.addRangeCheckedParam<Real>(
      "max_ionic_strength", 3.0, "max_ionic_strength >= 0.0", "Maximum value of ionic strength");
  params.addParam<unsigned>("extra_iterations_to_make_consistent",
                            0,
                            "Extra iterations to make the molalities, activities, etc, consistent "
                            "before commencing Newton iterations");
  params.addParam<Real>("temperature", 25, "Temperature of the aqueous solution");
  params.addRangeCheckedParam<Real>(
      "stoichiometry_tolerance",
      1E-6,
      "stoichiometry_tolerance >= 0.0",
      "Swapping involves inverting matrices via a singular value decomposition. During this "
      "process: (1) if abs(singular value) < stoi_tol * L1norm(singular values), then the "
      "matrix is deemed singular (so the basis swap is deemed invalid); (2) if abs(any "
      "stoichiometric coefficient) < stoi_tol then it is set to zero.");
  params.addRequiredParam<std::string>(
      "charge_balance_species",
      "Charge balance will be enforced on this basis species.  After swaps have been performed, "
      "this must be in the basis");
  params.addParam<std::vector<std::string>>(
      "prevent_precipitation",
      "Mineral species in this list will be prevented from precipitating, irrespective of their "
      "saturation index, unless they are in the basis");
  params.addParam<Real>("abs_tol",
                        1E-10,
                        "If the residual of the algebraic system (measured in mol) is lower than "
                        "this value, the Newton process is deemed to have converged");
  params.addParam<Real>("rel_tol",
                        1E-200,
                        "If the residual of the algebraic system is lower than this value times "
                        "the initial residual, the Newton process is deemed to have converged");
  params.addRangeCheckedParam<Real>(
      "min_initial_molality",
      1E-20,
      "min_initial_molality > 0.0",
      "Minimum value of the initial-guess molality used in the Newton process");
  params.addParam<unsigned>(
      "max_iter",
      100,
      "Maximum number of iterations allowed to solve one round of the algebraic system");
  params.addParam<bool>("verbose", false, "Print verbose information");
  params.addParam<Real>(
      "max_initial_residual",
      1E3,
      "Attempt to alter the initial-guess molalities so that the initial residual "
      "for the Newton process is less than this number of moles");
  params.addRangeCheckedParam<Real>(
      "swap_threshold",
      0.1,
      "swap_threshold >= 0.0",
      "If the molality of a basis species in the algebraic system falls below swap_threshold * "
      "abs_tol then it is swapped out of the basis.  The dimensions of swap_threshold are "
      "1/kg(solvent water)");
  params.addParam<unsigned>("ramp_max_ionic_strength",
                            20,
                            "The number of iterations over which to progressively increase the "
                            "maximum ionic strength (from zero to max_ionic_strength).  Increasing "
                            "this can help in convergence of the Newton process");
  params.addClassDescription("User object for solving geochemical reaction systems");

  return params;
}

ERS1::ERS1(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _egs(getUserObject<GeochemicalModelDefinition>("model_definition").getDatabase(),
         getParam<Real>("stoichiometry_tolerance"),
         getParam<std::vector<std::string>>("swap_out_of_basis"),
         getParam<std::vector<std::string>>("swap_into_basis"),
         getParam<std::string>("charge_balance_species"),
         getParam<std::vector<std::string>>("constraint_species"),
         getParam<std::vector<Real>>("constraint_value"),
         getParam<MultiMooseEnum>("constraint_meaning"),
         getParam<Real>("temperature"),
         getParam<Real>("max_ionic_strength") /
             (1.0 + getParam<unsigned>("ramp_max_ionic_strength")),
         getParam<unsigned>("extra_iterations_to_make_consistent"),
         getParam<Real>("min_initial_molality")),
    _num_basis(_egs.getNumInBasis()),
    _num_eqm(_egs.getNumInEquilibrium()),
    _mgd(_egs.getModelGeochemicalDatabase()),
    _num_in_algebraic_system(_egs.getNumInAlgebraicSystem()),
    _residual(_num_in_algebraic_system),
    _abs_residual(0.0),
    _jacobian(_num_in_algebraic_system, _num_in_algebraic_system),
    _new_mol(_num_in_algebraic_system),
    _precision(getParam<unsigned int>("precision")),
    _output_species(getParam<std::string>("output_species")),
    _prevent_precipitation(getParam<std::vector<std::string>>("prevent_precipitation")),
    _abs_tol(getParam<Real>("abs_tol")),
    _rel_tol(getParam<Real>("rel_tol")),
    _res0_times_rel(0.0),
    _max_iter(getParam<unsigned>("max_iter")),
    _verbose(getParam<bool>("verbose")),
    _max_initial_residual(getParam<Real>("max_initial_residual")),
    _swap_threshold(getParam<Real>("swap_threshold") * _abs_tol),
    _max_ionic_strength(getParam<Real>("max_ionic_strength")),
    _ramp_max_ionic_strength(getParam<unsigned>("ramp_max_ionic_strength")),
    _tot_iter(0)
{
  if (_ramp_max_ionic_strength > _max_iter)
    paramError("ramp_max_ionic_strength", "must be less than max_iter");
}

void
ERS1::initialize()
{
}

void
ERS1::execute()
{
  bool still_swapping_minerals = true;
  while (still_swapping_minerals)
  {
    _num_in_algebraic_system = _egs.getNumInAlgebraicSystem();
    _residual = DenseVector<Real>(_num_in_algebraic_system);
    _jacobian = DenseMatrix<Real>(_num_in_algebraic_system, _num_in_algebraic_system);
    _new_mol = DenseVector<Real>(_num_in_algebraic_system);

    unsigned iter = 0;
    _abs_residual = computeResidual(_residual);
    bool reducing_initial_molalities = (_abs_residual > _max_initial_residual);
    while (reducing_initial_molalities)
      reducing_initial_molalities = reduceInitialResidual();

    if (_verbose)
      _console << std::setprecision(_precision) << "iter = " << iter << " |R| = " << _abs_residual
               << std::endl;
    _res0_times_rel = _abs_residual * _rel_tol;
    while ((_abs_residual >= _res0_times_rel && _abs_residual >= _abs_tol && iter < _max_iter) ||
           iter < _ramp_max_ionic_strength)
    {
      iter += 1;
      _tot_iter += 1;
      _egs.computeJacobian(_residual, _jacobian);
      solveAndUnderrelax(_jacobian, _new_mol);
      _egs.setMaxIonicStrength(std::min(1.0, (iter + 1.0) / (_ramp_max_ionic_strength + 1.0)) *
                               _max_ionic_strength);
      _egs.setAlgebraicVariables(_new_mol);
      if (_egs.alterChargeBalanceSpecies(_swap_threshold) && _verbose)
        _console << "Changed change balance species to "
                 << _mgd.basis_species_name[_egs.getChargeBalanceBasisIndex()] << std::endl;
      _abs_residual = computeResidual(_residual);
      if (_verbose)
        _console << "iter = " << iter << " |R| = " << _abs_residual << std::endl;
    }

    if (iter >= _max_iter)
      mooseWarning("Number of iterations exceeds ", _max_iter, "\n");

    unsigned swap_out_of_basis = 0;
    unsigned swap_into_basis = 0;
    try
    {
      still_swapping_minerals = mineralSwapNeeded(swap_out_of_basis, swap_into_basis);
    }
    catch (const MooseException & e)
    {
      mooseError(e.what());
    }
    if (still_swapping_minerals)
    {
      // need to do a swap and re-solve
      try
      {
        _egs.performSwap(swap_out_of_basis, swap_into_basis);
      }
      catch (const MooseException & e)
      {
        mooseError(e.what());
      }
    }
  }
}

void
ERS1::finalize()
{
}

std::vector<std::string>
ERS1::speciesOfInterest() const
{
  return {""};
  // TODO
}

Real
ERS1::computeResidual(DenseVector<Real> & residual) const
{
  for (unsigned a_ind = 0; a_ind < _num_in_algebraic_system; ++a_ind)
    residual(a_ind) = _egs.getResidualComponent(a_ind);
  return residual.l1_norm();
}

void
ERS1::solveAndUnderrelax(DenseMatrix<Real> & jacobian, DenseVector<Real> & new_mol) const
{
  jacobian.lu_solve(_residual, new_mol);

  // at this point we want to do molality = molality - new_mol

  // Bethke recommends underrelaxation - probably want to do PETSc variational bounds in the
  // future
  Real one_over_delta = 1.0;
  const std::vector<Real> current_molality = _egs.getAlgebraicVariableValues();
  for (unsigned a_ind = 0; a_ind < _num_in_algebraic_system; ++a_ind)
    one_over_delta = std::max(one_over_delta, new_mol(a_ind) * 2.0 / current_molality[a_ind]);
  for (unsigned a_ind = 0; a_ind < _num_in_algebraic_system; ++a_ind)
    new_mol(a_ind) = current_molality[a_ind] - new_mol(a_ind) / one_over_delta;
}

bool
ERS1::mineralSwapNeeded(unsigned & swap_out_of_basis, unsigned & swap_into_basis) const
{
  bool swap_needed = false;
  unsigned ind = 0;

  // check if any basis minerals have negative free number of moles
  const std::vector<Real> basis_molality = _egs.getSolventMassAndFreeMolalityAndMineralMoles();
  std::vector<unsigned> molality_order(_num_basis);
  ind = 0;
  std::iota(molality_order.begin(), molality_order.end(), ind++);
  std::sort(molality_order.begin(), molality_order.end(), [&](int i, int j) {
    return basis_molality[i] < basis_molality[j];
  });

  // if the Newton did not converge then check for small molalities in the non-minerals
  if (_abs_residual > _res0_times_rel && _abs_residual > _abs_tol)
  {
    for (const auto & i : molality_order)
    {
      if (basis_molality[i] >= _swap_threshold)
      {
        // since we're going through basis_molality in ascending order, as soon as we get >=
        // _swap_threshold, we're safe
        break;
      }
      else if (!_mgd.basis_species_mineral[i] && _egs.getInAlgebraicSystem()[i] &&
               _egs.getChargeBalanceBasisIndex() != i)
      {
        // a non-mineral in the algebraic system has super low molality: try to find a legitimate
        // swap
        swap_out_of_basis = i;
        swap_into_basis = 0;

        bool legitimate_swap_found = false;
        Real best_stoi = 0.0;
        // the first non-gas is the best possible so far
        for (unsigned j = 0; j < _num_eqm; ++j)
        {
          if (_mgd.eqm_species_gas[j] || _mgd.eqm_stoichiometry(j, i) == 0.0)
            continue;
          best_stoi = std::abs(_mgd.eqm_stoichiometry(j, i)) * _egs.getEquilibriumMolality(j);
          swap_into_basis = j;
          legitimate_swap_found = true;
        }
        // now go through the remainder, trying to find a better swap
        for (unsigned j = swap_into_basis; j < _num_eqm; ++j)
        {
          if (_mgd.eqm_species_gas[j] || _mgd.eqm_stoichiometry(j, i) == 0.0)
            continue;
          const Real stoi = std::abs(_mgd.eqm_stoichiometry(j, i)) * _egs.getEquilibriumMolality(j);
          if (stoi > best_stoi)
          {
            best_stoi = stoi;
            swap_into_basis = j;
          }
        }
        if (legitimate_swap_found)
        {
          if (_verbose)
            _console << "Basis species " << _mgd.basis_species_name[swap_out_of_basis]
                     << " has very low molality of " << basis_molality[i] << " compared to "
                     << _swap_threshold << ".  Swapping it with equilibrium species "
                     << _mgd.eqm_species_name[swap_into_basis] << std::endl;
          swap_needed = true;
          break;
        }
        // if no legitimate swap is found, then loop around to the next basis species
      }
    }
  }
  if (swap_needed)
    return swap_needed;

  // now look through the molalities for minerals that are consumed
  for (const auto & i : molality_order)
  {
    if (basis_molality[i] >= 0.0)
    {
      // since we're going through basis_molality in ascending order, as soon as we get >=0, we're
      // safe
      break;
    }
    else if (_mgd.basis_species_mineral[i])
    {
      swap_needed = true;
      swap_out_of_basis = i;
      swap_into_basis = 0;

      bool legitimate_swap_found = false;
      Real best_stoi = 0.0;
      // the first non-mineral and non-gas is the best possible so far
      for (unsigned j = 0; j < _num_eqm; ++j)
      {
        if (_mgd.eqm_species_mineral[j] || _mgd.eqm_species_gas[j] ||
            _mgd.eqm_stoichiometry(j, i) == 0.0)
          continue;
        best_stoi = std::abs(_mgd.eqm_stoichiometry(j, i)) * _egs.getEquilibriumMolality(j);
        swap_into_basis = j;
        legitimate_swap_found = true;
      }
      // now go through the remainder, trying to find a better swap
      for (unsigned j = swap_into_basis; j < _num_eqm; ++j)
      {
        if (_mgd.eqm_species_mineral[j] || _mgd.eqm_species_gas[j] ||
            _mgd.eqm_stoichiometry(j, i) == 0.0)
          continue;
        const Real stoi = std::abs(_mgd.eqm_stoichiometry(j, i)) * _egs.getEquilibriumMolality(j);
        if (stoi > best_stoi)
        {
          best_stoi = stoi;
          swap_into_basis = j;
        }
      }
      if (!legitimate_swap_found)
        mooseException("Cannot find a legitimate swap for mineral ",
                       _mgd.basis_species_name[swap_out_of_basis]);
      if (_verbose)
        _console << "Mineral " << _mgd.basis_species_name[swap_out_of_basis]
                 << " consumed.  Swapping it with equilibrium species "
                 << _mgd.eqm_species_name[swap_into_basis] << std::endl;
      break;
    }
  }
  if (swap_needed)
    return swap_needed;

  // check maximum saturation index is not positive
  const std::vector<Real> eqm_SI = _egs.getSaturationIndices();
  std::vector<unsigned> mineral_SI_order(_num_eqm);
  ind = 0;
  std::iota(mineral_SI_order.begin(), mineral_SI_order.end(), ind++);
  std::sort(mineral_SI_order.begin(), mineral_SI_order.end(), [&](int i, int j) {
    return eqm_SI[i] > eqm_SI[j];
  });

  for (const auto & j : mineral_SI_order)
    if (eqm_SI[j] > 0.0)
      if (_mgd.eqm_species_mineral[j])
        if (std::find(_prevent_precipitation.begin(),
                      _prevent_precipitation.end(),
                      _mgd.eqm_species_name[j]) == _prevent_precipitation.end())
        {
          // mineral has positive saturation index and user is not preventing its precipitation.
          // determine the basis species to swap out
          swap_needed = true;
          swap_into_basis = j;
          bool legitimate_swap_found = false;
          swap_out_of_basis = 0;
          Real best_stoi = 0.0;
          for (unsigned i = 1; i < _num_basis; ++i) // never swap water (i=0)
          {
            if (basis_molality[i] > 0.0 && i != _egs.getChargeBalanceBasisIndex() &&
                !_mgd.basis_species_gas[i] && !_egs.getBasisActivityKnown()[i])
            {
              // don't want to swap out the charge-balance species or any gases of fixed fugacity
              // or any species with fixed activity
              const Real stoi = std::abs(_mgd.eqm_stoichiometry(j, i)) / basis_molality[i];
              if (stoi > best_stoi)
              {
                legitimate_swap_found = true;
                best_stoi = stoi;
                swap_out_of_basis = i;
              }
            }
          }
          if (!legitimate_swap_found)
            mooseException(
                "Cannot find a legitimate swap for the supersaturated equilibrium species ",
                _mgd.eqm_species_name[j]);
          if (_verbose)
            _console << "Mineral " << _mgd.eqm_species_name[j]
                     << " supersaturated.  Swapping it with basis species "
                     << _mgd.basis_species_name[swap_out_of_basis] << std::endl;
          break;
        }
  return swap_needed;
}

bool
ERS1::reduceInitialResidual()
{
  const Real initial_r = _abs_residual;
  const std::vector<Real> original_molality = _egs.getAlgebraicVariableValues();
  unsigned ind = 0;

  // to get an indication of whether we should increase or decrease molalities in the algorithm
  // below, find the median of the original molalities
  ind = 0;
  std::vector<unsigned> mol_order(_num_in_algebraic_system);
  std::iota(mol_order.begin(), mol_order.end(), ind++);
  std::sort(mol_order.begin(), mol_order.end(), [&](int i, int j) {
    return (original_molality[i] > original_molality[j]);
  });
  const Real median_molality = original_molality[mol_order[_num_in_algebraic_system / 2]];

  // get the index order of the residual vector (largest first)
  ind = 0;
  std::vector<unsigned> res_order(_num_in_algebraic_system);
  std::iota(res_order.begin(), res_order.end(), ind++);
  std::sort(res_order.begin(), res_order.end(), [&](int i, int j) {
    return std::abs(_residual(i)) > std::abs(_residual(j));
  });

  DenseVector<Real> new_molality(original_molality);
  for (const auto & a : res_order)
  {
    if (std::abs(_residual(a)) < _max_initial_residual)
      return false; // haven't managed to find a suitable new molality, and all remaining residual
                    // components are less than the cutoff, so cannot appreciably reduce from now
                    // on

    const Real multiplier = (original_molality[a] > median_molality) ? 0.5 : 2.0;
    // try using the multiplier
    new_molality(a) = multiplier * original_molality[a];
    _egs.setAlgebraicVariables(new_molality);
    _abs_residual = computeResidual(_residual);
    if (_abs_residual < initial_r)
      return true; // success: found a new molality that decreases the initial |R|

    // the above approach did not decrease |R|, so try using 1/multiplier
    new_molality(a) = original_molality[a] / multiplier;
    _egs.setAlgebraicVariables(new_molality);
    _abs_residual = computeResidual(_residual);
    if (_abs_residual < initial_r)
      return true; // success: found a new molality that decreases the initial |R|

    // the new molalities did not decrease |R|, so revert to the original molality, and move to
    // next-biggest residual component
    new_molality(a) = original_molality[a];
    _egs.setAlgebraicVariables(new_molality);
    _abs_residual = computeResidual(_residual);
  }
  return false;
}
