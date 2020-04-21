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

#include "EquilibriumReactionSolver.h"
#include "GeochemistryConstants.h"
#include "GeochemistryIonicStrength.h"
#include "GeochemistryActivity.h"

constexpr Real mol_cutoff = 40.0; // molalities are bounded below by 10^(-mol_cutoff)
constexpr Real mol_multiplier =
    1E100; // eqm_molalities are never greater than 1000 * sum(|basis_molalities|)

registerMooseObject("GeochemistryApp", EquilibriumReactionSolver);

defineLegacyParams(EquilibriumReactionSolver);

InputParameters
EquilibriumReactionSolver::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addRequiredParam<UserObjectName>(
      "model_definition",
      "The name of the GeochemicalModelDefinition user object.  Only equilibrium reactions are "
      "solved by EquilibriumReactionSolver, so the model_definition must not contain any kinetic "
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
  params.addParam<std::vector<std::string>>(
      "activity_species",
      "Names of the species that have their values fixed to constraint_value with meaning "
      "constraint_meaning.  All basis species (after swap_into_basis and swap_out_of_basis) must "
      "be provided with exactly one constraint");
  params.addParam<std::vector<Real>>("activity_values",
                                     "Numerical value of the containts on constraint_species");
  params.addParam<unsigned int>(
      "precision",
      4,
      "Precision for printing values.  Also, if the absolute value of a stoichiometric coefficient "
      "is less than 10^(-precision) then it is set to zero");
  params.addParam<std::string>("output_species",
                               "",
                               "Only output results for this species.  If not provided, results "
                               "for all species will be outputted");
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
  params.addParam<std::vector<std::string>>("bulk_moles_species", "asdf");
  params.addParam<std::vector<std::string>>("free_molality_species", "asdf");
  params.addParam<std::vector<std::string>>("free_moles_mineral_species", "asdf");
  params.addParam<std::vector<Real>>("bulk_moles_values", "asdf");
  params.addParam<std::vector<Real>>("free_molality_values", "asdf");
  params.addParam<std::vector<Real>>("free_moles_mineral_values", "asdf");
  params.addParam<Real>("moles_bulk_water", "asdf");
  params.addParam<Real>("mass_solvent_water", "asdf");
  params.addParam<unsigned>(
      "max_iter",
      100,
      "Maximum number of iterations allowed to solve one round of the algebraic system");
  params.addParam<bool>("verbose", false, "Print verbose information");
  params.addClassDescription("User object for solving geochemical reaction systems");

  return params;
}

EquilibriumReactionSolver::EquilibriumReactionSolver(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _mgd(getUserObject<GeochemicalModelDefinition>("model_definition").getDatabase()),
    _num_basis(_mgd.basis_species_index.size()),
    _num_eqm(_mgd.eqm_species_index.size()),
    _swapper(_num_basis, getParam<Real>("stoichiometry_tolerance")),
    _swap_out(getParam<std::vector<std::string>>("swap_out_of_basis")),
    _swap_in(getParam<std::vector<std::string>>("swap_into_basis")),
    _num_in_algebraic_system(0),
    _in_algebraic_system(_num_basis),
    _algebraic_index(_num_basis),
    _basis_index(_num_basis), // typically do not need all the slots here: just the first
                              // _num_in_algebraic_system are filled with meaningful data
    _nw(0.0),
    _basis_molality(_num_basis),
    _bulk_moles(_num_basis),
    _eqm_molality(_num_eqm),
    _residual(0),
    _abs_residual(0.0),
    _jacobian(0),
    _dmj_dmi(_num_eqm, _num_basis),
    _dmol(0),
    _precision(getParam<unsigned int>("precision")),
    _temperature(getParam<Real>("temperature")),
    _activity_species(getParam<std::vector<std::string>>("activity_species")),
    _activity_values(getParam<std::vector<Real>>("activity_values")),
    _basis_activity_known(_num_basis),
    _basis_activity(_num_basis),
    _eqm_activity(_num_eqm),
    _eqm_log10K(_num_eqm),
    _eqm_SI(_num_eqm),
    _dhA(0.0),
    _dhB(0.0),
    _dhBdot(0.0),
    _dha(0.0),
    _dhb(0.0),
    _dhc(0.0),
    _dhd(0.0),
    _dhatilde(0.0),
    _dhbtilde(0.0),
    _dhctilde(0.0),
    _dhdtilde(0.0),
    _output_species(getParam<std::string>("output_species")),
    _charge_balance_species(getParam<std::string>("charge_balance_species")),
    _charge_balance_index(0),
    _mass_solvent_water_provided(isParamValid("mass_solvent_water")),
    _mass_solvent_water(isParamValid("mass_solvent_water") ? getParam<Real>("mass_solvent_water")
                                                           : 0.0),
    _moles_bulk_water_provided(isParamValid("moles_bulk_water")),
    _moles_bulk_water(isParamValid("moles_bulk_water") ? getParam<Real>("moles_bulk_water") : 0.0),
    _bulk_moles_species(getParam<std::vector<std::string>>("bulk_moles_species")),
    _bulk_moles_values(getParam<std::vector<Real>>("bulk_moles_values")),
    _free_molality_species(getParam<std::vector<std::string>>("free_molality_species")),
    _free_molality_values(getParam<std::vector<Real>>("free_molality_values")),
    _free_moles_mineral_species(getParam<std::vector<std::string>>("free_moles_mineral_species")),
    _free_moles_mineral_values(getParam<std::vector<Real>>("free_moles_mineral_values")),
    _prevent_precipitation(getParam<std::vector<std::string>>("prevent_precipitation")),
    _abs_tol(getParam<Real>("abs_tol")),
    _rel_tol(getParam<Real>("rel_tol")),
    _res0_times_rel(0.0),
    _max_iter(getParam<unsigned>("max_iter")),
    _verbose(getParam<bool>("verbose")),
    _tot_iter(0)
{
  if (_swap_out.size() != _swap_in.size())
    paramError("swap_out_of_basis must have same length as swap_into_basis");
  if ((_mass_solvent_water_provided == false && _moles_bulk_water_provided == false) ||
      (_mass_solvent_water_provided == true && _moles_bulk_water_provided == true))
    paramError("mass_solvent_water or moles_bulk_water must be provided, but not both");
  if (_activity_species.size() != _activity_values.size())
    paramError("activity_species must have same length as activity_values");
  for (const auto & val : _activity_values)
    if (val <= 0.0)
      paramError("activity_values must all be positive");
  if (_bulk_moles_species.size() != _bulk_moles_values.size())
    paramError("bulk_moles_species must have same length as bulk_moles_values");
  if (_free_molality_species.size() != _free_molality_values.size())
    paramError("free_molality_species must have same length as free_molality_values");
  if (_free_moles_mineral_species.size() != _free_moles_mineral_values.size())
    paramError("free_moles_mineral_species must have same length as free_moles_mineral_values");
}

void
EquilibriumReactionSolver::initialize()
{
  // do swaps desired by user.  any exception here is an error
  for (unsigned i = 0; i < _swap_out.size(); ++i)
    try
    {
      if (_swap_out[i] == _charge_balance_species)
        mooseError("Cannot swap out ",
                   _charge_balance_species,
                   " because it is the charge-balance species\n");
      _swapper.performSwap(_mgd, _swap_out[i], _swap_in[i]);
    }
    catch (const MooseException & e)
    {
      mooseError(e.what());
    }

  // check charge_balance_species is a basis and has a charge
  _charge_balance_index = checkChargeBalanceSpecies();

  // check exactly the correct information has been provided regarding bulk number of moles, free
  // molality, fugacity, etc
  checkICs();
}

void
EquilibriumReactionSolver::execute()
{
  bool still_swapping_minerals = true;
  while (still_swapping_minerals)
  {
    buildTemperatureDependentQuantities();
    buildAlgebraicInfo(
        _in_algebraic_system, _num_in_algebraic_system, _algebraic_index, _basis_index);
    _residual = DenseVector<Real>(_num_in_algebraic_system);
    _dmol = DenseVector<Real>(_num_in_algebraic_system);
    _jacobian = DenseMatrix<Real>(_num_in_algebraic_system, _num_in_algebraic_system);
    initBulkAndFree(_nw, _bulk_moles, _basis_molality);
    buildKnownBasisActivities(_basis_activity_known, _basis_activity);
    for (unsigned j = 0; j < _num_eqm; ++j)
      _eqm_molality[j] = 0.0;

    unsigned iter = 0;
    buildActivities(_basis_activity, _eqm_activity);
    computeEqmMolalities(_eqm_molality);
    _abs_residual = computeResidual(_residual);
    if (_verbose)
      _console << std::setprecision(_precision) << "iter = " << iter << " |R| = " << _abs_residual
               << std::endl;
    _res0_times_rel = _abs_residual * _rel_tol;

    while (_abs_residual >= _res0_times_rel && _abs_residual >= _abs_tol && iter < _max_iter)
    {
      while (_abs_residual >= _res0_times_rel && _abs_residual >= _abs_tol && iter < _max_iter)
      {
        iter += 1;
        _tot_iter += 1;
        computeJacobian(_dmj_dmi, _jacobian);
        solveAndUnderrelax(_jacobian, _dmol);
        updateWithdMol(_nw, _basis_molality);
        computeBulkMoles(_bulk_moles);
        buildActivities(_basis_activity, _eqm_activity);
        computeEqmMolalities(_eqm_molality);
        _abs_residual = computeResidual(_residual);
        if (_verbose)
          _console << "iter = " << iter << " |R| = " << _abs_residual << std::endl;
      }
      enforceChargeNeutrality(_bulk_moles, _bulk_moles_values);
      _abs_residual = computeResidual(_residual);
      if (_verbose)
        _console << "after enforcing charge neutrality, |R| = " << _abs_residual << std::endl;
    }

    if (iter >= _max_iter)
      mooseWarning("Number of iterations exceeds ", _max_iter, "\n");

    // compute the free number of moles for minerals
    computeFreeMineralMoles(_basis_molality);
    // compute the mineral saturation indices
    saturationIndices(_eqm_SI);

    unsigned swap_out_of_basis = 0;
    unsigned swap_into_basis = 0;
    still_swapping_minerals = mineralSwapNeeded(swap_out_of_basis, swap_into_basis);
    if (still_swapping_minerals)
    {
      if (swap_out_of_basis == 0)
        mooseError("EquilibriumReactionSolver: attempting to swap out ",
                   _mgd.eqm_species_name[swap_into_basis],
                   " in favor of water!  This could be because the algorithm would like to "
                   "swap out the charge-balance species, in which case you should choose a "
                   "different charge-balance species\n");
      // perform the swap
      try
      {
        _swapper.performSwap(_mgd, _bulk_moles, swap_out_of_basis, swap_into_basis);
      }
      catch (const MooseException & e)
      {
        mooseError(e.what());
      }
      _charge_balance_index = checkChargeBalanceSpecies();
      recordBulkMoles(_moles_bulk_water, _bulk_moles_species, _bulk_moles_values);
    }
  }

  outputResults();
}

void
EquilibriumReactionSolver::finalize()
{
}

std::vector<std::string>
EquilibriumReactionSolver::speciesOfInterest() const
{
  if (_output_species == "")
  {
    std::vector<std::string> species = _mgd.basis_species_name;
    species.insert(species.end(), _mgd.eqm_species_name.begin(), _mgd.eqm_species_name.end());
    species.insert(species.end(), _mgd.kin_species_name.begin(), _mgd.kin_species_name.end());
    return species;
  }
  else if (_mgd.basis_species_index.count(_output_species) == 1 ||
           _mgd.eqm_species_index.count(_output_species) == 1 ||
           _mgd.kin_species_index.count(_output_species) == 1)
    return {_output_species};
  return {};
}

void
EquilibriumReactionSolver::buildKnownBasisActivities(std::vector<bool> & basis_activity_known,
                                                     std::vector<Real> & basis_activity) const
{
  basis_activity_known = std::vector<bool>(_num_basis, false);
  basis_activity.resize(_num_basis);
  const unsigned num_activities_provided = _activity_species.size();
  for (unsigned i = 0; i < num_activities_provided; ++i)
  {
    const std::string sp = _activity_species[i];
    if (_mgd.basis_species_index.count(sp) == 1)
    {
      // this must definitely be true for gases by the constructor check and the fact we never
      // swap gases out of the basis
      const unsigned basis_i = _mgd.basis_species_index.at(sp);
      basis_activity_known[basis_i] = true;
      basis_activity[basis_i] = _activity_values[i];
    }
  }

  // this will nicely overwrite a poor user choice of activity_mineral != 1.0, which is forbidden
  // in the geochemsitry module
  for (unsigned basis_i = 0; basis_i < _num_basis; ++basis_i)
    if (_mgd.basis_species_mineral[basis_i])
    {
      basis_activity_known[basis_i] = true;
      basis_activity[basis_i] = 1.0;
    }
}

void
EquilibriumReactionSolver::buildActivities(std::vector<Real> & basis_activity,
                                           std::vector<Real> & eqm_activity) const
{
  basis_activity.resize(_num_basis);
  eqm_activity.resize(_num_eqm);

  // TODO: make 3.0 a user-defined parameter
  const Real ionic_str = std::min(
      3.0, GeochemistryIonicStrength::ionicStrength(_mgd, _basis_molality, _eqm_molality, {}));
  const Real stoi_ionic_str = std::min(3.0,
                                       GeochemistryIonicStrength::stoichiometricIonicStrength(
                                           _mgd, _basis_molality, _eqm_molality, {}));

  for (unsigned basis_i = 0; basis_i < _num_basis; ++basis_i)
    if (_basis_activity_known[basis_i])
      continue;
    else if (basis_i == 0)
    {
      basis_activity[basis_i] = std::exp(GeochemistryActivity::lnActivityDHBdotWater(
          stoi_ionic_str, _dhA, _dhatilde, _dhbtilde, _dhctilde, _dhdtilde));
    }
    else if (_mgd.basis_species_radius[basis_i] == -0.5)
    {
      basis_activity[basis_i] = std::pow(
          10.0,
          GeochemistryActivity::log10ActCoeffDHBdotNeutral(ionic_str, _dha, _dhb, _dhc, _dhd));
    }
    else if (_mgd.basis_species_radius[basis_i] == -1.0)
    {
      basis_activity[basis_i] =
          std::pow(10.0, GeochemistryActivity::log10ActCoeffDHBdotAlternative(ionic_str, _dhBdot));
    }
    else if (_mgd.basis_species_charge[basis_i] == 0.0)
    {
      basis_activity[basis_i] = 1.0;
    }
    else
    {
      basis_activity[basis_i] =
          std::pow(10.0,
                   GeochemistryActivity::log10ActCoeffDHBdot(_mgd.basis_species_charge[basis_i],
                                                             _mgd.basis_species_radius[basis_i],
                                                             std::sqrt(ionic_str),
                                                             _dhA,
                                                             _dhB,
                                                             _dhBdot));
    }

  for (unsigned eqm_i = 0; eqm_i < _num_eqm; ++eqm_i)
    if (_mgd.eqm_species_mineral[eqm_i])
      eqm_activity[eqm_i] = 1.0;
    else if (_mgd.eqm_species_gas[eqm_i])
      mooseError("Equilibrium species cannot be gases in EquilibriumReactionSolver");
    else if (_mgd.eqm_species_radius[eqm_i] == -0.5)
    {
      eqm_activity[eqm_i] = std::pow(
          10.0,
          GeochemistryActivity::log10ActCoeffDHBdotNeutral(ionic_str, _dha, _dhb, _dhc, _dhd));
    }
    else if (_mgd.eqm_species_radius[eqm_i] == -1.0)
    {
      eqm_activity[eqm_i] =
          std::pow(10.0, GeochemistryActivity::log10ActCoeffDHBdotAlternative(ionic_str, _dhBdot));
    }
    else if (_mgd.eqm_species_charge[eqm_i] == 0.0)
    {
      eqm_activity[eqm_i] = 1.0;
    }
    else
    {
      eqm_activity[eqm_i] =
          std::pow(10.0,
                   GeochemistryActivity::log10ActCoeffDHBdot(_mgd.eqm_species_charge[eqm_i],
                                                             _mgd.eqm_species_radius[eqm_i],
                                                             std::sqrt(ionic_str),
                                                             _dhA,
                                                             _dhB,
                                                             _dhBdot));
    }
}

void
EquilibriumReactionSolver::buildAlgebraicInfo(std::vector<bool> & in_algebraic_system,
                                              unsigned & num_in_algebraic_system,
                                              std::vector<unsigned> & algebraic_index,
                                              std::vector<unsigned> & basis_index) const
{
  in_algebraic_system.resize(_num_basis);
  algebraic_index.resize(_num_basis);
  basis_index.resize(_num_basis);

  // build in_algebraic_system
  for (const auto & name_index : _mgd.basis_species_index)
  {
    const std::string name = name_index.first;
    const unsigned ind = name_index.second;
    if (name == "H2O")
      in_algebraic_system[ind] = _moles_bulk_water_provided;
    else if (_mgd.basis_species_mineral[ind])
      in_algebraic_system[ind] = false;
    else if (_mgd.basis_species_gas[ind])
      in_algebraic_system[ind] = false;
    else
      in_algebraic_system[ind] =
          (std::find(_free_molality_species.begin(), _free_molality_species.end(), name) ==
           _free_molality_species.end()); // if it has not been prescribed a free molality it is in
                                          // the algebraic system
  }

  // build algebraic_index and basis_index
  num_in_algebraic_system = 0;
  for (unsigned basis_i = 0; basis_i < _num_basis; ++basis_i)
    if (in_algebraic_system[basis_i])
    {
      algebraic_index[basis_i] = _num_in_algebraic_system;
      basis_index[_num_in_algebraic_system] = basis_i;
      num_in_algebraic_system += 1;
    }
}

void
EquilibriumReactionSolver::initBulkAndFree(Real & nw,
                                           DenseVector<Real> & bulk_moles,
                                           std::vector<Real> & basis_molality) const
{
  bulk_moles.resize(_num_basis);
  basis_molality.resize(_num_basis);

  if (_mass_solvent_water_provided)
  {
    nw = _mass_solvent_water; // fixed by user
    bulk_moles(0) =
        0.0; // actual value is determined after the solution to the algebraic system is found
  }
  else
  {
    bulk_moles(0) = _moles_bulk_water;                   // fixed by user
    nw = 0.999 * _moles_bulk_water / moles_per_kg_water; // initial guess
  }
  for (unsigned i = 0; i < _free_molality_species.size(); ++i)
  {
    const std::string name = _free_molality_species[i];
    if (_mgd.basis_species_index.count(name) == 0)
      continue;
    const unsigned basis_ind = _mgd.basis_species_index.at(name);
    basis_molality[basis_ind] = _free_molality_values[i]; // fixed by user
    bulk_moles(basis_ind) =
        _free_molality_values[i]; // actual value determined after algebraic solution found
  }
  for (unsigned i = 0; i < _bulk_moles_species.size(); ++i)
  {
    const std::string name = _bulk_moles_species[i];
    if (_mgd.basis_species_index.count(name) == 0)
      continue;
    const unsigned basis_ind = _mgd.basis_species_index.at(name);
    bulk_moles(basis_ind) = _bulk_moles_values[i];                // fixed by user
    basis_molality[basis_ind] = 0.9 * _bulk_moles_values[i] / nw; // initial guess
  }
}

void
EquilibriumReactionSolver::computeEqmMolalities(std::vector<Real> & eqm_molality) const
{
  eqm_molality.resize(_num_eqm);
  Real max_allowed = 0.0;
  for (unsigned i = 0; i < _num_basis; ++i)
    max_allowed += std::abs(_basis_molality[i]);
  max_allowed *= mol_multiplier;
  max_allowed = std::log10(max_allowed);

  for (unsigned j = 0; j < _num_eqm; ++j)
  {
    if (_mgd.eqm_species_mineral[j])
    {
      eqm_molality[j] = 0.0;
      continue;
    }
    // form log10(eqm_molality) first, to avoid overflows and underflows
    eqm_molality[j] = -std::log10(_eqm_activity[j]) - _eqm_log10K[j];
    const Real ap = activityProduct(j);
    if (ap == 0.0)
    {
      eqm_molality[j] = 0.0;
    }
    else
    {
      eqm_molality[j] += std::log10(ap);
      eqm_molality[j] = std::min(std::max(-mol_cutoff, eqm_molality[j]), max_allowed);
      eqm_molality[j] = std::pow(10.0, eqm_molality[j]);
    }
  }
}

Real
EquilibriumReactionSolver::computeResidual(DenseVector<Real> & residual) const
{
  for (unsigned a_ind = 0; a_ind < _num_in_algebraic_system; ++a_ind)
  {
    const unsigned basis_ind = _basis_index[a_ind];
    residual(a_ind) = -_bulk_moles(basis_ind);
    if (basis_ind == 0)
      residual(a_ind) += _nw * moles_per_kg_water;
    else
      residual(a_ind) += _nw * _basis_molality[basis_ind];
    for (unsigned eqm_i = 0; eqm_i < _num_eqm; ++eqm_i)
      residual(a_ind) += _nw * _mgd.eqm_stoichiometry(eqm_i, basis_ind) * _eqm_molality[eqm_i];
  }
  return residual.l1_norm();
}

void
EquilibriumReactionSolver::computeJacobian(DenseMatrix<Real> & dmj_dmi,
                                           DenseMatrix<Real> & jacobian) const
{
  dmj_dmi.zero();
  for (unsigned j = 0; j < _num_eqm; ++j)
    for (unsigned i = 0; i < _num_basis; ++i)
    {
      if (i == 0)
        dmj_dmi(j, i) = 0.0;
      else if (_basis_activity_known[i])
        dmj_dmi(j, i) = 0.0;
      else if (_eqm_molality[j] == 0.0)
        dmj_dmi(j, i) = 0.0; // incorrect if stoichiometry=1 and basis_molality=0, which is
                             // surely unlikely. This case also covers the mineral case where we
                             // are demanding _eqm_molality = 0
      else if (_mgd.eqm_stoichiometry(j, i) == 0.0)
        dmj_dmi(j, i) = 0.0;
      else
        dmj_dmi(j, i) = _mgd.eqm_stoichiometry(j, i) * _eqm_molality[j] /
                        _basis_molality[i]; // denominator cannot be zero because
                                            // otherwise _eqm_molality[j] would be zero
    }

  jacobian.zero();

  // derivative wrt _nw
  if (_basis_index[0] == 0)
  {
    // _nw is a degree of freedom
    for (unsigned a_ind = 0; a_ind < _num_in_algebraic_system; ++a_ind)
      jacobian(a_ind, 0) += _residual(a_ind) / _nw;
  }

  // derivative wrt m_i
  for (unsigned a_ind = 0; a_ind < _num_in_algebraic_system; ++a_ind)
  {
    if (_basis_index[a_ind] == 0)
    {
      // water
      continue;
    }
    else
      jacobian(a_ind, a_ind) += _nw;
  }

  // derivative wrt m_j
  for (unsigned a_ind = 0; a_ind < _num_in_algebraic_system; ++a_ind)
  {
    const unsigned basis_i = _basis_index[a_ind];
    for (unsigned b_ind = 0; b_ind < _num_in_algebraic_system; ++b_ind)
    {
      const unsigned basis_j = _basis_index[b_ind];
      for (unsigned eqm_ind = 0; eqm_ind < _num_eqm; ++eqm_ind)
        jacobian(a_ind, b_ind) +=
            _nw * _mgd.eqm_stoichiometry(eqm_ind, basis_i) * dmj_dmi(eqm_ind, basis_j);
    }
  }
}

void
EquilibriumReactionSolver::solveAndUnderrelax(DenseMatrix<Real> & jacobian,
                                              DenseVector<Real> & dmol) const
{
  jacobian.lu_solve(_residual, dmol);

  // at this point we want to do molality = molality - dmol

  // Bethke recommends underrelaxation
  Real one_over_delta = 1.0;
  for (unsigned a_ind = 0; a_ind < _num_in_algebraic_system; ++a_ind)
  {
    const unsigned basis_i = _basis_index[a_ind];
    if (basis_i == 0)
      one_over_delta = std::max(one_over_delta, dmol(a_ind) * 2.0 / _nw);
    else
      one_over_delta = std::max(one_over_delta, dmol(a_ind) * 2.0 / _basis_molality[basis_i]);
  }
  for (unsigned a_ind = 0; a_ind < _num_in_algebraic_system; ++a_ind)
    dmol(a_ind) = dmol(a_ind) / one_over_delta;
}

void
EquilibriumReactionSolver::updateWithdMol(Real & nw, std::vector<Real> & basis_molality) const
{
  for (unsigned a_ind = 0; a_ind < _num_in_algebraic_system; ++a_ind)
  {
    const unsigned basis_i = _basis_index[a_ind];
    if (basis_i == 0)
      nw -= _dmol(a_ind);
    else
      basis_molality[basis_i] -= _dmol(a_ind);
  }
}

void
EquilibriumReactionSolver::buildTemperatureDependentQuantities()
{
  // at 25degC currently
  for (unsigned i = 0; i < _num_eqm; ++i)
    _eqm_log10K[i] = _mgd.eqm_log10K(i, 1);
  // charged species Debye-Huckel at 25degC
  _dhA = 0.5092;
  _dhB = 0.3283;
  _dhBdot = 0.035;
  // water Debye-Huckel at 25degC
  _dhatilde = 1.45397;
  _dhbtilde = 0.022357;
  _dhctilde = 0.0093804;
  _dhdtilde = -0.0005262;
  // neutral species Debye-Huckel at 25degC
  _dha = 0.1127;
  _dhb = -0.01049;
  _dhc = 0.001545;
  _dhd = 0.0;
}

void
EquilibriumReactionSolver::computeBulkMoles(DenseVector<Real> & bulk_moles) const
{
  // water component
  if (_mass_solvent_water_provided)
  {
    bulk_moles(0) = moles_per_kg_water;
    for (unsigned j = 0; j < _num_eqm; ++j)
      bulk_moles(0) += _mgd.eqm_stoichiometry(j, 0) * _eqm_molality[j];
    bulk_moles(0) *= _nw;
  }
  else
  {
    bulk_moles(0) = _moles_bulk_water; // fixed by user
  }

  // all aqueous components that are provided with a free molality value (this can never be minerals
  // or gases or H2O)
  for (unsigned i = 0; i < _free_molality_species.size(); ++i)
  {
    const std::string name = _free_molality_species[i];
    if (_mgd.basis_species_index.count(name) == 0)
      continue;
    const unsigned basis_ind = _mgd.basis_species_index.at(name);
    bulk_moles(basis_ind) = _basis_molality[basis_ind];
    for (unsigned j = 0; j < _num_eqm; ++j)
      bulk_moles(basis_ind) += _mgd.eqm_stoichiometry(j, basis_ind) * _eqm_molality[j];
    bulk_moles(basis_ind) *= _nw;
  }

  // minerals that are provided with free_moles
  for (unsigned i = 0; i < _free_moles_mineral_species.size(); ++i)
  {
    const std::string name = _free_moles_mineral_species[i];
    if (_mgd.basis_species_index.count(name) == 0)
      continue;
    const unsigned basis_ind = _mgd.basis_species_index.at(name);
    if (!_mgd.basis_species_mineral[basis_ind])
      continue;
    bulk_moles(basis_ind) = _free_moles_mineral_values[i];
    for (unsigned i = 0; i < _num_eqm; ++i)
      bulk_moles(basis_ind) += _nw * _mgd.eqm_stoichiometry(i, basis_ind) * _eqm_molality[i];
  }

  // gas species
  for (unsigned i = 0; i < _num_basis; ++i)
  {
    if (!_mgd.basis_species_gas[i])
      continue;
    bulk_moles(i) = 0.0;
    for (unsigned j = 0; j < _num_eqm; ++j)
      bulk_moles(i) += _nw * _mgd.eqm_stoichiometry(j, i) * _eqm_molality[j];
  }

  // all components (including minerals and gases) that are provided with a bulk_moles.  Basis swaps
  // due to mineral dissolution/precipitation will have potentially modified _bulk_moles_species
  // (for those cases where the user has not provided the mass of solvent water, free molality (for
  // aqueous basis components) or free number of moles (for minerals)), so the following will
  // potentially overwrite the water component (if it was provided a bulk value by the user) and gas
  // values
  for (unsigned i = 0; i < _bulk_moles_species.size(); ++i)
  {
    const std::string name = _bulk_moles_species[i];
    if (_mgd.basis_species_index.count(name) == 0)
      continue;
    const unsigned basis_ind = _mgd.basis_species_index.at(name);
    bulk_moles(basis_ind) = _bulk_moles_values[i]; // fixed by user
  }
}

void
EquilibriumReactionSolver::enforceChargeNeutrality(DenseVector<Real> & bulk_moles,
                                                   std::vector<Real> & bulk_moles_values) const
{
  Real tot_charge = 0.0;
  for (unsigned i = 0; i < _num_basis; ++i)
  {
    if (i == _charge_balance_index)
      continue;
    tot_charge += bulk_moles(i) * _mgd.basis_species_charge[i];
  }
  bulk_moles(_charge_balance_index) =
      -tot_charge / _mgd.basis_species_charge[_charge_balance_index];
  // put this info into _bulk_moles_values so the algorithm does not resort to the incorrect
  // bulk_moles
  for (unsigned i = 0; i < _bulk_moles_species.size(); ++i)
    if (_bulk_moles_species[i] == _charge_balance_species)
      bulk_moles_values[i] = bulk_moles(_charge_balance_index);
}

void
EquilibriumReactionSolver::outputResults() const
{
  unsigned ind = 0;
  _console << std::setprecision(_precision);

  _console << "\nSummary:\n";

  _console << "Total number of iterations required = " << _tot_iter << "\n";
  _console << "Error in calculation = " << _abs_residual << "mol\n";

  Real charge = 0.0;
  for (unsigned i = 1; i < _num_basis; ++i) // do not loop over water
    charge += _mgd.basis_species_charge[i] * _bulk_moles(i);
  _console << "Charge of solution = " << charge << "mol\n";

  _console << "Mass of solvent water = " << _nw << "kg\n";

  // Output total aqueous mass
  Real mass = _nw;
  for (unsigned i = 1; i < _num_basis; ++i) // do not loop over water
    mass += _bulk_moles(i) * _mgd.basis_species_molecular_weight[i] / 1000.0;
  _console << "Mass of aqueous solution = " << mass << "kg\n";

  // Output the aqueous solution pH, if relevant
  if (_mgd.basis_species_index.count("H+"))
    _console << "pH = "
             << -std::log10(_basis_molality[_mgd.basis_species_index.at("H+")] *
                            _basis_activity[_mgd.basis_species_index.at("H+")])
             << "\n";
  if (_mgd.eqm_species_index.count("H+"))
    _console << "pH = "
             << -std::log10(_eqm_molality[_mgd.eqm_species_index.at("H+")] *
                            _eqm_activity[_mgd.eqm_species_index.at("H+")])
             << "\n";

  // Output the aqueous solution pe, if relevant
  if (_mgd.basis_species_index.count("e-"))
    _console << "pe = "
             << -std::log10(_basis_molality[_mgd.basis_species_index.at("e-")] *
                            _basis_activity[_mgd.basis_species_index.at("e-")])
             << "\n";
  if (_mgd.eqm_species_index.count("e-"))
    _console << "pe = "
             << -std::log10(_eqm_molality[_mgd.eqm_species_index.at("e-")] *
                            _eqm_activity[_mgd.eqm_species_index.at("e-")])
             << "\n";

  // Output ionic strengths
  _console << "Ionic strength = "
           << GeochemistryIonicStrength::ionicStrength(_mgd, _basis_molality, _eqm_molality, {})
           << "mol/kg(solvent water)\n";
  _console << "Stoichiometric ionic strength = "
           << GeochemistryIonicStrength::stoichiometricIonicStrength(
                  _mgd, _basis_molality, _eqm_molality, {})
           << "mol/kg(solvent water)\n";

  // Output activity of water
  _console << "Activity of water = " << _basis_activity[0] << "\n";

  // Output the basis species information, sorted by molality
  std::vector<unsigned> basis_order(_num_basis);
  ind = 0;
  std::iota(basis_order.begin(), basis_order.end(), ind++);
  std::sort(basis_order.begin(), basis_order.end(), [&](int i, int j) {
    return _basis_molality[i] > _basis_molality[j];
  });
  _console << "\nBasis Species:\n";
  for (const auto & i : basis_order)
    if (i == 0 || _mgd.basis_species_gas[i])
      continue;
    else
    {
      _console << _mgd.basis_species_name[i] << ";  bulk_moles = " << _bulk_moles(i)
               << "mol;  bulk_conc = "
               << _bulk_moles(i) * _mgd.basis_species_molecular_weight[i] * 1000.0 / mass
               << "mg/kg(solution);  molality = " << _basis_molality[i]
               << "mol/kg(solvent water);  free_conc = "
               << _basis_molality[i] * _mgd.basis_species_molecular_weight[i] * 1000.0
               << "mg/kg(solvent water)";
      if (_mgd.basis_species_mineral[i])
        _console << "\n";
      else
        _console << ";  act_coeff = " << _basis_activity[i]
                 << ";  log10(a) = " << std::log10(_basis_molality[i] * _basis_activity[i]) << "\n";
    }
  for (unsigned i = 0; i < _num_basis; ++i)
    if (_mgd.basis_species_gas[i])
      _console << _mgd.basis_species_name[i] << ";  fugacity = " << _basis_activity[i] << "\n";

  // Output the equilibrium species info, sorted by molality
  std::vector<unsigned> eqm_order(_num_eqm);
  ind = 0;
  std::iota(eqm_order.begin(), eqm_order.end(), ind++);
  std::sort(eqm_order.begin(), eqm_order.end(), [&](int i, int j) {
    return _eqm_molality[i] > _eqm_molality[j];
  });
  _console << "\nEquilibrium Species:\n";
  for (const auto & i : eqm_order)
    if (_eqm_molality[i] <= std::pow(10.0, -mol_cutoff))
      break;
    else if (_mgd.eqm_species_gas[i])
      continue;
    else
      _console << _mgd.eqm_species_name[i] << ";  molality = " << _eqm_molality[i]
               << "mol/kg(solvent water);  free_conc = "
               << _eqm_molality[i] * _mgd.eqm_species_molecular_weight[i] * 1000.0
               << "mg/kg(solvent water);  act_coeff = " << _eqm_activity[i]
               << ";  log10(a) = " << std::log10(_eqm_molality[i] * _eqm_activity[i])
               << ";  log10K = " << _eqm_log10K[i] << "\n";
  for (unsigned i = 0; i < _num_eqm; ++i)
    if (_mgd.eqm_species_gas[i])
      _console << _mgd.eqm_species_name[i] << ";  fugacity = " << _eqm_activity[i] << "\n";

  // Output the mineral info, sorted by saturation indices
  std::vector<unsigned> mineral_order(_num_eqm);
  ind = 0;
  std::iota(mineral_order.begin(), mineral_order.end(), ind++);
  std::sort(mineral_order.begin(), mineral_order.end(), [&](int i, int j) {
    return _eqm_SI[i] > _eqm_SI[j];
  });
  _console << "\nMinerals:\n";
  for (const auto & i : mineral_order)
    if (_mgd.eqm_species_mineral[i])
      _console << _mgd.eqm_species_name[i] << ";  SI = " << _eqm_SI[i] << "\n";
}

void
EquilibriumReactionSolver::saturationIndices(std::vector<Real> & si) const
{
  si.resize(_num_eqm);
  for (unsigned i = 0; i < _num_eqm; ++i)
    if (_mgd.eqm_species_mineral[i])
      si[i] = std::log10(activityProduct(i)) - _eqm_log10K[i];
    else
      si[i] = 0.0;
}

Real
EquilibriumReactionSolver::activityProduct(unsigned eqm_index) const
{
  Real ap = 1.0;
  for (unsigned i = 0; i < _num_basis; ++i)
  {
    if (_basis_activity_known[i] ||
        i == 0) // could be a mineral, or a gas of known fugacity, or a basis species with known
                // activity set by the user, or water.  In all these cases, molality is ignored
    {
      ap *= std::pow(_basis_activity[i], _mgd.eqm_stoichiometry(eqm_index, i));
    }
    else // activity not known
    {
      const Real a = _basis_activity[i] * _basis_molality[i];
      if (a <= 0.0) // Newton has strayed into unphysical regime
      {
        // note that if these checks are modified, similar checks in the jacobian must be modified
        if (_mgd.eqm_stoichiometry(eqm_index, i) != 0.0)
        {
          ap = 0.0;
          break;
        }
        else
        {
          // although activity <= 0, stoichiometric coefficient == 0, so no contribution
        }
      }
      else // basis_activity > 0
      {
        ap *= std::pow(a, _mgd.eqm_stoichiometry(eqm_index, i));
      }
    }
  }
  return ap;
}

unsigned
EquilibriumReactionSolver::checkChargeBalanceSpecies() const
{
  unsigned cbi = 0;
  if (_mgd.basis_species_index.count(_charge_balance_species) == 0)
    mooseError("Cannot enforce charge balance using ",
               _charge_balance_species,
               " because it is not in the basis");
  cbi = _mgd.basis_species_index.at(_charge_balance_species);
  if (_mgd.basis_species_charge[cbi] == 0.0)
    mooseError("Cannot enforce charge balance using ",
               _charge_balance_species,
               " because it has zero charge");
  return cbi;
}

void
EquilibriumReactionSolver::checkICs() const
{
  for (const auto & name_index : _mgd.basis_species_index)
  {
    const std::string name = name_index.first;
    if (name == "H2O")
    {
      if (std::find(_bulk_moles_species.begin(), _bulk_moles_species.end(), name) !=
              _bulk_moles_species.end() ||
          std::find(_free_molality_species.begin(), _free_molality_species.end(), name) !=
              _free_molality_species.end() ||
          std::find(_free_moles_mineral_species.begin(), _free_moles_mineral_species.end(), name) !=
              _free_moles_mineral_species.end())
        mooseError("H2O cannot be provided with a free_molality, free_moles_mineral or bulk_moles");
      continue;
    }
    const unsigned ind = name_index.second;
    if (_mgd.basis_species_mineral[ind])
    {
      if (std::find(_bulk_moles_species.begin(), _bulk_moles_species.end(), name) ==
              _bulk_moles_species.end() &&
          std::find(_free_moles_mineral_species.begin(), _free_moles_mineral_species.end(), name) ==
              _free_moles_mineral_species.end())
        mooseError(
            "The mineral ",
            name,
            " is in the basis so it must be provided with either a bulk or free number of moles");
      if (std::find(_bulk_moles_species.begin(), _bulk_moles_species.end(), name) !=
              _bulk_moles_species.end() &&
          std::find(_free_moles_mineral_species.begin(), _free_moles_mineral_species.end(), name) !=
              _free_moles_mineral_species.end())
        mooseError("The mineral ",
                   name,
                   " cannot be provided with both a bulk composition "
                   "(bulk_moles_species) and free number of moles (free_moles_mineral_species)");
      if (std::find(_free_molality_species.begin(), _free_molality_species.end(), name) !=
          _free_molality_species.end())
        mooseError("The mineral ",
                   name,
                   " cannot be provided with a free_molality because it is a mineral (use "
                   "free_moles_mineral_species instead)");
    }
    else if (_mgd.basis_species_gas[ind])
    {
      if (std::find(_activity_species.begin(), _activity_species.end(), name) ==
          _activity_species.end())
        mooseError("The gas ",
                   name,
                   " is in the basis so must have its fugacity fixed using activity_species and "
                   "activity_values");
      if (std::find(_bulk_moles_species.begin(), _bulk_moles_species.end(), name) !=
              _bulk_moles_species.end() ||
          std::find(_free_molality_species.begin(), _free_molality_species.end(), name) !=
              _free_molality_species.end() ||
          std::find(_free_moles_mineral_species.begin(), _free_moles_mineral_species.end(), name) !=
              _free_moles_mineral_species.end())
        mooseError("The gas ",
                   name,
                   " cannot be provided with a free_molality, free_moles_mineral or bulk_moles "
                   "because it is a gas");
    }
    else
    {
      // species is an aqueous primary
      if (std::find(_bulk_moles_species.begin(), _bulk_moles_species.end(), name) ==
              _bulk_moles_species.end() &&
          std::find(_free_molality_species.begin(), _free_molality_species.end(), name) ==
              _free_molality_species.end())
        mooseError("The species ",
                   name,
                   " is in the basis so it must be provided with either a bulk composition "
                   "(bulk_moles_species) or free molality (free_molality_species)");
      if (std::find(_bulk_moles_species.begin(), _bulk_moles_species.end(), name) !=
              _bulk_moles_species.end() &&
          std::find(_free_molality_species.begin(), _free_molality_species.end(), name) !=
              _free_molality_species.end())
        mooseError("The species ",
                   name,
                   " is in the basis.  It cannot be provided with both a bulk composition "
                   "(bulk_moles_species) and free molality (free_molality_species)");
    }
  }
}

void
EquilibriumReactionSolver::computeFreeMineralMoles(std::vector<Real> & basis_molality) const
{
  for (unsigned i = 0; i < _num_basis; ++i)
    if (_mgd.basis_species_mineral[i])
    {
      basis_molality[i] = _bulk_moles(i);
      for (unsigned j = 0; j < _num_eqm; ++j)
        basis_molality[i] -= _nw * _mgd.eqm_stoichiometry(j, i) * _eqm_molality[j];
    }
}

bool
EquilibriumReactionSolver::mineralSwapNeeded(unsigned & swap_out_of_basis,
                                             unsigned & swap_into_basis) const
{
  bool swap_needed = false;
  unsigned ind = 0;

  // check if any basis minerals have negative free number of moles
  std::vector<unsigned> mineral_n_order(_num_basis);
  ind = 0;
  std::iota(mineral_n_order.begin(), mineral_n_order.end(), ind++);
  std::sort(mineral_n_order.begin(), mineral_n_order.end(), [&](int i, int j) {
    return _basis_molality[i] < _basis_molality[j];
  });

  for (const auto & i : mineral_n_order)
  {
    if (_basis_molality[i] >= 0.0)
    {
      // since we're going through _basis_molality in ascending order, as soon as we get >=0, we're
      // safe
      break;
    }
    else if (_mgd.basis_species_mineral[i])
    {
      swap_needed = true;
      swap_out_of_basis = i;
      swap_into_basis = 0;
      Real best_stoi = 0.0;
      for (unsigned j = 0; j < _num_eqm; ++j)
      {
        const Real stoi = std::abs(_mgd.eqm_stoichiometry(j, i)) * _eqm_molality[j];
        if (stoi > best_stoi)
        {
          best_stoi = stoi;
          swap_into_basis = j;
        }
      }
      if (_verbose)
        _console << "Mineral " << _mgd.basis_species_name[swap_out_of_basis]
                 << " consumed.  Swapping it with equilibrium species "
                 << _mgd.eqm_species_name[swap_into_basis] << std::endl;
      break;
    }
  }
  if (swap_needed)
    return true;

  // check maximum saturation index is not positive
  std::vector<unsigned> mineral_SI_order(_num_eqm);
  ind = 0;
  std::iota(mineral_SI_order.begin(), mineral_SI_order.end(), ind++);
  std::sort(mineral_SI_order.begin(), mineral_SI_order.end(), [&](int i, int j) {
    return _eqm_SI[i] > _eqm_SI[j];
  });

  for (const auto & i : mineral_SI_order)
    if (_eqm_SI[i] > 0.0)
      if (_mgd.eqm_species_mineral[i])
        if (std::find(_prevent_precipitation.begin(),
                      _prevent_precipitation.end(),
                      _mgd.eqm_species_name[i]) == _prevent_precipitation.end())
        {
          // mineral has positive saturation index and user is not preventing its precipitation.
          // determine the basis species to swap out
          swap_needed = true;
          swap_into_basis = i;
          swap_out_of_basis = 0;
          Real best_stoi = 0.0;
          for (unsigned j = 1; j < _num_basis; ++j)
          {
            if (_basis_molality[j] > 0.0 && j != _charge_balance_index &&
                !_mgd.basis_species_gas[j]) // && !_mgd.basis_species_mineral[j])
            {
              // don't want to swap out the charge-balance species or any gases of fixed fugacity
              const Real stoi = std::abs(_mgd.eqm_stoichiometry(i, j)) / _basis_molality[j];
              if (stoi > best_stoi)
              {
                best_stoi = stoi;
                swap_out_of_basis = j;
              }
            }
          }
          if (_verbose)
            _console << "Mineral " << _mgd.eqm_species_name[i]
                     << " supersaturated.  Swapping it with basis species "
                     << _mgd.basis_species_name[swap_out_of_basis] << std::endl;
          break;
        }
  return swap_needed;
}

void
EquilibriumReactionSolver::recordBulkMoles(Real & moles_bulk_water,
                                           std::vector<std::string> & bulk_moles_species,
                                           std::vector<Real> & bulk_moles_values) const
{
  if (_moles_bulk_water_provided)
    moles_bulk_water = _bulk_moles(0);

  bulk_moles_species = std::vector<std::string>(0);
  bulk_moles_values = std::vector<Real>(0);
  for (unsigned i = 1; i < _num_basis; ++i) // do not loop over water
  {
    const std::string name = _mgd.basis_species_name[i];
    if (std::find(_free_molality_species.begin(), _free_molality_species.end(), name) !=
        _free_molality_species.end())
      continue;
    if (std::find(_free_moles_mineral_species.begin(), _free_moles_mineral_species.end(), name) !=
        _free_moles_mineral_species.end())
      continue;
    // free values not provided, so we need to record the bulk values
    bulk_moles_species.push_back(name);
    bulk_moles_values.push_back(_bulk_moles(i));
  }
}
