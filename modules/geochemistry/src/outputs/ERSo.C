//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ERSo.h"
#include "GeochemistryConstants.h"
#include "GeochemistryIonicStrength.h"
#include "GeochemistryFormattedOutput.h"

registerMooseObject("GeochemistryApp", ERSo);

InputParameters
ERSo::validParams()
{
  InputParameters params = Output::validParams();
  params.addRequiredParam<UserObjectName>(
      "model_definition",
      "The name of the GeochemicalModelDefinition user object.  Only equilibrium reactions are "
      "solved by ERSo, so the model_definition must not contain any kinetic "
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
      "nernst_swap_out_of_basis",
      "Before outputting Nernst potentials (after solving the system) these species are swapped "
      "out of the basis.  Often this is identical to swap_into_basis, so that the Nernst "
      "potentials are defined in terms of the original model definition.  There must be the same "
      "number of species in nernst_swap_out_of_basis and nernst_swap_into_basis.  If this list "
      "contains more than one species, the swapping is performed one-by-one, starting with the "
      "first pair (nernst_swap_out_of_basis[0] and nernst_swap_into_basis[0]) then the next pair, "
      "etc");
  params.addParam<std::vector<std::string>>("nernst_swap_into_basis",
                                            "Before outputting Nernst potentials (after solving "
                                            "the system) these species are swapped into the basis");
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
  params.addParam<Real>(
      "mol_cutoff",
      1E-40,
      "Information regarding species with molalities less than this amount will not be outputted");
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
  params.addParam<bool>(
      "ionic_str_using_basis_only",
      false,
      "If set to true, ionic strength and stoichiometric ionic strength will be computed using "
      "only the basis molalities, ignoring molalities of equilibrium species.  Since basis "
      "molality is usually greater than equilibrium molality, and the whole Debye-Huckel concept "
      "of activitiy coefficients depending on ionic strength is only approximate in practice, "
      "setting this parameter true often results in a reasonable approximation.  It can aid in "
      "convergence since it eliminates problems associated with unphysical huge equilibrium "
      "molalities that can occur during Newton iteration to the solution");
  params.addClassDescription("Solves a geochemical reaction system");

  params.set<ExecFlagEnum>("execute_on") = {EXEC_FINAL};
  return params;
}

ERSo::ERSo(const InputParameters & parameters)
  : Output(parameters),
    UserObjectInterface(this),
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
         getParam<Real>("min_initial_molality"),
         getParam<bool>("ionic_str_using_basis_only")),
    _original_redox_lhs(
        getUserObject<GeochemicalModelDefinition>("model_definition").getDatabase().redox_lhs),
    _num_basis(_egs.getNumInBasis()),
    _num_eqm(_egs.getNumInEquilibrium()),
    _mgd(_egs.getModelGeochemicalDatabase()),
    _num_in_algebraic_system(_egs.getNumInAlgebraicSystem()),
    _residual(_num_in_algebraic_system),
    _abs_residual(0.0),
    _jacobian(_num_in_algebraic_system, _num_in_algebraic_system),
    _new_mol(_num_in_algebraic_system),
    _precision(getParam<unsigned int>("precision")),
    _prevent_precipitation(getParam<std::vector<std::string>>("prevent_precipitation")),
    _abs_tol(getParam<Real>("abs_tol")),
    _rel_tol(getParam<Real>("rel_tol")),
    _res0_times_rel(0.0),
    _stoi_tol(getParam<Real>("stoichiometry_tolerance")),
    _max_iter(getParam<unsigned>("max_iter")),
    _verbose(getParam<bool>("verbose")),
    _max_initial_residual(getParam<Real>("max_initial_residual")),
    _swap_threshold(getParam<Real>("swap_threshold") * _abs_tol),
    _max_ionic_strength(getParam<Real>("max_ionic_strength")),
    _ramp_max_ionic_strength(getParam<unsigned>("ramp_max_ionic_strength")),
    _mol_cutoff(getParam<Real>("mol_cutoff")),
    _temperature(getParam<Real>("temperature")),
    _nernst_swap_out_of_basis(getParam<std::vector<std::string>>("nernst_swap_out_of_basis")),
    _nernst_swap_into_basis(getParam<std::vector<std::string>>("nernst_swap_into_basis")),
    _tot_iter(0)
{
  if (_ramp_max_ionic_strength > _max_iter)
    paramError("ramp_max_ionic_strength", "must be less than max_iter");
  if (_nernst_swap_out_of_basis.size() != _nernst_swap_into_basis.size())
    paramError("nernst_swap_out_of_basis", "must be of same size as nernst_swap_into_basis");
}

void
ERSo::output(const ExecFlagType & type)
{
  if (!shouldOutput(type))
    return;

  // do the work: solve the geochemical system
  solveSystem();

  // retrieve information
  const std::vector<Real> basis_molality = _egs.getSolventMassAndFreeMolalityAndMineralMoles();
  std::vector<Real> basis_activity;
  for (unsigned i = 0; i < _num_basis; ++i)
    basis_activity.push_back(_egs.getBasisActivity(i));
  std::vector<Real> basis_act_coef;
  for (unsigned i = 0; i < _num_basis; ++i)
    basis_act_coef.push_back(_egs.getBasisActivityCoefficient(i));
  const std::vector<Real> bulk_moles = _egs.getBulkMoles();
  std::vector<Real> eqm_molality;
  for (unsigned j = 0; j < _num_eqm; ++j)
    eqm_molality.push_back(_egs.getEquilibriumMolality(j));
  std::vector<Real> eqm_act_coef;
  for (unsigned j = 0; j < _num_eqm; ++j)
    eqm_act_coef.push_back(_egs.getEquilibriumActivityCoefficient(j));
  const std::vector<Real> eqm_SI = _egs.getSaturationIndices();

  unsigned ind = 0;
  _console << std::setprecision(_precision);

  _console << "\nSummary:\n";

  _console << "Total number of iterations required = " << _tot_iter << "\n";
  _console << "Error in calculation = " << _abs_residual << "mol\n";
  _console << "Charge of solution = " << _egs.getTotalCharge() << "mol";
  _console << " (charge-balance species = "
           << _mgd.basis_species_name[_egs.getChargeBalanceBasisIndex()] << ")\n";

  _console << "Mass of solvent water = " << basis_molality[0] << "kg\n";

  Real mass = basis_molality[0];
  for (unsigned i = 1; i < _num_basis; ++i) // do not loop over water
    mass += bulk_moles[i] * _mgd.basis_species_molecular_weight[i] / 1000.0;
  _console << "Mass of aqueous solution = " << mass << "kg\n";

  // Output the aqueous solution pH, if relevant
  if (_mgd.basis_species_index.count("H+"))
    _console << "pH = " << -std::log10(basis_activity[_mgd.basis_species_index.at("H+")]) << "\n";
  if (_mgd.eqm_species_index.count("H+"))
    _console << "pH = "
             << -std::log10(eqm_molality[_mgd.eqm_species_index.at("H+")] *
                            eqm_act_coef[_mgd.eqm_species_index.at("H+")])
             << "\n";

  // Output the aqueous solution pe, if relevant
  Real pe = 0.0;
  const bool pe_defined =
      (_mgd.basis_species_index.count("e-") == 1) || (_mgd.eqm_species_index.count("e-") == 1);
  if (_mgd.basis_species_index.count("e-"))
    pe = -std::log10(basis_activity[_mgd.basis_species_index.at("e-")]);
  if (_mgd.eqm_species_index.count("e-"))
    pe = -std::log10(eqm_molality[_mgd.eqm_species_index.at("e-")] *
                     eqm_act_coef[_mgd.eqm_species_index.at("e-")]);
  if (pe_defined)
    _console << "pe = " << pe << "\n";

  // Output ionic strengths
  _console << "Ionic strength = " << _egs.getIonicStrength() << "mol/kg(solvent water)\n";
  _console << "Stoichiometric ionic strength = " << _egs.getStoichiometricIonicStrength()
           << "mol/kg(solvent water)\n";

  // Output activity of water
  _console << "Activity of water = " << basis_activity[0] << "\n";

  // Output the basis species information, sorted by molality
  std::vector<unsigned> basis_order(_num_basis);
  ind = 0;
  std::iota(basis_order.begin(), basis_order.end(), ind++);
  std::sort(basis_order.begin(), basis_order.end(), [&](int i, int j) {
    return basis_molality[i] > basis_molality[j];
  });
  _console << "\nBasis Species:\n";
  for (const auto & i : basis_order)
    if (i == 0 || _mgd.basis_species_gas[i])
      continue;
    else
    {
      _console << _mgd.basis_species_name[i] << ";  bulk_moles = " << bulk_moles[i]
               << "mol;  bulk_conc = "
               << bulk_moles[i] * _mgd.basis_species_molecular_weight[i] * 1000.0 / mass
               << "mg/kg(solution);  molality = " << basis_molality[i]
               << "mol/kg(solvent water);  free_conc = "
               << basis_molality[i] * _mgd.basis_species_molecular_weight[i] * 1000.0
               << "mg/kg(solvent water)";
      if (_mgd.basis_species_mineral[i])
        _console << "\n";
      else
        _console << ";  act_coeff = " << basis_act_coef[i]
                 << ";  log10(a) = " << std::log10(basis_activity[i]) << "\n";
    }
  for (unsigned i = 0; i < _num_basis; ++i)
    if (_mgd.basis_species_gas[i])
      _console << _mgd.basis_species_name[i] << ";  fugacity = " << basis_activity[i] << "\n";

  // Output the equilibrium species info, sorted by molality
  std::vector<unsigned> eqm_order(_num_eqm);
  ind = 0;
  std::iota(eqm_order.begin(), eqm_order.end(), ind++);
  std::sort(eqm_order.begin(), eqm_order.end(), [&](int i, int j) {
    return eqm_molality[i] > eqm_molality[j];
  });
  _console << "\nEquilibrium Species:\n";
  for (const auto & i : eqm_order)
    if (eqm_molality[i] <= _mol_cutoff)
      break;
    else if (_mgd.eqm_species_gas[i])
      continue;
    else
      _console << _mgd.eqm_species_name[i] << ";  molality = " << eqm_molality[i]
               << "mol/kg(solvent water);  free_conc = "
               << eqm_molality[i] * _mgd.eqm_species_molecular_weight[i] * 1000.0
               << "mg/kg(solvent water);  act_coeff = " << eqm_act_coef[i]
               << ";  log10(a) = " << std::log10(eqm_molality[i] * eqm_act_coef[i]) << ";  "
               << _mgd.eqm_species_name[i] << " = "
               << GeochemistryFormattedOutput::reaction(
                      _mgd.eqm_stoichiometry, i, _mgd.basis_species_name, _stoi_tol, _precision)
               << ";  log10K = " << _egs.getLog10K(i) << "\n";
  for (unsigned i = 0; i < _num_eqm; ++i)
    if (_mgd.eqm_species_gas[i])
      _console << _mgd.eqm_species_name[i]
               << ";  act_coeff = " << _egs.getEquilibriumActivityCoefficient(i) << ";  "
               << _mgd.eqm_species_name[i] << " = "
               << GeochemistryFormattedOutput::reaction(
                      _mgd.eqm_stoichiometry, i, _mgd.basis_species_name, _stoi_tol, _precision)
               << ";  log10K = " << _egs.getLog10K(i) << "\n";

  // Output the mineral info, sorted by saturation indices
  std::vector<unsigned> mineral_order(_num_eqm);
  ind = 0;
  std::iota(mineral_order.begin(), mineral_order.end(), ind++);
  std::sort(mineral_order.begin(), mineral_order.end(), [&](int i, int j) {
    return eqm_SI[i] > eqm_SI[j];
  });
  _console << "\nMinerals:\n";
  for (const auto & i : mineral_order)
    if (_mgd.eqm_species_mineral[i])
      _console << _mgd.eqm_species_name[i] << " = "
               << GeochemistryFormattedOutput::reaction(
                      _mgd.eqm_stoichiometry, i, _mgd.basis_species_name, _stoi_tol, _precision)
               << ";  log10K = " << _egs.getLog10K(i) << ";  SI = " << eqm_SI[i] << "\n";

  // Output the Nernst potentials, if relevant
  const Real prefactor = -GeochemistryConstants::LOGTEN * GeochemistryConstants::GAS_CONSTANT *
                         (_temperature + GeochemistryConstants::CELSIUS_TO_KELVIN) /
                         GeochemistryConstants::FARADAY;
  _console << "\nNernst potentials:\n";
  if (pe_defined)
    _console << "e- = 0.5*H20 - 0.25*O2(aq) - 1*H+;  Eh = " << -prefactor * pe << "V\n";
  performNernstSwaps();
  // activities might have changed due to the swaps
  basis_activity.resize(0);
  for (unsigned i = 0; i < _num_basis; ++i)
    basis_activity.push_back(_egs.getBasisActivity(i));
  if (_mgd.redox_lhs == _original_redox_lhs)
    for (unsigned red = 0; red < _mgd.redox_stoichiometry.m(); ++red)
      _console << _mgd.redox_lhs << " = "
               << GeochemistryFormattedOutput::reaction(
                      _mgd.redox_stoichiometry, red, _mgd.basis_species_name, _stoi_tol, _precision)
               << ";  Eh = "
               << prefactor * (_egs.log10RedoxActivityProduct(red) - _egs.getRedoxLog10K(red))
               << "V\n";
}

void
ERSo::performNernstSwaps()
{
  // attempt the swaps specified by the user.  It is not an error if these are impossible
  for (unsigned sw = 0; sw < _nernst_swap_out_of_basis.size(); ++sw)
    if (_mgd.basis_species_index.count(_nernst_swap_out_of_basis[sw]) == 1 &&
        _mgd.eqm_species_index.count(_nernst_swap_into_basis[sw]) == 1)
      try
      {
        _egs.performSwap(_mgd.basis_species_index.at(_nernst_swap_out_of_basis[sw]),
                         _mgd.eqm_species_index.at(_nernst_swap_into_basis[sw]));
      }
      catch (const MooseException & e)
      {
        // it is not an error to be unable to make the swap, so just continue
        mooseWarning("Swapping ",
                     _nernst_swap_out_of_basis[sw],
                     " and ",
                     _nernst_swap_into_basis[sw],
                     ": ",
                     e.what());
      }
  if (_mgd.redox_lhs != _original_redox_lhs &&
      _mgd.basis_species_index.count(_original_redox_lhs) == 1 &&
      _mgd.eqm_species_index.count(_mgd.redox_lhs) == 1)
    try
    {
      // at some stage the original redox left-hand side must have been put into the basis, so let's
      // try to take it out and replace with _mgd.redox_lhs
      _egs.performSwap(_mgd.basis_species_index.at(_original_redox_lhs),
                       _mgd.eqm_species_index.at(_mgd.redox_lhs));
    }
    catch (const MooseException & e)
    {
      // it is not an error to be unable to make the swap, so just continue
    }
}

void
ERSo::solveSystem()
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
      if (_egs.alterChargeBalanceSpecies(_swap_threshold))
        _console << "Changed change balance species to "
                 << _mgd.basis_species_name[_egs.getChargeBalanceBasisIndex()] << std::endl;
      _abs_residual = computeResidual(_residual);
      if (_verbose)
        _console << "iter = " << iter << " |R| = " << _abs_residual << std::endl;
    }

    _egs.enforceChargeBalance(); // just to get _bulk_moles of the charge-balance species correct.
                                 // This is not used within EquilibriumGeochemicalSystem, but is
                                 // used in the output
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

Real
ERSo::computeResidual(DenseVector<Real> & residual) const
{
  for (unsigned a_ind = 0; a_ind < _num_in_algebraic_system; ++a_ind)
    residual(a_ind) = _egs.getResidualComponent(a_ind);
  return residual.l1_norm();
}

void
ERSo::solveAndUnderrelax(DenseMatrix<Real> & jacobian, DenseVector<Real> & new_mol) const
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
ERSo::mineralSwapNeeded(unsigned & swap_out_of_basis, unsigned & swap_into_basis) const
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
          _console << "Mineral " << _mgd.eqm_species_name[j]
                   << " supersaturated.  Swapping it with basis species "
                   << _mgd.basis_species_name[swap_out_of_basis] << std::endl;
          break;
        }
  return swap_needed;
}

bool
ERSo::reduceInitialResidual()
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
