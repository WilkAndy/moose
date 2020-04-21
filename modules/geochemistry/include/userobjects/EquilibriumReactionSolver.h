//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeochemicalModelDefinition.h"
#include "GeochemistrySpeciesSwapper.h"
#include "EquilibriumGeochemicalSystem.h"
#include "DenseMatrix.h"
#include <libmesh/dense_vector.h>

class EquilibriumReactionSolver;

template <>
InputParameters validParams<EquilibriumReactionSolver>();

/**
 * User object for solving reaction systems
 */
class EquilibriumReactionSolver : public GeneralUserObject
{
public:
  static InputParameters validParams();

  EquilibriumReactionSolver(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;

protected:
  /// the database of the user-defined model
  ModelGeochemicalDatabase _mgd;
  /// number of basis species
  const unsigned _num_basis;
  /// number of equilibrium species
  const unsigned _num_eqm;
  /// swapper that can swap species into and out of the basis
  GeochemistrySpeciesSwapper _swapper;
  const std::vector<std::string> _swap_out;
  const std::vector<std::string> _swap_in;
  /// number of unknowns in the algebraic system
  unsigned _num_in_algebraic_system;
  /// if _in_algebraic_system[i] == true then we need to solve for the i^th basis-species molality
  std::vector<bool> _in_algebraic_system;
  /// _algebraic_index[i] = index in the algebraic system of the basis species i.  _basis_index[_algebraic_index[i]] = i
  std::vector<unsigned> _algebraic_index;
  /// _basis_index[i] = index in the basis of the algebraic system species i, for i<num_in_algebraic_system.  _basis_index[_algebraic_index[i]] == i
  std::vector<unsigned> _basis_index;
  /// mass of solvent water
  Real _nw;
  /// _basis_molality[i] = molality of basis species i, except _molality[0] (water) which is meaningless
  std::vector<Real> _basis_molality;
  /// _bulk_moles(i) = number of moles (bulk composition) of basis species i
  DenseVector<Real> _bulk_moles;
  /// _eqm_molality[i] = molality of the equilibrium species i
  std::vector<Real> _eqm_molality;
  /// the residual of the algebraic system
  DenseVector<Real> _residual;
  /// L1 norm of residual
  Real _abs_residual;
  /// jacobian of the algebraic system
  DenseMatrix<Real> _jacobian;
  /// _dmj_dmi(j, i) = d(eqm_molality[j])/d(basis_molality[i])
  DenseMatrix<Real> _dmj_dmi;
  /// solution of _jacobian * _dmol = _residual
  DenseVector<Real> _dmol;
  const unsigned _precision;
  const Real _temperature;
  const std::vector<std::string> _activity_species;
  const std::vector<Real> _activity_values;
  /// whether activity[basis_species_i] is known
  std::vector<bool> _basis_activity_known;
  /// values of activity coefficient or activity (for water) or fugacity (for gases)
  std::vector<Real> _basis_activity;
  /// values of activity coefficient for equilibrium species
  std::vector<Real> _eqm_activity;
  /// equilibrium constant of the equilibrium species
  std::vector<Real> _eqm_log10K;
  /// values of the saturation index - for non-minerals this is zero
  std::vector<Real> _eqm_SI;
  /// Debye-Huckel parameter
  Real _dhA;
  /// Debye-Huckel parameter
  Real _dhB;
  /// Debye-Huckel parameter
  Real _dhBdot;
  /// Debye-Huckel parameter
  Real _dha;
  /// Debye-Huckel parameter
  Real _dhb;
  /// Debye-Huckel parameter
  Real _dhc;
  /// Debye-Huckel parameter
  Real _dhd;
  /// Debye-Huckel parameter
  Real _dhatilde;
  /// Debye-Huckel parameter
  Real _dhbtilde;
  /// Debye-Huckel parameter
  Real _dhctilde;
  /// Debye-Huckel parameter
  Real _dhdtilde;
  const std::string _output_species;
  const std::string _charge_balance_species;
  unsigned _charge_balance_index;
  const bool _mass_solvent_water_provided;
  const Real _mass_solvent_water;
  const bool _moles_bulk_water_provided;
  /// Moles of bulk water, if provided by the user.  This potentially gets modified after swaps due to mineral dissolution/precipitation
  Real _moles_bulk_water;
  /// bulk number of moles provided by user.  This potentially gets modified after swaps due to mineral dissolution/precipitation
  std::vector<std::string> _bulk_moles_species;
  /// bulk number of moles provided by user.  This potentially gets modified after swaps due to mineral dissolution/precipitation
  std::vector<Real> _bulk_moles_values;
  const std::vector<std::string> _free_molality_species;
  const std::vector<Real> _free_molality_values;
  const std::vector<std::string> _free_moles_mineral_species;
  const std::vector<Real> _free_moles_mineral_values;
  const std::vector<std::string> _prevent_precipitation;
  /// If the residual of the algebraic system falls below this value, the Newton process has converged
  const Real _abs_tol;
  /// If the residual of the algebraic system falls below this value times the initial residual, the Newton process has converged
  const Real _rel_tol;
  /// _res0_times_rel = _rel_tol * initial residual
  Real _res0_times_rel;
  /// Maximum number of Newton iterations allowed to solve algebraic system
  const unsigned _max_iter;
  /// Whether to print iteration residuals, swap information, etc
  const bool _verbose;
  /// Total number of Newton iterations used (including when minerals dissolve or precipitate)
  unsigned _tot_iter;

private:
  /// provide a list of the species of interest
  std::vector<std::string> speciesOfInterest() const;

  /**
   * Using the provided value of temperature, build:
   * _eqm_log10K for each eqm species
   * _dhA, _dhB, _dhBdot
   * _dha, _dhb, _dhc, _dhd
   * _dhatilde, _dhbtilde, _dhctilde, _dhdtilde
   */
  void buildTemperatureDependentQuantities();

  /// Builds in_algebraic_system, algebraic_index and basis_index, and sets num_in_algebraic_system appropriately
  void buildAlgebraicInfo(std::vector<bool> & in_algebraic_system,
                          unsigned & num_in_algebraic_system,
                          std::vector<unsigned> & algebraic_index,
                          std::vector<unsigned> & basis_index) const;

  /**
   * Initialise the bulk moles and molality for the basis after user-defined swaps
   * _nw
   * _bulk_moles (fully populated)
   * _basis_molality (fully populated)
   */
  void initBulkAndFree(Real & nw,
                       DenseVector<Real> & bulk_moles,
                       std::vector<Real> & basis_molality) const;

  /**
   * Fully populates basis_activity_known, which is true if activity has been set by the user, or
   * the basis species is a mineral.
   * Populates only those slots in _basis_activity for which basis_activity_known == true
   */
  void buildKnownBasisActivities(std::vector<bool> & basis_activity_known,
                                 std::vector<Real> & basis_activity) const;

  /**
   * Build the activites for the current iteration of the solve.
   * Does not modify any activities for which _basis_activity_known == true.
   * Using _basis_molality, _eqm_molality, compute the ionic strength
   * Then, using the Debye-Huckel parameters (_dhA, _dhB, etc) compute _basis_activity and
   * _eqm_activity
   */
  void buildActivities(std::vector<Real> & basis_activity, std::vector<Real> & eqm_activity) const;

  /**
   * Computes the molalities of all species at equilibrium with the basis
   * Using _basis_activity, eqm_activity, _eqm_log10K, compute _eqm_molality
   */
  void computeEqmMolalities(std::vector<Real> & eqm_molality) const;

  /**
   * Builds the residual of the algebraic system
   * @return the L1 norm of residual
   */
  Real computeResidual(DenseVector<Real> & residual) const;

  /// Builds the jacobian of the algebraic system
  void computeJacobian(DenseMatrix<Real> & dmj_dmi, DenseMatrix<Real> & jacobian) const;

  /// Solves _jacobian * _dmol = _residual for _dmol, then performs an underrelaxation
  void solveAndUnderrelax(DenseMatrix<Real> & jacobian, DenseVector<Real> & dmol) const;

  /// Update the degrees of freedom using _dmol
  void updateWithdMol(Real & nw, std::vector<Real> & basis_molality) const;

  /// Compute bulk number of moles for all basis species
  void computeBulkMoles(DenseVector<Real> & bulk_moles) const;

  /// Change the bulk concentration of the charge-balance species to enforce electrical neutrality, and record the information in bulk_moles_values
  void enforceChargeNeutrality(DenseVector<Real> & bulk_moles,
                               std::vector<Real> & bulk_moles_values) const;

  /// output the results to the console
  void outputResults() const;

  /// calculate mineral saturation indices and put them in the si vector
  void saturationIndices(std::vector<Real> & si) const;

  /**
   * Calculate the activity product for the given equilibrium species
   * @param eqm_index the index in _mgd for which the activity product is required
   * @return the activity product
   */
  Real activityProduct(unsigned eqm_index) const;

  /**
   * Check that the charge-balance species is in the basis and that it has nonzero charge
   * @return the basis index of the charge-balance species in _mgd
   */
  unsigned checkChargeBalanceSpecies() const;

  /// check enough information has been provided regarding bulk number of moles, free molality, fugacity, etc
  void checkICs() const;

  /// computes the number of free moles of basis minerals
  void computeFreeMineralMoles(std::vector<Real> & basis_molality) const;

  /**
   * Check if a basis swap is needed.  It is needed if:
   * - the free number of moles of a basis mineral is negative
   * - the saturation index of an equilibrium mineral is positive (and it is not in the
   * prevent_precipitation list)
   * @param swap_out_of_basis the index of the species in the basis that will be removed from the
   * basis
   * @param swap_into_basis the index of the equilibrium mineneral that will be added to the basis
   * @return true if a swap is needed, false otherwise
   */
  bool mineralSwapNeeded(unsigned & swap_out_of_basis, unsigned & swap_into_basis) const;

  /**
   * Record the values of _bulk_moles into moles_bulk_water, bulk_moles_species and
   * bulk_moles_values iff the free molality or free number of moles has not been provided for the
   * basis components
   * @param moles_bulk_water the number of moles of bulk water
   * @param bulk_moles_species the names of the basis species which have bulk number of moles in
   * bulk_moles_values
   * @param bulk_moles_values the number of bulk moles
   */
  void recordBulkMoles(Real & moles_bulk_water,
                       std::vector<std::string> & bulk_moles_species,
                       std::vector<Real> & bulk_moles_values) const;
};
