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
#include "EquilibriumGeochemicalSystem.h"
#include "EquilibriumGeochemicalSolver.h"
#include "Output.h"
#include "UserObjectInterface.h"
#include "DenseMatrix.h"
#include <libmesh/dense_vector.h>

/**
 * Solves a reaction system for molalities and outputs results
 */
class ERSo : public Output, public UserObjectInterface
{
public:
  static InputParameters validParams();

  ERSo(const InputParameters & parameters);

protected:
  virtual void output(const ExecFlagType & type) override;

  /// The equilibrium geochemical system that holds all the molalities, activities, etc
  EquilibriumGeochemicalSystem _egs;
  /// The solver
  EquilibriumGeochemicalSolver _solver;
  /// The left-hand side specified in the original model definition for redox half reactions
  const std::string _original_redox_lhs;
  /// number of basis species
  const unsigned _num_basis;
  /// number of equilibrium species
  const unsigned _num_eqm;
  /// reference to the underlying ModelGeochemicalDatabase
  const ModelGeochemicalDatabase & _mgd;
  /// number of molalities in the algebraic system
  unsigned _num_basis_in_algebraic_system;
  /// number of unknowns in the algebraic system
  unsigned _num_in_algebraic_system;
  /// residual of the algebraic system we wish to solve
  DenseVector<Real> _residual;
  /// L1 norm of residual
  Real _abs_residual;
  /// jacobian of the algebraic system
  DenseMatrix<Real> _jacobian;
  /// the new molality after finding the solution of _jacobian * neg_change_mol = _residual
  DenseVector<Real> _new_mol;
  /// precision of output
  const unsigned _precision;
  /// prevent precipitation of these minerals
  const std::vector<std::string> _prevent_precipitation;
  /// If the residual of the algebraic system falls below this value, the Newton process has converged
  const Real _abs_tol;
  /// If the residual of the algebraic system falls below this value times the initial residual, the Newton process has converged
  const Real _rel_tol;
  /// _res0_times_rel = _rel_tol * initial residual
  Real _res0_times_rel;
  /// Tolerance on stoichiometric coefficients before they are deemed to be zero
  const Real _stoi_tol;
  /// Maximum number of Newton iterations allowed to solve algebraic system
  const unsigned _max_iter;
  /// Whether to print iteration residuals, swap information, etc
  const bool _verbose;
  /// Iteratively reduce the initial-guess molalities so that the initial residual for the Newton process is less than:
  const Real _max_initial_residual;
  /// If a basis molality < swap_threshold, we attempt to swap it out of the basis
  const Real _swap_threshold;
  /// Maximum ionic strength
  const Real _max_ionic_strength;
  /// Number of iterations over which to increase the maximum ionic strength to _max_ionic_strength
  const unsigned _ramp_max_ionic_strength;
  /// Species with molalities less than mol_cutoff will not be outputted
  const Real _mol_cutoff;
  /// Temperature specified by user
  const Real _temperature;
  /// Species to swap out of basis prior to outputting the Nernst potentials
  const std::vector<std::string> _nernst_swap_out_of_basis;
  /// Species to swap into basis prior to outputting the Nernst potentials
  const std::vector<std::string> _nernst_swap_into_basis;
  /// Total number of Newton iterations used (including when minerals dissolve or precipitate)
  unsigned _tot_iter;

private:
  /**
   * Perform the swaps specified by the user, prior to outputting Nernst redox information.  In
   * addition, use swaps to attempt to make the left-hand side of the redox half reaction equal to
   * _original_redox_lhs
   */
  void performNernstSwaps();
};
