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
#include "DenseMatrix.h"
#include <libmesh/dense_vector.h>

class ERS1;

template <>
InputParameters validParams<ERS1>();

/**
 * User object for solving reaction systems
 */
class ERS1 : public GeneralUserObject
{
public:
  static InputParameters validParams();

  ERS1(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;

protected:
  /// The equilibrium geochemical system that holds all the molalities, activities, etc
  EquilibriumGeochemicalSystem _egs;
  /// number of basis species
  const unsigned _num_basis;
  /// number of equilibrium species
  const unsigned _num_eqm;
  /// reference to the underlying ModelGeochemicalDatabase
  const ModelGeochemicalDatabase & _mgd;
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
  /// focus output on this species
  const std::string _output_species;
  /// prevent precipitation of these minerals
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
  /// Iteratively reduce the initial-guess molalities so that the initial residual for the Newton process is less than:
  const Real _max_initial_residual;
  /// If a basis molality < swap_threshold, we attempt to swap it out of the basis
  const Real _swap_threshold;
  /// Maximum ionic strength
  const Real _max_ionic_strength;
  /// Number of iterations over which to increase the maximum ionic strength to _max_ionic_strength
  const unsigned _ramp_max_ionic_strength;
  /// Total number of Newton iterations used (including when minerals dissolve or precipitate)
  unsigned _tot_iter;

private:
  /// provide a list of the species of interest
  std::vector<std::string> speciesOfInterest() const;

  /**
   * Builds the residual of the algebraic system
   * @return the L1 norm of residual
   */
  Real computeResidual(DenseVector<Real> & residual) const;

  /**
   * Solves _jacobian * neg_change_mol = _residual for neg_change_mol, then performs an
   * underrelaxation to get new_mol
   * @param jacobian the jacobian of the system
   * @param new_mol upon exit, this will be the new molality values according to the underrelaxed
   * Newton process
   */
  void solveAndUnderrelax(DenseMatrix<Real> & jacobian, DenseVector<Real> & new_mol) const;

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
   * Progressively reduce the initial-guess molalities for the algebraic system to attempt to reduce
   * the residual
   */
  bool reduceInitialResidual();
};
