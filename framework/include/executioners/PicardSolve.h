//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "IterativeMultiAppSolve.h"

// System includes
#include <string>

class PicardSolve;

template <>
InputParameters validParams<PicardSolve>();

class PicardSolve : public IterativeMultiAppSolve
{
public:
  PicardSolve(Executioner * ex);

  static InputParameters validParams();

protected:
  /// Save the previous variable and postprocessor values as a SubApp
  virtual void savePreviousValuesAsSubApp() override final;

  /// Whether to use the coupling algorithm (relaxed Picard, Secant, ...) instead of Picard
  virtual bool useCouplingUpdateAlgorithm() override final;

  /// Save the previous variables and postprocessors as the main application
  virtual void savePreviousValuesAsMainApp() override final;

  /// Compute the new value of the coupling postprocessors based on the coupling algorithm selected
  virtual void updatePostprocessorsAsMainApp() override final;

  /// Compute the new variable values based on the coupling algorithm selected
  virtual void updateVariablesAsMainApp(const std::set<dof_id_type> & transformed_dofs) override final;

  /// Update variables and postprocessors as a SubApp
  virtual void updateAsSubApp(const std::set<dof_id_type> & secondary_transformed_dofs) override final;

  /// Print the convergence history of the coupling, at every coupling iteration
  virtual void printCouplingConvergenceHistory() override final;
};
