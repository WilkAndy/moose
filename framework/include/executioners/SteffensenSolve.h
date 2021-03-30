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

class SteffensenSolve : public IterativeMultiAppSolve
{
public:
  SteffensenSolve(Executioner * ex);

  static InputParameters validParams();

  /// Set relaxation postprocessors for the current solve as a SubApp
  virtual void setMultiAppRelaxationPostprocessors(const std::vector<std::string> & pps)
  {
    _secondary_transformed_pps = pps;
    _old_secondary_transformed_pps_values.resize(pps.size());
    _older_secondary_transformed_pps_values.resize(pps.size());
  }

private:
  /// Save the variable and postprocessor values as a SubApp
  virtual void savePreviousValuesAsSubApp();

  /// Whether to use the coupling algorithm (relaxed Picard, Steffensen, ...) instead of Picard
  virtual bool useCouplingUpdateAlgorithm();

  /// Save the previous values for the variables
  virtual void savePreviousValuesAsMainApp();

  /// Compute the new value of the coupling postprocessors based on the coupling algorithm selected
  virtual void updatePostprocessorsAsMainApp();

  /// Compute the new value variable values based on the coupling algorithm selected
  virtual void updateVariablesAsMainApp(const std::set<dof_id_type> & transformed_dofs);

  /// Update variables and postprocessors as a SubApp
  virtual void updateAsSubApp(const std::set<dof_id_type> & secondary_transformed_dofs);

  /// Print the convergence history of the coupling, at every coupling iteration
  virtual void printCouplingConvergenceHistory();

  /// Values of the relaxed postprocessors from two iterations prior
  std::vector<PostprocessorValue> _older_transformed_pps_values;

  /// Values of the postprocessors to relax outside of coupling iterations from two iterations prior (used as a subapp)
  std::vector<PostprocessorValue> _older_secondary_transformed_pps_values;
};
