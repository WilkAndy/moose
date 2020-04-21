//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PertinentGeochemicalSystem.h"

/**
 * Calculators to compute ionic strength and stoichiometric ionic strength
 */

namespace GeochemistryIonicStrength
{
/// Ionic strength
Real ionicStrength(const ModelGeochemicalDatabase & mgd,
                   const std::vector<Real> & basis_species_molality,
                   const std::vector<Real> & eqm_species_molality,
                   const std::vector<Real> & kin_species_molality);

/// Stoichiometric ionic strength
Real stoichiometricIonicStrength(const ModelGeochemicalDatabase & mgd,
                                 const std::vector<Real> & basis_species_molality,
                                 const std::vector<Real> & eqm_species_molality,
                                 const std::vector<Real> & kin_species_molality);
}
