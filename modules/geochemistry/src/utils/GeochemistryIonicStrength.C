//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GeochemistryIonicStrength.h"
#include "libmesh/utility.h"

namespace GeochemistryIonicStrength
{
Real
ionicStrength(const ModelGeochemicalDatabase & mgd,
              const std::vector<Real> & basis_species_molality,
              const std::vector<Real> & eqm_species_molality,
              const std::vector<Real> & kin_species_molality)
{
  const unsigned num_basis = mgd.basis_species_charge.size();
  const unsigned num_eqm = mgd.eqm_species_charge.size();
  const unsigned num_kin = mgd.kin_species_charge.size();
  mooseAssert(num_basis == basis_species_molality.size(),
              "Ionic strength calculation: Number of basis species in mgd not equal to the size of "
              "basis_species_molality");
  mooseAssert(
      num_eqm == eqm_species_molality.size(),
      "Ionic strength calculation: Number of equilibrium species in mgd not equal to the size of "
      "eqm_species_molality");
  mooseAssert(
      num_kin == kin_species_molality.size(),
      "Ionic strength calculation: Number of kinetic species in mgd not equal to the size of "
      "kin_species_molality");

  Real ionic_strength = 0.0;
  for (unsigned i = 0; i < num_basis; ++i)
    ionic_strength += Utility::pow<2>(mgd.basis_species_charge[i]) * basis_species_molality[i];
  for (unsigned i = 0; i < num_eqm; ++i)
    ionic_strength += Utility::pow<2>(mgd.eqm_species_charge[i]) * eqm_species_molality[i];
  for (unsigned i = 0; i < num_kin; ++i)
    ionic_strength += Utility::pow<2>(mgd.kin_species_charge[i]) * kin_species_molality[i];

  return 0.5 * ionic_strength;
}

Real
stoichiometricIonicStrength(const ModelGeochemicalDatabase & mgd,
                            const std::vector<Real> & basis_species_molality,
                            const std::vector<Real> & eqm_species_molality,
                            const std::vector<Real> & kin_species_molality)
{
  const unsigned num_basis = mgd.basis_species_charge.size();
  const unsigned num_eqm = mgd.eqm_species_charge.size();
  const unsigned num_kin = mgd.kin_species_charge.size();
  mooseAssert(num_basis == basis_species_molality.size(),
              "Stiochiometric ionic strength calculation: Number of basis species in mgd not equal "
              "to the size of basis_species_molality");
  mooseAssert(num_eqm == eqm_species_molality.size(),
              "Stiochiometric ionic strength calculation: Number of equilibrium species in mgd not "
              "equal to the size of eqm_species_molality");
  mooseAssert(num_kin == kin_species_molality.size(),
              "Stoichiometric ionic strength calculation: Number of kinetic species in mgd not "
              "equal to the size of kin_species_molality");

  Real ionic_strength = 0.0;
  for (unsigned i = 0; i < num_basis; ++i)
    ionic_strength += Utility::pow<2>(mgd.basis_species_charge[i]) * basis_species_molality[i];
  for (unsigned i = 0; i < num_eqm; ++i)
    if (mgd.eqm_species_charge[i] != 0.0)
      ionic_strength += Utility::pow<2>(mgd.eqm_species_charge[i]) * eqm_species_molality[i];
    else
    {
      for (unsigned j = 0; j < num_basis; ++j)
        ionic_strength += Utility::pow<2>(mgd.basis_species_charge[j]) * eqm_species_molality[i] *
                          mgd.eqm_stoichiometry(i, j);
    }
  for (unsigned i = 0; i < num_kin; ++i)
    if (mgd.kin_species_charge[i] != 0.0)
      ionic_strength += Utility::pow<2>(mgd.kin_species_charge[i]) * kin_species_molality[i];
    else
    {
      for (unsigned j = 0; j < num_basis; ++j)
        ionic_strength += Utility::pow<2>(mgd.basis_species_charge[j]) * kin_species_molality[i] *
                          mgd.kin_stoichiometry(i, j);
    }

  return 0.5 * ionic_strength;
}
} // namespace GeochemistryIonicStrength
