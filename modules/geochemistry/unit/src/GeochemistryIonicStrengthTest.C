//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "gtest/gtest.h"

#include "GeochemistryIonicStrength.h"

TEST(GeochemistryIonicStrengthTest, sizeExceptions)
{
  GeochemicalDatabaseReader database("database/moose_testdb.json");

  // this should pass: CH4(g) depends on the redox couple CH4(aq) which in turn depends on the basis
  // species provided
  PertinentGeochemicalSystem model(
      database, {"H2O", "Ca++", "H+", "HCO3-", "O2(aq)"}, {"Calcite"}, {"CH4(g)"}, {}, {}, {});

  ModelGeochemicalDatabase mgd = model.modelGeochemicalDatabase();
}

TEST(GeochemistryIonicStrengthTest, ionicStrength)
{
  GeochemicalDatabaseReader database("database/moose_testdb.json");

  // The following system has secondary species: CO2(aq), CO3--, CaCO3, CaOH+, OH-, (O-phth)--,
  // >(s)FeO-, e-
  PertinentGeochemicalSystem model(database,
                                   {"H2O", "H+", "HCO3-", "O2(aq)", "Ca++", ">(s)FeOH"},
                                   {"Calcite"},
                                   {},
                                   {"Calcite_asdf"},
                                   {"CH4(aq)"},
                                   {">(s)FeOCa+"});

  ModelGeochemicalDatabase mgd = model.modelGeochemicalDatabase();

  Real gold_ionic_str = 0.0;

  std::vector<Real> basis_m(6);
  basis_m[mgd.basis_species_index["H2O"]] = 1.0;
  basis_m[mgd.basis_species_index["H+"]] = 2.0;
  gold_ionic_str += 2.0;
  basis_m[mgd.basis_species_index["HCO3-"]] = 3.0;
  gold_ionic_str += 3.0;
  basis_m[mgd.basis_species_index["O2(aq)"]] = 4.0;
  basis_m[mgd.basis_species_index["Ca++"]] = 5.0;
  gold_ionic_str += 4 * 5.0;
  basis_m[mgd.basis_species_index[">(s)FeOH"]] = 6.0;

  std::vector<Real> eqm_m(9);
  eqm_m[mgd.eqm_species_index["CO2(aq)"]] = 1.1;
  eqm_m[mgd.eqm_species_index["CO3--"]] = 2.2;
  gold_ionic_str += 4 * 2.2;
  eqm_m[mgd.eqm_species_index["CaCO3"]] = 3.3;
  eqm_m[mgd.eqm_species_index["CaOH+"]] = 4.4;
  gold_ionic_str += 4.4;
  eqm_m[mgd.eqm_species_index["OH-"]] = 5.5;
  gold_ionic_str += 5.5;
  eqm_m[mgd.eqm_species_index["(O-phth)--"]] = 6.6;
  gold_ionic_str += 4 * 6.6;
  eqm_m[mgd.eqm_species_index[">(s)FeO-"]] = 7.7;
  gold_ionic_str += 7.7;
  eqm_m[mgd.eqm_species_index["e-"]] = 8.8;
  gold_ionic_str += 8.8;
  eqm_m[mgd.eqm_species_index["Calcite"]] = 9.9;

  std::vector<Real> kin_m(3);
  kin_m[mgd.kin_species_index["Calcite_asdf"]] = -1.1;
  kin_m[mgd.kin_species_index["CH4(aq)"]] = -2.2;
  kin_m[mgd.kin_species_index[">(s)FeOCa+"]] = -3.3;
  gold_ionic_str += -3.3;

  gold_ionic_str *= 0.5;
  const Real ionic_str = GeochemistryIonicStrength::ionicStrength(mgd, basis_m, eqm_m, kin_m);

  ASSERT_NEAR(ionic_str, gold_ionic_str, 1E-9);
}

TEST(GeochemistryIonicStrengthTest, stoichiometricIonicStrength)
{
  GeochemicalDatabaseReader database("database/moose_testdb.json");

  // The following system has secondary species: CO2(aq), CO3--, CaCO3, CaOH+, OH-, (O-phth)--,
  // >(s)FeO-, e-
  PertinentGeochemicalSystem model(database,
                                   {"H2O", "H+", "HCO3-", "O2(aq)", "Ca++", ">(s)FeOH"},
                                   {"Calcite"},
                                   {},
                                   {"Calcite_asdf"},
                                   {"CH4(aq)"},
                                   {">(s)FeOCa+"});

  ModelGeochemicalDatabase mgd = model.modelGeochemicalDatabase();

  Real gold_ionic_str = 0.0;

  std::vector<Real> basis_m(6);
  basis_m[mgd.basis_species_index["H2O"]] = 1.0;
  basis_m[mgd.basis_species_index["H+"]] = 2.0;
  gold_ionic_str += 2.0;
  basis_m[mgd.basis_species_index["HCO3-"]] = 3.0;
  gold_ionic_str += 3.0;
  basis_m[mgd.basis_species_index["O2(aq)"]] = 4.0;
  basis_m[mgd.basis_species_index["Ca++"]] = 5.0;
  gold_ionic_str += 5.0 * 4;
  basis_m[mgd.basis_species_index[">(s)FeOH"]] = 6.0;

  std::vector<Real> eqm_m(9);
  eqm_m[mgd.eqm_species_index["CO2(aq)"]] = 1.1;
  gold_ionic_str += 1.1 * (1 + 1);
  eqm_m[mgd.eqm_species_index["CO3--"]] = 2.2;
  gold_ionic_str += 2.2 * 4;
  eqm_m[mgd.eqm_species_index["CaCO3"]] = 3.3;
  gold_ionic_str += 3.3 * (4 + 1 - 1);
  eqm_m[mgd.eqm_species_index["CaOH+"]] = 4.4;
  gold_ionic_str += 4.4 * (1);
  eqm_m[mgd.eqm_species_index["OH-"]] = 5.5;
  gold_ionic_str += 5.5 * (1);
  eqm_m[mgd.eqm_species_index["(O-phth)--"]] = 6.6;
  gold_ionic_str += 6.6 * (4);
  eqm_m[mgd.eqm_species_index[">(s)FeO-"]] = 7.7;
  gold_ionic_str += 7.7 * (1);
  eqm_m[mgd.eqm_species_index["e-"]] = 8.8;
  gold_ionic_str += 8.8 * (1);
  eqm_m[mgd.eqm_species_index["Calcite"]] = 9.9;
  gold_ionic_str += 9.9 * (4 + 1 - 1);

  std::vector<Real> kin_m(3);
  kin_m[mgd.kin_species_index["Calcite_asdf"]] = -1.1;
  gold_ionic_str += -1.1 * (2 * 4 + 1 - 1);
  kin_m[mgd.kin_species_index["CH4(aq)"]] = -2.2;
  gold_ionic_str += -2.2 * (1 + 1);
  kin_m[mgd.kin_species_index[">(s)FeOCa+"]] = -3.3;
  gold_ionic_str += -3.3 * (1);

  gold_ionic_str *= 0.5;
  const Real ionic_str =
      GeochemistryIonicStrength::stoichiometricIonicStrength(mgd, basis_m, eqm_m, kin_m);

  ASSERT_NEAR(ionic_str, gold_ionic_str, 1E-9);
}
