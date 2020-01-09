# Design of PorousFlow and ChemicalReactions: January 2020

Enhancements are planned for both `PorousFlow` and `ChemicalReactions` during 2020.  A refactor of both modules may be desirable to best incorporate these additions.  This document describes the current architecture of both modules, some of the desired enhancements, and proposes new designs.

## Current `PorousFlow`

### Underlying, unique requirements

There are three pecularities, which are largely unique to `PorousFlow`, that fundamentally impact its current design.

1. `PorousFlow` is [solving](governing_equations.md) versions of the advection equation: $\dot{u} + \nabla\cdot({\mathbf{v}} w) = 0$.  In this equation, $u$, $\mathbf{v}$ and $w$ are all functions of the independent variables, which may be pressure and temperature, for example.  The equations of `PorousFlow` also contain other terms, such as couplings with solid mechanics and sources, but because advection is part of the physics, [numerical stabilization](stabilization.md) is necessary.  The types of stabilization used require that $u$, $\mathbf{v}$, etc, be evaluated at the nodes, not just the quadpoints.

2. Most of `PorousFlow` works for any choice of independent `Variables`.  For instance, the Kernels for advection, such as [PorousFlowAdvectiveFlux.md] are the same for (pressure, temperature) `Variables` as they are for (pressure, saturation, enthalpy, primary concentration) `Variables`, etc.  The reason for this design is that:

- `PorousFlow` is very flexible and can simulate many different scenarios - with arbitrary number of phases, fluid components, isothermal or temperature-dependent, etc - and it makes no sense to create separate classes for every conceivable scenario;
- the choice of `Variables` can greatly impact the speed of numerical convergence in certain simulations.

3. `PorousFlow` computes the exact Jacobian in every case.

### Nodal Materials

For better or worse, `PorousFlow` currently uses "nodal Materials", based on `PorousFlowMaterial` to satisfy the requirements described above.

- The usual `Material` data structures are used, with properties indexed by `_qp`.  However, `_qp` actually refers to a node number.  For example, the code for [PorousFlowSingleComponentFluid.md] contains `_qp` throughout the code, but note that if `_nodal_material==true` the `_porepressure` is actually the `MaterialProperty` called "PorousFlow_porepressure_nodal" that was created by [PorousFlow1PhaseP.md] (for instance) that used `coupledDofValues("porepressure")` and *not* `coupledValue("porepressure")`.

- In nearly all cases, the choice of whether a nodal or (normal) quadpoint material is made automatically by PorousFlow, so that the end user does not need to set the 'at_nodes' parameter for each material. When the input file is parsed, the `PorousFlowAddMaterialAction` checks which type of PorousFlow material (nodal, quadpoint or even both) is required by the other objects (Kernels, BCs, Postprocessors etc), and ensures that the `at_nodes` parameter is set appropriately. This action determines which form of each material is required by checking `PorousFlowDependencies`. This class contains a static list of dependencies for each object that must be manually maintained. 

- Because `Materials` are used, a single class can store multiple properties.  Usually the properties stored are a physical property, such as permeability, which need not be a single real value (permeability is a 2-tensor) as well as its derivatives with respect to the arbitrary `Variables`.

- Usually `_qp=i` corresponds to nodenumber=i, which is fine even if nodenumber=0 is not the nearest node to `_qp=0` (for instance) because this holds consistently throughout the code.

As an aside, there are cases where it is appropriate to demand that the *nearest* quadpoint to each node is used to store properties at that node:

- [PorousFlowPorosity.md] can depend on volumetric strain.  The volumetric strain will be evaluated correctly at the quadpoints (it is not a nodal Material), so the Porosity Material should use the correct strain.  This is mostly important for finite elements where the number of quadpoints do not equal the number of nodes (3D hexahedra on boundaries have 8 nodes, but their boundaries (used in BCs) only have 4 quadpoints; 3D hexahedra involved in DiracKernels have 8 nodes, but can contain any number of quadpoints).
- This impacts Kernels, such as [PorousFlowMassTimeDerivative.md] that has to compute derivatives with respect to strain (gradients of displacement).
- The DiracKernels [PorousFlowPeacemanBorehole.md] and [PorousFlowPolyLineSink.md] require `strain_at_nearest_qp=true` to be set if Porosity depend on volumetric strain, to ensure the nodal Porosity Material uses volumetric strain at the Dirac quadpoints, and can therefore be computed.

### Implementation of PorousFlow Materials

- The input file contains a set of "start" Materials.  These are special because they work for only *one* choice of Variables.  Their purpose is to compute porepressure, temperature, saturation and mass fractions given the `Variables`.  Most commonly, because the `Variables` are indeed porepressure, temperature, etc, they simply take the `Variables` and store their values (nodal values for nodal Materials, quadpoint values for ordinary Materials) in the Material data structures at `_qp`, along with their (rather trivial) derivatives with repect to the `Variables`.  A good example of this is [PorousFlowTemperature.md].  Less trivial cases include [PorousFlow1PhaseP.md] which computes a saturation (and its derivatives) along with passing through the user-supplied porepressure.

- The remaining Materials (relative permeability, fluid density, etc) assume that porepressure, etc, have been computed and recorded in Properties such as `PorousFlow_porepressure_nodal`.  Derivatives with respect to the `Variables` are computed using the chain rule and the derivatives computed by the "start" Materials.

- There is one complication that is invisible to the user.  Because `PorousFlow` works for arbitrary number of phases (and fluid species), all the `Materials`, `Kernels`, `BCs`, etc, must operate on arbitrarily-sized `std::vectors` Material Properties of porepressure, saturation and mass fractions.  They obtain the number of phases and fluid components from the [PorousFlowDictator.md] that must be present in every simulation, and can thus check for input-file inconsistencies.

### PorousFlow Joiners

Properties such as fluid density, fluid viscosity and relative permeability are defined for each phase in the input file.  However, all `Materials`, `Kernels`, etc, operate on `std::vectors` of Material properties corresponding to each phase.

For instance, the user might define relative permeability for phase 0 using [PorousFlowRelativePermeabilityVG.md] that will result in Property `PorousFlow_relative_permeability_nodal0`, and for phase 1 using [PorousFlowRelativePermeabilityConst.md] that will result in Property `PorousFlow_relative_permeability_nodal1`.  However, the [PorousFlowAdvectiveFlux.md] Kernel asks for `getMaterialProperty<std::vector<Real>>("PorousFlow_relative_permeability_nodal")`, so the two phase relative permeabilities have to be joined into a `std::vector`.

The individual phase properties are automatically combined into a `std::vector` material property by `PorousFlowJoiner` materials that are added by the action system. This happens in the background, so that the user does not have to worry about this step. When this input file is parsed, the `PorousFlowAddMaterialJoiner` action checks if
any phase-specific materials (such as fluid properties or relative permeability) are present. If they are, then this action adds the appropriate `PorousFlowJoiner` material to ensure that the correctly sized `std::vector` of properties is formed for use by other objects.

### Kernels and BCs

All the `Kernels` and `BCs` retrieve MaterialProperties corresponding to permeability, fluid density, viscosity, porepressure, temeratpure, etc, as well as all of the derivative information, and assemble the relevant residual and Jacobian entries.  The fully-upwinded and mass-lumped `Kernels` typically use nodal-Material information, while the others use standard quad-point Material information.  The `BCs` are all fully-upwinded, so they always retrieve nodal-Material information.

To do [mass lumping](mass_lumping.md), porosity is lumped to the nodes.  [porosity.md] is tricky in many regards:

- its in-situ value is a standard Material value, evaluated at the quadpoints of elements, because it is a property of the porous skeleton
- it can depend on volumetric strain, which could in principal vary throughout each element
- it can depend on porepressure and temperature, which are nodal
- it can have a cyclic dependence on precipitate volume fraction.

Matrix internal energy could be similarly complicated, but at present it is simply a standard Material value of heat capacity multiplied by a nodal value of temperature, which is lumped to the nodes.

Permeability and thermal conductivity are never upwinded, even if they depend on porepressure, temperature, etc, so only their quadpoint values are ever computed and used.

The TVD advection Kernels use the AdvectiveFluxCalculator `UserObjects`.

### AuxKernels

[PorousFlowPropertyAux.md] retrieves only usual "quadpoint" Material properties and places them in constant, monomial AuxVariables.

[PorousFlowDarcyVelocityComponent.md] retrieves only usual "quadpoint" Material properties (gradient of porepressure, fluid density, viscosity, relative permeability, permeability) to construct the Darcy velocity.

### DiracKernels

These fully upwind, so use "nodal" Material values of internal energy, enthalpy, mass fractions, relative permeability, viscosity, density, etc.  However, [PorousFlowPeacemanBorehole.md] and [PorousFlowPolyLineSink.md] use the quadpoint values of temperature and porepressure to initially define the flows (before potential modification with the fully-upwinded values of relative permeability, etc), since their values at the current Dirac quadpoint dictate the flows.

### Postprocessors

These are fully-upwinded in the same way as `Kernels` and `BCs` to avoid discrepancies: they use "nodal" Material values of fluid density, saturation, mass fraction, matrix internal energy and fluid internal energy.

### UserObjects

These are mostly fairly trivial.  Relevant to the current discussion, the TVD [advective flux calculators](kt.md) use "nodal" Material values of fluid density, viscosity, enthalpy, mass fraction and relative permeability, and, to define the basic fluid velocity, the usual quadpoint Material values of permeability, fluid density and gradient of porepressure.

### Chemical reactions currently in PorousFlow

PorousFlow can currently model aqueous equilibrium reactions as well as precipitation-dissolution reactions ([tutorial_07.md] and [tutorial_13.md]) .  This functionality is similar to the current `ChemicalReactions` module.  There main differences are:

1. The molar volumes must be specified in `PorousFlow` (in the [PorousFlowAqueousPreDisChemistry.md] Material).  This is so that the concentrations may be measured in $m^{3}/m^{3}$ rather than mol.m$^{-3}$.

2. Users of `PorousFlow` must specify the stoichiometric coefficients (in the [PorousFlowMassFractionAqueousEquilibriumChemistry.md] and [PorousFlowAqueousPreDisChemistry.md] Materials and the [PorousFlowPreDis.md] Kernel)

3. The density of the mineral species must be specified (in [PorousFlowPreDis.md]) since PorousFlow is based on mass fluxes, rather than volume fluxes.

4. Mass lumping for the precipitation-dissolution is used: porosity, saturation and mineral reaction rate are all retrieved from "nodal" Materials.

5. The advection of fluid species is typically performed by a [PorousFlowAdvectiveFlux.md] Kernel, which employs full upwinding.  It is assumed that mass fraction is equal to volume fraction, because only aqueous reactions are currently considered.

The method used in PorousFlow is:

- The concentrations of primary species are `Variables`.

- A [PorousFlowMassFractionAqueousEquilibriumChemistry.md] Material computes the `PorousFlow` mass-fractions.  Physically, these are total concentrations of the primary species of the aqueous equilibrium chemical reactions, along with a final mass fraction of the final component, which is assumed to be pure water.

- Advection of the `PorousFlow` mass-fractions (total concentrations) is performed by a [PorousFlowAdvectiveFlux.md] Kernel.  Dispeersion is performed by a [PorousFlowDispersiveFlux.md] Kernel.

- A [PorousFlowAqueousPreDisChemistry.md] Material computes the reaction rate for each precipitation/dissolution reaction.

- The precipitation/dissolution reaction rate feeds into a [PorousFlowPreDis.md] Kernel to ensure primary species are removed-from or added-to the aqueous solution appropriately.

- The precipitation/dissolution reaction rate also feeds to a [PorousFlowAqueousPreDisMineral.md] Material to compute the concentrations of mineral (so-called secondary concentrations).  This can then be fed to a [PorousFlowPorosity.md] Material to alter porosity.

All derivatives of total concentrations, precipitation/dissolution reaction rates, and mineral concentrations with respect to the `Variables` (primary species and temperature) are computed, and fed through the usual chain-rule procedure to compute the Jacobian entries.

## Current `ChemicalReactions`

The `ChemicalReactions` module is documented [here](chemical_reactions/index.md).  Most of the internal workings of `ChemicalReactions` is hidden from the user because it relies heavily on Actions that add the required Kernels, etc.  Here we describe the internal workings.

`ChemicalReactions` describes aqueous equilibrium and precipitation/dissolution chemistry.  Temperature is assumed to be prescribed (either by some dynamics, or being constant, or an AuxVariable, etc).  Parameters such as equilibrium constants may depend on temperature, since they can be `AuxVariables` (for instance, computed using [EquilibriumConstantAux.md]), but no derivatives of residuals with respect to temperature are computed.  The kinetic rate constant for the precipitation/dissolution is assumed to follow a temperature-dependent Arrhenius relation (and its derivatives are not used).

In terms of the concentration of primary species, $C_{j}$, and the concentration of the secondary species, $C_{i}$, the total concentration of the $j^{\mathrm{th}}$ primary species, $\Psi_{j}$, is
\begin{equation}
\Psi_{j} = C_{j} + \sum_{i}\nu_{ji}C_{i} \ .
\end{equation}
(Note here the conventional, but mathematically unusualy, notation where subscript $j$ indicates "primary" and subscript $i$ indicates "secondary".  That is, subscripts give special meaning to the quantity they are subscripting.)  The concentration of secondary species, $C_{i}$ is a complicated function of all the $C_{j}$.

In `ChemicalReactions` the concentrations of primary species are:

- measured in mol.litre$^{-1}$ (question: or is it mol.m$^{-3}$?);
- always `Variables`.

The continuity equation reads:
\begin{equation}
0 = \phi \frac{\partial}{\partial t}\Psi_{j} + \nabla\cdot ({\mathbf{v}}\Psi_{j}) - \nabla\cdot (\phi D\nabla \Psi_{j}) + \mathrm{P}_{j} - Q_{j} \ .
\end{equation}
The terms are as follows.

1. $\phi \frac{\partial}{\partial t}\Psi_{j}$ is computed by [PrimaryTimeDerivative.md] (the time-derivative of $C_{j}$) and a number of [CoupledBEEquilibriumSub.md] (corresponding to the time derivatives of all the $C_{i}$).  The porosity, $\phi$, is assumed to be constant and no lumping is performed.  $C_{i}$ are computed inside the Kernel.  Derivatives of the residual with respect to all $C_{j}$ are worked out.

2. $\nabla\cdot ({\mathbf{v}}\Psi_{j})$ is computed by [PrimaryConvection.md] (the $C_{j}$ part of $\Psi_{j}$) and [CoupledConvectionReactionSub.md] (the $C_{i}$ part of $\Psi_{j}$).  Here ${\mathbf{v}}$ is the Darcy velocity, which is defined in terms of a coupled pressure `Variable` or `AuxVariable`.  The hydraulic conductivity in the Darcy velocity is assumed to be constant, and no numerical stabilization is performed.  $C_{i}$ is computed inside the Kernel.  Derivatives of the residual with respect to all $C_{j}$ and the pressure are computed.

3. $\nabla\cdot (\phi D\nabla \Psi_{j})$ is computed by [PrimaryDiffusion.md] (the $C_{j}$ part of $\Psi_{j}$) and [CoupledDiffusionReactionSub.md] (the $C_{i}$ part of $\Psi_{j}$).  The product $\phi D$ is assumed to be contant and no numerical stabilization is performed.  $C_{i}$ is computed inside the Kernel.  Derivatives of the residual with respect to all $C_{j}$ are computed.  [!style color=red](*IS $\phi$ REALLY HERE?!*)

4. $\mathrm{P}_{j}$ is the rate of mineral precipitation.  It is included in the continuity equation by the [CoupledBEKinetic.md] Kernel.  This requires a coupled `AuxVariable` that is the mineral concentration, computed by [KineticDisPreConcAux.md] (and the time rate-of-change of this is used in the Kernel's residual).  The mineral concentration is performed at the nodes, and the [CoupledBEKinertic.md] uses quadpoint values.  No derivatives of the residual are computed.

5. $Q_{j}$ is a source.

If pressure is a `Variable` then it is solved in a quasi-static fashion (solving $\nabla\cdot(\nabla P - \rho {\mathbf{g}}) = 0$) using the [DarcyFluxPressure.md] Kernel.

`ChemicalReactions` contains a number of `AuxKernels`. to compute secondary concentrations, mineral concentrations and mineralization rates.  [EquilibriumConstantAux.md] can compute temperature-dependent equilibrium constants that can then be used in the Kernels.  [KineticDisPreConcAux.md] has already been mentioned: it computes the mineral concentrations.  The following are used for visualisation purposes only: [AqueousEquilibriumRxnAux.md] computes the concentration of secondary species; [KineticDisPreRateAux.md] computes the precipitation rate; [PHAux.md] computes the pH; [TotalConcentrationAux.md] computes the total concentration.

The [ChemicalOutFlowBC.md] sets the normal derivative of a variable, which can be used to remove primary species from a boundary to which they are flowing.




## Enhancements

- Effort should be put into simplifying the input file syntax for the end users. Suggestions from users would be most welcome to help guide this work.
- Ability to decouple the physics, via operator splitting and implemented via MultiApps, to solve reactions and then transport, or other approaches.  For example, TOUGHREACT solves fluid+heat equations to get temperature distributions, phase saturations and flow velocities.  It then uses these to transport the components.  Then it solves the chemistry.  The transport-reaction process can be iterated until convergence.  Finally, it updates porosity, permeability and capillary pressure due to changes in minerals; and water/gas liberated act as sources for the fluid+heat equations in the next timestep.  Need to look at PFLOTRAN in detail for another example.
- ChemicalReactions database reader.  One fundamental thing Andy doesn't understand is whether $K=K(P,T)$, or $K=K(T)$, where $K$ is the equilibrium constant.  Andy thinks the databases only have $K=K(T)$, but somewhere in the ToughReact manual it mentions $P$ dependence.
- Different possibilities for reaction rate of mineralisation that depends on a different functional form of $\Omega$ (see ToughReact Eqns B.9)
- Different possibilities for reaction rate-"constant" of mineralisation, which depends on pH or other species instead of just Arrhenius (see ToughReact Eqns B.10-B.13).
- Chemical reactions of GeoTES (not sure what this is)
- ChemicalReactions gas chemistry (see ToughReact section B.5 - i'm not sure whether this is actually quite simple or hard!)
- ChemicalReactions cation exchange and surface complexation (See ToughReact Appendix B).  Do we need these?  I think this just modifies the total concentration equations.
- Include numerical stabilization in ChemicalReactions.
- Include the PorousFlow checks on negative mineral concentration (etc) into ChemicalReactions.
- Presumably use PorousFlow's porosity and permeability classes.
- Reactive surface area a function of mineral volume fraction, water saturation and density (ToughReact Appendix G).
- Activity coefficients functions of molality and temperature (ToughReact Appendix H).
- PorousFlow hysteresis.
- Perhaps put all the transport of ChemicalReactions into PorousFlow, so ChemicalReactions only does ODEs, in which case PorousFlow needs to be able to easily solve the current type of transport that's in ChemicalReactions (fullly saturated, etc).

## Impact minimization

A refactor of `ChemicalReactions` will probably have limited impact:

- We are enormously increasing its capabilities anyway, which is likely to be welcomed by users;
- We believe it has few current users;
- Most input files contain `Actions`, so the user does not know what goes on "behind-the-scenes" anyway;

A refactor of `PorousFlow` will have to be handled carefully, with backwards compatibility.

## Proposed new design 1

1. Need to be able to send porosity, dispersion/diffusion and Darcy velocity from PorousFlow to ChemicalReactions in a MultiApp situation.  Similarly, need to be able to send total concentrations from ChemicalReactions to PorousFlow in order that things like viscosity can depend on mass fraction.

2. Use a `NodalUserObject` to compute the total concentrations (mol.litre$^{-1}$) at the nodes, given

- primary species concentrations
- stoichiometric coefficients (assumed constant)
- activity coefficients (assumed constant?)
- equilibrium constants (assumed to depend on $T$ only)

Derivatives of the total concentrations with respect to the primary concentrations and temperature are also computed.  [!style color=red](Will this ever be different for different blocks?)

3. Also need to compute $\nabla\Psi_{j}$, at the quadpoints, for use in diffusion/dispersion.  Presumably this can be a Material, or perhaps we can compute it in an `ElementUserObject`, based on the `NodalUserObject` information (assuming no block-dependence).  It is super annoying that we need both nodal and quadpoint total concentrations!

4. Write a `NodalUserObject` that computes the mineral reaction rates at the nodes, given

- primary species concentrations
- stoichiometric coefficients (assumed constant)
- activity coefficients (assumed constant?)
- equilibrium constants (assumed to depend on $T$ only)
- rate constants (assumed to depend on $T$, and maybe pH and primary species concentrations)
- reactive surface areas (assumed constant?)

Derivatives of the reaction rate with respect to the primary concentrations and temperature are also computed.  [!style color=red](Will this ever be different for different blocks?)

4. Write a `NodalUserObject` to take these mineral reaction rates, and their derivatives, and compute mineral concentrations (and their derivatives).

5. Also need the mineral concentrations (and derivatives) at the quadpoints, to compute porosity at the quadpoints, which is needed in PorousFlow (but possibly not in ChemicalReactions: see [!style color=red](red) text above).  An `ElementUserObject` or `Material` or something else?

6. The mineral concentrations need to feed into [PorousFlowPorosity.md] to compute porosity both at the nodes and quadpoints in a PorousFlow-only simulation, a MultiApp simulation, or a strongly coupled simulation.

7. The mineral concentration needs to feed into a new ChemicalReactions porosity object.  This needs to be a "nodal" Material, because it needs to be evaluated at the nodes (assuming we're doing mass lumping), but in-situ porosity is block dependent, so the nodal porosity depends on which element we're in.  If $\phi$ doesn't appear in the diffusion term (above - see [!style color=red](red) text) then this porosity only needs to be evaluated at the nodes (assuming we're doing lumping to the nodes).  In a standalone ChemicalReactions situation, $\phi$, and its derivatives with respect to mineral concentrations (and hence temperature and primary concentrations) can be computed directly using the initial porosity and the mineral concentrations.  In MultiApp simulations, $\phi$ will need to come from PorousFlow at the start of the timestep, and be modified by the mineralization.  In fully-coupled simulations, we just need to use [PorousFlowPorosity.md], and not this new ChemicalReactions porosity object.

8. Write new Kernels for ChemicalReactions that solve for the primary concentrations, by evolving the total concentrations (computed using the NodalUserObjects), using prescribed values of porosity, dispersion/diffusion and Darcy velocity, and also using the mineral reaction rates computed by the `NodalUserObject`.  Porosity comes from the new ChemicalReactions porosity class when using MultiApps or in a standalone ChemicalReactions simulation (in a strongly-coupled situation we use the PorousFlow Kernels instead).  Porosity can change during this step in the solution problem.  Almost definitely mass lumping and full upwinding should be used.  This transport can be used stand-alone (without PorousFlow) or as part of a MultiApp situation.  It's a bit of a bummer that we can't use the PorousFlow Kernels instead, since they essentially do the same thing, but that would preclude us from running ChemicalReaction-only simulations.

9. Use the NodalUserObject for total concentration in computations in [PorousFlowMassFractionAqueousEquilibriumChemistry.md] instead of the current code.

10. Use the NodalUserObject for reaction rate in the [PorousFlowPreDis.md] Kernel instead of the current code.

11. Use the NodalUserObject for total concentration in [TotalConcentrationAux.md]



## Proposed new design 2
