# Transport

Author: Andy Wilkins

Only rudimentary transport is available as part of the `geochemistry` module.  Coupling with the `PorousFlow` module allows users to access more advanced features such as:

- pressure and temperature that tightly couple with fluid flows, instead of being specified by uncoupled dynamics
- densities, viscosities, etc, that depend on solute concentrations, temperature and pressure, instead of being constant
- porosity and permeability that change with precipitation and dissolution
- multiphase flow
- coupling with geomechanics
- sophisticated numerical stabilization

Notation and definitions are described in [nomenclature.md].

## Molalities, volumes and concentrations

When just considering geochemistry without considering spatial transport, here is no reason to think of finite elements at all.  When considering transport, the volume of the finite element, or volume associated with a node in the mesh, is important.  All the formulae presented in the chemical reactions (eg [equilibrium.md]) make no reference to this volume, but obviously they could be divided by a constant volume without any change.

## Transporting total concentrations

Because density is assumed to be constant, transport acts on the [basis](basis.md) species' total concentrations.  Recall that the basis is
\begin{equation}
\mathrm{basis} = (A_{w}, A_{i}, A_{k}, A_{m}, A_{p}) \ ,
\end{equation}
Only the water, the primary aqueous species, the minerals that are mobile and the dissolved gases, are transported.  Minerals and sorbed species are generally assumed to be immobile.  Equilibrium  provides [equations](equilibrium.md) for the bulk composition of water, primary aqueous species, etc.  It is simple to remove the contributions from the immobile mineral and sorbed species from those formulae to provide the concentrations involved in transport:
\begin{equation}
\begin{aligned}
C_{w} & = \frac{1}{V}\left( M_{w} - n_{w}\sum_{q}\nu_{wq}m_{q}\right) \\
C_{i} & = \frac{1}{V}\left( M_{i} - n_{w}\sum_{q}\nu_{iq}m_{q}\right) \\
C_{k} & = \frac{1}{V}\left( M_{k} - n_{k} - n_{w}\sum_{q}\nu_{kq}m_{q}\right) \\
C_{m} & = \frac{1}{V}\left( M_{m} - n_{w}\sum_{q}\nu_{mq}m_{q}\right) \ .
\end{aligned}
\end{equation}
In these equations, $V$ is the volume occupied by the components.  In the third equation, the immobile $n_{k}$ has been removed from the bulk composition $M_{k}$.

## Transport equations

Denote mobile concentrations collectively $(C_{w}, C_{i}, C_{k}, C_{m})$ by $C_{r}$.  Mass conservation (the continuity equation) reads
\begin{equation}
0 = \frac{\partial}{\partial t}(\phi C_{r}) + \nabla\cdot (\mathbf{q} - \phi D\nabla) C_{r} + \frac{1}{V}\sum_{\bar{k}}\tilde{\nu}_{r\bar{k}}r_{\bar{k}} \ .
\end{equation}
Here:

- $t$ \[s\] is time
- $\nabla$ is the vector spatial derivatives
- $\phi$ \[dimensionless\] is the rock porosity, which is assumed to be given to the `geomechanics` module by an external agent
- $\mathbf{q}$ \[m.s$^{-1}$\] is the Darcy flux vector.  usually $q_{i} = \frac{k_{ij}}{\mu}(\nabla_{j}P - \rho g_{j})$, where $k_{ij}$ \[m$^{2}$\] is the permeability, $\mu$ is the fluid viscosity, $P$ is its porepressure, $\rho$ is its density, and $g_{j}$ is the acceleration due to gravity.  (In these equations, $i$ and $j$ indicate spatial directions, not primary and secondary species!)  In `geomechanics` $\mathbf{q}$ is assumed to be given by an external agent, which could be a user-defined set of functions, or a field that varies spatially and temporally by solving Darcy's equation in conjunction with the `geomechanics` equations.
- $D$ \[m$^{2}$.s$^{-1}$\] is the hydrodynamic dispersion tensor, which is assumed to be given to the `geomechanics` module.
- $\tilde{\nu}$ are stoichiometric coefficients
- $r_{\bar{k}}$ \[mol.s$^{-1}$\] are kinetic reaction rates

The stoichiometric coefficients, $\tilde{\nu}$ could be complicated.  The [kinetic](kinetics.md) equations provide the rates of change for $(M_{w}, M_{i}, M_{p}, n_{\bar{k}})$ which then must be used to derive the rates of change for the $C_{r}$.  Similar computations to the "Redox kinetics" case [here](kinetics.md) must be performed.  In cases where there is no sorption and the kinetic reactions do not involve $A_{k}$, this is actually straightforward, since $\tilde{\nu}_{\ast\bar{k}} = \nu_{\ast\bar{k}}$.  In practice, the evaluation of $\tilde{\nu}$ is is not necessary, due to the operator-splitting method employed in the `geomechanics` module.

## Operator splitting

The equations are solved by splitting the spatial transport from the chemistry.  In particular:

- the transport equations $0 = \frac{\partial}{\partial t}\phi C_{r} + \nabla\cdot (\mathbf{q} - \phi D\nabla) C_{r}$ are solved to provide updates for $(M_{w}, M_{i}, M_{k}, M_{m})$.
- These updates are used as source terms in the equilibrium or kinetic equations.




!bibtex bibliography
