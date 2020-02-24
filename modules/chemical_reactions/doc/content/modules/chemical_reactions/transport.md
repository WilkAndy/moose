# Transport

Author: Andy Wilkins

The `geochemistry` module is designed to solve reactive transport in geochemical systems.  The `geochemistry` module is currently in development.  The documentation describes the module's functionality, but does not currently include any description of the MOOSE implementation.

The functionality of the `geochemistry` module is currently focussed on simulation of GeoTES systems, which operate between around $25^{\circ}$C and $300^{\circ}$C.  This means, for instance, that there are no plans to implement the virial Pitzer/HMW models for activity.  The `geochemistry` module's functionality is a subset of that described in [!cite](bethke_2007).

I have tried to write thorough documentation that is understandable to someone with no background in (geo)chemistry.  The `geochemistry` module consists of two broad parts: the reactions and the transport.

!bibtex bibliography
