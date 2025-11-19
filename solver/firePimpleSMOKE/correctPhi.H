#include "fvc.H"

surfaceScalarField phi
(
    IOobject("phi", runTime.name(), mesh),
    fvc::interpolate(rho * U) & mesh.Sf()
);

U.correctBoundaryConditions();

// phi.correctBoundaryConditions();

adjustPhi(phi, U, p);

// Make the fluxes relative to the mesh-motion
fvc::makeRelative(phi, rho, U);