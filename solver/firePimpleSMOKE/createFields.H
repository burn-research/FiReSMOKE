/*——————————————————————————————————————————————————————————————————————-*\
|                              _____ __  __  ____  _  ________            |
|            ___              / ____|  \/  |/ __ \| |/ /  ____|           |
|           /  _| _  ___  ___| (___ | \  / | |  | | ' /| |__              |
|           | |_ | ||  _|/ _ \\___ \| |\/| | |  | |  < |  __|             |
|           |  _|| || | |  __/ ___) | |  | | |__| | . \| |____.           |
|           |_|  |_||_|  \___|_____/|_|  |_|\____/|_|\_\______|           |
|                                                                         |
|   Authors: A. Cuoci, R. Amaduzzi, A. Péquin, A. Parente                 |
|                                                                         |
|   Contacts: Alberto Cuoci                                               |
|   email: alberto.cuoci@polimi.it                                        |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano (Italy)                      |
|                                                                         |
|   Contacts: Ruggero Amaduzzi, Arthur Péquin, Alessandro Parente         |
|	email: alessandro.parente@ulb.be			                          |
|   Aero-Thermo-Mechanical Department                                     |
|   Université Libre de Bruxelles                                         |
|   Avenue F. D. Roosevelt 50, 1050 Bruxelles (Belgium)                   |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of fireSMOKE solver.                                |
|                                                                         |
|       License                                                           |
|                                                                         |
|   Copyright(C) 2017-2014 A. Cuoci, A. Parente                           |
|   fireSMOKE is free software: you can redistribute it and/or modify     |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   fireSMOKE is distributed in the hope that it will be useful,          |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with fireSMOKE. If not, see <http://www.gnu.org/licenses/>.     |
|                                                                         |
\*------------------------------------------------------------------------*/

// #include "createRDeltaT.H"

Info<< nl << "Reading thermophysicalProperties" << endl;

autoPtr<fluidMulticomponentThermo> thermoPtr //v12
(
    fluidMulticomponentThermo::New(mesh)
);
fluidMulticomponentThermo& thermo = thermoPtr(); //v12
thermo.validate(args.executable(), "h", "e");

PtrList<volScalarField>& Y = const_cast<PtrList<volScalarField>&>(thermo.Y()); //v12

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.name(), //v12
        mesh
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.name(), //v12
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField& p = thermo.p();
const volScalarField& psi = thermo.psi();
      volScalarField& T = thermo.T();

#include "compressibleCreatePhi.H"

pressureReference pressureReference(p, pimple.dict(), false);	// OFv10
mesh.schemes().setFluxRequired(p.name());

dimensionedScalar rhoMax
(
    dimensionedScalar::lookupOrDefault
    (
	"rhoMax",
	pimple.dict(),
	dimDensity,
	GREAT
    )
);

dimensionedScalar rhoMin
(
    dimensionedScalar::lookupOrDefault
    (
	"rhoMin",
	pimple.dict(),
	dimDensity,
	0
    )
);

Info << "Creating turbulence model.\n" << nl;

autoPtr<compressible::momentumTransportModel> momentumTransport //v12
(
    compressible::momentumTransportModel::New(rho, U, phi, thermo)
);

Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidMulticomponentThermophysicalTransportModel> thermophysicalTransportPtr
(
    fluidMulticomponentThermophysicalTransportModel::New
    (
        momentumTransport(),
        thermo
    )
);

fluidMulticomponentThermophysicalTransportModel& thermophysicalTransport = thermophysicalTransportPtr();

Info<< "Creating reaction model\n" << endl;	// OFv10
autoPtr<combustionModel> reaction //v12
(
    combustionModel::New(thermo, *momentumTransport)
);

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.name(), //v12
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

dimensionedScalar initialMass = fvc::domainIntegrate(rho);	// OFv10

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

#include "createMRF.H"
#include "createFvModels.H"	// OFv10
#include "createFvConstraints.H"	// OFv10

fileName path = args.rootPath()+"/"+args.globalCaseName();
OFstream Pav(path+"/p.Vessel");


Info<< "Creating field number of active species\n" << endl;

volScalarField drg_important_species
(
    IOobject
    (
        "drg::important_species",
        runTime.name(), //v12
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("drg::important_species", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

volScalarField drg_important_reactions
(
    IOobject
    (
        "drg::important_reactions",
        runTime.name(), //v12
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("drg::important_reactions", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

volScalarField cpu_reaction_times
(
    IOobject
    (
        "cpu::reaction_times",
        runTime.name(), //v12
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("cpu::reaction_times", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

#if SPARC==1
    #include "extensions/sparc/SPARC_create_fields.H"
#endif
