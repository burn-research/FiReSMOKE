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

	Info<< "Creating reaction model\n" << endl;

	#if OPENFOAM_VERSION >=60
	autoPtr<psiReactionThermo> pThermo
	(
		psiReactionThermo::New(mesh)
	);
	psiReactionThermo& thermo = pThermo();
	#else
	autoPtr<combustionModels::psiCombustionModel> reaction
	(
	    combustionModels::psiCombustionModel::New(mesh)
	);
	psiReactionThermo& thermo = reaction->thermo();
	#endif

	thermo.validate(args.executable(), "h", "e");

	basicMultiComponentMixture& composition = thermo.composition();
	PtrList<volScalarField>& Y = composition.Y();

	word inertSpecie(thermo.lookup("inertSpecie"));



    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();
    volScalarField& T = thermo.T();

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "compressibleCreatePhi.H"


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, simple.dict(), pRefCell, pRefValue);

    #if OPENFOAM_VERSION >= 40

	#if DEVVERSION == 1
	pressureControl pressureControl(p, rho, simple.dict());
	#endif

	mesh.setFluxRequired(p.name());

	Info<< "Creating turbulence model\n" << endl;
	autoPtr<compressible::turbulenceModel> turbulence
	(
		compressible::turbulenceModel::New
		(
			rho,
			U,
			phi,
			thermo
		)
	);

	#include "createMRF.H"

    #else

	dimensionedScalar rhoMax(simple.dict().lookup("rhoMax"));
    	dimensionedScalar rhoMin(simple.dict().lookup("rhoMin"));

	Info<< "Creating turbulence model\n" << endl;
	autoPtr<compressible::turbulenceModel> turbulence
	(
		compressible::turbulenceModel::New
		(
			rho,
			U,
			phi,
			thermo
		)
	);

    #endif

	#if OPENFOAM_VERSION >=60
	Info<< "Creating reaction model\n" << endl;
	autoPtr<CombustionModel<psiReactionThermo>> reaction
	(
	    CombustionModel<psiReactionThermo>::New(thermo, turbulence())
	);
	#else
	// Set the turbulence into the reaction model
	reaction->setTurbulence(turbulence());
	#endif

    dimensionedScalar initialMass = fvc::domainIntegrate(rho);

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());


Info<< "Creating field number of active species\n" << endl;

volScalarField drg_important_species
(
    IOobject
    (
        "drg::important_species",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("drg::important_species", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

volScalarField drg_important_reactions
(
    IOobject
    (
        "drg::important_reactions",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("drg::important_reactions", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

volScalarField cpu_reaction_times
(
    IOobject
    (
        "cpu::reaction_times",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("cpu::reaction_times", dimensionSet(0, 0, 0, 0, 0), 0.0)
);


#include "extensions/sparc/SPARC_create_fields.H"