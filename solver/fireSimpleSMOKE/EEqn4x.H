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

if (energyEquation == true)
{
	T.max(minT);
    	T.min(maxT);

	if (combustionModel == "EDC")
	{
		TStar_[0].max(minT);
    		TStar_[0].min(maxT);
	}
	volScalarField& he = thermo.he();

	if (laminarDiffusion == true)
	{
		dimensionedScalar zero("zero", dimensionSet(1,-1,-3,0,0,0,0), scalar(0.)); 
		massDiffusionInEnergyEquation = zero;

		forAll(Y, i)
		{ 
			volScalarField& Yi = Y[i];
			volScalarField& Dmixi = Dmix_[i];
			volScalarField& HeSpeciesi = HeSpecies_[i];

			massDiffusionInEnergyEquation += HeSpeciesi* (fvc::laplacian(rho*Dmixi, Yi)); 
		}
	}

	// Energy equation
	{
		
		fvScalarMatrix EEqn
		(
			fvm::div(phi, he)
			+ 	(
					he.name() == "e"
					? fvc::div(phi, volScalarField("Ekp", 0.5*magSqr(U) + p/rho))
					: fvc::div(phi, volScalarField("K", 0.5*magSqr(U)))
			)
			- fvm::laplacian(turbulence->mut()/Prt, he)
			#if DEVVERSION == 1
			- radiation->Sh(thermo, he)
			#else
			- radiation->Sh(thermo)
			#endif
                        + massDiffusionInEnergyEquation
			==
			Qdot +
			fvOptions(rho, he)
		);

    		EEqn.relax();

    		fvOptions.constrain(EEqn);

   		EEqn.solve();

    		fvOptions.correct(he);

    		thermo.correct();
    		radiation->correct();
	}

	T.max(minT);
    	T.min(maxT);
	
	Info<< "min/max(T) = " << min(T).value() << ", " << max(T).value() << endl;

	#include "sparkPatch.H"
}
