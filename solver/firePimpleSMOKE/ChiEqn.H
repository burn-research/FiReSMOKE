/*-----------------------------------------------------------------------*\
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
|   email: alessandro.parente@ulb.be                                      |
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
\*-----------------------------------------------------------------------*/

{
	if (ChiEquation == true)
	{
	// Solving equation of Scalar Dissipation Rate 
	volScalarField D1 = Cd1*rho*sqr(Chi)/(varf+SMALL);
	volScalarField D2 = Cd2*rho*turbulence->epsilon()/(turbulence->k()+ksmall_)*Chi;
	volScalarField P1 = 2.00*Cp1*turbulence->epsilon()/(turbulence->k()+ksmall_)*turbulence->mut()/Sct*magSqr(fvc::grad(f));
	volScalarField P2 = Cp2*turbulence->mut()*Chi/(turbulence->k()+ksmall_)*(fvc::grad(U) && dev(twoSymm(fvc::grad(U))));

	volScalarField S_chi = P1 + P2 - D1 - D2;	

	fvScalarMatrix ChiEqn
	(
		fvm::ddt(rho, Chi)
		+ fvm::div(phi, Chi)
		- fvm::laplacian(thermo.alpha()+turbulence->mut()/Sct, Chi)
		==
		S_chi
	);

	ChiEqn.relax();
	ChiEqn.solve();
	Chi.max(0.00000001);
	Chi.min(maxChi);

    	Info << "Scalar dissipation rate transport equation solved" << endl;
    	Info<< "min/max(Chi) = " << min(Chi).value() << ", " << max(Chi).value() << endl;
	}
	else
	{
		Info << "Algebraic equation for Chi" << endl;
		scalarField& ChiCells = Chi.ref();
		//ChiCells = Cx/rho*(thermo.alpha()+turbulence->mut()/Sct)*magSqr(fvc::grad(f));
		ChiCells = Cx*turbulence->epsilon()/(turbulence->k()+ksmall_)*varf;
		Info<< "min/max(Chi) = " << min(Chi).value() << ", " << max(Chi).value() << endl;
	}
}

