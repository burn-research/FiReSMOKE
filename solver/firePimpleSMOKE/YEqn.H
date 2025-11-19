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

tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.schemes().div("div(phi,Yi)")
    )
);


{
	reaction->correct();
	#include "chemistry.H"

	label inertIndex = -1;
	volScalarField Yt(0.0*Y[0]);

	const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

	if (laminarDiffusion == true)
	{
		forAll(Y, i)
		{
			if (Y[i].name() != inertSpecie)
			{
				volScalarField& Yi = Y[i];
				volScalarField& Dmixi = Dmix_[i];

				fvScalarMatrix YiEqn
				(
					fvm::ddt(rho, Yi)
					+ mvConvection->fvmDiv(phi, Yi)
					- fvm::laplacian(rho*Dmixi + rho*momentumTransport->nut()/Sct, Yi) //v12
					==
					R_[i]
					+ fvModels.source(rho, Yi)
				);

				YiEqn.relax();

				fvConstraints.constrain(YiEqn);
                                
				YiEqn.solve("Yi");
                                
				fvConstraints.constrain(Yi);

				Yi.max(0.0);
				Yt += Yi;
			}
			else
			{
				inertIndex = i;
			}
		}
	}
	else
	{

		forAll(Y, i)
		{
			if (Y[i].name() != inertSpecie)
			{
				volScalarField& Yi = Y[i];

				fvScalarMatrix YiEqn
				(
					fvm::ddt(rho, Yi)
					+ mvConvection->fvmDiv(phi, Yi)
					- fvm::laplacian(rho*momentumTransport->nuEff()/Sct, Yi) //v12
					==
					R_[i]
					+ fvModels.source(rho, Yi)
				);

				YiEqn.relax();

				fvConstraints.constrain(YiEqn);
                                
				YiEqn.solve("Yi");

				fvConstraints.constrain(Yi);

				Yi.max(0.0);
				Yt += Yi;
			}
			else
			{
				inertIndex = i;
			}
		}
    }

    Y[inertIndex] = scalar(1.0) - Yt;
    Y[inertIndex].max(0.0);

    const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

    Info << "CPU time for species transport (overall): " << tEnd - tStart << " s " << "(" << (tEnd - tStart)/double(mesh.C().size())*1e3 << " ms per cell)" << endl;
}



