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

if (homogeneousReactions == true)
{
    Info << " ************************************************************** " << endl;
    Info << " * Using " << combustionModel << " model (if T > " << minTemperature4Chemistry << " [K])" << endl;
	const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

	if (combustionModel == "ED" || combustionModel == "EDFR")
	{
		#include "EDFR.H"
	}
	else if (combustionModel == "EDC")
	{
		#if SPARC==1
			if (SPARCswitch == true)
			{
				#include "extensions/sparc/EDC_SPARC.H"
			}
			else if (SPARCswitch == false)
			{
				#include "EDC.H"
			}
		#else
				#include "EDC.H"
		#endif
	}
	else if (combustionModel == "PaSR")
	{
		#if SPARC==1
			if (SPARCswitch == true)
			{
				#include "extensions/sparc/PaSR_SPARC.H"
			}
			else if (SPARCswitch == false)
			{
				#include "PaSR.H"
			}
		#else
				#include "PaSR.H"
		#endif
	}	
	else if (combustionModel == "QLFR")
	{
		#if SPARC==1
			if (SPARCswitch == true)
			{
				#include "extensions/sparc/QLFR_SPARC.H"
			}
			else if (SPARCswitch == false)
			{
				#include "QLFR.H"
			}
		#else
				#include "QLFR.H"
		#endif
	}
	else if (combustionModel == "PaSR_LFR") // SPARC not currently available for PaSR_LFR model
	{
		#include "PaSR_LFR.H"
	}
	else if (combustionModel == "LFR") // SPARC not currently available for LFR model
	{
		#include "LFR.H"
	}
	else if (combustionModel == "mPaSR")
	{
		#if SPARC==1
			if (SPARCswitch == true)
			{
				#include "extensions/sparc/mPaSR_SPARC.H"
			}
			else if (SPARCswitch == false)
			{
				#include "mPaSR.H"
			}
		#else
				#include "mPaSR.H"
		#endif
	}

    const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
    CPUtime_comb_model_part = CPUtime_comb_model_part + (tEnd - tStart);
    Info << " * Time for combustion model [s] = " << tEnd - tStart << endl;
    Info << " ************************************************************** " << endl;
}

