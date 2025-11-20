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

{
	if (counterToUpdateChemistry == numberIterationsToUpdateChemistry)
	{
		counterToUpdateChemistry = 0;

		const double threshold_tau_c = 0.10;

		bool debug = false;
		Info << "Fine structures integration (SPARC)..." << endl;
		#include "propertiesFineStructures.H"

		const int unsigned nsTot = Y.size();
		const int unsigned neTot = ne_total;

		OpenSMOKE::OpenSMOKEVectorDouble h_species(nsTot);

		const scalarField& pCells = thermo.p().internalField();
		const scalarField& TCells = thermo.T().internalField();
		const scalarField& kappaCells = kappa.internalField();

		scalarField& QdotCells = Qdot.ref();
		scalarField& dQCells = dQ.ref();
		scalarField& cpuReactionsCells = cpu_reaction_times.ref();
		const scalarField& clusterIndexCells = clusterIndex.ref();
	
		// Statistics for PSR/NLS
		double nls_norm2_cumulative_ = 0.;
		double nls_norm2_max_ = 0.;
		int nls_failures_ = 0;
		int nls_number_iterations_cumulative_ = 0;

		// Statistics
		unsigned int counter = 0;
		double cpuTimeAllCells_ODE = 0.;

		const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
		forAll(TCells, celli)
		{	
			const scalar T_cell     = TCells[celli];

			if (T_cell > minTemperature4Chemistry)
			{
				activation_map[celli] = 1.;
				const unsigned int iCluster = clusterIndexCells[celli];
				const unsigned int NClocal  = vector_thermodynamicsMapXML[iCluster]->NumberOfSpecies();
				const unsigned int NRlocal  = vector_kineticsMapXML[iCluster]->NumberOfReactions(); 
			
				const int unsigned ns = NClocal;
				const int unsigned ne = ns+3;

				Eigen::VectorXd yStar0(ne);
				Eigen::VectorXd yStarf(ne);
				Eigen::VectorXd yMin(ne);
				Eigen::VectorXd yMax(ne);
	
				for(unsigned int i=0;i<ns;i++)	yMin(i) = 0.; yMin(ne-3) = minT;	yMin(ne-2) = 0.;	yMin(ne-1) = 0.;
				for(unsigned int i=0;i<ns;i++)	yMax(i) = 1.; yMax(ne-3) = maxT;	yMax(ne-2) = 1.e16;	yMax(ne-1) = 1.e16;

				// Integration time
				const double t0 = 0.;
				double tf = 0.;

				// Extract the fine strucutres mass fractions
				if (reactorType == REACTORTYPE_PSR)
				{
					Info << "PSR not available" << endl;
					abort();
				}
				else if (reactorType == REACTORTYPE_PFR)
				{
					if (iLES == true) // LES tauStar
					{
						if ( tauStarType == DT_LES)
						{
							tf = runTime.deltaTValue();
						}
						else if ( tauStarType == MIN_DTLES_TAUMIX)
						{
							tf = min(tauMix[celli], runTime.deltaTValue());
						}
					}
					else // RANS tauStar
					{
						tf = min(tauMix[celli], tauChem[celli]);
					}

					for(unsigned int i=0;i<ns;i++)
					{
						const unsigned int j = classifier_indices[iCluster](i);
						yStar0(i) = Y[j].internalField()[celli];
					}
					yStar0(ne-3) = TCells[celli];
					yStar0(ne-2) = tf;
					yStar0(ne-1) = 0.;

					vector_fine_structure_ode_pfr[iCluster].Set(pCells[celli]);
				}			

				// Solve the ODE system
				{
					const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

					if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
					{
						if (reactorType == REACTORTYPE_PSR)
						{
							Info << "PSR not available" << endl;
							abort();
						}
						else if (reactorType == REACTORTYPE_PFR)
						{
							#include "extensions/sparc/solvePFR_SPARC.H"
						}
					}
			
					const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

					cpuTimeAllCells_ODE += (tEnd-tStart);
				}
	
				{		
					for(unsigned int i=0;i<ns;i++)
					{
						const unsigned int j = classifier_indices[iCluster](i);
						YStar_[j].ref()[celli] = yStarf(i);
					}
					TStar_[0].ref()[celli] = yStarf(ne-3);
				
					for(unsigned int i=0;i<ns;i++)
					{
						const unsigned int j = classifier_indices[iCluster](i);
						R_[j].ref()[celli] = rho[celli]*kappa[celli]*(YStar_[j].internalField()[celli]-Y[j].internalField()[celli])/tf;
					}
				}

				thermodynamicsMapXML->SetTemperature(TCells[celli]);
				thermodynamicsMapXML->SetPressure(pCells[celli]);
				thermodynamicsMapXML->hMolar_Species(h_species.GetHandle());

				double HeatRelease = 0.;
				for(unsigned int i=0;i<ns;i++)
				{
					const unsigned int j = classifier_indices[iCluster](i);
					HeatRelease += -h_species[j+1]/thermodynamicsMapXML->MW(j) * R_[j].internalField()[celli];
				}
				QdotCells[celli] = HeatRelease;
				dQCells[celli] = HeatRelease * mesh.V()[celli];

				if (laminarDiffusion == true)
				{
					for(unsigned int i=0;i<ns;i++)
					{
						const unsigned int j = classifier_indices[iCluster](i);
						HeSpecies_[j].ref()[celli] = -h_species[j+1]/thermodynamicsMapXML->MW(j);
					}
				}

				if (debug == true)
				{
					Info << celli << " " << TCells[celli] << " " << gammaStarNum[celli] << " " << mDotStar[celli] << " " << HeatRelease << endl;
					for(unsigned int i=0;i<ns;i++)
					{
						const unsigned int j = classifier_indices[iCluster](i);
						Info << Y[j].name() << " " << yStar0[j] << " " << yStarf[j] << " " << R_[j].internalField()[celli] << endl;
					}
					getchar();
				}
			}
			else
			{
				activation_map[celli] = 0.;
				for(unsigned int i=0;i<nsTot;i++)
					R_[i].ref()[celli] *= 0.;//[kg.m-3.s-1]
			}


			if (counter%(int(0.20*mesh.nCells())+1) == 0)
			{
				if (solverType == SOLVERTYPE_ODE)
					Info << "   Accomplished: " << counter << "/" << mesh.nCells() << endl;
				else
				{
					Info 	<< "  Accomplished: " 		<< counter << "/" << mesh.nCells() << " ";
					Info 	<< "  Norm2(avg): " 		<< nls_norm2_cumulative_/static_cast<double>(counter+1) 	
						<< "  Norm2(max): " 		<< nls_norm2_max_ 
						<< "  Iterations(avg): " 	<< nls_number_iterations_cumulative_/static_cast<double>(counter+1) 
						<< "  Ode: " 			<< nls_failures_ << endl;
				}
			}

			counter++;

		}

		const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			
		Info << "CPU time for chemistry (PaSR): " << tEnd - tStart << " s " << "(" << (tEnd - tStart)/double(TCells.size())*1e3 << " ms per cell)" << endl;
		Info << "CPU time for ODE solution:    " << cpuTimeAllCells_ODE << " s " << "(" << (cpuTimeAllCells_ODE)/double(TCells.size())*1e3 << " ms per cell)" << endl;

	
		// Correct boundary conditions
		Qdot.correctBoundaryConditions();
		dQ.correctBoundaryConditions();
	}

	counterToUpdateChemistry++;
}

