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
#include "propertiesFineStructures.H"

	const int unsigned ns = Y.size();

	Eigen::VectorXd omega(ns);
	Eigen::VectorXd R(ns);
	Eigen::VectorXd r(kineticsMapXML->NumberOfReactions());
	OpenSMOKE::OpenSMOKEVectorDouble h_species(ns);

	// These are used only for the EDFR model
	OpenSMOKE::OpenSMOKEVectorDouble omega_(ns);
	OpenSMOKE::OpenSMOKEVectorDouble c_(ns);

	// OpenFOAM scalar fields
	const scalarField &pCells = thermo.p().internalField();
	const scalarField &rhoCells = rho.internalField();
	const scalarField &TCells = thermo.T().internalField();

	scalarField &QdotCells = Qdot.ref();
	scalarField &dQCells = dQ.ref();

	forAll(TCells, celli)
	{
		const scalar T_cell = TCells[celli];

		if (T_cell > minTemperature4Chemistry)
		{
			activation_map[celli] = 1.;
			// Calculating constant coefficient
			const double rho_epsilon_over_k = rhoCells[celli] * one_over_tauStar.internalField()[celli];

			// Vector containing mass fractions and mole formation rates
			for (unsigned int i = 0; i < ns; i++)
				omega(i) = Y[i].ref()[celli];
			R.setZero();

			// Set spark
			if ((spark == true) &&
				(runTime.value() >= ignitionTime) && (runTime.value() <= (ignitionTime + ignitionDuration)))
			{
				for (unsigned int i = 0; i < ns; i++)
					omega(i) = (omega(i) + ignitionAdditionalMassFraction) / (1. + ns * ignitionAdditionalMassFraction);
			}

			// Calculating the reaction rates
			if (combustionModel == "ED")
			{
				for (unsigned int k = 0; k < kineticsMapXML->NumberOfReactions(); k++)
				{
					double sumOmegaProducts = 0.;
					for (unsigned int j = 0; j < ed_product_index_[k].size(); j++)
						sumOmegaProducts += omega(ed_product_index_[k][j]);

					r[k] = constA * constB * sumOmegaProducts / ed_denominator_products_[k];
					for (unsigned int j = 0; j < ed_reactant_index_[k].size(); j++)
						r[k] = std::min(r[k], constA * omega(ed_reactant_index_[k][j]) / ed_denominator_reactants_[k][j]);
					r[k] *= rho_epsilon_over_k;
				}
			}
			else if (combustionModel == "EDFR")
			{
				double MW;
				for (unsigned int i = 0; i < ns; i++)
					omega_[i + 1] = omega(i);

				// Calculating the concentrations
				double cTot = pCells[celli] / (PhysicalConstants::R_J_kmol * TCells[celli]);
				thermodynamicsMapXML->MoleFractions_From_MassFractions(c_.GetHandle(), MW, omega_.GetHandle());
				Product(cTot, &c_);

				// Calculates kinetics
				thermodynamicsMapXML->SetTemperature(TCells[celli]);
				thermodynamicsMapXML->SetPressure(pCells[celli]);
				kineticsMapXML->SetTemperature(TCells[celli]);
				kineticsMapXML->SetPressure(pCells[celli]);
				kineticsMapXML->ReactionEnthalpiesAndEntropies();
				kineticsMapXML->KineticConstants();
				kineticsMapXML->ReactionRates(c_.GetHandle());

				for (unsigned int k = 0; k < kineticsMapXML->NumberOfReactions(); k++)
				{
					double sumOmegaProducts = 0.;
					for (unsigned int j = 0; j < ed_product_index_[k].size(); j++)
						sumOmegaProducts += omega(ed_product_index_[k][j]);

					r[k] = constA * constB * sumOmegaProducts / ed_denominator_products_[k];
					for (unsigned int j = 0; j < ed_reactant_index_[k].size(); j++)
						r[k] = std::min(r[k], constA * omega(ed_reactant_index_[k][j]) / ed_denominator_reactants_[k][j]);

					// Constant coefficient
					r[k] *= rho_epsilon_over_k;

					// Comparison with chemistry
					r[k] = std::min(r[k], kineticsMapXML->GiveMeReactionRates()[k]);
				}
			}

			// Calculating the formation rates [kmol/m3/s]
			for (unsigned int k = 0; k < kineticsMapXML->NumberOfReactions(); k++)
			{
				for (unsigned int j = 0; j < ed_reactant_index_[k].size(); j++)
					R[ed_reactant_index_[k][j]] -= ed_reactant_nu_[k][j] * r[k];
				for (unsigned int j = 0; j < ed_product_index_[k].size(); j++)
					R[ed_product_index_[k][j]] += ed_product_nu_[k][j] * r[k];
			}

			// Formation rates [kg/m3/s]
			for (unsigned int j = 0; j < ns; j++)
				R_[j].ref()[celli] = R[j] * thermodynamicsMapXML->MW(j);
		}
		else
		{
			activation_map[celli] = 0.;
			for (unsigned int i = 0; i < ns; i++)
				R_[i].ref()[celli] *= 0.; //[kg.m-3.s-1]
		}

		// Calculating the heat release
		{
			thermodynamicsMapXML->SetTemperature(TCells[celli]);
			thermodynamicsMapXML->SetPressure(pCells[celli]);
			thermodynamicsMapXML->hMolar_Species(h_species.GetHandle());

			double HeatRelease = 0.;
			for (unsigned int i = 0; i < ns; i++)
				HeatRelease += -h_species[i + 1] * R[i];
			QdotCells[celli] = HeatRelease;
			dQCells[celli] = HeatRelease * mesh.V()[celli];
		}
	}

	// Correct boundary conditions
	Qdot.correctBoundaryConditions();
	dQ.correctBoundaryConditions();
	Info << " * min/max(Q) = " << min(Qdot).value() << ", " << max(Qdot).value() << endl;
}
