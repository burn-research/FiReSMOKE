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

//#if SPARC==1
//#include "extensions/sparc/SPARC_setup.H"

Eigen::MatrixXd CharacteristicChemicalTimes::SPARC_computeNumericalJacobian(const double T, const double P, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermoMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap)
{
	// ------------------------------------------------------------------------------------------ //
	// mPaSR function to return the Numerical Jacobian to mPaSR model
	// @A. Péquin, Université Libre de Bruxelles

	// T : temperature
	// P : pressure
	// y : vector of mass fractions (vec_y)

	// Jac_i_j 		= d y_dot_i / d y_j
	// d y_dot_i 	= ( y_dot_i(vec_y with perturbation on y_j) - y_dot_i(reference vec_y) )
	// d y_j 		= perturbation applied on y_j
	// ------------------------------------------------------------------------------------------ //
	// std::vector<double> testVector(100, 1.0);  // Allocate 100 doubles
	// Info << "testVector allocated successfully." << endl;

	// Definitions
	int N = y.Size();

	//Info << "Size of y computeJc: " << y.Size() << endl;  // Check the size of the input y
	//Info << "N before y_plus init computeJc: " << N << endl; // Crucial check!
	//Info << "Size of nsCluster computeJc: " << thermoMap.NumberOfSpecies() << endl;

	OpenSMOKE::OpenSMOKEVectorDouble y_plus(N);
	OpenSMOKE::OpenSMOKEVectorDouble y_dot_reference(N);		// rates of reference
	OpenSMOKE::OpenSMOKEVectorDouble y_dot_plus(N);				// rates from y_plus
	OpenSMOKE::OpenSMOKEVectorDouble omega(N);					// mass fractions modified indices
	OpenSMOKE::OpenSMOKEVectorDouble omega_plus(N);				// mass fractions modified indices holding a perturbation
	OpenSMOKE::OpenSMOKEVectorDouble x(N); 						// mole fractions
	OpenSMOKE::OpenSMOKEVectorDouble x_plus(N); 				// mole fractions holding a perturbation
	OpenSMOKE::OpenSMOKEVectorDouble c(N); 						// molar concentrations
	OpenSMOKE::OpenSMOKEVectorDouble c_plus(N); 				// molar concentrations holding a perturbation
	OpenSMOKE::OpenSMOKEVectorDouble dc_over_dt(N);				// molar rates [kmol.m-3.s-1]
	OpenSMOKE::OpenSMOKEVectorDouble dc_over_dt_plus(N);		// molar rates [kmol.m-3.s-1] holding a perturbation

	Eigen::MatrixXd J(N, N);									// Jacobian matrix

	// Perturbation
	const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
	const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);			
	const double TOLR = 100. * OpenSMOKE::OPENSMOKE_MACH_EPS_FLOAT;
	const double TOLA = 1.e-8;

	// Setting up: temperature and pressure
	T_ = T;
	P_ = P;

	//Info << "check kineticsMap " << kineticsMap.NumberOfReactions() << endl;
	thermoMap.SetTemperature(T);
	thermoMap.SetPressure(P);
	kineticsMap.SetTemperature(T);
	kineticsMap.SetPressure(P);
	kineticsMap.KineticConstants();

	// Calculations of the reference rates: y_dot_reference
	for(unsigned int i=0;i<N;i++)
		omega[i+1] = y[i];												  // filling up mass fraction vector with modified indices

	double MW;
	thermoMap.MoleFractions_From_MassFractions(x.GetHandle(), MW, omega.GetHandle()); 
	const double cTot = P/(PhysicalConstants::R_J_kmol * T);
	const double rho = cTot*MW;

	//Info << "Check compute cTot "  << cTot << endl;
	//Info << "Check compute rho "  << rho << endl;

	Product(cTot, x, &c);																								// getting molar concentrations c

	//Equations(c, dc_over_dt);																							// getting molar rates
	kineticsMap.ReactionRates(c.GetHandle());
	kineticsMap.FormationRates(dc_over_dt.GetHandle());


	// Info << "Check compute Jc before perturbation"  << endl;
	// for (unsigned int i = 0; i < thermoMap.NumberOfSpecies();  i++)
	// 		{
	// 			Info << thermoMap.MW(i) << " ";
	// 		}
	// Info << endl;

	// for (unsigned int i = 0; i < dc_over_dt.Size();  i++)
	// 		{
	// 			Info << dc_over_dt[i] << " ";
	// 		}
	// Info << endl;

	for(unsigned int i=0;i<N;i++)
		y_dot_reference[i] = dc_over_dt[i+1]*thermoMap.MW(i)/rho;
	
	// for (unsigned int i = 0; i < y_dot_reference.Size();  i++)
	// {
	// 	Info << y_dot_reference[i] << " ";
	// }
	// Info << endl;
	// Compute perturbed y_dot and fill the Jacobian matrix
	for(unsigned int j=0;j<N;j++)
	{
		for(unsigned int i=0;i<N;i++)
			y_plus[i] = y[i];

		double hf = 1.e0;
		double error_weight = 1./(TOLA+TOLR*std::fabs(y[j]));
		double hJ = ETA2 * std::fabs(std::max(y[j], 1./error_weight));
		double hJf = hf/error_weight;
		hJ = std::max(hJ, hJf);
		hJ = std::max(hJ, ZERO_DER);
		double dy = std::min(hJ, 1.e-3 + 1e-3*std::fabs(y[j]));

		y_plus[j] = y[j] + dy;																							// perturbed mass fractions

		for(unsigned int i=0;i<N;i++)
			omega_plus[i+1] = y_plus[i];																				// modified indices

		double MW_plus;
		thermoMap.MoleFractions_From_MassFractions(x_plus.GetHandle(), MW_plus, omega_plus.GetHandle());
		Product(cTot, x_plus, &c_plus);																					// perturbed molar concentrations

		//Equations(c_plus, dc_over_dt_plus);																				// perturbed molar rates
		kineticsMap.ReactionRates(c_plus.GetHandle());
		kineticsMap.FormationRates(dc_over_dt_plus.GetHandle());

		for(unsigned int i=0;i<N;i++)
			y_dot_plus[i] = dc_over_dt_plus[i+1]*thermoMap.MW(i)/rho;										// perturbed rates [s-1]

		// for (unsigned int i = 0; i < y_dot_plus.Size();  i++)
		// {
		// 	Info << y_dot_plus[i] << " ";
		// }
		// Info << endl;
		// Populate Jacobian matrix
		for(unsigned int i=0;i<N;i++)
			J(i,j) = (y_dot_plus[i]-y_dot_reference[i])/dy;

		// Info << "Check compute Jc "  << J.size() << endl;
		// for (unsigned int i = 0; i < y.Size();  i++)
		// 	for (unsigned int j = 0; i < y.Size();  i++)
        //             {
        //                 Info << J(i,j) << " ";
        //             }
        // Info << endl;
	}

	return J;
}

//#endif

