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

#ifndef CharacteristicChemicalTimes_H
#define CharacteristicChemicalTimes_H

class CharacteristicChemicalTimes
{
	public:

	enum CharacteristicChemicalTimesType { CHEMICAL_TIMES_FORMATION_RATES, CHEMICAL_TIMES_REACTION_RATES, CHEMICAL_TIMES_EIGENVALUES };

	public:

	CharacteristicChemicalTimes(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& 		thermodynamicsMapXML, 
					OpenSMOKE::KineticsMap_CHEMKIN& 		kineticsMapXML);

	void SetThresholdChemicalTime(const double threshold);
	void SetThresholdTemperatureChemicalTime(const double threshold);
	void SetType(const CharacteristicChemicalTimesType type);
	void SetLapackMode(const bool flag);
	double CalculateCharacteristicChemicalTime(const double T, const double P, const OpenSMOKE::OpenSMOKEVectorDouble& omega);

	private:

	double FromFormationRates(const double T, const double P, const OpenSMOKE::OpenSMOKEVectorDouble& omega);
	double FromReactionRates(const double T, const double P, const OpenSMOKE::OpenSMOKEVectorDouble& omega);
	void Equations(const OpenSMOKE::OpenSMOKEVectorDouble& c, OpenSMOKE::OpenSMOKEVectorDouble& dc_over_dt);

	public:

	double FromEigenValueAnalysis(const double T, const double P, const OpenSMOKE::OpenSMOKEVectorDouble& omega);

	void NumericalJacobian(const OpenSMOKE::OpenSMOKEVectorDouble& y, Eigen::MatrixXd& J);
	void SearchConservativeModes(const std::vector<double>& lambda_real, const std::vector<double>& lambda_mod, std::vector<double>& lambda_real_cleaned);

	Eigen::MatrixXd mPaSR_computeNumericalJacobian(const double T, const double P, const OpenSMOKE::OpenSMOKEVectorDouble& y);
	Eigen::MatrixXd mPaSR_computeNumericalJacobian_const(const double T, const double P, const OpenSMOKE::OpenSMOKEVectorDouble& y);

	#if SPARC==1
		Eigen::MatrixXd SPARC_computeNumericalJacobian(const double T, const double P, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermoMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap);
	#endif

	private:

	OpenSMOKE::ThermodynamicsMap_CHEMKIN& 	thermodynamicsMapXML_;
	OpenSMOKE::KineticsMap_CHEMKIN& 	kineticsMapXML_;

	unsigned int ne_;
	unsigned int ns_;
	unsigned int nr_;
	unsigned int n_;

	bool lapack_mode_;

	double T_;
	double P_;

	Eigen::MatrixXd Jc_;
	Eigen::VectorXd sum_nu_;

	Eigen::VectorXd vr_;
	Eigen::VectorXd vl_;

	OpenSMOKE::OpenSMOKEVectorDouble x_;
	OpenSMOKE::OpenSMOKEVectorDouble c_;
	OpenSMOKE::OpenSMOKEVectorDouble R_;
	OpenSMOKE::OpenSMOKEVectorDouble rf_;

	// Local variables
	std::vector<double> lambda_real_;
	std::vector<double> lambda_imag_;
	std::vector<double> lambda_mod_;
	std::vector<double> lambda_real_cleaned_;

	double tauChem_threshold_;
	double tauChem_temperature_threshold_;

	CharacteristicChemicalTimesType type_;

	// OpenSMOKE::OpenSMOKEVectorDouble y_plus_;					// vec_y holding a perturbation
	// OpenSMOKE::OpenSMOKEVectorDouble y_dot_reference_;		// rates of reference
	// OpenSMOKE::OpenSMOKEVectorDouble y_dot_plus_;				// rates from y_plus

	// OpenSMOKE::OpenSMOKEVectorDouble omega_;					// mass fractions modified indices
	// OpenSMOKE::OpenSMOKEVectorDouble omega_plus_;				// mass fractions modified indices holding a perturbation
	// //OpenSMOKE::OpenSMOKEVectorDouble x_; 						// mole fractions
	// OpenSMOKE::OpenSMOKEVectorDouble x_plus_; 				// mole fractions holding a perturbation
	// //OpenSMOKE::OpenSMOKEVectorDouble c_; 						// molar concentrations
	// OpenSMOKE::OpenSMOKEVectorDouble c_plus_; 				// molar concentrations holding a perturbation
	// OpenSMOKE::OpenSMOKEVectorDouble dc_over_dt_;				// molar rates [kmol.m-3.s-1]
	// OpenSMOKE::OpenSMOKEVectorDouble dc_over_dt_plus_;	

};

#include "CharacteristicChemicalTimes.C"

#endif // CharacteristicChemicalTimes_H
