/*----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                               |
|                                                                         |
|   Copyright(C) 2014, 2013, 2012  Alberto Cuoci                          |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef OpenSMOKE_ThermodynamicsMap_H
#define OpenSMOKE_ThermodynamicsMap_H

namespace OpenSMOKE
{
	//!  A virtual class to provide a common interface to thermodynamic maps
	/*!
			This virtual class provides a common interface to thermodynamic maps
	*/

	class ThermodynamicsMap
	{
	
	public:

		/**
		* @brief returns the number of species
		*/
		unsigned int NumberOfSpecies() const { return nspecies_; }

		/**
		* @brief returns the names of the species
		*/
		const std::vector<std::string>& NamesOfSpecies() const { return names_; }

		/**
		* @brief returns the molecular weights of a species [kg/kmol]
		*/
		const std::vector<double>& MWs() const { return MW__; }

		/**
		* @brief returns the molecular weights of a species [kg/kmol]
		*/
		const double MW(const unsigned int i) const { return MW__[i];}

		/**
		* @brief Sets the temperature (in K)
		*/
		virtual void SetTemperature(const double& T) = 0;
	
		/**
		* @brief Sets the pressure (in Pa)
		*/
		virtual void SetPressure(const double& P) = 0;

		/**
		* @brief Calculates the molecular weight of a mixture from the mole fractions
		*/
		virtual double MolecularWeight_From_MoleFractions(const double* x) = 0;

		/**
		* @brief Calculates the molecular weight of a mixture from the mass fractions
		*/
		virtual double MolecularWeight_From_MassFractions(const double* y) = 0;

		/**
		* @brief Calculates the mass fractions from the mole fractions
		*/
		virtual void MassFractions_From_MoleFractions(double* y, double& MW, const double* x) = 0;

		/**
		* @brief Calculates the mole fractions from the mass fractions
		*/
		virtual void MoleFractions_From_MassFractions(double* x, double& MW, const double* y) = 0;

		/**
		* @brief Calculates the molar specific heat of the mixture from the mole fractions
		*/
		virtual double cpMolar_Mixture_From_MoleFractions(const double* x) = 0;

		/**
		* @brief Calculates the molar enthalpy of the mixture from the mole fractions
		*/
		virtual double hMolar_Mixture_From_MoleFractions(const double* x) = 0;

		/**
		* @brief Calculates the molar entropy of the mixture from the mole fractions
		*/
		virtual double sMolar_Mixture_From_MoleFractions(const double* x) = 0;

		/**
		* @brief Calculates the molar internal energy of the mixture from the mole fractions
		*/
		virtual double uMolar_Mixture_From_MoleFractions(const double* x) = 0;

		/**
		* @brief Calculates the molar free Gibbs' energy of the mixture from the mole fractions
		*/
		virtual double gMolar_Mixture_From_MoleFractions(const double* x) = 0;

		/**
		* @brief Calculates the molar free Helmoltz' energy of the mixture from the mole fractions
		*/
		virtual double aMolar_Mixture_From_MoleFractions(const double* x) = 0;

		/**
		* @brief Calculates the molar specific heats of the species from the mole fractions
		*/
		virtual void cpMolar_Species(double* cp_species) = 0;

		/**
		* @brief Calculates the molar enthalpies of the species from the mole fractions
		*/
		virtual void hMolar_Species(double* h_species) = 0;

		/**
		* @brief Calculates the molar entropies of the species from the mole fractions
		*/
		virtual void sMolar_Species(double* s_species) = 0;

		/**
		* @brief Calculates the molar internal energies of the species from the mole fractions
		*/
		virtual void uMolar_Species(double* u_species) = 0;

		/**
		* @brief Calculates the molar Gibbs' free energies of the species from the mole fractions
		*/
		virtual void gMolar_Species(double* g_species) = 0;

		/**
		* @brief Calculates the molar Helmoltz' free energies of the species from the mole fractions
		*/
		virtual void aMolar_Species(double* a_species) = 0;

		/**
		* @brief Calculates the mixture averaged molar entropies of the species from the mole fractions
		*/
		virtual void sMolar_Species_MixtureAveraged_From_MoleFractions(double* s_species, const double* x) = 0;

		/**
		* @brief Calculates the mixture averaged free Gibbs' energies of the species from the mole fractions
		*/
		virtual void gMolar_Species_MixtureAveraged_From_MoleFractions(double* g_species, const double* x) = 0;
		
		/**
		* @brief Calculates the mixture averaged molar Helmoltz's free energies of the species from the mole fractions
		*/		
		virtual void aMolar_Species_MixtureAveraged_From_MoleFractions(double* a_species, const double* x) = 0;

		/**
		* @brief returns the normalized reaction enthalpies for the species
		*/	
		virtual const std::vector<double>& Species_H_over_RT()  = 0;

		/**
		* @brief returns the normalized reaction entropies for the species
		*/
		virtual const std::vector<double>& Species_S_over_R() = 0;

		/**
		* @brief returns the names of the species
		*/
		const std::vector<std::string>& names() const { return names_; }

		/**
		* @brief returns the index (1-based) of the species with the requested name
		         if the species is not present, a fatal error will be reported
		*/
		unsigned int IndexOfSpecies(const std::string name) const;

		/**
		* @brief returns the index (1-based) of the species with the requested name
		         if the species is not present, the returned index is equal to 0
		*/
		unsigned int IndexOfSpeciesWithoutError(const std::string name) const;

		/**
		* @brief returns the index (1-based) of the element with the requested name
		         if the element is not present, a fatal error will be reported
		*/
		unsigned int IndexOfElement(const std::string name) const;

		/**
		* @brief returns the index (1-based) of the element with the requested name
		         if the element is not present, the returned index is equal to 0
		*/
		unsigned int IndexOfElementWithoutError(const std::string name) const;


		/**
		*@brief Import elements from a file in XML format
		*/
		void ImportElementsFromXMLFile(rapidxml::xml_document<>& doc);

		double GetMixtureFractionFromMassFractions(
		const std::vector<double>& mass,
		const std::vector<std::string>& fuel_names, 	const std::vector<double>& mass_fuel,
		const std::vector<std::string>& oxidizer_names,	const std::vector<double>& mass_oxidizer);

		std::vector<double> GetMoleFractionsFromEquivalenceRatio( const double equivalence_ratio,
						const std::vector<std::string>& fuel_names, const std::vector<double>& moles_fuel,
						const std::vector<std::string>& oxidizer_names,	const std::vector<double>& moles_oxidizer );

		double GetLocalEquivalenceRatio( 	const std::vector<double>& moles, const std::vector<double>& moles_st,
							const std::vector<std::string>& fuel_names);

		// Provisional
		//virtual void Test(const int nLoops, const double& T, int* index) = 0;

	protected:

		unsigned int nspecies_;					//!< number of species

		std::vector<std::string> names_;		//!< names of the species
		std::vector<double> MW__;				//!< molecular weights of the species

		double T_;								//!< temperature [K]
		double P_;								//!< pressure [Pa]

		Eigen::MatrixXd atomic_composition_;	//!< atomic composition
		std::vector<std::string> elements_;		//!< names of elements
		std::vector<double> mw_elements__;		//!< molecular weights of elements
	};
}

#include "ThermodynamicsMap.hpp"

#endif // OpenSMOKE_ThermodynamicsMap_H
