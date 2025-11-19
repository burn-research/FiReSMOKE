/*-----------------------------------------------------------------------*\
|                                                                         |
|   ╭╮╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╭━━━┳━╮╭━┳━━━┳╮╭━┳━━━╮                               |
|   ┃┃╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱┃╭━╮┃┃╰╯┃┃╭━╮┃┃┃╭┫╭━━╯                               |
|   ┃┃╭━━┳╮╭┳┳━╮╭━━┳━┫╰━━┫╭╮╭╮┃┃╱┃┃╰╯╯┃╰━━┳╮╱╭╮                           |
|   ┃┃┃╭╮┃╰╯┣┫╭╮┫╭╮┃╭┻━━╮┃┃┃┃┃┃┃╱┃┃╭╮┃┃╭━┳╯╰┳╯╰╮                          |
|   ┃╰┫╭╮┃┃┃┃┃┃┃┃╭╮┃┃┃╰━╯┃┃┃┃┃┃╰━╯┃┃┃╰┫╰━┻╮╭┻╮╭╯                          |
|   ╰━┻╯╰┻┻┻┻┻╯╰┻╯╰┻╯╰━━━┻╯╰╯╰┻━━━┻╯╰━┻━━━┻╯╱╰╯                           |
|                                                                         |
|   Authors: Alberto Cuoci                                                |
|                                                                         |
|   Contacts: Alberto Cuoci                                               |
|   email: alberto.cuoci@polimi.it                                        |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano (Italy)                      |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of laminarSMOKE++ solver.                           |
|                                                                         |
|   License                                                               |
|                                                                         |
|   Copyright(C) 2020, 2021 Alberto Cuoci                                 |
|   laminarSMOKE++ is free software: you can redistribute it and/or       |
|   modify it under the terms of the GNU General Public License           |
|   as published by the Free Software Foundation, either version 3 of     |
|   the License, or (at your option) any later version.                   |
|                                                                         |
|   laminarSMOKE++ is distributed in the hope that it will be useful,     |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with laminarSMOKE++.                                            |
|   If not, see <http://www.gnu.org/licenses/>.                           |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef BatchReactorHomogeneousConstantPressure_H
#define	BatchReactorHomogeneousConstantPressure_H

namespace OpenSMOKE
{
	class ThermodynamicsMap_CHEMKIN;
	class KineticsMap_CHEMKIN;
}

class BatchReactorHomogeneousConstantPressure
{
public:

	//- Default constructor
	BatchReactorHomogeneousConstantPressure(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& 		thermodynamicsMapXML, 
			OpenSMOKE::KineticsMap_CHEMKIN& 		kineticsMapXML);

	//- Set reactor pressure (in Pa)
	//void SetReactor( const double P0 );
	void SetReactor( const double P0 )
	{
		P0_ 		= P0;
	}

	//- Set reactor mass fractions
	void SetMassFractions( const OpenSMOKE::OpenSMOKEVectorDouble& omega );
	
	//- Enable the solution of energy equation
	void SetEnergyEquation(const bool flag) { energyEquation_ = flag; }
	
	//- Enable mass fraction checkings (default: false)
	void SetCheckMassFractions( const bool flag ) { checkMassFractions_ = flag; };

	//- Returns the reactor temperature
	double T() const;

	//- Returns the vector of formation rates of species (in kmol/m3/s)
	const OpenSMOKE::OpenSMOKEVectorDouble& R() const { return R_; }

	//- Returns the reaction heat (in W/m3)
	double QR() const { return QR_; }

	//- Return the total number of equations
	//unsigned int NumberOfEquations() const;
	unsigned int NumberOfEquations() const { return number_of_equations_; }

	// OFv10 Set debug
	void SetDebug(	const double flag)
	{
		debug_ = flag;
	}

	//- Returns the differential equations
	int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

	//- Print function
	int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);


private:

	OpenSMOKE::ThermodynamicsMap_CHEMKIN& 	thermodynamicsMapXML_;
	OpenSMOKE::KineticsMap_CHEMKIN& 	kineticsMapXML_;

	double T_;			//!< temperature (in K)
	double P0_;			//!< pressure (in Pa)
	
	double rho_;			//!< density (in kg/m3)
	double cTot_;			//!< total concentration (in kmol/m3)
	double MW_;			//!< molecular weight (in kg/kmol)
	double CpMixMass_;		//!< constant pressure specific heat (in J/kg/K)
	double QR_;			//!< heat release (in W/m3)

	unsigned int number_of_gas_species_;		//!< number of species
	unsigned int number_of_reactions_;		//!< number of reactions
	unsigned int number_of_equations_;		//!< numer of equations

	OpenSMOKE::OpenSMOKEVectorDouble omega_;		//!< mass fractions
	OpenSMOKE::OpenSMOKEVectorDouble x_;		//!< mole fractions
	OpenSMOKE::OpenSMOKEVectorDouble c_;		//!< concentrations (in kmol/m3)
	OpenSMOKE::OpenSMOKEVectorDouble R_;		//!< formation rates (in kmol/m3/s)

	bool checkMassFractions_;	//!< if true, mass fractions are normalized
	bool energyEquation_;		//!< if true, the energy equation is solved
	bool debug_;			//!< if true, additional debug (time-consuming) operations are carried out
};

#include "BatchReactorHomogeneousConstantPressure.C"

#endif // BatchReactorHomogeneousConstantPressure_H