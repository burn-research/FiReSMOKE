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

// OpenSMOKE
#include "OpenSMOKE_Definitions.h"
#include <string>
#include <iostream>
#include <numeric>
#include <Eigen/Dense>

// Base classes
#include "kernel/thermo/ThermoPolicy_CHEMKIN.h"
#include "kernel/kinetics/ReactionPolicy_CHEMKIN.h"
#include "math/PhysicalConstants.h"
#include "math/OpenSMOKEUtilities.h"

// DRG
#include "DRG.H"

// Maps
#include "maps/ThermodynamicsMap_CHEMKIN.h"
#include "maps/KineticsMap_CHEMKIN.h"


BatchReactorHomogeneousConstantPressure::BatchReactorHomogeneousConstantPressure
(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMapXML, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMapXML ):
	thermodynamicsMapXML_(thermodynamicsMapXML), 
	kineticsMapXML_(kineticsMapXML)
{
	number_of_gas_species_ = thermodynamicsMapXML_.NumberOfSpecies();
	number_of_reactions_ = kineticsMapXML_.NumberOfReactions();
	number_of_equations_ = number_of_gas_species_ + 1 ;	// species and temperature
	QR_ = 0.;

	ChangeDimensions(number_of_gas_species_, &omega_, true);
	ChangeDimensions(number_of_gas_species_, &x_, true);
	ChangeDimensions(number_of_gas_species_, &c_, true);
	ChangeDimensions(number_of_gas_species_, &R_, true);
	//ChangeDimensions(number_of_reactions_, 	 &rStar_, true);

	checkMassFractions_ = false;
	energyEquation_ = true;
	debug_ = false;
}

void BatchReactorHomogeneousConstantPressure::SetMassFractions( const OpenSMOKE::OpenSMOKEVectorDouble& omega )
{
	omega_ = omega;
}

double BatchReactorHomogeneousConstantPressure::T() const
{
	return T_;
}

int BatchReactorHomogeneousConstantPressure::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
{
	// Recover mass fractions
	if (checkMassFractions_ == true)
	{	
		for(unsigned int i=0;i<number_of_gas_species_;i++)
			omega_(i) = std::max(y[i], 0.);
	}
	else
	{
		for(unsigned int i=0;i<number_of_gas_species_;i++)
			omega_(i) = y[i];
	}

	// Recover temperature
	T_ = y[number_of_gas_species_]; 

	// Calculates the pressure and the concentrations of species
	thermodynamicsMapXML_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());
	cTot_ = P0_/PhysicalConstants::R_J_kmol/T_;
	for (unsigned int i = 0; i < x_.Size(); ++i)
	{
    c_(i) = cTot_ * x_(i); //	c_ = cTot_*x_;
	}
    	rho_ = cTot_*MW_;

	// Calculates thermodynamic properties
	thermodynamicsMapXML_.SetTemperature(T_);
	thermodynamicsMapXML_.SetPressure(P0_);

	// Calculates kinetics
	kineticsMapXML_.SetTemperature(T_);
	kineticsMapXML_.SetPressure(P0_);
	kineticsMapXML_.KineticConstants();
	kineticsMapXML_.ReactionRates(c_.GetHandle());
	kineticsMapXML_.FormationRates(R_.GetHandle());

	// Species equations
	for (unsigned int i=0;i<number_of_gas_species_;++i)	
		dy(i) = thermodynamicsMapXML_.MW(i)*R_(i)/rho_;

    	// Energy equation
    	dy(number_of_gas_species_) = 0.;     
    	if (energyEquation_ == true)
    	{
		CpMixMass_ = thermodynamicsMapXML_.cpMolar_Mixture_From_MoleFractions(x_.GetHandle()) / MW_;
		QR_ = kineticsMapXML_.HeatRelease(R_.GetHandle());
	
		dy(number_of_gas_species_)  = QR_ / (rho_*CpMixMass_);
	}

	return 0;
}

int BatchReactorHomogeneousConstantPressure::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
	return 0;
}