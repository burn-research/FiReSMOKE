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

#include <iostream>

{
    // Thermo-chemical map
    CharacteristicChemicalTimes* chemicalTimes = new CharacteristicChemicalTimes(*thermodynamicsMapXML, *kineticsMapXML);

    int ns      = thermodynamicsMapXML->NumberOfSpecies();// nb of species
    int nelem   = thermodynamicsMapXML->elements().size();// nb of elements

    const std::vector<std::string> species_names = thermodynamicsMapXML->NamesOfSpecies();
    const scalarField& TCells = thermo.T().internalField();
    const scalarField& pCells = thermo.p().internalField();
    const scalarField& rhoCells = rho.internalField();

    // Laminar
    Eigen::VectorXd R_base(ns);
    Eigen::VectorXd R_laminar(ns);
    OpenSMOKE::OpenSMOKEVectorDouble vec_Y_laminar(ns);
    OpenSMOKE::OpenSMOKEVectorDouble x(ns);
    OpenSMOKE::OpenSMOKEVectorDouble c(ns);
    OpenSMOKE::OpenSMOKEVectorDouble dc_over_dt_laminar(ns);

    // Heat release
    OpenSMOKE::OpenSMOKEVectorDouble h_species(ns);
    scalarField& QdotCells = Qdot.internalFieldRef(); //v12
    scalarField& dQCells = dQ.internalFieldRef(); //v12

    // Getting kappa
    #include "propertiesFineStructures.H"
    const scalarField& kappaCells = kappa.internalField();

    //--------------------- Looping over cells ------------------------//
    forAll(TCells, celli)
    {
        const scalar T_cell     = TCells[celli];
        const scalar p_cell     = pCells[celli];
        const scalar rho_cell   = rhoCells[celli];

		if (T_cell > minTemperature4Chemistry)
		{
            activation_map[celli] = 1.;
            
			//--------------- Laminar rates ------------------//
			double MW;
			for(unsigned int i=0;i<ns;i++)
				vec_Y_laminar[i+1]      = Y[i].internalField()[celli];
			thermodynamicsMapXML->MoleFractions_From_MassFractions(x.GetHandle(), MW, vec_Y_laminar.GetHandle());
			const double cTot = p_cell/(PhysicalConstants::R_J_kmol * T_cell);
			Product(cTot, x, &c);

			thermodynamicsMapXML->SetTemperature(T_cell);
			thermodynamicsMapXML->SetPressure(p_cell);
			kineticsMapXML->SetTemperature(T_cell);
			kineticsMapXML->SetPressure(p_cell);
			kineticsMapXML->KineticConstants();
			kineticsMapXML->ReactionRates(c.GetHandle());
			kineticsMapXML->FormationRates(dc_over_dt_laminar.GetHandle());

			for(unsigned int i=0;i<ns;i++)
				R_laminar[i] = dc_over_dt_laminar[i+1]*thermodynamicsMapXML->MW(i);//[kg.m-3.s-1]

            for(unsigned int i=0;i<ns;i++)
                R_[i].internalFieldRef()[celli] = kappaCells[celli]*R_laminar[i];//[kg.m-3.s-1] //v12
		}
		else
		{
            activation_map[celli] = 0.;
			for(unsigned int i=0;i<ns;i++)
				R_[i].internalFieldRef()[celli] *= 0.;//[kg.m-3.s-1] //v12
		}

        //--------------- Heat release ------------------//
        thermodynamicsMapXML->SetTemperature(T_cell);
        thermodynamicsMapXML->SetPressure(p_cell);
        thermodynamicsMapXML->hMolar_Species(h_species.GetHandle());// [J.kmol-1]

        double HeatRelease = 0.;
        for(unsigned int i=0;i<ns;i++)
            HeatRelease += -h_species[i+1]/thermodynamicsMapXML->MW(i) * R_[i].internalField()[celli];// [W.m-3]

        QdotCells[celli] = HeatRelease;
        dQCells[celli] = HeatRelease * mesh.V()[celli];

		//--------------- Laminar Diffusion ------------------//
		if (laminarDiffusion == true)
		{		
			for(unsigned int i=0;i<ns;i++)
				HeSpecies_[i].internalFieldRef()[celli] = -h_species[i+1]/thermodynamicsMapXML->MW(i); //v12
		}
	}

    // Correct boundary conditions
    Qdot.correctBoundaryConditions();
    dQ.correctBoundaryConditions();
    Info << " * min/max(Q) = " << min(Qdot).value() << ", " << max(Qdot).value() << endl;
}