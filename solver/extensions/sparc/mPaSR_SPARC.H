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
#include <iostream>

{
	if (counterToUpdateChemistry == numberIterationsToUpdateChemistry)
	{
		Info << " * swich_ODE = "           << switch_ODE       << endl;// bool switch_ODE(false);// default true = ODE rates, false = laminar rates
        Info << " * switch_LAPACK = "       << switch_LAPACK    << endl;// bool switch_LAPACK(true);   // default true: true = LAPACK solver, false = EIGEN solver
        Info << " * switch_Jconst = "       << switch_Jconst    << endl;// bool switch_Jconst(false); 	// default false: true = constrained Jacobian, false = kinetic Jacobian
        Info << " * switch_Ekeep = "        << switch_Ekeep     << endl;// bool switch_Ekeep(true); // default true: true = keep conserved E subspace, false = set E modal fractions to 0
        Info << " * activation tolerance = " << activation_threshold << endl;// double activation_threshold = 1.e-8;// default= 1.e-32, as a rate threshold to perform CSP decomposition

        // Thermo-chemical map
        CharacteristicChemicalTimes* chemicalTimes = new CharacteristicChemicalTimes(*thermodynamicsMapXML, *kineticsMapXML);
        int ns      = thermodynamicsMapXML->NumberOfSpecies();// nb of species
        int nelem   = thermodynamicsMapXML->elements().size();// nb of elements

        //OpenSMOKE::OpenSMOKEVectorDouble h_species(ns);

        const std::vector<std::string> species_names = thermodynamicsMapXML->NamesOfSpecies();
        const scalarField& TCells = thermo.T().internalField();
        const scalarField& pCells = thermo.p().internalField();
        const scalarField& rhoCells = rho.internalField();
        //OpenSMOKE::OpenSMOKEVectorDouble vec_Y(ns);


        // Fractions
        #include "propertiesFineStructures.H"
        const double tauc_max   = max(threshold_chemical_time,1e+3);

        // // Heat release
        OpenSMOKE::OpenSMOKEVectorDouble h_species(ns);
        scalarField& QdotCells = Qdot.internalFieldRef(); //v12
        scalarField& dQCells = dQ.internalFieldRef(); //v12

        const scalarField& clusterIndexCells = clusterIndex.internalFieldRef(); //v12

        Info << "Check initializataion"  << endl;

        //--------------------- Looping over cells ------------------------//
        forAll(TCells, celli)
        {
            //Info << "start loop cells"   << endl;
            const scalar T_cell     = TCells[celli];
            const scalar p_cell     = pCells[celli];
            const scalar rho_cell   = rhoCells[celli];

            const unsigned int iCluster = clusterIndexCells[celli];
            const unsigned int NClocal  = vector_thermodynamicsMapXML[iCluster]->NumberOfSpecies();
            const unsigned int NRlocal  = vector_kineticsMapXML[iCluster]->NumberOfReactions(); 
            
            const int unsigned nsCluster = NClocal;
            const int unsigned ne = nsCluster+3;
            
            //CharacteristicChemicalTimes* chemicalTimes_SPARC = new CharacteristicChemicalTimes(*vector_thermodynamicsMapXML[iCluster], *vector_kineticsMapXML[iCluster]);
            // OpenSMOKE::ThermodynamicsMap_CHEMKIN thermoMap = *vector_thermodynamicsMapXML[iCluster];
            // OpenSMOKE::KineticsMap_CHEMKIN kineticsMap = *vector_kineticsMapXML[iCluster];

            // Info << "check thermoMap " << thermoMap.NumberOfSpecies() << endl;
            // Info << "check kineticsMap " << kineticsMap.NumberOfReactions() << endl;

            // Jacobian decomposition
            Eigen::MatrixXd Jc(nsCluster, nsCluster);
            std::vector<double> lambda_real(nsCluster);
            std::vector<double> lambda_imag(nsCluster);
            std::vector<double> lambda_mod(nsCluster);
            std::vector<double> lambda_mod_sorted(nsCluster);
            Eigen::MatrixXd vr(nsCluster, nsCluster);
            Eigen::MatrixXd vr_sorted(nsCluster, nsCluster);
            Eigen::MatrixXd vl(nsCluster, nsCluster);


            Eigen::VectorXd tauc(nsCluster);
            Eigen::VectorXd kappa(nsCluster);

            // Rates
            Eigen::VectorXd R_base(nsCluster);
            Eigen::VectorXd R_mPaSR(nsCluster);
            Eigen::VectorXd f_vector(nsCluster);
            Eigen::MatrixXd kappaVR(nsCluster,nsCluster);
            // ODE
            Eigen::VectorXd R_ODE(nsCluster);
            Eigen::VectorXd yStar0(ne);
            Eigen::VectorXd yStarf(ne);
            Eigen::VectorXd yMin(ne);
            Eigen::VectorXd yMax(ne);

            for(unsigned int i=0;i<nsCluster;i++)	yMin(i) = 0.; yMin(ne-3) = 200.;	yMin(ne-2) = 0.;	yMin(ne-1) = 0.;
            for(unsigned int i=0;i<nsCluster;i++)	yMax(i) = 1.; yMax(ne-3) = 5000.;	yMax(ne-2) = 1.e16;	yMax(ne-1) = 1.e16;

            // Laminar
            Eigen::VectorXd R_laminar(nsCluster);
            OpenSMOKE::OpenSMOKEVectorDouble vec_Y_laminar(nsCluster);
            OpenSMOKE::OpenSMOKEVectorDouble x(nsCluster);
            OpenSMOKE::OpenSMOKEVectorDouble c(nsCluster);
            OpenSMOKE::OpenSMOKEVectorDouble dc_over_dt_laminar(nsCluster);

            OpenSMOKE::OpenSMOKEVectorDouble vec_Y(nsCluster);

            // Heat release
            //OpenSMOKE::OpenSMOKEVectorDouble h_species(nsCluster);
            
            //Info << "Check declarations"  << endl;

            if (TCells[celli] > minTemperature4Chemistry)
            {
                activation_map[celli] = 1.;

                //--------------- Base rates ------------------//
                if (switch_ODE == true)
                {
                    //--------------- Residence Time ------------------//
                    double tf       = 0.;

                    if (iLES == true) // LES
                    {
                        if (tauStarType == DT_LES)
                        {
                            tf = runTime.deltaTValue();
                        }
                        else if (tauStarType == MIN_DTLES_TAUMIX)
                        {
                            tf = min(tauMix[celli], runTime.deltaTValue());
                        }
                        else if (tauStarType == TAUMIX)
                        {
                            tf = tauMix[celli];
                        }
                    }
                    else // RANS
                        tf = tauMix[celli];

                    //--------------- ODE rates ------------------//
                    for(unsigned int i=0;i<nsCluster;i++)
                    {
                        const unsigned int j = classifier_indices[iCluster](i);
                        yStar0(i) = Y[j].internalField()[celli];
                    }
                    yStar0(ne-3) = T_cell;
                    yStar0(ne-2) = tf;
                    yStar0(ne-1) = 0.;

                    vector_fine_structure_ode_pfr[iCluster].Set(p_cell);

                    const double t0 = 0.;

                    // Solve the ODE system
                    {
                        if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
                        {
                            #include "extensions/sparc/solvePFR_SPARC.H"
                        }
                    }

                    for(unsigned int i=0;i<nsCluster;i++)
                    {
                        const unsigned int j = classifier_indices[iCluster](i);
                        YStar_[j].internalFieldRef()[celli] = yStarf(i); //v12
                    }
                    TStar_[0].internalFieldRef()[celli] = yStarf(ne-3); //v12

                    for(unsigned int i=0;i<nsCluster;i++)
                    {
                        const unsigned int j = classifier_indices[iCluster](i);
                        R_ODE[i] = rho_cell*(YStar_[j].internalField()[celli]-Y[j].internalField()[celli])/tf;//[kg.m-3.s-1]
                    }
                    R_base = R_ODE;
                }
                else
                {
                    //--------------- Laminar rates ------------------//
                    double MW;
                    //Info << "Check laminar 1 start"  << endl;
                    for(unsigned int i=0;i<nsCluster;i++)
                    {
                        const unsigned int j = classifier_indices[iCluster](i);
                        vec_Y_laminar[i+1]      = Y[j].internalField()[celli];
                    }

                    //thermodynamicsMapXML->MoleFractions_From_MassFractions(x.GetHandle(), MW, vec_Y_laminar.GetHandle());
                    //thermoMap.MoleFractions_From_MassFractions(x.GetHandle(), MW, vec_Y_laminar.GetHandle());
                    vector_thermodynamicsMapXML[iCluster]->MoleFractions_From_MassFractions(x.GetHandle(), MW, vec_Y_laminar.GetHandle());
                    const double cTot = p_cell/(PhysicalConstants::R_J_kmol * T_cell);
                    Product(cTot, x, &c);
                    //Info << "Check laminar 2 start"  << endl;
                    vector_thermodynamicsMapXML[iCluster]->SetTemperature(T_cell);
                    //Info << "Check laminar setT"  << endl;
                    vector_thermodynamicsMapXML[iCluster]->SetPressure(p_cell);
                    //Info << "Check laminar setP"  << endl;
                    vector_kineticsMapXML[iCluster]->SetTemperature(T_cell);
                    //Info << "Check laminar kinsetT"  << endl;
                    vector_kineticsMapXML[iCluster]->SetPressure(p_cell);
                    //Info << "Check laminar kinsetP"  << endl;
                    vector_kineticsMapXML[iCluster]->KineticConstants();
                    //Info << "Check laminar kinConst"  << endl;
                    vector_kineticsMapXML[iCluster]->ReactionRates(c.GetHandle());
                    //Info << "Check laminar kinC"  << endl;
	                vector_kineticsMapXML[iCluster]->FormationRates(dc_over_dt_laminar.GetHandle());
                    // thermodynamicsMapXML->SetTemperature(T_cell);
                    // thermodynamicsMapXML->SetPressure(p_cell);
                    // kineticsMapXML->SetTemperature(T_cell);
                    // kineticsMapXML->SetPressure(p_cell);
                    // kineticsMapXML->KineticConstants();
                    // kineticsMapXML->ReactionRates(c.GetHandle());
                    // kineticsMapXML->FormationRates(dc_over_dt_laminar.GetHandle());
                    //Info << "Check laminar start"  << endl;

                    for(unsigned int i=0;i<nsCluster;i++)
                    {
                        //const unsigned int j = classifier_indices[iCluster](i);
                        //R_laminar[i] = dc_over_dt_laminar[j+1]*thermodynamicsMapXML->MW(j);//[kg.m-3.s-1]
                        R_laminar[i] = dc_over_dt_laminar[i+1]*vector_thermodynamicsMapXML[iCluster]->MW(i);//[kg.m-3.s-1]
                    }
                    R_base = R_laminar;
                }
                //Info << "Check laminar ends"  << endl;
                double R_base_max_abs = max(abs(R_base.minCoeff()),abs(R_base.maxCoeff()));


                if (R_base_max_abs > activation_threshold)// Condition to perform modal decomposition
                {
                    activation_map[celli] = 2.;

                    //--------------- Jacobian ------------------//
                    // Kinetic Jacobian         -> mPaSR_computeNumericalJacobian
                    // Constrained Jacobian     -> mPaSR_computeNumericalJacobian_const
                    //Info << "iCluster: " << iCluster << " valid? " << (iCluster < classifier_indices.size()) << endl;
                    Jc.setZero();
                    //OpenSMOKE::ThermodynamicsMap_CHEMKIN thermoMap = *vector_thermodynamicsMapXML[iCluster];

                    for(unsigned int i=0;i<nsCluster;i++)
                    {
                        const unsigned int j = classifier_indices[iCluster](i);
                        vec_Y[i]      = Y[j].internalField()[celli];
                    }

                    // Info << "Check jacobian decomposition start"  << endl;
                    // Info << "Size of nsCluster: " << nsCluster << endl;
                    // Info << "Size of Jc: " << Jc.size() << endl;
                    // Info << "Size of thermoMap: " << vector_thermodynamicsMapXML[iCluster]->NumberOfSpecies() << endl;
                    // Info << "Size of thermoMap: " << thermoMap.NumberOfSpecies() << endl;
                    // Info << "Size of vec_Y: " << vec_Y.Size() << endl;
                    // for (unsigned int i = 0; i < vec_Y.Size();  i++)
                    // {
                    //     Info << vec_Y[i] << " ";
                    // }
                    // Info << endl;
                    if (switch_Jconst == true)
                    {
                        Jc = chemicalTimes->mPaSR_computeNumericalJacobian_const(T_cell, p_cell, vec_Y);
                    }
                    else
                    {
                        Jc = chemicalTimes->SPARC_computeNumericalJacobian(T_cell, p_cell, vec_Y, *vector_thermodynamicsMapXML[iCluster], *vector_kineticsMapXML[iCluster]);
                        //Jc = chemicalTimes_SPARC->mPaSR_computeNumericalJacobian(T_cell, p_cell, vec_Y);
                    }

                    //Info << "Size of Jc after calculation: " << Jc.size() << endl;
                    //Info << "Check jacobian decomposition end"  << endl;

                    if (switch_LAPACK == true)
                    {
                        //--------------- LAPACK solver ------------------//
                        int info = LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'V', nsCluster, Jc.data(), nsCluster, lambda_real.data(), lambda_imag.data(), vl.data(), nsCluster, vr.data(), nsCluster);
                        //Info << "Check lapack solver in" << info << endl;

                        if (info != 0)
                        {
                            Info << " *  OnTheFlyCEMA::Calculate: LAPACKE_dgeev failure -> status:" << info << ", using EIGEN solver instead!" << endl;
                            Eigen::EigenSolver<Eigen::MatrixXd> eigenJc(Jc);
                            //Info << "Check eigen lapack solver " << eigenJc.eigenvalues().size() << endl;

                            for(unsigned int m=0;m<nsCluster;m++)
                            {
                                lambda_real[m] = eigenJc.eigenvalues()[m].real();
                                lambda_imag[m] = eigenJc.eigenvalues()[m].imag();
                            }
                        }
                    }
                    else
                    {
                        //---------------- EIGEN solver ------------------//
                        Eigen::EigenSolver<Eigen::MatrixXd> eigenJc(Jc);

                        for(unsigned int m=0;m<nsCluster;m++)
                        {
                            lambda_real[m] = eigenJc.eigenvalues()[m].real();
                            lambda_imag[m] = eigenJc.eigenvalues()[m].imag();
                        }
                    }
                    //Info << "Check eigen lapack solver " << lambda_real.size() << endl;
                    //Info << "Check eigen lapack solver " << lambda_imag.size() << endl;

                    for(unsigned int m=0;m<nsCluster;m++)
                    {
                        lambda_mod[m] = std::sqrt(lambda_real[m] * lambda_real[m] + lambda_imag[m] * lambda_imag[m]);
                    }
                    // Info << "Check unsorted lambdas " << lambda_mod.size() << endl;
                    // for (unsigned int i = 0; i < lambda_mod.size();  i++)
                    // {
                    //     Info << lambda_mod[i] << " ";
                    // }
                    // Info << endl;
                    std::vector<size_t> sorted_indices(nsCluster);
                    sorted_indices = OpenSMOKE::SortAndTrackIndicesDecreasing(lambda_mod);

                    //Info << "Check sorting end"  << endl;

                    for(unsigned int m=0;m<nsCluster;m++)
                    {
                        lambda_mod_sorted[m] = lambda_mod[sorted_indices[m]];
                    }
                    lambda_mod = lambda_mod_sorted;// returning the sorted version

                    vr_sorted = vr;
                    for(unsigned int m=0;m<nsCluster;m++)
                    {
                        for(unsigned int i=0;i<nsCluster;i++)
                        {
                            vr_sorted(i,m) = vr(i,sorted_indices[m]);
                        }
                    }
                    vr = vr_sorted;// returning the sorted version

                    if (vr.determinant() != 0)
                        vl = vr.inverse();
                    else
                    {
                        vr.setZero();
                        vl.setZero();
                    }
                    
                    //Info << "Check sorting vr end"  << endl;

                    //--------------- Fractions ------------------//
                    for(unsigned int m=0;m<nsCluster;m++)
                    {
                        tauc[m] = min(1./(lambda_mod[m] + SMALL),tauc_max);
                        kappa[m] = tauc[m] / (tauc[m] + tauMix[celli]);
                        if (switch_Ekeep == false)// Setting to 0 the modal fractions of the conserved space E
                        {
                            if (m > (nsCluster-ne-1))
                                kappa[m] = 0.;
                        }
                        modal_tauc_[m].internalFieldRef()[celli] = tauc[m]; //v12
                        modal_kappa_[m].internalFieldRef()[celli] = kappa[m]; //v12
                    }
                    //Info << "Check fractions end"  << endl;
                    //--------------- mPaSR rates ------------------//
                    kappaVR     = vr * kappa.asDiagonal();  // Modes
                    f_vector    = vl * R_base;              // Amplitudes
                    R_mPaSR     = kappaVR * f_vector;       // mPaSR rates

                    for(unsigned int i=0;i<nsCluster;i++)
                    {
                        const unsigned int j = classifier_indices[iCluster](i);
                        R_[j].internalFieldRef()[celli] = R_mPaSR[i];//[kg.m-3.s-1] //v12
                    }
                    //Info << "Check mpasr rates end"  << endl;
                    //Info << "Check R_mPaSR size " << R_mPaSR.size() << endl;

                }
                else
                {
                    R_mPaSR = R_base;// laminar

                    for(unsigned int i=0;i<nsCluster;i++)
                    {
                        const unsigned int j = classifier_indices[iCluster](i);
                        R_[j].internalFieldRef()[celli] = R_mPaSR[i];//[kg.m-3.s-1] //v12
                    }
                }
            }
            else
            {
                activation_map[celli] = 0.;
                // for(unsigned int i=0;i<nsCluster;i++)
                // {
                //     const unsigned int j = classifier_indices[iCluster](i);
                //     R_[j].internalFieldRef()[celli] *= 0.;//[kg.m-3.s-1] //v12
                // }
                for(unsigned int i=0;i<ns;i++)
                {
                    //const unsigned int j = classifier_indices[iCluster](i);
                    R_[i].internalFieldRef()[celli] *= 0.;//[kg.m-3.s-1] //v12
                }
            }
            //Info << "Check R size " << R_.size()  << endl;  

            //--------------- Heat release ------------------//
            // thermodynamicsMapXML->SetTemperature(T_cell);
            // thermodynamicsMapXML->SetPressure(p_cell);
            // thermodynamicsMapXML->hMolar_Species(h_species.GetHandle());// [J.kmol-1]
            // thermoMap.SetTemperature(T_cell);
            // thermoMap.SetPressure(p_cell);
            // thermoMap.hMolar_Species(h_species.GetHandle()); // [J.kmol-1]

            vector_thermodynamicsMapXML[iCluster]->SetTemperature(T_cell);
            vector_thermodynamicsMapXML[iCluster]->SetPressure(p_cell);
            vector_thermodynamicsMapXML[iCluster]->hMolar_Species(h_species.GetHandle());
            //Info << "Check heat release maps set"  << endl;
            
            double HeatRelease = 0.;
            for(unsigned int i=0;i<nsCluster;i++)
            {
                const unsigned int j = classifier_indices[iCluster](i);
                HeatRelease += -h_species[j+1]/thermodynamicsMapXML->MW(j) * R_[j].internalField()[celli];// [W.m-3]
            }
            // Info << "Check heat release set"  << endl;
            // Info << "Check heat release size " << HeatRelease  << endl;
            // Info << "Check heat release type " << typeid(HeatRelease).name()  << endl;
            // Info << "celli: " << celli << ", QdotCells size: " << QdotCells.size() << ", dQCells size: " << dQCells.size() << ", mesh.V size: " << mesh.V().size() << endl;
            // Info << "Check QdotCells type " << typeid(QdotCells).name()  << endl;

            QdotCells[celli] = HeatRelease;
            dQCells[celli] = HeatRelease * mesh.V()[celli];

            //Info << "Check heat release calculated"   << endl;
        }
        
        //Info << "Check out of cell loop"   << endl;
        // Correct boundary conditions
        Qdot.correctBoundaryConditions();
        dQ.correctBoundaryConditions();
        Info << " * min/max(Q) = " << min(Qdot).value() << ", " << max(Qdot).value() << endl;
    }
}