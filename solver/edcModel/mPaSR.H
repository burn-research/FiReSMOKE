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
    Info << " * swich_ODE = "           << switch_ODE       << endl;// bool switch_ODE(false);// default true = ODE rates, false = laminar rates
    Info << " * switch_LAPACK = "       << switch_LAPACK    << endl;// bool switch_LAPACK(true);   // default true: true = LAPACK solver, false = EIGEN solver
    Info << " * switch_Jconst = "       << switch_Jconst    << endl;// bool switch_Jconst(false); 	// default false: true = constrained Jacobian, false = kinetic Jacobian
    Info << " * switch_Ekeep = "        << switch_Ekeep     << endl;// bool switch_Ekeep(true); // default true: true = keep conserved E subspace, false = set E modal fractions to 0
    Info << " * activation tolerance = " << activation_threshold << endl;// double activation_threshold = 1.e-8;// default= 1.e-32, as a rate threshold to perform CSP decomposition

    // Thermo-chemical map
    CharacteristicChemicalTimes* chemicalTimes = new CharacteristicChemicalTimes(*thermodynamicsMapXML, *kineticsMapXML);

    int ns      = thermodynamicsMapXML->NumberOfSpecies();// nb of species
    int nelem   = thermodynamicsMapXML->elements().size();// nb of elements

    const std::vector<std::string> species_names = thermodynamicsMapXML->NamesOfSpecies();
    const scalarField& TCells = thermo.T().internalField();
    const scalarField& pCells = thermo.p().internalField();
    const scalarField& rhoCells = rho.internalField();
    OpenSMOKE::OpenSMOKEVectorDouble vec_Y(ns);
 
    // Jacobian decomposition
    Eigen::MatrixXd Jc(ns, ns);
    std::vector<double> lambda_real(ns);
    std::vector<double> lambda_imag(ns);
    std::vector<double> lambda_mod(ns);
    std::vector<double> lambda_mod_sorted(ns);
    Eigen::MatrixXd vr(ns, ns);
    Eigen::MatrixXd vr_sorted(ns, ns);
    Eigen::MatrixXd vl(ns, ns);

    // Fractions
    #include "propertiesFineStructures.H"
    const double tauc_max   = max(threshold_chemical_time,1e+3);
    Eigen::VectorXd tauc(ns);
    Eigen::VectorXd kappa(ns);
    
    // Rates
    Eigen::VectorXd R_base(ns);
    Eigen::VectorXd R_mPaSR(ns);
    Eigen::VectorXd f_vector(ns);
    Eigen::MatrixXd kappaVR(ns,ns);

    // ODE
    Eigen::VectorXd R_ODE(ns);
    const int unsigned ne = ne_total;// number of transported equations
    Eigen::VectorXd yStar0(ne);
    Eigen::VectorXd yStarf(ne);
    Eigen::VectorXd yMin(ne);
    Eigen::VectorXd yMax(ne);

    for(unsigned int i=0;i<ns;i++)	yMin(i) = 0.; yMin(ne-3) = 200.;	yMin(ne-2) = 0.;	yMin(ne-1) = 0.;
    for(unsigned int i=0;i<ns;i++)	yMax(i) = 1.; yMax(ne-3) = 5000.;	yMax(ne-2) = 1.e16;	yMax(ne-1) = 1.e16;

    // Laminar
    Eigen::VectorXd R_laminar(ns);
    OpenSMOKE::OpenSMOKEVectorDouble vec_Y_laminar(ns);
    OpenSMOKE::OpenSMOKEVectorDouble x(ns);
    OpenSMOKE::OpenSMOKEVectorDouble c(ns);
    OpenSMOKE::OpenSMOKEVectorDouble dc_over_dt_laminar(ns);

    // Heat release
    OpenSMOKE::OpenSMOKEVectorDouble h_species(ns);
    scalarField& QdotCells = Qdot.internalFieldRef();
    scalarField& dQCells = dQ.internalFieldRef();

    //--------------------- Looping over cells ------------------------//
    forAll(TCells, celli)
    {
        const scalar T_cell     = TCells[celli];
        const scalar p_cell     = pCells[celli];
        const scalar rho_cell   = rhoCells[celli];

		if (T_cell > minTemperature4Chemistry)
		{
            activation_map[celli] = 1.;

            //--------------- Base rates ------------------//
            if (switch_ODE == true)
            {
                //--------------- Residence Time ------------------//
                scalar tf       = 0.;

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
                for(unsigned int i=0;i<ns;i++)
                    yStar0(i) = Y[i].internalField()[celli];
                yStar0(ne-3) = T_cell;
                yStar0(ne-2) = tf;
                yStar0(ne-1) = 0.;

                fine_structure_ode_pfr.Set(p_cell);

                const double t0 = 0.;
                // Solve the ODE system
                {
                    if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
                    {
                        #include "solvePFR_DI.H"
                    }
                    #if OPENSMOKE_USE_DVODE == 1
                    else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
                    {
                        if (reactorType == REACTORTYPE_PFR)
                        {
                            SolverDVODEPFR->SetInitialValues(t0, yStar0.data());
                            SolverDVODEPFR->Solve(tf);
                            SolverDVODEPFR->Solution(yStarf.data());
                        }
                    }
                    #endif
                    #if OPENSMOKE_USE_ODEPACK == 1
                    else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
                    {
                        if (reactorType == REACTORTYPE_PFR)
                        {
                            SolverDLSODEPFR->SetInitialValues(t0, yStar0.data());
                            SolverDLSODEPFR->Solve(tf);
                            SolverDLSODEPFR->Solution(yStarf.data());
                        }
                    }
                    #endif
                    #if OPENSMOKE_USE_ODEPACK == 1
                    else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA)
                    {
                        if (reactorType == REACTORTYPE_PFR)
                        {
                            SolverDLSODAPFR->SetInitialValues(t0, yStar0.data());
                            SolverDLSODAPFR->Solve(tf);
                            SolverDLSODAPFR->Solution(yStarf.data());
                        }
                    }
                    #endif
                    #if OPENSMOKE_USE_DASPK == 1
                    else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
                    {
                        if (reactorType == REACTORTYPE_PFR)
                        {
                            SolverDASPKPFR->SetInitialValues(t0, yStar0.data());
                            SolverDASPKPFR->Solve(tf);
                            SolverDASPKPFR->Solution(yStarf.data());
                        }
                    }
                    #endif
                    #if OPENSMOKE_USE_RADAU == 1
                    else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
                    {
                        if (reactorType == REACTORTYPE_PFR)
                        {
                            SolverRADAU5PFR->SetInitialValues(t0, yStar0.data());
                            SolverRADAU5PFR->Solve(tf);
                            SolverRADAU5PFR->Solution(yStarf.data());
                        }
                    }
                    #endif
                    #if OPENSMOKE_USE_MEBDF == 1
                    else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
                    {
                        if (reactorType == REACTORTYPE_PFR)
                        {
                            SolverMEBDFPFR->SetInitialValues(t0, yStar0.data());
                            SolverMEBDFPFR->Solve(tf);
                            SolverMEBDFPFR->Solution(yStarf.data());
                        }
                    }
                    #endif
                    #if OPENSMOKE_USE_SUNDIALS == 1
                    else if (odeParametersFineStructures.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
                    {
                        if (reactorType == REACTORTYPE_PFR)
                        {
                            SolverCVODEPFR->SetInitialValues(t0, yStar0.data());
                            SolverCVODEPFR->Solve(tf);
                            SolverCVODEPFR->Solution(yStarf.data());
                        }
                    }
                    #endif
                }

                for(unsigned int i=0;i<ns;i++)
                    YStar_[i].internalFieldRef()[celli] = yStarf(i);
                TStar_[0].internalFieldRef()[celli] = yStarf(ne-3);

                for(unsigned int i=0;i<ns;i++)
                    R_ODE[i] = rho_cell*(YStar_[i].internalField()[celli]-Y[i].internalField()[celli])/tf;//[kg.m-3.s-1]

                R_base = R_ODE;
            }
            else
            {
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

                R_base = R_laminar;
            }

            double R_base_max_abs = max(abs(R_base.minCoeff()),abs(R_base.maxCoeff()));

            if (R_base_max_abs > activation_threshold)// Condition to perform modal decomposition
            {
                activation_map[celli] = 2.;

                //--------------- Jacobian ------------------//
                // Kinetic Jacobian         -> mPaSR_computeNumericalJacobian
                // Constrained Jacobian     -> mPaSR_computeNumericalJacobian_const

                Jc.setZero();
                for(unsigned int i=0;i<ns;i++)
                    vec_Y[i]      = Y[i].internalField()[celli];

                if (switch_Jconst == true)
                {
                    Jc = chemicalTimes->mPaSR_computeNumericalJacobian_const(T_cell, p_cell, vec_Y);
                }
                else
                {
                    Jc = chemicalTimes->mPaSR_computeNumericalJacobian(T_cell, p_cell, vec_Y);
                }

                if (switch_LAPACK == true)
                {
                //#if OPENSMOKE_USE_MKL == 1 
                    //--------------- LAPACK solver ------------------//
                    int info = LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'V', ns, Jc.data(), ns, lambda_real.data(), lambda_imag.data(), vl.data(), ns, vr.data(), ns);

                    if (info != 0)
                    {
                        Info << " *  OnTheFlyCEMA::Calculate: LAPACKE_dgeev failure -> status:" << info << ", using EIGEN solver instead!" << endl;
                        Eigen::EigenSolver<Eigen::MatrixXd> eigenJc(Jc);

                        for(unsigned int m=0;m<ns;m++)
                        {
                            lambda_real[m] = eigenJc.eigenvalues()[m].real();
                            lambda_imag[m] = eigenJc.eigenvalues()[m].imag();
                        }
                    }
                //#endif
                }
                else
                {
                    //---------------- EIGEN solver ------------------//
                    Eigen::EigenSolver<Eigen::MatrixXd> eigenJc(Jc);

                    for(unsigned int m=0;m<ns;m++)
                    {
                        lambda_real[m] = eigenJc.eigenvalues()[m].real();
                        lambda_imag[m] = eigenJc.eigenvalues()[m].imag();
                    }
                }

                for(unsigned int m=0;m<ns;m++)
                    lambda_mod[m] = std::sqrt(lambda_real[m] * lambda_real[m] + lambda_imag[m] * lambda_imag[m]);

                std::vector<size_t> sorted_indices(ns);
                sorted_indices = OpenSMOKE::SortAndTrackIndicesDecreasing(lambda_mod);

                for(unsigned int m=0;m<ns;m++)
                    lambda_mod_sorted[m] = lambda_mod[sorted_indices[m]];

                lambda_mod = lambda_mod_sorted;// returning the sorted version

                vr_sorted = vr;
                for(unsigned int m=0;m<ns;m++)
                {
                    for(unsigned int i=0;i<ns;i++)
                        vr_sorted(i,m) = vr(i,sorted_indices[m]);
                }
                vr = vr_sorted;// returning the sorted version

                if (vr.determinant() != 0)
                    vl = vr.inverse();
                else
                {
                    vr.setZero();
                    vl.setZero();
                }

                //--------------- Fractions ------------------//
                for(unsigned int m=0;m<ns;m++)
                {
                    tauc[m] = min(1./(lambda_mod[m] + SMALL),tauc_max);
                    kappa[m] = tauc[m] / (tauc[m] + tauMix[celli]);
                    if (switch_Ekeep == false)// Setting to 0 the modal fractions of the conserved space E
                    {
                        if (m > (ns-nelem-1))
                            kappa[m] = 0.;
                    }
                    modal_tauc_[m].internalFieldRef()[celli] = tauc[m];
                    modal_kappa_[m].internalFieldRef()[celli] = kappa[m];
                }

                //--------------- mPaSR rates ------------------//
                kappaVR     = vr * kappa.asDiagonal();  // Modes
                f_vector    = vl * R_base;              // Amplitudes
                R_mPaSR     = kappaVR * f_vector;       // mPaSR rates

                for(unsigned int i=0;i<ns;i++)
                    R_[i].internalFieldRef()[celli] = R_mPaSR[i];//[kg.m-3.s-1]
            }
            else
            {
                R_mPaSR = R_base;// laminar

                for(unsigned int i=0;i<ns;i++)
                    R_[i].internalFieldRef()[celli] = R_mPaSR[i];//[kg.m-3.s-1]
            }
		}
		else
		{
            activation_map[celli] = 0.;
			for(unsigned int i=0;i<ns;i++)
				R_[i].internalFieldRef()[celli] *= 0.;//[kg.m-3.s-1]
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
    }

    // Correct boundary conditions
    Qdot.correctBoundaryConditions();
    dQ.correctBoundaryConditions();
    Info << " * min/max(Q) = " << min(Qdot).value() << ", " << max(Qdot).value() << endl;
}