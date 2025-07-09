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

{	 
	const double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

	tmp<volScalarField> tepsilon(turbulence->epsilon());
	const volScalarField& epsilon = tepsilon();
	
	tmp<volScalarField> tk(turbulence->k());
	const volScalarField& k = tk();
	
	tmp<volScalarField> tmu(thermo.mu());
	const volScalarField& mu = tmu();
	
	tmp<volScalarField> trho(thermo.rho());
	const volScalarField& rho = trho();

	if (combustionModel == "ED" || combustionModel == "EDFR")
	{
		if (iLES == true) // LES
		{
			tauStar = Ctau * sqrt( (mu+musmall_)/rho / (epsilon+epsilonsmall_) );
			one_over_tauStar = 1./tauStar;
		}
		else // RANS
		{
			tauStar = k/(epsilon+epsilonsmall_);
			one_over_tauStar = 1./tauStar;
		}
		
		// Info<< " * min/max(epsilon)   = " << min(epsilon).value() 	<< ", " << max(epsilon).value() << endl;
		// Info<< " * min/max(k)         = " << min(k).value()		<< ", " << max(k).value() 	<< endl;
		Info<< " * min/max(tauStar)   = " << min(tauStar).value() 	<< ", " << max(tauStar).value() << endl;
	}

	if (combustionModel == "EDC")
	{
		if (edc_version == "v2016")
		{
			// Calculation of characteristic times
			{
				CharacteristicChemicalTimes* chemicalTimes = new CharacteristicChemicalTimes(*thermodynamicsMapXML, *kineticsMapXML);
				chemicalTimes->SetThresholdChemicalTime(threshold_chemical_time);
				chemicalTimes->SetType(characteristic_chemical_times_type);
				chemicalTimes->SetThresholdTemperatureChemicalTime(threshold_temperature_chemical_time);

				const unsigned int ns = thermodynamicsMapXML->NumberOfSpecies();
				OpenSMOKE::OpenSMOKEVectorDouble vec_Y_for_tauc(ns);

				const scalarField& TCells = T.internalField();
				const scalarField& pCells = p.internalField(); 
				
				scalarField& tauChemCells = tauChem.ref(); 
				scalarField& DaCells = Da.ref(); 
				scalarField& ReTCells = ReT.ref();
				
				forAll(tauStar, i)
				{
					const double Ti = TCells[i];
					const double pi = pCells[i];

					for(unsigned int j=0;j<ns;j++)
						vec_Y_for_tauc[j+1] = Y[j].internalField()[i];

					// Chemical time
					const double tauLocal = chemicalTimes->CalculateCharacteristicChemicalTime(Ti,pi,vec_Y_for_tauc);
					
					// Coefficients
					const scalar nu = (mu[i]+SMALL)/(rho[i] + SMALL);
					const scalar DaLocal = max(min(Foam::sqrt(nu/(epsilon[i] + SMALL))/tauLocal, maximum_Da), minimum_Da);
					const scalar ReTLocal = CRe_coefficient*sqr(k[i])/(nu*epsilon[i] + SMALL);

					const scalar C1 = Foam::sqrt(Ctau_coefficient/3.);
					const scalar CtauLocal = max( min(C1/(DaLocal*Foam::sqrt(ReTLocal + 1)), 5), 0.40825);
					const scalar CgammaLocal = max( min( Cgamma_coefficient *Foam::sqrt(DaLocal*(ReTLocal+1.)), 2.1377), 0.5);
					
					// Fill fields
					tauChemCells[i] = tauLocal;
					DaCells[i] = DaLocal;
					ReTCells[i] = ReTLocal;
				
					// Fine structure properties
					tauStar[i] = CtauLocal * Foam::sqrt( nu / (epsilon[i]+SMALL) );
					gammaL[i]  = CgammaLocal * pow025(nu*epsilon[i]/(sqr(k[i]) + SMALL));
				}

				// Correct boundary conditions
				tauChem.correctBoundaryConditions();
				Da.correctBoundaryConditions();
				ReT.correctBoundaryConditions();
			}

			dimensionedScalar tauStarSMALL("tauStarSMALL", dimensionSet(0,0,1,0,0,0,0), scalar(SMALL));
			// Fine structures: characteristic time [s]
			mDotStar = 1./(tauStar+tauStarSMALL);

			// Fine structure properties: volumetric fraction
			gammaL = min(maxGammaL, gammaL);
			gammaStarNum = pow(gammaL,edc_exp_1);
			gammaStarDen = pow(gammaL,edc_exp_2);
		}
		else
		{
			// Fine structures: characteristic time [s]
			tauStar = Ctau * sqrt( (mu+musmall_)/(rho+rhosmall_) / (epsilon+epsilonsmall_) );
			mDotStar = 1./tauStar;

			// Fine structure properties: volumetric fraction
			gammaL = Ccsi * pow025( (mu+musmall_)/(rho+rhosmall_) * epsilon / sqr(k+ksmall_));	
			gammaL = min(maxGammaL, gammaL);
			gammaStarNum = pow(gammaL,edc_exp_1);
			gammaStarDen = pow(gammaL,edc_exp_2);
		}

		// EDC coefficient
		edc_coefficient = rho * mDotStar * gammaStarNum / (1.-gammaStarDen);
		tfStar = (1.-gammaStarDen)/(mDotStar);		

		// Info
		Info<< " * min/max(epsilon)   = " << min(epsilon).value() 	<< ", " << max(epsilon).value() 	<< endl;
		Info<< " * min/max(k)         = " << min(k).value()		<< ", " << max(k).value()		<< endl;
		Info<< " * min/max(gammaL)    = " << min(gammaL).value()	<< ", " << max(gammaL).value()		<< endl;
		Info<< " * min/max(gammaStar) = " << min(gammaStarDen)		<< ", " << max(gammaStarDen) 		<< endl;
		Info<< " * min/max(tauStar)   = " << min(tauStar).value()	<< ", " << max(tauStar).value() 	<< endl;
		Info<< " * min/max(tfStar)    = " << min(tfStar)	 	<< ", " << max(tfStar)	 		<< endl;

		if (edc_version == "v2016")
		{
			Info<< " * min/max(tc)   = " << min(tauChem).value() 	<< ", " << max(tauChem).value()	<< endl;
			Info<< " * min/max(Da)   = " << min(Da).value()		<< ", " << max(Da).value()	<< endl;
			Info<< " * min/max(ReT)  = " << min(ReT).value()	<< ", " << max(ReT).value()	<< endl;
		}
	}

	if (combustionModel == "PaSR")
	{
		#if STEADYSTATE == 0
		if (iLES == true)
		{
			if ( tauMixType == SCALE_MEAN)
			{
				const Foam::compressible::LESModel& lesModel = U.db().lookupObject<Foam::compressible::LESModel>("turbulenceProperties");
				epsilonSGS = 0.54433*pow((k+ksmall_), 1.5)/lesModel.delta();
				Info << " * using geometrical mean mixing time scale..." << endl;
				tauMix =  sqrt( (1.2247*lesModel.delta()/sqrt(k+ksmall_))*sqrt(mu/rho/epsilonSGS) ); //geometrical mean of Kolmogorov scale and sub-grid velocity stretch time (filter size)
			}
			else if ( tauMixType == SCALE_KOLMOGOROV)
			{
				const Foam::compressible::LESModel& lesModel = U.db().lookupObject<Foam::compressible::LESModel>("turbulenceProperties");
				epsilonSGS = 0.54433*pow((k+ksmall_), 1.5)/lesModel.delta();
				Info << " * using Kolmogorov mixing time scale..." << endl;
				tauMix = sqrt(mu/rho/epsilonSGS);
			}
			else if ( tauMixType == SCALE_STRETCH)
			{
				const Foam::compressible::LESModel& lesModel = U.db().lookupObject<Foam::compressible::LESModel>("turbulenceProperties");
				Info << " * using sub-grid velocityr stretch time scale..." << endl;
				tauMix = 1.2247*lesModel.delta()/sqrt(k+ksmall_);
			}
			else if ( tauMixType == SCALE_DYNAMIC )
			{	
				const Foam::compressible::LESModel& lesModel = U.db().lookupObject<Foam::compressible::LESModel>("turbulenceProperties");
				epsilonSGS = 0.54433*pow((k+ksmall_), 1.5)/lesModel.delta();
				Info << " * using dynamic mixing time scale..." << endl;
				tauMix = min(max (varf/(Chi+chismall_), tauMixsmall_), ((k+ksmall_)/(epsilonSGS+epsilonsmall_)));  // dynamic scale
			}
			else if ( tauMixType == SCALE_DT_LES )
			{
				Info << " * using deltaT as mixing time scale..." << endl;
				tauMix = runTime.deltaT();  // deltaT scale
			}
			//Info<< " * min/max(epsilon(ssg)) = " << min(epsilonSGS).value() 	<< ", " << max(epsilonSGS).value() 		<< endl;
			//Info<< " * min/max(epsilon(les)) = " << min(lesModel.epsilon()).value() << ", " << max(lesModel.epsilon()).value() 	<< endl;
			//Info<< " * min/max(mut)          = " << min(turbulence->mut()).value() 	<< ", " << max(turbulence->mut()).value() 	<< endl;
			Info<< " * min/max(taumix)       = " << min(tauMix).value() 		<< ", " << max(tauMix).value()	 		<< endl;
		}
		else
		#endif	// RANS
		{
			if ( tauMixType == SCALE_KOLMOGOROV)
			{
				Info << " * using Kolmogorov mixing time scale..." << endl;
				tauMix =  sqrt( (mu+musmall_)/rho/(epsilon+epsilonsmall_) );	// Kolmogorov mixing scale
				Cmix = tauMix/((k+ksmall_)/(epsilon+epsilonsmall_));  //global scale
			}
			else if ( tauMixType == SCALE_MEAN)
			{
				Info << " * using geometrical mean mixing time scale..." << endl;
				tauMix =  sqrt( ((k+ksmall_)/(epsilon+epsilonsmall_))*sqrt((mu+musmall_)/rho/(epsilon+epsilonsmall_)) ); //geometrical mean of integral scale and Kolmogorov scale
				Cmix = tauMix/((k+ksmall_)/(epsilon+epsilonsmall_));  //global scale
			}
			else if ( tauMixType == SCALE_GLOBAL)
			{
				Info << " * using global mixing time scale..." << endl;
				tauMix = Cmixglobal*(k+ksmall_)/(epsilon+epsilonsmall_);  //global scale
			}
			
			else if ( tauMixType == SCALE_DYNAMIC )
			{
				Info << " * using dynamic mixing time scale..." << endl;
				tauMix = min(max (varf/(Chi+chismall_), tauMixsmall_), ((k+ksmall_)/(epsilon+epsilonsmall_)));  // dynamic scale
				Cmix = tauMix/((k+ksmall_)/(epsilon+epsilonsmall_));  // equivalent Cmix
			}

			// Info<< " * min/max(epsilon)  = " << min(epsilon).value() 		<< ", " << max(epsilon).value() << endl;
			// Info<< " * min/max(mut)      = " << min(turbulence->mut()).value()	<< ", " << max(turbulence->mut()).value() << endl;
			Info<< " * min/max(taumix)   = " << min(tauMix).value() 		<< ", " << max(tauMix).value() << endl;
		}

		// Calculation of characteristic times
		{
			CharacteristicChemicalTimes* chemicalTimes = new CharacteristicChemicalTimes(*thermodynamicsMapXML, *kineticsMapXML);
			chemicalTimes->SetThresholdChemicalTime(threshold_chemical_time);
			chemicalTimes->SetType(characteristic_chemical_times_type);
			chemicalTimes->SetThresholdTemperatureChemicalTime(threshold_temperature_chemical_time);

			const unsigned int ns = thermodynamicsMapXML->NumberOfSpecies();
			OpenSMOKE::OpenSMOKEVectorDouble vec_Y_for_tauc(ns);

			const scalarField& TCells = T.internalField();
			const scalarField& pCells = p.internalField(); 
			const scalarField& tauMixCells = tauMix.internalField(); 
			
			scalarField& tauChemCells = tauChem.ref(); 
			scalarField& DaCells = Da.ref(); 
			scalarField& kappaCells = kappa.ref(); 
			
			forAll(tauStar, i)
			{
				const double Ti = TCells[i];
				const double pi = pCells[i];

				for(unsigned int j=0;j<ns;j++)
					vec_Y_for_tauc[j+1] = Y[j].internalField()[i];

				// Chemical time
				const double tauLocal = chemicalTimes->CalculateCharacteristicChemicalTime(Ti,pi,vec_Y_for_tauc);

				// Mixing time				
				const scalar tauMixLocal = tauMixCells[i];

				// PASR coefficients
				const scalar kappaLocal = tauLocal/(tauLocal+tauMixLocal);
				const scalar DaLocal = max( min(tauMixLocal/tauLocal, maximum_Da), minimum_Da);

				// Fill fields
				tauChemCells[i] = tauLocal;
				DaCells[i] = DaLocal;
				kappaCells[i] = kappaLocal;
			}

			// Correct boundary conditions
			tauChem.correctBoundaryConditions();
			Da.correctBoundaryConditions();
			kappa.correctBoundaryConditions();

			Info<< " * min/max(tc)    = " << min(tauChem).value() 	<< ", " << max(tauChem).value() 	<< endl;
			Info<< " * min/max(Da)    = " << min(Da).value()	<< ", " << max(Da).value()		<< endl;
			Info<< " * min/max(kappa) = " << min(kappa).value()	<< ", " << max(kappa).value()		<< endl;
		}
	}

	if (combustionModel == "QLFR")
	{
		#if STEADYSTATE == 0
		if (iLES == true)
		{
			if ( tauMixType == SCALE_MEAN)
			{
				const Foam::compressible::LESModel& lesModel = U.db().lookupObject<Foam::compressible::LESModel>("turbulenceProperties");
				epsilonSGS = 0.54433*pow((k+ksmall_), 1.5)/lesModel.delta();
				Info << " * using geometrical mean mixing time scale..." << endl;
				tauMix =  sqrt( (1.2247*lesModel.delta()/sqrt(k+ksmall_))*sqrt(mu/rho/epsilonSGS) ); //geometrical mean of Kolmogorov scale and sub-grid velocity stretch time (filter size)
			}
			else if ( tauMixType == SCALE_KOLMOGOROV)
			{
				const Foam::compressible::LESModel& lesModel = U.db().lookupObject<Foam::compressible::LESModel>("turbulenceProperties");
				epsilonSGS = 0.54433*pow((k+ksmall_), 1.5)/lesModel.delta();
				Info << " * using Kolmogorov mixing time scale..." << endl;
				tauMix = sqrt(mu/rho/epsilonSGS);
			}
			else if ( tauMixType == SCALE_STRETCH)
			{
				const Foam::compressible::LESModel& lesModel = U.db().lookupObject<Foam::compressible::LESModel>("turbulenceProperties");
				Info << " * using sub-grid velocityr stretch time scale..." << endl;
				tauMix = 1.2247*lesModel.delta()/sqrt(k+ksmall_);
			}
			else if ( tauMixType == SCALE_DYNAMIC )
			{	
				const Foam::compressible::LESModel& lesModel = U.db().lookupObject<Foam::compressible::LESModel>("turbulenceProperties");
				epsilonSGS = 0.54433*pow((k+ksmall_), 1.5)/lesModel.delta();
				Info << " * using dynamic mixing time scale..." << endl;
				tauMix = min(max (varf/(Chi+chismall_), tauMixsmall_), ((k+ksmall_)/(epsilonSGS+epsilonsmall_)));  // dynamic scale
			}
			else if ( tauMixType == SCALE_DT_LES )
			{
				Info << " * using deltaT as mixing time scale..." << endl;
				tauMix = runTime.deltaT();  // deltaT scale
			}

			//Info<< " * min/max(epsilon(ssg)) = " << min(epsilonSGS).value() 	    << ", " << max(epsilonSGS).value() 		<< endl;
			//Info<< " * min/max(epsilon(les)) = " << min(lesModel.epsilon()).value() << ", " << max(lesModel.epsilon()).value() 	<< endl;
			//Info<< " * min/max(mut)          = " << min(turbulence->mut()).value() 	<< ", " << max(turbulence->mut()).value() 	<< endl;
			Info<< " * min/max(taumix)       = " << min(tauMix).value() 		<< ", " << max(tauMix).value()	 		<< endl;
		}
		else
		#endif	// RANS
		{
			Info << "QLFR model to be used only for LES" << endl;
			abort();
		}

		// Calculation of characteristic times
		{
			const unsigned int ns = thermodynamicsMapXML->NumberOfSpecies();
			OpenSMOKE::OpenSMOKEVectorDouble yStar(ns);

			const scalarField& TCells = T.internalField();
			const scalarField& pCells = p.internalField(); 
			
			scalarField& kappaCells = kappa.ref(); 

			forAll(tauStar, i)
			{
				const double Ti = TCells[i];
				const double pi = pCells[i];

				for(unsigned int j=0;j<ns;j++)
					yStar[j+1] = YStar_[j].internalField()[i];

				const scalar kappaLocal = 1;
			}

			// Correct boundary conditions
			kappa.correctBoundaryConditions();
		}
	}

	if (combustionModel == "PaSR_LFR")
	{
		#if STEADYSTATE == 0
		if (iLES == true)
		{
			if ( tauMixType == SCALE_MEAN)
			{
				const Foam::compressible::LESModel& lesModel = U.db().lookupObject<Foam::compressible::LESModel>("turbulenceProperties");
				epsilonSGS = 0.54433*pow((k+ksmall_), 1.5)/lesModel.delta();
				Info << " * using geometrical mean mixing time scale..." << endl;
				tauMix =  sqrt( (1.2247*lesModel.delta()/sqrt(k+ksmall_))*sqrt(mu/rho/epsilonSGS) ); //geometrical mean of Kolmogorov scale and sub-grid velocity stretch time (filter size)
			}
			else if ( tauMixType == SCALE_KOLMOGOROV)
			{
				const Foam::compressible::LESModel& lesModel = U.db().lookupObject<Foam::compressible::LESModel>("turbulenceProperties");
				epsilonSGS = 0.54433*pow((k+ksmall_), 1.5)/lesModel.delta();
				Info << " * using Kolmogorov mixing time scale..." << endl;
				tauMix = sqrt(mu/rho/epsilonSGS);
			}
			else if ( tauMixType == SCALE_STRETCH)
			{
				const Foam::compressible::LESModel& lesModel = U.db().lookupObject<Foam::compressible::LESModel>("turbulenceProperties");
				Info << " * using sub-grid velocityr stretch time scale..." << endl;
				tauMix = 1.2247*lesModel.delta()/sqrt(k+ksmall_);
			}
			else if ( tauMixType == SCALE_DYNAMIC )
			{	
				const Foam::compressible::LESModel& lesModel = U.db().lookupObject<Foam::compressible::LESModel>("turbulenceProperties");
				epsilonSGS = 0.54433*pow((k+ksmall_), 1.5)/lesModel.delta();
				Info << " * using dynamic mixing time scale..." << endl;
				tauMix = min(max (varf/(Chi+chismall_), tauMixsmall_), ((k+ksmall_)/(epsilonSGS+epsilonsmall_)));  // dynamic scale
			}
			else if ( tauMixType == SCALE_DT_LES )
			{
				Info << " * using deltaT as mixing time scale..." << endl;
				tauMix = runTime.deltaT();  // deltaT scale
			}
			//Info<< " * min/max(epsilon(ssg)) = " << min(epsilonSGS).value() 	<< ", " << max(epsilonSGS).value() 		<< endl;
			//Info<< " * min/max(epsilon(les)) = " << min(lesModel.epsilon()).value() << ", " << max(lesModel.epsilon()).value() 	<< endl;
			//Info<< " * min/max(mut)          = " << min(turbulence->mut()).value() 	<< ", " << max(turbulence->mut()).value() 	<< endl;
			Info<< " * min/max(taumix)       = " << min(tauMix).value() 		<< ", " << max(tauMix).value()	 		<< endl;
		}
		else
		#endif	// RANS
		{
			if ( tauMixType == SCALE_KOLMOGOROV)
			{
				Info << " * using Kolmogorov mixing time scale..." << endl;
				tauMix =  sqrt( (mu+musmall_)/rho/(epsilon+epsilonsmall_) );	// Kolmogorov mixing scale
				Cmix = tauMix/((k+ksmall_)/(epsilon+epsilonsmall_));  //global scale
			}
			else if ( tauMixType == SCALE_MEAN)
			{
				Info << " * using geometrical mean mixing time scale..." << endl;
				tauMix =  sqrt( ((k+ksmall_)/(epsilon+epsilonsmall_))*sqrt((mu+musmall_)/rho/(epsilon+epsilonsmall_)) ); //geometrical mean of integral scale and Kolmogorov scale
				Cmix = tauMix/((k+ksmall_)/(epsilon+epsilonsmall_));  //global scale
			}
			else if ( tauMixType == SCALE_GLOBAL)
			{
				Info << " * using global mixing time scale..." << endl;
				tauMix = Cmixglobal*(k+ksmall_)/(epsilon+epsilonsmall_);  //global scale
			}
			
			else if ( tauMixType == SCALE_DYNAMIC )
			{
				Info << " * using dynamic mixing time scale..." << endl;
				tauMix = min(max (varf/(Chi+chismall_), tauMixsmall_), ((k+ksmall_)/(epsilon+epsilonsmall_)));  // dynamic scale
				Cmix = tauMix/((k+ksmall_)/(epsilon+epsilonsmall_));  // equivalent Cmix
			}

			// Info<< " * min/max(epsilon)  = " << min(epsilon).value() 		<< ", " << max(epsilon).value() << endl;
			// Info<< " * min/max(mut)      = " << min(turbulence->mut()).value()	<< ", " << max(turbulence->mut()).value() << endl;
			Info<< " * min/max(taumix)   = " << min(tauMix).value() 		<< ", " << max(tauMix).value() << endl;
		}

		// Calculation of characteristic times
		{
			CharacteristicChemicalTimes* chemicalTimes = new CharacteristicChemicalTimes(*thermodynamicsMapXML, *kineticsMapXML);
			chemicalTimes->SetThresholdChemicalTime(threshold_chemical_time);
			chemicalTimes->SetType(characteristic_chemical_times_type);
			chemicalTimes->SetThresholdTemperatureChemicalTime(threshold_temperature_chemical_time);

			const unsigned int ns = thermodynamicsMapXML->NumberOfSpecies();
			OpenSMOKE::OpenSMOKEVectorDouble vec_Y_for_tauc(ns);

			const scalarField& TCells = T.internalField();
			const scalarField& pCells = p.internalField(); 
			const scalarField& tauMixCells = tauMix.internalField(); 
			
			scalarField& tauChemCells = tauChem.ref(); 
			scalarField& DaCells = Da.ref(); 
			scalarField& kappaCells = kappa.ref(); 
			
			forAll(tauStar, i)
			{
				const double Ti = TCells[i];
				const double pi = pCells[i];

				for(unsigned int j=0;j<ns;j++)
					vec_Y_for_tauc[j+1] = Y[j].internalField()[i];

				// Chemical time
				const double tauLocal = chemicalTimes->CalculateCharacteristicChemicalTime(Ti,pi,vec_Y_for_tauc);

				// Mixing time				
				const scalar tauMixLocal = tauMixCells[i];

				// PaSR_LFR coefficients
				const scalar kappaLocal = tauLocal/(tauLocal+tauMixLocal);
				const scalar DaLocal = max( min(tauMixLocal/tauLocal, maximum_Da), minimum_Da);

				// Fill fields
				tauChemCells[i] = tauLocal;
				DaCells[i] = DaLocal;
				kappaCells[i] = kappaLocal;
			}

			// Correct boundary conditions
			tauChem.correctBoundaryConditions();
			Da.correctBoundaryConditions();
			kappa.correctBoundaryConditions();

			Info<< " * min/max(tc)    = " << min(tauChem).value() 	<< ", " << max(tauChem).value() 	<< endl;
			Info<< " * min/max(Da)    = " << min(Da).value()	<< ", " << max(Da).value()		<< endl;
			Info<< " * min/max(kappa) = " << min(kappa).value()	<< ", " << max(kappa).value()		<< endl;
		}
	}

	if (combustionModel == "mPaSR")
	{
		#if STEADYSTATE == 0
		if (iLES == true)
		{
			if ( tauMixType == SCALE_MEAN)
			{
				const Foam::compressible::LESModel& lesModel = U.db().lookupObject<Foam::compressible::LESModel>("turbulenceProperties");
				epsilonSGS = 0.54433*pow((k+ksmall_), 1.5)/lesModel.delta();
				Info << " * using geometrical mean mixing time scale..." << endl;
				tauMix =  sqrt( (1.2247*lesModel.delta()/sqrt(k+ksmall_))*sqrt(mu/rho/epsilonSGS) ); //geometrical mean of Kolmogorov scale and sub-grid velocity stretch time (filter size)
			}
			else if ( tauMixType == SCALE_KOLMOGOROV)
			{
				const Foam::compressible::LESModel& lesModel = U.db().lookupObject<Foam::compressible::LESModel>("turbulenceProperties");
				epsilonSGS = 0.54433*pow((k+ksmall_), 1.5)/lesModel.delta();
				Info << " * using Kolmogorov mixing time scale..." << endl;
				tauMix = sqrt(mu/rho/epsilonSGS);
			}
			else if ( tauMixType == SCALE_STRETCH)
			{
				const Foam::compressible::LESModel& lesModel = U.db().lookupObject<Foam::compressible::LESModel>("turbulenceProperties");
				Info << " * using sub-grid velocityr stretch time scale..." << endl;
				tauMix = 1.2247*lesModel.delta()/sqrt(k+ksmall_);
			}
			else if ( tauMixType == SCALE_DYNAMIC )
			{	
				const Foam::compressible::LESModel& lesModel = U.db().lookupObject<Foam::compressible::LESModel>("turbulenceProperties");
				epsilonSGS = 0.54433*pow((k+ksmall_), 1.5)/lesModel.delta();
				Info << " * using dynamic mixing time scale..." << endl;
				tauMix = min(max (varf/(Chi+chismall_), tauMixsmall_), ((k+ksmall_)/(epsilonSGS+epsilonsmall_)));  // dynamic scale
			}
			else if ( tauMixType == SCALE_DT_LES )
			{
				Info << " * using deltaT as mixing time scale..." << endl;
				tauMix = runTime.deltaT();  // deltaT scale
			}

			//Info<< " * min/max(epsilon(ssg)) = " << min(epsilonSGS).value() 	<< ", " << max(epsilonSGS).value() 		<< endl;
			//Info<< " * min/max(epsilon(les)) = " << min(lesModel.epsilon()).value() << ", " << max(lesModel.epsilon()).value() 	<< endl;
			//Info<< " * min/max(mut)          = " << min(turbulence->mut()).value() 	<< ", " << max(turbulence->mut()).value() 	<< endl;
			Info<< " * min/max(taumix)       = " << min(tauMix).value() 		<< ", " << max(tauMix).value()	 		<< endl;
		}
		else
		#endif	// RANS
		{
			if ( tauMixType == SCALE_KOLMOGOROV)
			{
				Info << " * using Kolmogorov mixing time scale..." << endl;
				tauMix =  sqrt( (mu+musmall_)/rho/(epsilon+epsilonsmall_) );	// Kolmogorov mixing scale
				Cmix = tauMix/((k+ksmall_)/(epsilon+epsilonsmall_));  //global scale
			}
			else if ( tauMixType == SCALE_MEAN)
			{
				Info << " * using geometrical mean mixing time scale..." << endl;
				tauMix =  sqrt( ((k+ksmall_)/(epsilon+epsilonsmall_))*sqrt((mu+musmall_)/rho/(epsilon+epsilonsmall_)) ); //geometrical mean of integral scale and Kolmogorov scale
				Cmix = tauMix/((k+ksmall_)/(epsilon+epsilonsmall_));  //global scale
			}
			else if ( tauMixType == SCALE_GLOBAL)
			{
				Info << " * using global mixing time scale..." << endl;
				tauMix = Cmixglobal*(k+ksmall_)/(epsilon+epsilonsmall_);  //global scale
			}
			
			else if ( tauMixType == SCALE_DYNAMIC)
			{
				Info << " * using dynamic mixing time scale..." << endl;
				tauMix = min(max (varf/(Chi+chismall_), tauMixsmall_), ((k+ksmall_)/(epsilon+epsilonsmall_)));  // dynamic scale
				Cmix = tauMix/((k+ksmall_)/(epsilon+epsilonsmall_));  // equivalent Cmix
			}

			// Info<< " * min/max(epsilon)  = " << min(epsilon).value() 		<< ", " << max(epsilon).value() << endl;
			// Info<< " * min/max(mut)      = " << min(turbulence->mut()).value()	<< ", " << max(turbulence->mut()).value() << endl;
			Info<< " * min/max(taumix)   = " << min(tauMix).value() 		<< ", " << max(tauMix).value() << endl;
		}
	}

	const double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
	Info<< " * Time in propFineStructures [s] = " << tEnd - tStart << endl;
}
    

