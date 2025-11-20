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

// This is a unsteady simulation
#define STEADYSTATE 0

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// OpenSMOKE++ Dictionaries
#include "dictionary/OpenSMOKE_Dictionary"

// ODE solvers
#include "math/native-ode-solvers/MultiValueSolver"
#include "math/external-ode-solvers/ODE_Parameters.h"

// NLS solvers
#include "math/native-nls-solvers/NonLinearSystemSolver"
#include "math/native-nls-solvers/parameters/NonLinearSolver_Parameters.h"

// OpenFOAM
#include "fvCFD.H"
#include "fluidReactionThermo.H"	
#include "combustionModel.H"
#include "compressibleMomentumTransportModels.H"
#include "fluidReactionThermophysicalTransportModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "CorrectPhi.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "radiationModel.H"

// Utilities
#include "Utilities.H"

// DRG
#include "DRG.H"

// ODE systems
#include "ODE_PSR.H"
#include "ODE_PSR_Interface.H"
#include "ODE_PFR.H"
#include "ODE_PFR_Interface.H"
#include "BatchReactorHomogeneousConstantPressure.H"
#include "BatchReactorHomogeneousConstantPressure_ODE_Interface.H"

// NLS Systems
#include "NLS_PSR.H"
#include "NLS_PSR_Interface.H"

// Characteristic chemical times
#include "CharacteristicChemicalTimes.H"

// ISAT
#if FIRESMOKE_USE_ISAT == 1
    #include "ISAT.h"
    #include "numericalJacobian4ISAT.H"
    #include "mappingGradients/mappingGradient4OpenFOAM.h"
#endif

// SPARC
#if SPARC==1
	#include "myNeuralNetworkBasedOnPCA.h"
	#include "classifyPoint.h"
	#include "classifyPoint_initialize.h"
	#include "classifyPoint_terminate.h"
	#include "extensions/sparc/SPARC_classifier_VQ2.H"
	#include "extensions/sparc/SPARC_classifier_SOFTMAX.H"
	#include "extensions/sparc/ODE_PFR_SPARC.H"
	#include "extensions/sparc/ODE_PFR_SPARC_Interface.H"
	#include "extensions/sparc/SPARC_classifier_NEURAL.H"
	#include "extensions/sparc/SPARC_predictor_NEURAL.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    unsigned int runTimeStep = 0;

	#include "postProcess.H"

	#include "setRootCaseLists.H"
	#include "createTime.H"
	#include "createMesh.H"
	#include "readGravitationalAcceleration.H"
	#include "createDyMControls.H"
	#include "initContinuityErrs.H"
	#include "createFields.H"
	#include "createOpenSMOKEFields.H"
	#include "createRhoUfIfPresent.H"
	#include "createRadiationModel.H"

	turbulence->validate();

	if (!LTS)
	{
		#include "compressibleCourantNo.H"
		#include "setInitialDeltaT.H"
	}

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    double CPUtime_comb_model = 0.;
    double CPUtime_comb_model_part = 0.;

    while (runTime.run())
    {
        #include "readDyMControls.H"

	// Store divrhoU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence
        autoPtr<volScalarField> divrhoU;
        if (correctPhi)
        {
            divrhoU = new volScalarField
            (
                "divrhoU",
                fvc::div(fvc::absolute(phi, rho, U))
            );
        }

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        fvModels.preUpdateMesh();

        // Store momentum to set rhoUf for introduced faces.
        autoPtr<volVectorField> rhoU;
        if (rhoUf.valid())
        {
            rhoU = new volVectorField("rhoU", rho*U);
        }

        // Update the mesh for topology change, mesh to mesh mapping
        mesh.update();

        runTime++;
		runTimeStep++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

	// --- Pressure-velocity PIMPLE corrector loop
	while (pimple.loop())
	{
		if (!pimple.flow())
		{
			if (pimple.models())
			{
				fvModels.correct();
			}

			if (pimple.thermophysics())
			{
				#include "properties.H"
				const double t0 = runTime.value() - runTime.deltaT().value();
				const double tf = runTime.value();
                		#include "YEqn.H"
                		#include "EEqn.H"
						#include "fEqn.H"
						#include "varfEqn.H"
			}
		}
		else
		{
			if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
                	{
                    		// Move the mesh
                    		mesh.move();

                   		if (mesh.changing())
                    		{
                        		MRF.update();

                        		if (correctPhi)
                        		{
                            		#include "correctPhi.H"
                       			}

                        		if (checkMeshCourantNo)
                        		{
                            		#include "meshCourantNo.H"
                        		}
                    		}
                	}

                	if (pimple.firstPimpleIter() && !pimple.simpleRho())
                	{
                    		#include "rhoEqn.H"
                	}

                	if (pimple.models())
                	{
                    		fvModels.correct();
                	}

                	#include "UEqn.H"

                	if (pimple.thermophysics())
			{
				#include "properties.H"
				const double t0 = runTime.value() - runTime.deltaT().value();
				const double tf = runTime.value();
				#include "YEqn.H"
                #include "EEqn.H"
			}

		 	if (dynamicCmixEquations == true)
		 	{	
				#include "fEqn.H"
		    	#include "varfEqn.H"
		    	#include "ChiEqn.H"
		 	}

			// --- Pressure corrector loop
			while (pimple.correct())
			{
				#include "pEqn.H"
			}
		
			if (pimple.turbCorr())
			{
				turbulence->correct();
			}
		}
	}

	if (SPARCswitch == true)
	{
    		#include "extensions/sparc/SPARC_local_post_processing.H"
	}
	
	rho = thermo.rho();

	runTime.write();

	Pav << runTime.timeName() << "\t" << p.weightedAverage(mesh.V()).value() << endl;

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << " ClockTime = " << runTime.elapsedClockTime() << " s"
		<< "TotalTime combustion model = " << CPUtime_comb_model << " s,  "
	    << "TotalTime combustion model part = " << CPUtime_comb_model_part << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
