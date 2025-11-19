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
#include "fvMesh.H" //v12
#include "volFields.H" //v12
#include "surfaceFields.H" //v12
#include "fvc.H" //v12
#include "fvm.H" //v12
#include "compressibleMomentumTransportModels.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "radiationModel.H"
#include "constrainHbyA.H" //v12
#include "adjustPhi.H" //v12
#include "argList.H" //v12
#include "fluidMulticomponentThermo.H" //v12
#include "dimensionedTypes.H" //v12
#include "OFstream.H" //v12
#include "timeSelector.H" //v12
#include "functionObjectList.H" //v12
#include "combustionModel.H" //v12
#include "uniformDimensionedFields.H" //v12
#include "constrainPressure.H" //v12
#include "OFstream.H" //v12
#include "timeSelector.H" //v12
#include "fluidThermo.H" //v12
#include "rhoThermo.H" //v12
#include "zeroGradientFvPatchFields.H" //v12
#include "fvcFlux.H"
#include "fvcDdt.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"
#include "fvcVolumeIntegrate.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "hydrostaticInitialisation.H"
#include "IOMRFZoneList.H"
#include "localEulerDdtScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "linear.H"
#include "fvcDiv.H"
#include "fvcMeshPhi.H"
#include "fluidMulticomponentThermophysicalTransportModel.H"

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

using namespace Foam; //v12

int main(int argc, char *argv[])
{
    unsigned int runTimeStep = 0;

    argList::addNote("Solver for firePimpleSMOKE"); //v12
    argList args(argc, argv); //v12

	#include "createTime.H"
	#include "createMesh.H"
	#include "readGravitationalAcceleration.H"
	#include "createDyMControls.H"
	#include "initContinuityErrs.H"
	#include "createFields.H"
	#include "createOpenSMOKEFields.H"
	#include "createRadiationModel.H"

	momentumTransport->validate(); //v12

	#include "compressibleCourantNo.H"
	#include "setInitialDeltaT.H"

	#include "volFields.H"

    Info<< "\nStarting time loop\n" << endl;

    double CPUtime_comb_model = 0.;
    double CPUtime_comb_model_part = 0.;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        autoPtr<volScalarField> divrhoU;
        if (correctPhi)
        {
            divrhoU = new volScalarField
            (
                "divrhoU",
                fvc::div(fvc::absolute(phi, rho, U))
            );
        }
		#include "compressibleCourantNo.H"
		#include "setDeltaT.H"

        fvModels.preUpdateMesh();
		Foam::IOMRFZoneList MRF(mesh);
		autoPtr<surfaceVectorField> rhoUf
		(
			new surfaceVectorField
			(
				IOobject
				(
					"rhoUf",
					runTime.name(),
					mesh,
					IOobject::READ_IF_PRESENT,
					IOobject::AUTO_WRITE
				),
				fvc::interpolate(rho * U)
			)
		);

		if ((mesh.dynamic() || MRF.size()))
		{
			Info<< "Constructing face momentum rhoUf" << endl;

			// Ensure the U BCs are up-to-date before constructing Uf
			U.correctBoundaryConditions();

			rhoUf() = fvc::interpolate(rho*U);
		}

		if (!mesh.schemes().steady())
		{
			rho = thermo.rho();
			fvc::correctRhoUf(rhoUf, rho, U, phi, MRF);
		}

        // Update the mesh for topology change, mesh to mesh mapping
        mesh.update();

        runTime++;
		runTimeStep++;
        Info<< "Time = " << runTime.name() << nl << endl; //v12

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
			if (pimple.firstPimpleIter())
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
		
			if (pimple.correctTransport()) //v12
			{
				momentumTransport->correct(); //v12
			}
		}
	}

	if (SPARCswitch == true)
	{
    		#include "extensions/sparc/SPARC_local_post_processing.H"
	}
	
	rho = thermo.rho();

	runTime.write();

	Pav << runTime.name() << "\t" << p.weightedAverage(mesh.V()).value() << endl; //v12

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
