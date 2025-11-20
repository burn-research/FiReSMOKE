FiReSMOKE with SPARC and mPaSR model
========

Finite-rate combustion solvers (ED, ED/FR, EDC, QLFR, PaSR, mPaSR) for OpenFOAM based on the [OpenSMOKE++ framework][1].
Each branch of fireSMOKE requires one of the following OpenFOAM versions, respectively:
- OpenFOAM-12.x (main), OpenFOAM-10.x (of10), OpenFOAM-7.x (of7)

Docker installation provides a self-contained environment for OpenFOAM-12, OpenFOAM-10, OpenFOAM-7 with a pre-compiled version of FiReSMOKE.

If you use fireSMOKE for your publications, we kindly ask you to cite the following two papers:

> Parente, A., Malik, R.M., Contino, F., Cuoci, A., Dally, B., 
> Extension of the Eddy Dissipation Concept for turbulence/chemistry interactions to MILD combustion
> (2016) Fuel, 163, pp. 98-111, DOI: 10.1016/j.fuel.2015.09.020

> Cuoci, A., Frassoldati, A., Faravelli, T., Ranzi, E., 
> OpenSMOKE++: An object-oriented framework for the numerical modeling of reactive systems with detailed kinetic mechanisms 
> (2015) Computer Physics Communications, 192, pp. 237-264, DOI: 10.1016/j.cpc.2015.02.014

If you use the SPARC plugin for your publications, we kindly ask you to cite the following two papers:

> D'Alessio, G., Parente, A., Stagni, A., Cuoci, A., 
> Adaptive chemistry via pre-partitioning of composition space and mechanism reduction. 
> (2020) Combustion and Flame, Volume 211, pp. 68-82. DOI: 10.1016/j.combustflame.2019.09.010.

> Amaduzzi, R., D'Alessio, G., Pagani, P., Cuoci, A., Malpica Galassi, R., Parente, A., 
> Automated adaptive chemistry for Large Eddy Simulations of turbulent reacting flows. 
> (2024) Combustion and Flame, Volume 259, pp. 113-136. DOI: 10.1016/j.combustflame.2023.113136.

If you use the mPaSR model for your publications, we kindly ask you to cite the following two papers:

> Quadarella, E., Péquin, A., Stagni, A., Parente, A., Faravelli, T., Im, H. G.,
> A generalized partially stirred reactor model for turbulent closure.
> (2023) Proceedings of the Combustion Institute, Volume 39, Issue 4, pp 5329-5338. DOI: 10.1016/j.proci.2022.08.061

> Péquin, A., Quadarella, E., Malpica Galassi, R., Iavarone, S., Im, H. G., Parente, A.,
> A modal decomposition-based partially stirred reactor (mPaSR) model for turbulent combustion closure: Implementation details and a posteriori validation.
> (2025) Combustion and Flame, Volume 279, 2025, 114269. DOI: 10.1016/j.combustflame.2025.114269.

![](https://github.com/apequin/FiReSMOKE_mPaSR/blob/main/runs/1_SandiaFlameD_RANS/flameD_ignition_mPaSR.gif)

Compulsory libraries
--------------------
- OpenSMOKE++ (already included in fireSMOKE)
- Eigen (http://eigen.tuxfamily.org/index.php?title=Main_Page)
- Boost C++ (http://www.boost.org/)

Optional libraries
------------------
- Intel MKL (https://software.intel.com/en-us/intel-mkl)
- ODEPACK (http://computation.llnl.gov/casc/odepack/odepack_home.html)
- DVODE (http://computation.llnl.gov/casc/odepack/odepack_home.html)
- DASPK (http://www.engineering.ucsb.edu/~cse/software.html)
- Sundials (http://computation.llnl.gov/casc/sundials/main.html)
- MEBDF (http://wwwf.imperial.ac.uk/~jcash/IVP_software/readme.html)
- RADAU (http://www.unige.ch/~hairer/software.html)
- Armadillo (for SPARC compilation, https://arma.sourceforge.net/) 

Compilation
-----------
Three different options are available to compile the code, according to the level of support for the solution of ODE systems
1. Minimalist: no external, optional libraries are required. Only the native OpenSMOKE++ ODE solver can be used.
2. Minimalist + Intel MKL: only the native OpenSMOKE++ ODE solver can be used, but linear algebra operations are managed by the Intel MKL libraries
3. Complete: all the optional libraries are linked to the code, in order to have the possibility to work with different ODE solvers

<a/>

-----------------------------------------------------
1. Instructions to compile the Minimalist version

  - Open the `mybashrc.minimalist` and adjust the paths to the compulsory external libraries (in particular choose the OpenFOAM version you are working with)
  - Type: `source mybashrc.minimalist`
  - Compile the steady-state solver: from the `solver/fireSimpleSMOKE` folder type `wmake`
  - Compile the unsteady solver: from the `solver/firePimpleSMOKE` folder type `wmake`

-----------------------------------------------------
2. Instructions to compile the Minimalist+MKL version

  - Open the `mybashrc.minimalist.mkl` and adjust the paths to the compulsory external libraries and the paths to the Intel MKL library (in particular choose the OpenFOAM version you are working with)
  -  Type: `source mybashrc.minimalist.mkl`
  - Compile the steady-state solver: from the `solver/fireSimpleSMOKE` folder type `wmake`
  - Compile the unsteady solver: from the `solver/firePimpleSMOKE` folder type `wmake`

-----------------------------------------------------
3. Instructions to compile the Complete version
  - Open the `mybashrc.complete` and adjust the paths to the compulsory external libraries and the Intel MKL library (in particular choose the OpenFOAM version you are working with). You can choose the additional external libraries you want to add to edcSMOKE, by modifying the `EXTERNAL_ODE_SOLVERS` variable: in particular `1` means that the support is requested, while `0` means that no support is requested. Obviously, for each requested library, you need to provide the correct path.
  - Type: `source mybashrc.complete`
  - Compile the steady-state solver: from the `solver/fireSimpleSMOKE` folder type `wmake`
  - Compile the unsteady solver: from the `solver/firePimpleSMOKE` folder type `wmake`

-----------------------------------------------------
4. Instructions to compile with SPARC plugin

  - Adjust the paths in the desired `mybashrc` file to the armadillo libraries.
  -  Type: `source mybashrc.minimalist.mkl`
  - Compile the steady-state solver: from the `solver/fireSimpleSMOKE` folder type `wmake`
  - Compile the unsteady solver: from the `solver/firePimpleSMOKE` folder type `wmake`

<a/>

## Docker Installation 
Only tested for MacOS.

1. Install Docker for your OS [Docker Installation](https://docs.docker.com/engine/install/)
2. OpenFOAM7 + FiReSMOKE is launched from the script `firesmoke2-macos` in this repository. The script needs to be located somewhere on the user’s `PATH` for convenient execution. The following commands will then install in the system-wide /usr/local/bin directory and make the script executable:
  
   `sudo curl --create-dirs -o /usr/local/bin/firesmoke2-macos firesmoke2-macos`
   
   `sudo chmod 755 /usr/local/bin/firesmoke2-macos`
   
   + if you do not have permissions to open the application: `xattr -d com.apple.quarantine /usr/local/bin/firesmoke2-macos`
   
4. The Docker container mounts the user’s file system so that case files are stored permanently. The container mounts the directory from where `firesmoke2-macos` is launched by default, but the user can also specify the directory using the “-d” option.  Mounting the user’s $HOME directory is disallowed.  Where a case-sensitive volume has been created, the container mount directory would typically coincide with the mount directory (or sub-directory) of the volume.  For example, for a case-sensitive volume mounted in the default location, `$HOME/openfoam`:
   
   `cd $HOME/openfoam`
   
   `firesmoke2-macos`



[1]: https://www.opensmokepp.polimi.it/
