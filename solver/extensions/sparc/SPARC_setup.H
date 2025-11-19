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

#if SPARC==1

// Clustering
int nclusters = 0;
std::vector<Eigen::VectorXi> classifier_indices;
std::vector<OpenSMOKE::ThermodynamicsMap_CHEMKIN*> vector_thermodynamicsMapXML;
std::vector<OpenSMOKE::KineticsMap_CHEMKIN*> vector_kineticsMapXML;
std::vector<OpenSMOKE::TransportPropertiesMap_CHEMKIN*> vector_transportMapXML; 


// FITCTREE
const dictionary& fitctreeDictionary = combustionOptionsDictionary.subDict("FitCTree");
const Switch iFitCTree = Switch(fitctreeDictionary.lookup(word("fitctree")));
Eigen::VectorXd fitctree_mu;
Eigen::VectorXd fitctree_sigma;
Eigen::MatrixXd fitctree_pca_coeff;
int fitctree_npcs = 0;
int fitctree_npcs_retained = 0;
autoPtr<std::ofstream> fFitctree;
label fitctreeIntegralUpdate  = 10;
label fitctreeIntegralCounter = fitctreeIntegralUpdate;
if (SPARCswitch == true)
{
if (iFitCTree == true)
{
	Foam::string file_name = fitctreeDictionary.lookup("xml");	
	boost::filesystem::path path_xml_input_file = file_name;
	Foam::string folder_name = fitctreeDictionary.lookup("folder");	
	boost::filesystem::path path_kinetics_folder = folder_name;
	
	std::cout << "Opening XML file containing FitCTree data..." << std::endl;

	// Open the XML file
	boost::property_tree::ptree ptree; // OFv10
  	boost::property_tree::read_xml( (path_xml_input_file).c_str(), ptree ); //OFv10
	//rapidxml::xml_document<> doc;
	std::vector<char> xml_string;
	//OpenSMOKE::OpenInputFileXML(doc, xml_string, path_xml_input_file);
	//rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");
	boost::property_tree::ptree opensmoke_node = ptree.get_child("opensmoke");
	//rapidxml::xml_node<>* opensmoke_node = ptree.first_node("opensmoke");

	// Read number of clusters
	{
		//rapidxml::xml_node<>* nclusters_node = opensmoke_node->first_node("classes");
		boost::property_tree::ptree nclusters_node = opensmoke_node.get_child("classes");
		try
		{
			//nclusters = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(nclusters_node->value())));
			nclusters = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(nclusters_node.data())));
			std::cout << " * Number of clusters: " << nclusters << std::endl;
		}
		catch (...)
		{
			OpenSMOKE::ErrorMessage("Importing data from MATLAB(R) preprocessing", "Error in reading the number of clusters.");
		}
	}

	fitctree_mu.resize(thermodynamicsMapXML->NumberOfSpecies()+1);
	fitctree_sigma.resize(thermodynamicsMapXML->NumberOfSpecies()+1);

	// Read mean values of filtered variables
	{
		std::cout << " * Reading mean values (mu)..." << std::endl;
		//rapidxml::xml_node<>* mu_node = opensmoke_node->first_node("mu");
		boost::property_tree::ptree mu_node = opensmoke_node.get_child("mu");
		try
		{
			std::stringstream fInput;
			//fInput << mu_node->value();
			fInput << mu_node.data();

			for (unsigned int j = 0; j < thermodynamicsMapXML->NumberOfSpecies()+1; j++)
				fInput >> fitctree_mu(j);
		}
		catch (...)
		{
			OpenSMOKE::ErrorMessage("Importing data from MATLAB(R) preprocessing", "Error in reading means of filtered variables.");
		}
	}

	// Read std deviations of filtered variables
	{
		std::cout << " * Reading std deviation values (sigma)..." << std::endl;
		//rapidxml::xml_node<>* sigma_node = opensmoke_node->first_node("sigma");
		boost::property_tree::ptree sigma_node = opensmoke_node.get_child("sigma");
		try
		{
			std::stringstream fInput;
			//fInput << sigma_node->value();
			fInput << sigma_node.data();

			for (unsigned int j = 0; j <  thermodynamicsMapXML->NumberOfSpecies() + 1; j++)
				fInput >> fitctree_sigma(j);
		}
		catch (...)
		{
			OpenSMOKE::ErrorMessage("Importing data from MATLAB(R) preprocessing", "Error in reading std deviations of filtered variables.");
		}
	}

	// Read PC if any
	{
		//rapidxml::xml_node<>* npcs_node = opensmoke_node->first_node("npcs");
		boost::property_tree::ptree npcs_node = opensmoke_node.get_child("npcs");
		try
		{
			//fitctree_npcs = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(npcs_node->value())));
			fitctree_npcs = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(npcs_node.data())));
			std::cout << " * Number of PCs (total): " << fitctree_npcs << std::endl;
		}
		catch (...)
		{
			OpenSMOKE::ErrorMessage("Importing data from MATLAB(R) preprocessing", "Error in reading the number of PCs.");
		}

		if (fitctree_npcs != 0)
		{
			//rapidxml::xml_node<>* npcs_retained_node = opensmoke_node->first_node("npcs-retained");
			boost::property_tree::ptree npcs_retained_node = opensmoke_node.get_child("npcs-retained");
			try
			{
				//fitctree_npcs_retained = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(npcs_retained_node->value())));
				fitctree_npcs_retained = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(npcs_retained_node.data())));
				std::cout << " * Number of PCs (retained): " << fitctree_npcs_retained << std::endl;
			}
			catch (...)
			{
				OpenSMOKE::ErrorMessage("Importing data from MATLAB(R) preprocessing", "Error in reading the number of retained PCs.");
			}

			fitctree_pca_coeff.resize(fitctree_npcs, fitctree_npcs_retained);
			fitctree_pca_coeff.setZero();

			//rapidxml::xml_node<>* pca_coeff_node = opensmoke_node->first_node("pca_coeff");
			boost::property_tree::ptree pca_coeff_node = opensmoke_node.get_child("pca_coeff");
			try
			{
				std::stringstream fInput;
				//fInput << pca_coeff_node->value();
				fInput << pca_coeff_node.data();

				for (unsigned int j = 0; j < fitctree_npcs; j++)
					for (unsigned int i = 0; i < fitctree_npcs_retained; i++)
						fInput >> fitctree_pca_coeff(j,i);
			}
			catch (...)
			{
				OpenSMOKE::ErrorMessage("Importing data from MATLAB(R) preprocessing", "Error in reading weights of PCA.");
			}
		}
	}

	// Read thermodynamics and kinetics maps
	std::cout << " * Reading reduced kinetic mechanisms..." << std::endl;
	vector_thermodynamicsMapXML.resize(nclusters);
	vector_transportMapXML.resize(nclusters);
	vector_kineticsMapXML.resize(nclusters);
	classifier_indices.resize(nclusters);
	for (int j = 0; j < nclusters; j++)
	{
		std::cout << "   ** Reading reduced mechanism " << j << "..." << std::endl;

		std::stringstream label; label << j;
		const std::string folder_name = "kinetics." + label.str();

		boost::property_tree::read_xml( (path_kinetics_folder / folder_name / "kinetics.xml").c_str(), ptree ); //OFv10
		//rapidxml::xml_document<> doc;
		std::vector<char> xml_string;
		//OpenSMOKE::OpenInputFileXML(doc, xml_string, path_kinetics_folder / folder_name / "kinetics.xml");

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
		vector_thermodynamicsMapXML[j] = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree); 
		vector_transportMapXML[j] = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(ptree); 
		vector_kineticsMapXML[j] = new OpenSMOKE::KineticsMap_CHEMKIN(*vector_thermodynamicsMapXML[j], ptree); 
		//vector_thermodynamicsMapXML[j] = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(doc);
		//vector_transportMapXML[j] = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(doc);
		//vector_kineticsMapXML[j] = new OpenSMOKE::KineticsMap_CHEMKIN(*vector_thermodynamicsMapXML[j], doc);
		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
		std::cout << "Time to read XML file: " << tEnd - tStart << std::endl;

		// Mapping
		{
			unsigned int NC = thermodynamicsMapXML->NumberOfSpecies();
			unsigned int NClocal = vector_thermodynamicsMapXML[j]->NumberOfSpecies();
			classifier_indices[j].resize(NClocal);
			classifier_indices[j].setConstant(-1);
			for (unsigned int i=0;i<NClocal;i++)
				for (unsigned int k=0;k<NC;k++)
					if (thermodynamicsMapXML->NamesOfSpecies()[k] == vector_thermodynamicsMapXML[j]->NamesOfSpecies()[i])
						classifier_indices[j](i) = k;

			for (unsigned int i=0;i<NClocal;i++)
				if (classifier_indices[j](i) == -1)
				{
					Info << "Kinetic mechanism #" << i << endl;
					Info << "Species " << vector_thermodynamicsMapXML[j]->NamesOfSpecies()[i] << " is not included in the main kinetic mechanism" << endl;
					abort();
				}
		}
	}

	if ( !boost::filesystem::exists("fitctree.out") )
	{
		fFitctree.reset(new std::ofstream("fitctree.out", std::ios::out));
		fFitctree().setf(std::ios::scientific);

		fFitctree() << std::setw(20) << "time";
		fFitctree() << std::setw(20) << "mean_species";
		fFitctree() << std::setw(20) << "max_species";
		fFitctree() << std::setw(20) << "mean_reactions";
		fFitctree() << std::setw(20) << "max_reactions";
		fFitctree() << std::setw(20) << "meanCPUchem";
		fFitctree() << std::setw(20) << "maxCPUchem";
		fFitctree() << std::endl;
	}
	else
	{
		fFitctree.reset(new std::ofstream("fitctree.out", std::ios::app));
		fFitctree().setf(std::ios::scientific);
	}	
}
}

// VQ2
const dictionary& vq2Dictionary = combustionOptionsDictionary.subDict("VQ2");
const Switch iVQ2 = Switch(vq2Dictionary.lookup(word("vq2")));
const Switch iVQ2dubug = Switch(vq2Dictionary.lookup(word("debug")));
autoPtr<std::ofstream> fVQ2;
label vq2IntegralUpdate  = 10;
label vq2IntegralCounter = vq2IntegralUpdate;
ClassifierVQ2 VQ2;

if (SPARCswitch == true)
{
if (iVQ2 == true)
{
	// This folder contains a 'kinetics' and a 'data' subfolders
	Foam::string main_folder_name = vq2Dictionary.lookup("folder");	
	boost::filesystem::path path_main_folder = main_folder_name;
	boost::filesystem::path path_kinetics_folder = path_main_folder / "kinetics";
	boost::filesystem::path path_data_folder = path_main_folder / "data";	

	int scaling = 0;
	word scaling_type(vq2Dictionary.lookup("scaling"));
	if (scaling_type == "AUTO")		scaling = 0;
	else if (scaling_type == "VAST")	scaling = 1;
	else if (scaling_type == "PARETO")	scaling = 2;
	else
	{
		Info << "Wrong scaling option: AUTO || VAST || PARETO" << endl;
		abort();
	}

	label eigens  = vq2Dictionary.lookupOrDefault("eigens", 1);

	std::cout << "Opening files containing VQ2 data..." << std::endl;
	VQ2.Setup(path_data_folder.c_str(), eigens, scaling);
	nclusters = VQ2.nc();
	
	if (VQ2.nv() != (thermodynamicsMapXML->NumberOfSpecies()+1))
		OpenSMOKE::ErrorMessage("Importing data for VQ2", "Mismatch between the data and the kinetic mechanism");

	// Read thermodynamics and kinetics maps
	vector_thermodynamicsMapXML.resize(nclusters);
	vector_transportMapXML.resize(nclusters);
	vector_kineticsMapXML.resize(nclusters);
	classifier_indices.resize(nclusters);
	
	for (int j = 0; j < nclusters; j++)
	{
		std::stringstream label; label << j;
		const std::string folder_name = "kinetics." + label.str();

		Foam::string kinetics = combustionOptionsDictionary.lookup("folder"); //OFv10
		boost::filesystem::path path_kinetics = kinetics;
  		boost::property_tree::ptree ptree; // OFv10
  		boost::property_tree::read_xml( (path_kinetics_folder / folder_name / "kinetics.xml").c_str(), ptree ); //OFv10
		//rapidxml::xml_document<> doc;
		std::vector<char> xml_string;
		//OpenSMOKE::OpenInputFileXML(doc, xml_string, path_kinetics_folder / folder_name / "kinetics.xml");

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
		vector_thermodynamicsMapXML[j] = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree); 
		vector_transportMapXML[j] = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(ptree); 
		vector_kineticsMapXML[j] = new OpenSMOKE::KineticsMap_CHEMKIN(*vector_thermodynamicsMapXML[j], ptree);
		//vector_thermodynamicsMapXML[j] = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(doc);
		//vector_transportMapXML[j] = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(doc);
		//vector_kineticsMapXML[j] = new OpenSMOKE::KineticsMap_CHEMKIN(*vector_thermodynamicsMapXML[j], doc);
		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
		std::cout << "Time to read XML file: " << tEnd - tStart << std::endl;

		// Mapping
		{
			unsigned int NC = thermodynamicsMapXML->NumberOfSpecies();
			unsigned int NClocal = vector_thermodynamicsMapXML[j]->NumberOfSpecies();
			classifier_indices[j].resize(NClocal);
			classifier_indices[j].setConstant(-1);
			for (unsigned int i=0;i<NClocal;i++)
				for (unsigned int k=0;k<NC;k++)
					if (thermodynamicsMapXML->NamesOfSpecies()[k] == vector_thermodynamicsMapXML[j]->NamesOfSpecies()[i])
						classifier_indices[j](i) = k;

			for (unsigned int i=0;i<NClocal;i++)
				if (classifier_indices[j](i) == -1)
				{
					Info << "Kinetic mechanism #" << i << endl;
					Info << "Species " << vector_thermodynamicsMapXML[j]->NamesOfSpecies()[i] << " is not included in the main kinetic mechanism" << endl;
					abort();
				}
		}
	}

	if ( !boost::filesystem::exists("vq2.out") )
        {
                fVQ2.reset(new std::ofstream("vq2.out", std::ios::out));
                fVQ2().setf(std::ios::scientific);

                fVQ2() << std::setw(20) << "time";
                fVQ2() << std::setw(20) << "mean_species";
                fVQ2() << std::setw(20) << "max_species";
                fVQ2() << std::setw(20) << "mean_reactions";
                fVQ2() << std::setw(20) << "max_reactions";
                fVQ2() << std::setw(20) << "meanCPUchem";
                fVQ2() << std::setw(20) << "maxCPUchem";
                fVQ2() << std::endl;
        }
        else
        {
                fVQ2.reset(new std::ofstream("vq2.out", std::ios::app));
                fVQ2().setf(std::ios::scientific);
        }
}
}
// SOM
int npca = 0;
Eigen::MatrixXd pca_weights;
Eigen::VectorXd pca_mu;
Eigen::VectorXd pca_sigma;
int nf;
Eigen::VectorXi listretspecies;
const dictionary& somDictionary = combustionOptionsDictionary.subDict("SelfOrganizingMap");
const Switch iSOM = Switch(somDictionary.lookup(word("som")));
autoPtr<std::ofstream> fSom;
label somIntegralUpdate  = 10;
label somIntegralCounter = somIntegralUpdate;
if (SPARCswitch == true)
{
if (iSOM == true)
{
	Foam::string file_name = somDictionary.lookup("xml");	
	boost::filesystem::path path_xml_input_file = file_name;
	Foam::string folder_name = somDictionary.lookup("folder");	
	boost::filesystem::path path_kinetics_folder = folder_name;
	
	std::cout << "Opening XML file containing SOM data..." << std::endl;

	// Open the XML file
	boost::property_tree::ptree ptree; // OFv10
  	boost::property_tree::read_xml( (path_xml_input_file).c_str(), ptree ); //OFv10
	//rapidxml::xml_document<> doc;
	std::vector<char> xml_string;
	//OpenSMOKE::OpenInputFileXML(doc, xml_string, path_xml_input_file);
	boost::property_tree::ptree opensmoke_node = ptree.get_child("opensmoke");
	//rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");

	// Read number of clusters
	{
		//rapidxml::xml_node<>* nclusters_node = opensmoke_node->first_node("classes");
		boost::property_tree::ptree nclusters_node = opensmoke_node.get_child("classes");
		try
		{
			//nclusters = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(nclusters_node->value())));
			nclusters = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(nclusters_node.data())));
			std::cout << " * Number of clusters: " << nclusters << std::endl;
		}
		catch (...)
		{
			OpenSMOKE::ErrorMessage("Importing data from MATLAB(R) preprocessing", "Error in reading the number of clusters.");
		}
	}

	// Read number of principal components
	{
		//rapidxml::xml_node<>* npca_node = opensmoke_node->first_node("principal-components");
		boost::property_tree::ptree npca_node = opensmoke_node.get_child("principal-components");
		try
		{
			//npca = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(npca_node->value())));
			npca = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(npca_node.data())));
			std::cout << " * Number of principal components: " << npca << std::endl;
		}
		catch (...)
		{
			OpenSMOKE::ErrorMessage("Importing data from MATLAB(R) preprocessing", "Error in reading the number principal components.");
		}
	}

	// Read number retained species 
	{
		//rapidxml::xml_node<>* nretspecies_node = opensmoke_node->first_node("number-retained-species");
		boost::property_tree::ptree nretspecies_node = opensmoke_node.get_child("number-retained-species");
		try
		{
			//nf = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(nretspecies_node->value())));
			nf = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(nretspecies_node.data())));
			std::cout << " * Number of retained species: " << nf << std::endl;
		}
		catch (...)
		{
			OpenSMOKE::ErrorMessage("Importing data from MATLAB(R) preprocessing", "Error in reading the number of retained species.");
		}
	}

	// Read list retained species 
	if (nf != thermodynamicsMapXML->NumberOfSpecies())
	{
		//rapidxml::xml_node<>* listretspecies_node = opensmoke_node->first_node("list-retained-species");
		boost::property_tree::ptree listretspecies_node = opensmoke_node.get_child("list-retained-species");
		try
		{
			listretspecies.resize(nf);

			std::stringstream fInput;
			//fInput << listretspecies_node->value();
			fInput << listretspecies_node.data();

			for (unsigned int j = 0; j < nf; j++)
				fInput >> listretspecies(j);	// (1-index based)

			Info << "Retained species: " << nf << endl;
			for (unsigned int j = 0; j < nf; j++)
				Info << thermodynamicsMapXML->NamesOfSpecies()[listretspecies(j) - 1] << endl;
		}
		catch (...)
		{
			OpenSMOKE::ErrorMessage("Importing data from MATLAB(R) preprocessing", "Error in reading the list of retained species.");
		}
	}

	pca_weights.resize(nf+1,npca);
	pca_mu.resize(nf+1);
	pca_sigma.resize(nf+1);

	// Read mean values of filtered variables
	{
		//rapidxml::xml_node<>* mu_node = opensmoke_node->first_node("mu");
		boost::property_tree::ptree mu_node = opensmoke_node.get_child("mu");
		try
		{
			std::stringstream fInput;
			//fInput << mu_node->value();
			fInput << mu_node.data();

			for (unsigned int j = 0; j < nf+1; j++)
				fInput >> pca_mu(j);
		}
		catch (...)
		{
			OpenSMOKE::ErrorMessage("Importing data from MATLAB(R) preprocessing", "Error in reading means of filtered variables.");
		}
	}

	// Read std deviations of filtered variables
	{
		//rapidxml::xml_node<>* sigma_node = opensmoke_node->first_node("sigma");
		boost::property_tree::ptree sigma_node = opensmoke_node.get_child("sigma");
		try
		{
			std::stringstream fInput;
			//fInput << sigma_node->value();
			fInput << sigma_node.data();

			for (unsigned int j = 0; j < nf + 1; j++)
				fInput >> pca_sigma(j);
		}
		catch (...)
		{
			OpenSMOKE::ErrorMessage("Importing data from MATLAB(R) preprocessing", "Error in reading std deviations of filtered variables.");
		}
	}

	// Read PCA weights
	{
		//rapidxml::xml_node<>* w_node = opensmoke_node->first_node("weights");
		boost::property_tree::ptree w_node = opensmoke_node.get_child("weights");
		try
		{
			std::stringstream fInput;
			//fInput << w_node->value();
			fInput << w_node.data();

			for (unsigned int j = 0; j < nf + 1; j++)
				for (unsigned int i = 0; i < npca; i++)
					fInput >> pca_weights(j,i);
		}
		catch (...)
		{
			OpenSMOKE::ErrorMessage("Importing data from MATLAB(R) preprocessing", "Error in reading weights of PCA.");
		}
	}

	// Read thermodynamics and kinetics maps
	vector_thermodynamicsMapXML.resize(nclusters);
	vector_transportMapXML.resize(nclusters);
	vector_kineticsMapXML.resize(nclusters);
	classifier_indices.resize(nclusters);
	for (int j = 0; j < nclusters; j++)
	{
		std::stringstream label; label << j;
		const std::string folder_name = "kinetics." + label.str();

		boost::property_tree::read_xml( (path_kinetics_folder / folder_name / "kinetics.xml").c_str(), ptree ); //OFv10
		//rapidxml::xml_document<> doc;
		std::vector<char> xml_string;
		//OpenSMOKE::OpenInputFileXML(doc, xml_string, path_kinetics_folder / folder_name / "kinetics.xml");

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
		vector_thermodynamicsMapXML[j] = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree); 
		vector_transportMapXML[j] = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(ptree); 
		vector_kineticsMapXML[j] = new OpenSMOKE::KineticsMap_CHEMKIN(*vector_thermodynamicsMapXML[j], ptree); 
		//vector_thermodynamicsMapXML[j] = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(doc);
		//vector_transportMapXML[j] = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(doc);
		//vector_kineticsMapXML[j] = new OpenSMOKE::KineticsMap_CHEMKIN(*vector_thermodynamicsMapXML[j], doc);
		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
		std::cout << "Time to read XML file: " << tEnd - tStart << std::endl;

		// Mapping
		{
			unsigned int NC = thermodynamicsMapXML->NumberOfSpecies();
			unsigned int NClocal = vector_thermodynamicsMapXML[j]->NumberOfSpecies();
			classifier_indices[j].resize(NClocal);
			classifier_indices[j].setConstant(-1);
			for (unsigned int i=0;i<NClocal;i++)
				for (unsigned int k=0;k<NC;k++)
					if (thermodynamicsMapXML->NamesOfSpecies()[k] == vector_thermodynamicsMapXML[j]->NamesOfSpecies()[i])
						classifier_indices[j](i) = k;

			for (unsigned int i=0;i<NClocal;i++)
				if (classifier_indices[j](i) == -1)
				{
					Info << "Kinetic mechanism #" << i << endl;
					Info << "Species " << vector_thermodynamicsMapXML[j]->NamesOfSpecies()[i] << " is not included in the main kinetic mechanism" << endl;
					abort();
				}
		}
	}

	if ( !boost::filesystem::exists("som.out") )
                        {
                                fSom.reset(new std::ofstream("som.out", std::ios::out));
                                fSom().setf(std::ios::scientific);

                                fSom() << std::setw(20) << "time";
                                fSom() << std::setw(20) << "mean_species";
                                fSom() << std::setw(20) << "max_species";
                                fSom() << std::setw(20) << "mean_reactions";
                                fSom() << std::setw(20) << "max_reactions";
                                fSom() << std::setw(20) << "meanCPUchem";
                                fSom() << std::setw(20) << "maxCPUchem";
                                fSom() << std::endl;
                        }
                        else
                        {
                                fSom.reset(new std::ofstream("som.out", std::ios::app));
                                fSom().setf(std::ios::scientific);
                        }	
}
}

// SOFTMAX
const dictionary& softmaxDictionary = combustionOptionsDictionary.subDict("SOFTMAX");
const Switch iSOFTMAX = Switch(softmaxDictionary.lookup(word("SOFTMAX")));
autoPtr<std::ofstream> fSOFTMAX;
label softmaxIntegralUpdate  = 10;
label softmaxIntegralCounter = softmaxIntegralUpdate;

ClassifierSOFTMAX SOFTMAX;

if (SPARCswitch == true)
{
if (iSOFTMAX == true)
{
	// This folder contains a 'kinetics' and a 'data' subfolders
	Foam::string main_folder_name = softmaxDictionary.lookup("folder");	
	boost::filesystem::path path_main_folder = main_folder_name;
	boost::filesystem::path path_kinetics_folder = path_main_folder / "kinetics";
	boost::filesystem::path path_data_folder = path_main_folder / "data";

	std::cout << "Opening files containing SOFTMAX data..." << std::endl;
	SOFTMAX.Setup(path_data_folder.c_str());
	nclusters = SOFTMAX.nc();
	
	if (SOFTMAX.nv() != (thermodynamicsMapXML->NumberOfSpecies()+1))
		OpenSMOKE::ErrorMessage("Importing data for SOFTMAX", "Mismatch between the data and the kinetic mechanism");

	// Read thermodynamics and kinetics maps
	vector_thermodynamicsMapXML.resize(nclusters);
	vector_transportMapXML.resize(nclusters);
	vector_kineticsMapXML.resize(nclusters);
	classifier_indices.resize(nclusters);
	
	for (int j = 0; j < nclusters; j++)
	{
		std::stringstream label; label << j;
		const std::string folder_name = "kinetics." + label.str();

		boost::property_tree::ptree ptree; // OFv10
		boost::property_tree::read_xml( (path_kinetics_folder / folder_name / "kinetics.xml").c_str(), ptree ); //OFv10
		//rapidxml::xml_document<> doc;
		std::vector<char> xml_string;
		//OpenSMOKE::OpenInputFileXML(doc, xml_string, path_kinetics_folder / folder_name / "kinetics.xml");

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
		vector_thermodynamicsMapXML[j] = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree); 
		vector_transportMapXML[j] = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(ptree); 
		vector_kineticsMapXML[j] = new OpenSMOKE::KineticsMap_CHEMKIN(*vector_thermodynamicsMapXML[j], ptree); 
		//vector_thermodynamicsMapXML[j] = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(doc);
		//vector_transportMapXML[j] = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(doc);
		//vector_kineticsMapXML[j] = new OpenSMOKE::KineticsMap_CHEMKIN(*vector_thermodynamicsMapXML[j], doc);
		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
		std::cout << "Time to read XML file: " << tEnd - tStart << std::endl;
		// Mapping
		{
			unsigned int NC = thermodynamicsMapXML->NumberOfSpecies();
			unsigned int NClocal = vector_thermodynamicsMapXML[j]->NumberOfSpecies();
			classifier_indices[j].resize(NClocal);
			classifier_indices[j].setConstant(-1);
			for (unsigned int i=0;i<NClocal;i++)
				for (unsigned int k=0;k<NC;k++)
					if (thermodynamicsMapXML->NamesOfSpecies()[k] == vector_thermodynamicsMapXML[j]->NamesOfSpecies()[i])
						classifier_indices[j](i) = k;

			for (unsigned int i=0;i<NClocal;i++)
				if (classifier_indices[j](i) == -1)
				{
					Info << "Kinetic mechanism #" << i << endl;
					Info << "Species " << vector_thermodynamicsMapXML[j]->NamesOfSpecies()[i] << " is not included in the main kinetic mechanism" << endl;
					abort();
				}
		}
	}

	if ( !boost::filesystem::exists("softmax.out") )
        {
                fSOFTMAX.reset(new std::ofstream("softmax.out", std::ios::out));
                fSOFTMAX().setf(std::ios::scientific);

                fSOFTMAX() << std::setw(20) << "time";
                fSOFTMAX() << std::setw(20) << "mean_species";
                fSOFTMAX() << std::setw(20) << "max_species";
                fSOFTMAX() << std::setw(20) << "mean_reactions";
                fSOFTMAX() << std::setw(20) << "max_reactions";
                fSOFTMAX() << std::setw(20) << "meanCPUchem";
                fSOFTMAX() << std::setw(20) << "maxCPUchem";
                fSOFTMAX() << std::endl;
        }
        else
        {
                fSOFTMAX.reset(new std::ofstream("softmax.out", std::ios::app));
                fSOFTMAX().setf(std::ios::scientific);
        }
}
}

// NEURAL
const dictionary& neuralDictionary = combustionOptionsDictionary.subDict("NEURAL");
const Switch iNEURAL = Switch(neuralDictionary.lookup(word("NEURAL")));
autoPtr<std::ofstream> fNEURAL;
label neuralIntegralUpdate  = 10;
label neuralIntegralCounter = neuralIntegralUpdate;

ClassifierNEURAL NEURAL;

if (SPARCswitch == true)
{
if (iNEURAL == true)
{
        // This folder contains a 'kinetics' and a 'data' subfolders
        Foam::string main_folder_name = neuralDictionary.lookup("folder");
        boost::filesystem::path path_main_folder = main_folder_name;
        boost::filesystem::path path_kinetics_folder = path_main_folder / "kinetics";
        boost::filesystem::path path_data_folder = path_main_folder / "data";

        std::cout << "Opening files containing NEURAL data..." << std::endl;
        NEURAL.Setup(path_data_folder.c_str());
        nclusters = NEURAL.nc();

        if (NEURAL.nv() != (thermodynamicsMapXML->NumberOfSpecies()+1))
                OpenSMOKE::ErrorMessage("Importing data for NEURAL", "Mismatch between the data and the kinetic mechanism");

        // Read thermodynamics and kinetics maps
        vector_thermodynamicsMapXML.resize(nclusters);
        vector_transportMapXML.resize(nclusters);
        vector_kineticsMapXML.resize(nclusters);
        classifier_indices.resize(nclusters);

        for (int j = 0; j < nclusters; j++)
        {
                std::stringstream label; label << j;
                const std::string folder_name = "kinetics." + label.str();

				boost::property_tree::ptree ptree; // OFv10
                boost::property_tree::read_xml( (path_kinetics_folder / folder_name / "kinetics.xml").c_str(), ptree ); //OFv10
				//rapidxml::xml_document<> doc;
				std::vector<char> xml_string;
				//OpenSMOKE::OpenInputFileXML(doc, xml_string, path_kinetics_folder / folder_name / "kinetics.xml");

				double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
				vector_thermodynamicsMapXML[j] = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree); 
				vector_transportMapXML[j] = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(ptree); 
				vector_kineticsMapXML[j] = new OpenSMOKE::KineticsMap_CHEMKIN(*vector_thermodynamicsMapXML[j], ptree); 
				//vector_thermodynamicsMapXML[j] = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(doc);
				//vector_transportMapXML[j] = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(doc);
				//vector_kineticsMapXML[j] = new OpenSMOKE::KineticsMap_CHEMKIN(*vector_thermodynamicsMapXML[j], doc);
				double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
				std::cout << "Time to read XML file: " << tEnd - tStart << std::endl;

                // Mapping
                {
                        unsigned int NC = thermodynamicsMapXML->NumberOfSpecies();
                        unsigned int NClocal = vector_thermodynamicsMapXML[j]->NumberOfSpecies();
                        classifier_indices[j].resize(NClocal);
                        classifier_indices[j].setConstant(-1);
                        for (unsigned int i=0;i<NClocal;i++)
                                for (unsigned int k=0;k<NC;k++)
                                        if (thermodynamicsMapXML->NamesOfSpecies()[k] == vector_thermodynamicsMapXML[j]->NamesOfSpecies()[i])
                                                classifier_indices[j](i) = k;

                        for (unsigned int i=0;i<NClocal;i++)
                                if (classifier_indices[j](i) == -1)
                                {
                                        Info << "Kinetic mechanism #" << i << endl;
                                        Info << "Species " << vector_thermodynamicsMapXML[j]->NamesOfSpecies()[i] << " is not included in the main kinetic mechanism" << endl;
                                        abort();
                                }
                }
        }

        if ( !boost::filesystem::exists("neural.out") )
        {
                fNEURAL.reset(new std::ofstream("neural.out", std::ios::out));
                fNEURAL().setf(std::ios::scientific);

                fNEURAL() << std::setw(20) << "time";
                fNEURAL() << std::setw(20) << "mean_species";
                fNEURAL() << std::setw(20) << "max_species";
                fNEURAL() << std::setw(20) << "mean_reactions";
                fNEURAL() << std::setw(20) << "max_reactions";
                fNEURAL() << std::setw(20) << "meanCPUchem";
                fNEURAL() << std::setw(20) << "maxCPUchem";
                fNEURAL() << std::endl;
        }
        else
        {
                fNEURAL.reset(new std::ofstream("neural.out", std::ios::app));
                fNEURAL().setf(std::ios::scientific);
        }
}
}


//PREDICTIONS
const dictionary& predictionsDictionary = combustionOptionsDictionary.subDict("PREDICTIONS");
const Switch iPREDICTIONS = Switch(predictionsDictionary.lookup(word("PREDICTIONS")));
autoPtr<std::ofstream> fPREDICTIONS;
label predictionsIntegralUpdate  = 10;
label predictionsIntegralCounter = predictionsIntegralUpdate;

PredictorNEURAL PREDICTIONS;

if (SPARCswitch == true)
{
if (iPREDICTIONS == true)
{
        // This folder contains a 'kinetics' and a 'data' subfolders
        Foam::string main_folder_name = predictionsDictionary.lookup("folder");
        boost::filesystem::path path_main_folder = main_folder_name;
        boost::filesystem::path path_kinetics_folder = path_main_folder / "kinetics";
        boost::filesystem::path path_data_folder = path_main_folder / "data";

/*	double muOUT = 0.;
        muOUT = readScalar(predictionsDictionary.lookup("muOUT"));

 
        double sigmaOUT = 0.;
        sigmaOUT = readScalar(predictionsDictionary.lookup("sigmaOUT"));
*/
        std::cout << "Opening files containing PREDICTIONS data..." << std::endl;
        PREDICTIONS.Setup(path_data_folder.c_str());
        nclusters = PREDICTIONS.nc();

        if (PREDICTIONS.nv() != (thermodynamicsMapXML->NumberOfSpecies()+1))
                OpenSMOKE::ErrorMessage("Importing data for PREDICTIONS", "Mismatch between the data and the kinetic mechanism");

        // Read thermodynamics and kinetics maps
        vector_thermodynamicsMapXML.resize(nclusters);
        vector_transportMapXML.resize(nclusters);
        vector_kineticsMapXML.resize(nclusters);
        classifier_indices.resize(nclusters);

        for (int j = 0; j < nclusters; j++)
        {
                std::stringstream label; label << j;
                const std::string folder_name = "kinetics." + label.str();
				
				boost::property_tree::ptree ptree; // OFv10
                boost::property_tree::read_xml( (path_kinetics_folder / folder_name / "kinetics.xml").c_str(), ptree ); //OFv10
				//rapidxml::xml_document<> doc;
				std::vector<char> xml_string;
				//OpenSMOKE::OpenInputFileXML(doc, xml_string, path_kinetics_folder / folder_name / "kinetics.xml");

				double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
				vector_thermodynamicsMapXML[j] = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree); 
				vector_transportMapXML[j] = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(ptree); 
				vector_kineticsMapXML[j] = new OpenSMOKE::KineticsMap_CHEMKIN(*vector_thermodynamicsMapXML[j], ptree); 
				//vector_thermodynamicsMapXML[j] = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(doc);
				//vector_transportMapXML[j] = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(doc);
				//vector_kineticsMapXML[j] = new OpenSMOKE::KineticsMap_CHEMKIN(*vector_thermodynamicsMapXML[j], doc);
				double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
				std::cout << "Time to read XML file: " << tEnd - tStart << std::endl;

                // Mapping
                {
                        unsigned int NC = thermodynamicsMapXML->NumberOfSpecies();
                        unsigned int NClocal = vector_thermodynamicsMapXML[j]->NumberOfSpecies();
                        classifier_indices[j].resize(NClocal);
                        classifier_indices[j].setConstant(-1);
                        for (unsigned int i=0;i<NClocal;i++)
                                for (unsigned int k=0;k<NC;k++)
                                        if (thermodynamicsMapXML->NamesOfSpecies()[k] == vector_thermodynamicsMapXML[j]->NamesOfSpecies()[i])
                                                classifier_indices[j](i) = k;

                        for (unsigned int i=0;i<NClocal;i++)
                                if (classifier_indices[j](i) == -1)
                                {
                                        Info << "Kinetic mechanism #" << i << endl;
                                        Info << "Species " << vector_thermodynamicsMapXML[j]->NamesOfSpecies()[i] << " is not included in the main kinetic mechanism" << endl;
                                        abort();
                                }
                }
        }

        if ( !boost::filesystem::exists("predictions.out") )
        {
                fPREDICTIONS.reset(new std::ofstream("predictions.out", std::ios::out));
                fPREDICTIONS().setf(std::ios::scientific);

                fPREDICTIONS() << std::setw(20) << "time";
                fPREDICTIONS() << std::setw(20) << "mean_species";
                fPREDICTIONS() << std::setw(20) << "max_species";
                fPREDICTIONS() << std::setw(20) << "mean_reactions";
                fPREDICTIONS() << std::setw(20) << "max_reactions";
                fPREDICTIONS() << std::setw(20) << "meanCPUchem";
                fPREDICTIONS() << std::setw(20) << "maxCPUchem";
                fPREDICTIONS() << std::endl;
        }
        else
        {
                fPREDICTIONS.reset(new std::ofstream("predictions.out", std::ios::app));
                fPREDICTIONS().setf(std::ios::scientific);
        }
}
}

#endif
