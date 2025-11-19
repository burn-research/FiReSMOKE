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

#ifndef ClassifierNEURAL_HPP
#define ClassifierNEURAL_HPP

#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <math.h>
#include <numeric>
#include <armadillo>


arma::field<arma::Mat<double>> GetMuSigma(const arma::Mat<double>& matrix);
arma::Mat<double> softmax(arma::Mat<double> layer);
arma::Mat<double> relu_layer(arma::Mat<double> layer);



class ClassifierNEURAL
{
	public:
	
		void Setup(const std::string folder_name);
		void ClassifyNN(arma::Mat<double> Y, std::vector<int>& idx);

		unsigned int nc() const { return nc_; }
		unsigned int nv() const { return nv_; }
	
	private:

		void FatalErrorMessage(const std::string message);
	
		unsigned int nv_;
		unsigned int nc_;
		unsigned int nl_;

		arma::Mat<double> weightsHL1Matrix_;
                arma::Mat<double> biasHL1Matrix_;
                arma::Mat<double> weightsHL2Matrix_;
                arma::Mat<double> biasHL2Matrix_;
		arma::Mat<double> weightsCLMatrix_;
                arma::Mat<double> biasCLMatrix_;
		arma::Mat<double> Xtraining_;
		arma::field<arma::Mat<double>> MuSigma_;

};

void ClassifierNEURAL::Setup(const std::string folder_name)
{
	//DEEP LEARNING OPTION
	//
	//
	

	std::string filename_weightsHL1 = folder_name + "/" + "weightsHL1.txt";
        std::string filename_biasHL1 = folder_name + "/" + "biasHL1.txt";
	std::string filename_weightsHL2 = folder_name + "/" + "weightsHL2.txt";
        std::string filename_biasHL2 = folder_name + "/" + "biasHL2.txt";
        std::string filename_weightsCL = folder_name + "/" + "weightsCL.txt"; 
        std::string filename_biasCL = folder_name + "/" + "biasCL.txt";
        std::string filename_Xtrain = folder_name + "/" + "X.txt";        
//	END DEEP LEARNING MOD



	 std::cout << "Reading hidden layer n1 weights.." << std::endl;
        weightsHL1Matrix_ = GetMatrix<double>(filename_weightsHL1); 
        std::cout << "Hidden layer weights n1 loaded successfully" << std::endl;

        std::cout << "Reading hidden layer n1 biases.." << std::endl;
        biasHL1Matrix_ = GetMatrix<double>(filename_biasHL1); 
        std::cout << "Hidden layer n1 biases loaded successfully" << std::endl;

	std::cout << "Reading hidden layer n2 weights.." << std::endl;
        weightsHL2Matrix_ = GetMatrix<double>(filename_weightsHL2);
        std::cout << "Hidden layer weights n2 loaded successfully" << std::endl;

        std::cout << "Reading hidden layer n2 biases.." << std::endl;
        biasHL2Matrix_ = GetMatrix<double>(filename_biasHL2);
        std::cout << "Hidden layer n2 biases loaded successfully" << std::endl;

	

        std::cout << "Reading classification layer weights.." << std::endl;
        weightsCLMatrix_ = GetMatrix<double>(filename_weightsCL); 
        std::cout << "Classification layer weights loaded successfully" << std::endl;

        std::cout << "Reading classification layer biases.." << std::endl;
        biasCLMatrix_ = GetMatrix<double>(filename_biasCL); 
        std::cout << "Classification layer biases loaded successfully" << std::endl;


	std::cout << "The weights must be transposed (PYTHON is activated)" << std::endl;
	weightsHL1Matrix_ = weightsHL1Matrix_.t();
	weightsHL2Matrix_ = weightsHL2Matrix_.t();
	weightsCLMatrix_ = weightsCLMatrix_.t();
	std::cout << "Transposed." <<std::endl;

// END OF DEEP MODIFICATION N.2


	std::cout << "Reading training matrix" << std::endl;
        Xtraining_ = GetMatrix<double>(filename_Xtrain);
        std::cout << "Training matrix loaded successfully" << std::endl;

	nv_ = weightsHL1Matrix_.n_cols;
	nl_ = weightsHL1Matrix_.n_rows;
    	nc_ = weightsCLMatrix_.n_rows;


	std::cout << "Number of variables (nv): " << nv_ << std::endl;
	std::cout << "Number of neurons in the hidden layer of the net (nl):  " << nl_ << std::endl;
    	std::cout << "Number of clusters (nc): " << nc_ << std::endl;

	std::cout << "Calculating scaling factors from training matrix..." << std::endl;
    	MuSigma_ = GetMuSigma(Xtraining_);
}

void ClassifierNEURAL::ClassifyNN(arma::Mat<double> Y, std::vector<int>& idx)
{    

	arma::Row<double> mu = MuSigma_(0);
    	arma::Mat<double> sigma = MuSigma_(1);

	const unsigned int n_oss = Y.n_rows;
	const unsigned int n_var = Y.n_cols;

	for(int i=0; i<n_oss; i++){
        	Y.row(i) = Y.row(i) - mu;
	}


	for(int i =0; i<n_oss; i++){
        	for(int j =0; j<n_var; j++){
                        Y(i,j) /= sigma(0,j);
                }
        }

	int n_classes = weightsCLMatrix_.n_rows;
	

	Y = Y.t();

	// DEEP CLASS
	//
	
	arma::Mat<double> tmp = weightsHL1Matrix_*Y;
        arma::Mat<double> temp = tmp.each_col() + biasHL1Matrix_;
	arma::Mat<double> y1 = relu_layer(temp);
	arma::Mat<double> tmp2 = weightsHL2Matrix_*y1;
        arma::Mat<double> temp2 = tmp2.each_col() + biasHL2Matrix_;
        arma::Mat<double> y2 = relu_layer(temp2);
	arma::Mat<double> tmp3 = weightsCLMatrix_*y2;
	arma::Mat<double> temp3 = tmp3.each_col() + biasCLMatrix_;
        arma::Mat<double> output = softmax(temp3);

   	//END OF DEEP CLASS

        for(int i=0; i<output.n_cols; i++){
            idx[i] = output.col(i).index_max(); 
        }
   

}

void ClassifierNEURAL::FatalErrorMessage(const std::string message)
{
	std::cout << "ClassifierNEURAL Fatal error message" << std::endl;
	std::cout << message << std::endl;
	exit(-1);
}


arma::field<arma::Mat<double>> GetMuSigma(const arma::Mat<double>& matrix)
{
    const int n_oss = matrix.n_rows;
    const int n_variab = matrix.n_cols;
    
    arma::Mat<double> rowave(1,matrix.n_cols);
    arma::Mat<double> d(1,matrix.n_cols);

    for(int j=0; j<n_variab; j++){
	rowave.col(j) = arma::min(matrix.col(j));
    }
    
    for(int j=0; j<n_variab; j++){
    	d.col(j) = arma::max(matrix.col(j)) - arma::min(matrix.col(j)) + 1e-16; 
    } 


    arma::field<arma::Mat<double>> MuSigma(2,1);
    MuSigma(0,0) = rowave;
    MuSigma(1,0) = d;

    return MuSigma;
}

arma::Mat<double> softmax(arma::Mat<double> layer){
    double eps = 1e-8;
    arma::Mat<double> maxCols = arma::max(layer,0);
    for(int i = 0; i<layer.n_cols; i++){
        layer.col(i) = layer.col(i) - maxCols(i);
    }

    layer = arma::exp(layer);
    arma::Mat<double> matSum = arma::sum(layer,0);
    for(int i = 0; i<layer.n_cols; i++){
        layer.col(i) = layer.col(i)/matSum(i);
    }

    return layer;
}

arma::Mat<double> relu_layer(arma::Mat<double> layer){
	
	const bool leaky = false;
        const double alpha = 0.001;

    arma::Mat<double> activated_layer(layer.n_rows, layer.n_cols);

    for(int i = 0; i<layer.n_rows; i++){
        for(int j =0; j<layer.n_cols; j++){
            if(layer(i,j) < 0 && leaky == true){
                activated_layer(i,j) = alpha*layer(i,j);
            } else if(layer(i,j) < 0 && leaky == false) {
                activated_layer(i,j) = 0;
            } else if(layer(i,j) >= 0){
                activated_layer(i,j) = layer(i,j);
            }
        }
    }
 
    return activated_layer;
}


#endif // ClassifierNEURAL_HPP
