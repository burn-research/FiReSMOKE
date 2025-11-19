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



class PredictorNEURAL
{
        public:
		//double muOUT;
		//double sigmaOUT;
			
                void Setup(const std::string folder_name);
                void PredictNN(arma::Mat<double> Y, std::vector<int>& idx);

                unsigned int nc() const { return nc_; }
                unsigned int nv() const { return nv_; }

                double predictedTarget(int i) const{ return predictedTarget_[i]; }

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
                arma::Mat<double> RweightsHL1Matrix_;
                arma::Mat<double> RbiasHL1Matrix_;
                arma::Mat<double> RweightsHL2Matrix_;
                arma::Mat<double> RbiasHL2Matrix_;
                arma::Mat<double> RweightsOUTMatrix_;
                arma::Mat<double> RbiasOUTMatrix_;
                arma::Mat<double> Xtraining_;
                arma::field<arma::Mat<double>> MuSigma_;
                arma::Mat<double> predictedTarget_;

		arma::Mat<double> muOUT;
                arma::Mat<double> sigmaOUT;

};

void PredictorNEURAL::Setup(const std::string folder_name)
{
        std::string filename_Xtrain = folder_name + "/" + "X.txt";
	std::string filename_RweightsHL1 = folder_name + "/" + "w1.txt";
        std::string filename_RbiasHL1 = folder_name + "/" + "b1.txt";
        std::string filename_RweightsHL2 = folder_name + "/" + "w2.txt";
        std::string filename_RbiasHL2 = folder_name + "/" + "b2.txt";
        std::string filename_RweightsOUT = folder_name + "/" + "wOUT.txt";
        std::string filename_RbiasOUT = folder_name + "/" + "bOUT.txt";

	std::string filename_muOUT = folder_name + "/" + "muOUT.txt";
        std::string filename_sigmaOUT = folder_name + "/" + "sigmaOUT.txt";

	std::cout << "PATH: " << folder_name << std::endl;
        std::cout << "Nome1: " << filename_RweightsHL1 << std::endl;


	std::cout << "Reading REGRESSION hidden layer n1 weights.." << std::endl;
        RweightsHL1Matrix_ = GetMatrix<double>(filename_RweightsHL1);
        std::cout << "REGRESSION hidden layer weights n1 loaded successfully" << std::endl;

        std::cout << "Reading REGRESSION hidden layer n1 biases.." << std::endl;
        RbiasHL1Matrix_ = GetMatrix<double>(filename_RbiasHL1);
        std::cout << "REGRESSION hidden layer n1 biases loaded successfully" << std::endl;

        std::cout << "Reading REGRESSION hidden layer n2 weights.." << std::endl;
        RweightsHL2Matrix_ = GetMatrix<double>(filename_RweightsHL2);
        std::cout << "REGRESSION hidden layer weights n2 loaded successfully" << std::endl;

        std::cout << "Reading REGRESSION  hidden layer n2 biases.." << std::endl;
        RbiasHL2Matrix_ = GetMatrix<double>(filename_RbiasHL2);
        std::cout << "REGRESSION hidden layer n2 biases loaded successfully" << std::endl;

        std::cout << "Reading REGRESSION output layer weights.." << std::endl;
        RweightsOUTMatrix_ = GetMatrix<double>(filename_RweightsOUT);
        std::cout << "REGRESSION output layer weights loaded successfully" << std::endl;

        std::cout << "Reading REGRESSION output layer biases.." << std::endl;
        RbiasOUTMatrix_ = GetMatrix<double>(filename_RbiasOUT);
        std::cout << "REGRESSION output layer biases loaded successfully" << std::endl;
	
	muOUT = GetMatrix<double>(filename_muOUT);
        sigmaOUT = GetMatrix<double>(filename_sigmaOUT);


	std::cout << "The REGRESSION weights must be transposed (PYTHON is activated)" << std::endl;
        RweightsHL1Matrix_ = RweightsHL1Matrix_.t();
        RweightsHL2Matrix_ = RweightsHL2Matrix_.t();
        RweightsOUTMatrix_ = RweightsOUTMatrix_.t();
        std::cout << "Transposed." <<std::endl;
        


	std::cout << "Reading training matrix" << std::endl;
        Xtraining_ = GetMatrix<double>(filename_Xtrain);
        std::cout << "Training matrix loaded successfully" << std::endl;

        nv_ = RweightsHL1Matrix_.n_cols;
        nl_ = RweightsHL1Matrix_.n_rows;
        nc_ = 1;


        std::cout << "Number of variables (nv): " << nv_ << std::endl;
        std::cout << "Number of neurons in the hidden layer of the net (nl):  " << nl_ << std::endl;
	
	std::cout << "CHECK ON MU!!!!: " << muOUT << std::endl;
        std::cout << "CHECK ON SIGMA!!!!: " << sigmaOUT << std::endl;        

        std::cout << "Calculating scaling factors from training matrix..." << std::endl;
        MuSigma_ = GetMuSigma(Xtraining_);
}



void PredictorNEURAL::PredictNN(arma::Mat<double> Y, std::vector<int>& idx)
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


        Y = Y.t();

       for(int i=0; i<idx.size(); i++){
                idx[i] = 0;
        }


	arma::Mat<double> Rtmp = RweightsHL1Matrix_*Y;
        arma::Mat<double> Rtemp = Rtmp.each_col() + RbiasHL1Matrix_;
        arma::Mat<double> Ry1 = relu_layer(Rtemp);
        arma::Mat<double> Rtmp2 = RweightsHL2Matrix_*Ry1;
        arma::Mat<double> Rtemp2 = Rtmp2.each_col() + RbiasHL2Matrix_;
        arma::Mat<double> Ry2 = relu_layer(Rtemp2);
        arma::Mat<double> Rtmp3 = RweightsOUTMatrix_*Ry2;
        arma::Mat<double> Rtemp3 = Rtmp3.each_col() + RbiasOUTMatrix_;
        predictedTarget_ = Rtemp3;

	std::cout << "SECOND CHECK ON MU: " << muOUT(0,0) << std::endl;
        std::cout << "SECOND CHECK ON SIGMA: " << sigmaOUT(0,0) << std::endl;
    
//unscale the prediction
        for(int i=0; i<n_oss; i++){
                predictedTarget_[i] *= (sigmaOUT(0,0) + 1E-16);
        }



        //uncenter the prediction
        for(int i=0; i<n_oss; i++){
                predictedTarget_[i] += (muOUT(0,0));
        }

}



void PredictorNEURAL::FatalErrorMessage(const std::string message)
{
        std::cout << "ClassifierNEURAL Fatal error message" << std::endl;
        std::cout << message << std::endl;
        exit(-1);
}


#endif // PredictorNEURAL_HPP                                                                                                                                                                                                                               107,0-1       12%
