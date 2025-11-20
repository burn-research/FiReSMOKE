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

#ifndef NLSSystemVirtualClassWithOpenSMOKEVectors_H
#define	NLSSystemVirtualClassWithOpenSMOKEVectors_H

namespace OpenSMOKE
{
	
	class NLSSystemVirtualClassWithOpenSMOKEVectors
	{
	public:

		unsigned int NumberOfEquations() { return ne_; }
		virtual void GetEquations(const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy) = 0;
		virtual void PrintResults(const OpenSMOKE::OpenSMOKEVectorDouble& y) { };

	protected:

		unsigned int ne_;

		void MemoryAllocation()
		{
			ChangeDimensions(ne_, &y_, true);
			ChangeDimensions(ne_, &dy_, false);
		}

		void Equations(const Eigen::VectorXd &Y, Eigen::VectorXd &DY)
		{
			y_.CopyFrom(Y.data());
			GetEquations(y_, dy_);
			dy_.CopyTo(DY.data());
		}

		void Jacobian(const Eigen::VectorXd &Y, const double t, Eigen::MatrixXd &J) { };

		void Print(unsigned int number_iterations, const double t, const double phiW, const Eigen::VectorXd &Y, const Eigen::VectorXd &f)
		{
			y_.CopyFrom(Y.data());
			PrintResults(y_);
		}

	private:

		OpenSMOKE::OpenSMOKEVectorDouble  y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	};
}
	
#endif	// NLSSystemVirtualClassWithOpenSMOKEVectors
