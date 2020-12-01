// This file contains  some common helper functions used by other classes in this project.


// See readme on Github for references
#ifndef HELPER_H
#define HELPER_H

#include "constants.C"



double rad2degree(double radians){
	return radians*180/PI;
}

double degree2rad(double degrees){
	return degrees*PI/180;
}



// This method gives the 4-vectors of incident and outgoing particles for elastic scattering process.
// See definition of e,p,e_dash,p_dash,q below.
// These 4 momenta a setup on the basis on a caclulation made in my notes ( 8 Oct 2020 )
// Energy of incident electron in eV
// theta between 0 to pi ( in radians )
ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* 
					elastic4Vectors(double E, double theta){

	double c = LIGHT_SPEED; 
	double mp = MASS_PROTON; 
	double me = MASS_ELECTRON;
	double cos = COS(theta);


	double c4 = c*c*c*c;  // TODO: cleanup the variable c to imporove efficiency

	//double E = 10*1e10;  // 10 GeV   // Being give via arguement now.
	//double theta = TMath::Pi()/4;		// Being give via arguement now.

	// Solving equation of form Ax^2 + Bx + C = 0
	double A = 4*c*c*(E*E - me*me*c4 - TMath::Power((E + mp*c*c)/cos,2));
	double B = 8*c*TMath::Sqrt(E*E - me*me*c4)*(E*mp + me*me*c*c)*c*c;
	double C = TMath::Power(2*c*c*(E*mp + me*me*c*c),2) - 4*me*me*c4*TMath::Power(E + mp*c*c,2);

	double D = B*B - 4*A*C;  // Discriminant of quadratic euqation


	double zero_plus = (-B - TMath::Sqrt(D))/(2*A);
	double zero_minus = (-B + TMath::Sqrt(D))/(2*A);
	double pz = TMath::Max(zero_plus,zero_minus);
	if ( cos < 0 )
		pz = TMath::Min(zero_plus,zero_minus);

	if ( cos*pz < 0 ){
		cout<<" Momentum of scattered electron ( = pz ) has ambiguous sign.";
		exit(1);
	}
	// If theta from 0 to pi/2, scattering is forward => choose pz positive
	// If theta from pi/2 to pi, scattering is backward => choose pz positive


	// e = 4 vector of incident electron
	// e_dash = 4 vector of scattered electron
	// p = 4 vector of incident proton
	// p_dash = 4 vector of recoiled proton
	// q = e - e_dash = 4 momentum transfer
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > e,p,e_dash,p_dash,q;

	e.SetE(E);
	e.SetPz(TMath::Sqrt(E*E - me*me*c4)/c);

	p.SetE(mp*c*c);

	e_dash.SetPz(pz);
	e_dash.SetPy(TMath::Tan(theta)*pz);
	e_dash.SetE(TMath::Sqrt(me*me*c4 + c*c*pz*pz*TMath::Power(1/cos,2)));

	// TODO: compare with the actual formula for E' from (5.15) given in [1]
	// Ans: Tested. My formula is also correct.


	p_dash = e + p - e_dash;
	// TODO: do some sanity testing on calculated 4 vectors and exit in case of negative energiesfor example.
	// Ans: Done


	q = e - e_dash;  // TODO: does the order make a difference here ? Ans: Yes. It does.


	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > *arrayOf4Vectors = 
						new ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >[5];
	arrayOf4Vectors[0] = e;
	arrayOf4Vectors[1] = p;
	arrayOf4Vectors[2] = e_dash;
	arrayOf4Vectors[3] = p_dash;
	arrayOf4Vectors[4] = q;

	return arrayOf4Vectors;

}




// What should be the incidentElectronEnergy in order for the elastically scattered electron to have
// a given value of Q2 at theta ? 
// Use this method if you want to fix Q2 instead of fixing incident electron energy
double incidentEnergyGenerator( double Q2, double theta){
	double a = ( 1 - COS(theta))*Q2/MASS_PROTON;
	double b = POWER(SIN(theta/2),2);
	double c = SQRT(a*a + 16*b*Q2);
	double small = (a-c)/(8*b);
	double large = (a+c)/(8*b);
	return large;
}



// E = incidentElectronEnergy
double calculateScatteredElectronEnergy(double E, double theta){
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* vectors = elastic4Vectors(E, theta);
	return vectors[2].E();
}


// gives the probabily of an electron coming out with energy E through a radiator of length t
// and incident energy E0
// based upon A.1 from [4]
double energyThroughExternalRadiator( double E, double E0, double t){
	//for this function t MUST BE greater than of equal to ln(2)
	double ln2 = LN(2.0);

	double num = POWER(LN(E0/E),t/ln2 - 1) / E0;
	double den = TMath::Gamma(t/ln2);

	if ( num < 0 || den < 0){

		// I thought that numerator and denominator must be positive. But that's not the case.
		//cout<<"  Invalid argument: E="<<E<<" ,E0="<<E0<<" ,t="<<t<<" ,num="<<num<<" ,den="<<den<<endl;
   		//exit(1);
	}

	return num/den;
}

// see equation 58 in [3] for formula and definitions
double lambda_e(double incidentElectronEnergy){
	double ans = ALPHA/PI;
	ans *= LN(4* (incidentElectronEnergy*incidentElectronEnergy - MASS_ELECTRON*MASS_ELECTRON)/ (MASS_ELECTRON*MASS_ELECTRON) )-1;
	return ans;
}



// extra_lambda is defined as lambda_dash - lambda
// see equation 62 in [3] for formula and definitions of lambda and lambda_dash
double extra_lambda( double k, double k_dash, double theta){
	double ans = ALPHA/PI;
	double temp1 = 2*LN(k/k_dash);
	double temp2 = LN( (1 - COS(theta))/2 );
	return ans*(temp1+temp2);
}


double simulateIncomingRadiator(double & E0, TRandom *randomObject, double t, double & lostEnergy){
	double electron_upper_limit = E0 - MIN_ENERGY_PHOTON;
	double electron_lower_limit = MAX(0,E0 - MAX_ENERGY_PHOTON);

	// Only photon between MIN_ENERGY_PHOTON and MAX_ENERGY_PHOTON should be produced.
	// This puts a constraint on the energy E of electron emitted from radiator
	double E = randomObject->Uniform(electron_lower_limit, electron_upper_limit);
	double crossSection = energyThroughExternalRadiator(E,E0,t);

	// calculating the photon energy before changing E0 by reference
	lostEnergy = E0 - E;
	E0 = E;
	return crossSection*(electron_upper_limit - electron_lower_limit);

}

double simulateOutgoingRadiator(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* vectors, TRandom *randomObject, double t, double & lostEnergy){
	double E0 = vectors[2].E();

	double electron_upper_limit = E0 - MIN_ENERGY_PHOTON;
	double electron_lower_limit = MAX(0,E0 - MAX_ENERGY_PHOTON);

	// Only photon between MIN_ENERGY_PHOTON and MAX_ENERGY_PHOTON should be produced.
	// This puts a constraint on the energy E of electron emitted from radiator
	double E = randomObject->Uniform(electron_lower_limit, electron_upper_limit);
	double crossSection = energyThroughExternalRadiator(E,E0,t);

	lostEnergy = E0 - E;

	// Changing the energy of outgoing electron now. It will jeopardize the on-shell relation
	// Don't use the momenta of this vector after this point before fixing this. 
	// TODO: Fix the outgoing momenta below
	vectors[2].SetE(E);

	return crossSection*(electron_upper_limit - electron_lower_limit);

}


#endif