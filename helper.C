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

// Using 5.15 from [1]	
double calculateScatteredElectronEnergy2(double E, double theta){
	double M = MASS_PROTON;
	double numerator = E;
	double cos = COS(theta);
	double denominator = 1 + (E/M)*(1-cos);
	return numerator/denominator;

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


	double c4 = c*c*c*c; 
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
	double correct_e_dash = calculateScatteredElectronEnergy2(E,theta);
	e_dash.SetE(correct_e_dash);


	p_dash = e + p - e_dash;

	q = e - e_dash; 


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


// Using inverse of Rosenbluth formula, this functions tells what the 
// value of incident energy should be to get the particular value of 
// scattered energy at the given theta
// Using (5.15) from Particles and Nuclie book
double incidentEnergyFromScatteredEnergy( double E_scattered, double theta){
	//cout.precision(20);
	//cout<<"Need scattered energy: "<<E_scattered<<endl;
	double ans = E_scattered*MASS_PROTON;
	ans *= 1/(MASS_PROTON - E_scattered*(1-COS(theta)));
	//cout<<"Incident energy should be: "<<ans<<endl;
	return ans;
}

// this calculates q2 for an arbitrary scattering energy of lepton ( not necessarily for elastic scattering )
double calcluateq2arbitrary(double E0, double E1, double theta){
	double ans = 0.0;
	ans += (E0-E1)*(E0-E1);
	double p0 = SQRT(E0*E0 - MASS_ELECTRON*MASS_ELECTRON);
	double p1 = SQRT(E1*E1 - MASS_ELECTRON*MASS_ELECTRON);
	ans += -( p0 - p1*COS(theta))*( p0 - p1*COS(theta));
	ans += -(p1*SIN(theta))*(p1*SIN(theta));
	return ans;
}


double calculateq2(double E, double theta){
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* vectors = elastic4Vectors(E, theta);
	return vectors[4].Dot(vectors[4]);
}

// gives the probabily of an electron coming out with energy E through a radiator of length t
// and incident energy E0
// based upon A.1 from [4]
double energyThroughExternalRadiator2( double E, double E0, double t){
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

// based on (75) from [3]
double energyThroughExternalRadiator( double E, double E0, double t){
	double Eext = E0 - E;
	double bt = t;
	double k = SQRT(E0*E0 - MASS_ELECTRON*MASS_ELECTRON);
	double phi = 1 - Eext/k + (3.0/4)*(Eext/k)*(Eext/k);
	//double phi = 1 - Eext/k;
	double ans = 1.0/GAMMA(1+bt);
	ans *= bt/Eext;
	ans *= POWER(Eext/k,bt);
	ans *= phi;
	return ans;
}

// calculates the square of four vector x as define in (2.3) of [5]
double x2(double e1, double e3, double theta){
	double M = MASS_PROTON;
	double q2 = calcluateq2arbitrary(e1, e3, theta);
	double Q = SQRT( -q2 );
	return POWER(Q + SQRT(Q*Q + 4*M*M),2);
}


double lambda_incident(double incidentElectronEnergy, double scatteredElectronEnergy, double theta){
	double k = SQRT(incidentElectronEnergy*incidentElectronEnergy - MASS_ELECTRON*MASS_ELECTRON);
	double k_dash = SQRT(scatteredElectronEnergy*scatteredElectronEnergy - MASS_ELECTRON*MASS_ELECTRON);
	double m = MASS_ELECTRON;
	double ans = (ALPHA/PI)*(LN(4*k*k/(m*m))-1);
	double correction = (ALPHA/PI)*2*LN(k/k_dash);
	correction += (ALPHA/PI)*LN( (1- COS(theta))/2);

	ans += correction;
	//cout<<"Incoming thickness : "<<ans<<endl;
	//cout<<"With correction : "<<ans + correction<<endl;
	return ans;
}

double lambda_scattered(double incidentElectronEnergy, double scatteredElectronEnergy, double theta){
	double k = SQRT(incidentElectronEnergy*incidentElectronEnergy - MASS_ELECTRON*MASS_ELECTRON);
	double k_dash = SQRT(scatteredElectronEnergy*scatteredElectronEnergy - MASS_ELECTRON*MASS_ELECTRON);
	double m = MASS_ELECTRON;
	double ans = (ALPHA/PI)*(LN(4*k_dash*k_dash/(m*m))-1);
	double correction = (ALPHA/PI)*2*LN(k/k_dash);
	correction += (ALPHA/PI)*LN( (1- COS(theta))/2);
	ans += correction;
	//cout<<"Outgoing thickness : "<<ans<<endl;
	//cout<<"With correction : "<<ans + correction<<endl;
	return ans;
}


double simulateRadiator(double E0, TRandom *randomObject, double t, double photonEnergy){

	double E = E0 - photonEnergy;
	double crossSection = energyThroughExternalRadiator(E,E0,t*1.0);
	//cout<<"Probablity emission: "<<crossSection<<endl;
	//cout<<"Second Probablity emission: "<<energyThroughExternalRadiator2(E,E0,t*1.0)<<endl;
	if( crossSection > 1 ){
		cout<<"-------------------------------------------Probablity too large: "<<crossSection<<endl;
	}

	// calculating the photon energy before changing E0 by reference
	E0 = E;
	if( isinf(crossSection)){
		cout<<"WARNING: Returning probablity for "<<E0<<" to "<<E<<" : "<<crossSection<<endl;	
	}
	
	return crossSection;

}


//Using (29) from reference [3] to account for vacuum polarization
double vacuumPolarization(double q2){
	double ans = 1/(3*PI);
	ans *= -5.0/3.0 + LN(-q2/(MASS_ELECTRON*MASS_ELECTRON));
	return ans;
}


// Using (33) from reference [3] to account for corrections given in Figure 3 of reference [3]
double hardCorrections(double q2){
	double ans = 2*ALPHA;
	ans *= -3*LN(-q2/(MASS_ELECTRON*MASS_ELECTRON))/(4*PI) + 1/PI - vacuumPolarization(q2);
	//cout<<" Q2 is : "<<-q2<<endl;
	//cout<<"Hard correction is: "<<ans<<endl;
	return ans;
}


#endif