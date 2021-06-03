/* This file contains  some common helper functions used by other classes in this project.
*  See readme on Github for references
*/
#ifndef HELPER_H
#define HELPER_H

#include "constants.C"



double rad2degree(double radians){
	return radians*180/PI;
}

double degree2rad(double degrees){
	return degrees*PI/180;
}

// Using 5.15 from [3]	
double calculateScatteredElectronEnergy(double E, double theta){
	double M = MASS_PROTON;
	double numerator = E;
	double cos = COS(theta);
	double denominator = 1 + (E/M)*(1-cos);
	return numerator/denominator;

}


/*
*	What should be the incidentElectronEnergy in order for the elastically scattered electron to have
*	a given value of Q2 at theta ? 
*	Use this method to fix Q2 instead of fixing incident electron energy
*/
double incidentEnergyGenerator( double Q2, double theta){
	double a = ( 1 - COS(theta))*Q2/MASS_PROTON;
	double b = POWER(SIN(theta/2),2);
	double c = SQRT(a*a + 16*b*Q2);
	double small = (a-c)/(8*b);
	double large = (a+c)/(8*b);
	return large;
}


/*
*	Using inverse of Rosenbluth formula, this functions tells what the 
*	value of incident energy should be to get the particular value of 
*	scattered energy at the given theta
*	Using (5.15) from Particles and Nuclie book
*/
double incidentEnergyFromScatteredEnergy( double E_scattered, double theta){
	double ans = E_scattered*MASS_PROTON;
	ans *= 1/(MASS_PROTON - E_scattered*(1-COS(theta)));
	return ans;
}

// this calculates q2 for an arbitrary scattering energy of lepton ( not necessarily for elastic scattering )
double calculateq2(double E0, double E1, double theta){
	double ans = 0.0;
	ans += (E0-E1)*(E0-E1);
	double p0 = SQRT(E0*E0 - MASS_ELECTRON*MASS_ELECTRON);
	double p1 = SQRT(E1*E1 - MASS_ELECTRON*MASS_ELECTRON);
	ans += -( p0 - p1*COS(theta))*( p0 - p1*COS(theta));
	ans += -(p1*SIN(theta))*(p1*SIN(theta));
	return ans;
}

double calculateq2(double E, double theta){
	double E_scattered = calculateScatteredElectronEnergy(E,theta);
	return  calculateq2(E, E_scattered, theta);
}

/*
*	gives the probabily of an electron coming out with energy E through a radiator of length t and incident energy E0.
*	Based upon A.1 from [4]
*/
double energyThroughExternalRadiator( double E, double E0, double t){
	double ln2 = LN(2.0);

	double num = POWER(LN(E0/E),t/ln2 - 1) / E0;
	double den = TMath::Gamma(t/ln2);

	return num/den;
}

// based on (75) from [1]
double energyThroughExternalRadiator2( double E, double E0, double t){
	double Eext = E0 - E;
	double bt = t;
	double k = SQRT(E0*E0 - MASS_ELECTRON*MASS_ELECTRON);
	double phi = 1 - Eext/k + (3.0/4)*(Eext/k)*(Eext/k);
	double ans = 1.0/GAMMA(1+bt);
	ans *= bt/Eext;
	ans *= POWER(Eext/k,bt);
	ans *= phi;
	return ans;
}

// calculates the square of four vector x as define in (2.3) of [5]
double x2(double e1, double e3, double theta){
	double M = MASS_PROTON;
	double q2 = calculateq2(e1, e3, theta);
	double Q = SQRT( -q2 );
	return POWER(Q + SQRT(Q*Q + 4*M*M),2);
}


double lambda_incident(double incidentElectronEnergy, double scatteredElectronEnergy, double theta){
	double k = SQRT(incidentElectronEnergy*incidentElectronEnergy - MASS_ELECTRON*MASS_ELECTRON);
	double k_dash = SQRT(scatteredElectronEnergy*scatteredElectronEnergy - MASS_ELECTRON*MASS_ELECTRON);
	double m = MASS_ELECTRON;


	double deltaEP =  (ALPHA/PI)*LN(scatteredElectronEnergy/incidentElectronEnergy)*4*LN(k/k_dash);
	double triangle =  (ALPHA/PI)*LN(scatteredElectronEnergy/incidentElectronEnergy)*2*LN(((1-COS(theta))/2.0));

	double ans = (ALPHA/PI)*(LN(4*k*k/(m*m))-1);

	// Uncomment these lines to add the effects of (60) and (61) from [1]
	// double correction = (ALPHA/PI)*2*LN(k/k_dash);
	// correction += (ALPHA/PI)*LN( (1- COS(theta))/2);
	// ans += correction;

	return ans;
}

double lambda_scattered(double incidentElectronEnergy, double scatteredElectronEnergy, double theta){
	double k = SQRT(incidentElectronEnergy*incidentElectronEnergy - MASS_ELECTRON*MASS_ELECTRON);
	double k_dash = SQRT(scatteredElectronEnergy*scatteredElectronEnergy - MASS_ELECTRON*MASS_ELECTRON);
	double m = MASS_ELECTRON;

	double deltaEP =  (ALPHA/PI)*LN(scatteredElectronEnergy/incidentElectronEnergy)*4*LN(k/k_dash);
	double triangle =  (ALPHA/PI)*LN(scatteredElectronEnergy/incidentElectronEnergy)*2*LN(((1-COS(theta))/2.0));

	double ans = (ALPHA/PI)*(LN(4*k_dash*k_dash/(m*m))-1);
	
	// Uncomment these lines to add the effects of (60) and (61) from [1]
	// double correction = (ALPHA/PI)*2*LN(k/k_dash);
	// correction += (ALPHA/PI)*LN( (1- COS(theta))/2);
	// ans += correction;

	return ans;
}


double simulateRadiator(double E0, TRandom *randomObject, double t, double E){

	double photonEnergy = E0 - E;
	double crossSection = energyThroughExternalRadiator(E,E0,t*1.0);
	
	return crossSection;

}


//Using (29) from reference [1] to account for vacuum polarization
double vacuumPolarization(double q2){
	double ans = 1/(3*PI);
	ans *= -5.0/3.0 + LN(-q2/(MASS_ELECTRON*MASS_ELECTRON));
	return ans;
}


// Using (33) from reference [1] to account for corrections given in Figure 3 of reference [1]
double hardCorrections(double q2){
	double ans = 2*ALPHA;
	ans *= -3*LN(-q2/(MASS_ELECTRON*MASS_ELECTRON))/(4*PI) + 1/PI - vacuumPolarization(q2);
	return ans;
}


#endif