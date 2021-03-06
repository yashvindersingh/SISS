#ifndef CONSTANTS_H
#define CONSTANTS_H


// These contants are used throughout the program. Natural units are used here
#define ONE_DEGREE TMath::Pi()/180
#define EULER 2.71828
#define PI TMath::Pi()
#define LIGHT_SPEED 1
#define ATOMIC_NUMBER 1
#define HBAR 1
#define ALPHA 1.0/137
#define MASS_PROTON 938e6			// every mass and energy is in eV
#define MASS_ELECTRON 0.511e6
#define DETECTOR_THRESHOLD 0e6      // The mininum energy of of an electron detectable by the detector.
									// In TDIS, this threshold is 56 MeV/c ( line no. 870 )

#define GeV2mb 0.389379e18			// From wikipedia. Multiply by these factors to convert cross-section
#define GeV2ub 389.389379e18        // from natural units to milli or micro-barns.



#define LOG 	TMath::Log10
#define LN 	    TMath::Log
#define TAN 	TMath::Tan
#define COS 	TMath::Cos
#define ACOS 	TMath::ACos
#define SIN 	TMath::Sin
#define POWER	TMath::Power
#define SQRT	TMath::Sqrt
#define GAMMA	TMath::Gamma
#define DILOG	TMath::DiLog
#define ABS	TMath::Abs



// Below are the parameters for each run. Check readme page on Github for more information.
#define MINIMUM_DETECTOR_ANGLE 15.95*TMath::Pi()/180
#define MAXIMUM_DETECTOR_ANGLE 16.05*TMath::Pi()/180
#define EXACT_DETECTOR_ANGLE 16.0*TMath::Pi()/180
#define DELTAEMAX 500e6
#define DESTINATION_FILE "../testing/t102.root"


// Number of bins for E_g,E_l,E_p,theta_p histgrams
#define nbins_E_g 500
#define nbins_E_l 500
#define nbins_E_p 500
#define nbins_theta_p 120

// Lower value of x-axis range for E_g,E_l,E_p,theta_p histgrams ( in MeV and radians )
#define xmin_E_g 0
#define xmin_E_l 1800
#define xmin_E_p 0
#define xmin_theta_p 0

// Upper value of x-axis range for E_g,E_l,E_p,theta_p histgrams ( in MeV and radians )
#define xmax_E_g 30
#define xmax_E_l 2300
#define xmax_E_p 2500
#define xmax_theta_p 1.6

#endif