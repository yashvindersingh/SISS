// This ROOT macro uses (5.1) equation in [5] for scattering at a fixed angle EXACT_DETECTOR_ANGLE
#include "helper.C"
#include "detector.C"

/*
*		Returns Mott cross section of electron scattered at angle theta radians. 
*		Implementation of equation (6.1) in [3].
*/
double mott_formula(double theta, double incidentElectronEnergy){
	
	int Z = ATOMIC_NUMBER;
	double alpha = ALPHA;
	double c = LIGHT_SPEED; 
	double hbar = HBAR;

	double scatteredElectronEnergy = calculateScatteredElectronEnergy(incidentElectronEnergy, theta);
	double ans = Z*Z*alpha*alpha*hbar*hbar*c*c/4;
	ans = ans*TMath::Power(TMath::Cos(theta/2),2);
	ans = ans/TMath::Power(TMath::Sin(theta/2),4);
	ans = ans*scatteredElectronEnergy/TMath::Power(incidentElectronEnergy,3);
	
	return ABS(ans);
	
}

/*
*		Implementation of equation (6.12) in [3].
*/
double dipole_form_factor(double Q2){
	double temp = 1 + Q2*LIGHT_SPEED*LIGHT_SPEED/(0.71e18);
	return TMath::Power(temp,-2);
}


/*
*		Returns rosenbluth cross section of electron scattered at theta radians.
*		Implementation of equation (6.10) in [3].
*/
double rosenbluth_formula(double theta, double incidentElectronEnergy){
	
	double mott = mott_formula(theta,incidentElectronEnergy);
	double q2 = calculateq2(incidentElectronEnergy, theta);
	double Q2 = -1*q2;
	double tau = Q2/(4*MASS_PROTON*MASS_PROTON*LIGHT_SPEED*LIGHT_SPEED);
	double dipole_electric = dipole_form_factor(Q2);
	double dipole_magnetic = 2.79*dipole_form_factor(Q2);
	double ans = mott*(
						(dipole_electric*dipole_electric + tau*dipole_magnetic*dipole_magnetic)/(1+tau) +
	 					2*tau*TMath::Power(dipole_magnetic*TMath::Tan(theta/2),2));
	return ABS(ans);

}

// This formula calculates the coefficient of Z^0 term in (5.1) equation in [5]
double Z0(double q2,double e1, double e3, double theta, double deltaE){
	double m = MASS_ELECTRON;
	double eta = e1/e3;
	double lqm = LN(-1.0*q2/(m*m));
	double temp = 13.0 * lqm/6;
	temp += -28.0/9;
	temp += -(lqm - 1 )*LN(4*e1*e3/POWER(2*eta*deltaE,2));
	temp += -0.5*LN(eta)*LN(eta);
	temp += DILOG( POWER( COS(theta/2),2 ));
	temp += - PI*PI/6;
	return temp*ALPHA/PI;
}

// This formula calculates the coefficient of Z^1 term in (5.1) equation in [5]
double Z1(double q2,double e1, double e3, double theta, double deltaE){
	double m = MASS_ELECTRON;
	double M = MASS_PROTON;
	double eta = e1/e3;
	double Z = 1;
	double x = x2(e1,e3,theta)/(4*M*M);
	double temp = -1*LN(eta)*LN(-1*q2*x/POWER(2.0*eta*deltaE,2));
	temp += DILOG( 1 - eta/x);
	temp += -DILOG(1.0-1.0/(eta*x));
	return temp*(2*ALPHA*Z/PI);
}

// This formula calculates the coefficient of Z^2 term in (5.1) equation in [5]
double Z2(double q2,double e1, double e3, double theta, double deltaE){
	double m = MASS_ELECTRON;
	double M = MASS_PROTON;
	double eta = e1/e3;
	double Z = 1;
	double x = x2(e1,e3,theta)/(4*M*M);
	double e4 = e1 + M - e3;
	double p4 = SQRT( e4*e4 - M*M);
	double rho2 = -q2 + 4*M*M;
	double temp = (e4/p4) * ( -0.5*LN(x)*LN(x) - LN(x)*LN(rho2/(M*M)) + LN(x) );
	temp -= ( (e4/p4)*LN(x) - 1 )*LN(M*M/POWER(2.0*eta*deltaE,2));
	temp += 1;
	temp += (e4/p4)*( -DILOG(1-1/(x*x)) + 2*DILOG(-1/x) + PI*PI/6 );
	return temp*(ALPHA*Z*Z/PI);
}

// This formula calculates the derivative of the coefficient of Z^0 term ( with respect to delta_E ) in (5.1) equation in [5]
double Z0_differential(double q2,double e1, double e3, double theta, double deltaE){
	double m = MASS_ELECTRON;
	double eta = e1/e3;
	double lqm = LN(-1.0*q2/(m*m));
	double temp = 2*ALPHA/PI;
	temp *= lqm - 1 ;
	temp *= 1.0/deltaE;
	return temp;
}

// This formula calculates the derivative of the coefficient of Z^2 term ( with respect to delta_E ) in (5.1) equation in [5]
double Z1_differential(double q2,double e1, double e3, double theta, double deltaE){
	double m = MASS_ELECTRON;
	double M = MASS_PROTON;
	double eta = e1/e3;
	double Z = 1;
	double temp = 4*ALPHA*Z/PI;
	temp *= LN(eta)/deltaE;
	return temp;
}

// This formula calculates the derivative of the coefficient of Z^2 term ( with respect to delta_E ) in (5.1) equation in [5]
double Z2_differential(double q2,double e1, double e3, double theta, double deltaE){
	double m = MASS_ELECTRON;
	double M = MASS_PROTON;
	double eta = e1/e3;
	double Z = 1;
	double x = x2(e1,e3,theta)/(4*M*M);
	double e4 = e1 + M - e3;
	double p4 = SQRT( e4*e4 - M*M);
	double rho2 = -q2 + 4*M*M;
	double temp = 2*ALPHA*Z*Z/PI;
	temp *= (e4/p4)*LN(x) - 1;
	temp *= 1.0/deltaE;
	return temp;
}

// This formula calculates the value of delta in (5.1) equation in [5]
double delta_differential(double q2,double e1, double e3, double theta, double deltaE){

	double ans = Z0_differential(q2, e1, e3, theta, deltaE) 
					+ 
				Z1_differential(q2, e1, e3, theta, deltaE)
					+
				Z2_differential(q2, e1, e3, theta, deltaE);
	return ans;
}

// This formula calculates the derivative of delta ( with respect to delta_E ) in (5.1) equation in [5]
double delta(double q2,double e1, double e3, double theta, double deltaE){
	double ans = Z0(q2, e1, e3, theta, deltaE) 
					+ 
				Z1(q2, e1, e3, theta, deltaE)
					+
				Z2(q2, e1, e3, theta, deltaE);
	return ans;
}

void maximon(double incidentElectronEnergy, long count){

	TRandom *randomObject = new TRandom(0);
	double theta = -1.0;
	double phi = -1.0;
	double crossSection = -1.0;
	double rosenbluthCrossSection = -1.0;

	// Histogram to plot the energies of bremsstrahlung photons. See readme at Github.
	TH1D *photonEnergyHistogram = new TH1D("E_m","E_m",nbins_E_g,xmin_E_g,xmax_E_g);

	Detector *electronDetector = new Detector(rad2degree(MINIMUM_DETECTOR_ANGLE),rad2degree(MAXIMUM_DETECTOR_ANGLE),0,360,incidentElectronEnergy);

	double totalCrossSection = 0 ;
	double totalRosenbluthCrossSection = 0 ;

	for( long events_needed = count ; events_needed > 0 ; ){

		for ( long i = 10; i >= 1; --i){
			if ( events_needed == (i*count)/10){
				cout<<"Generated "<<(10-i)*10<<" percent events ... "<<endl;
				continue;
			}
		}


		// start with the normalized cross section per event and change it to microbarns. It will be 
		// modified by reference in equivalent radiators
		crossSection = (2*PI)/(count);
		crossSection *= GeV2ub;
		rosenbluthCrossSection = crossSection;

		theta = EXACT_DETECTOR_ANGLE;

		double e1 = incidentElectronEnergy;
		double e3el = calculateScatteredElectronEnergy(e1, theta);
		double e3 = randomObject->Uniform(e3el - DELTAEMAX,e3el);
		double deltaE = e3el - e3;
		double lostEnergy = deltaE;
		double q2 = calculateq2(e1,e3el,theta);
		double deltaCorrection = delta_differential(q2,e1, e3el, theta, deltaE);
		deltaCorrection *= POWER(EULER,delta(q2,e1, e3el, theta, deltaE));
		deltaCorrection *= DELTAEMAX;


		phi = randomObject->Uniform(0,2*PI);

		//----------   DETECTOR WORK BEGIN ---------

		// Detectors need angles in degrees. So do conversion here.
		double thetaDegrees = rad2degree(theta);
		double phiDegrees = rad2degree(phi);

		rosenbluthCrossSection *= rosenbluth_formula( theta, e1);
		crossSection = rosenbluthCrossSection*(deltaCorrection);
		if ( electronDetector->detect(thetaDegrees, phiDegrees, e3, crossSection) ) {
			totalCrossSection += crossSection;
			totalRosenbluthCrossSection += rosenbluthCrossSection;
			events_needed--;
			photonEnergyHistogram->Fill((lostEnergy)*1e-6, crossSection); 
		}

		//---------------DETECTOR WORK END----------------------------



	}


	// Save the histogram of photon energies.
    TFile *file = new TFile(DESTINATION_FILE,"UPDATE");
    photonEnergyHistogram->GetYaxis()->SetTitle("cross-section(ub)");
    photonEnergyHistogram->GetXaxis()->SetTitle("Total Energy(MeV)");
	photonEnergyHistogram->Write("E_m");
    file->Close();

    electronDetector->showHistogram();

	cout<<"Successfully generated "<<count<<" events."<<endl;
	cout<<"Total cross section was ( in ub ) "<<totalCrossSection<<endl;
	cout<<"Total Rosenbluth cross section was ( in ub ) "<<totalRosenbluthCrossSection<<endl;
	cout<<"Delta is : "<<(totalCrossSection/totalRosenbluthCrossSection)<<endl;
	cout<<"Histograms have been written in file: "<<DESTINATION_FILE<<endl;

}





