/*
*	This ROOT macro performs simulations at a fixed angle. The value of the angle 
*	is read from EXACT_DETECTOR_ANGLE ( defined in constants.C ). The results of the
*	simulation are saved in a file. The location of the file is read from the variable
*	DESTINATION_FILE ( defined in constants.C )
*
*/
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

/*
*		This method for called for each event. 'e1' denotes incident electron energy. 'e3' denotes the scattered electron energy. 
*		This method produces a bremmstrahulg photon using the equivalent radiator approximation.
*		Energy given to photon is stored in the variable 'photonEnergy' ( passed by reference ).
*		'crossSection' is also modified/multiplied by the probability of this event to happen. It is also passed by reference.
*/
void scatteringWithRadiators( double e1, double e3, double theta, TRandom*randomObject, double & crossSection, double & photonEnergy ){
 	
 	bool incomingRadiator = true;
 	double e3_elastic = calculateScatteredElectronEnergy(e1,theta);
 	double deltaE = e3_elastic - e3;

 	// next line decides with 50-50% chance if photon is emitted from initial state radiation or the final state radiation
 	if ( randomObject->Uniform(-1.0,1.0) > 0.0 ) incomingRadiator = false;

	if ( incomingRadiator ){
		double lambda_e_tilda = lambda_incident(e1, e3, theta);
		double e1_required = incidentEnergyFromScatteredEnergy(e3,theta);
		photonEnergy = e1 - e1_required;
		crossSection *= simulateRadiator(e1, randomObject, lambda_e_tilda, e1_required);
		crossSection *= rosenbluth_formula( theta, e1_required);
	}


	if ( !incomingRadiator ){
		crossSection *= rosenbluth_formula( theta, e1);
		double lambda_e_dash_tilda = lambda_scattered(e1, e3, theta);
		photonEnergy = e3_elastic - e3;
		crossSection *= simulateRadiator(e3_elastic, randomObject, lambda_e_dash_tilda, e3);

	}

	crossSection *= DELTAEMAX;

}

/*
* Entry point to the program.
* Call this ROOT macro to start simulation.
* 'count' denotes the number of events.
*/
void simulateScattering(double incidentElectronEnergy, long count){

	TRandom *randomObject = new TRandom(0);
	double theta = -1.0;
	double phi = -1.0;
	double crossSection = -1.0;
	double rosenbluthCrossSection = -1.0;

	// Histogram to plot the energies of bremsstrahlung photons. See readme at Github.
	// TH1D *photonEnergyHistogram = new TH1D("E_g","E_g",nbins_E_g,xmin_E_g,xmax_E_g);
	// TH1D *missingEnergyHistogram = new TH1D("E_m","E_m",nbins_E_g,xmin_E_g,xmax_E_g);

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
		


		// Start with the normalized cross section per event and change it to microbarns. 
		// 'crossSection' will be modified by reference in equivalent radiators.
		crossSection = (2*PI)/(count);
		crossSection *= GeV2ub;
		rosenbluthCrossSection = crossSection;
		theta = EXACT_DETECTOR_ANGLE;

		// Calculate the 'expected' energy of elastically scattered electron because we use it to calculate radiator width
		double e3_el = calculateScatteredElectronEnergy(incidentElectronEnergy, theta);
		double photonEnergy = -1;
		double e3 = randomObject->Uniform(e3_el - DELTAEMAX,e3_el);
		double deltaE = e3_el - e3;

		scatteringWithRadiators(incidentElectronEnergy, e3, theta, randomObject, crossSection, photonEnergy);
		double q2 = calculateq2(incidentElectronEnergy, e3_el, theta);
		crossSection = (1-hardCorrections(q2))*crossSection;


		if ( photonEnergy < 0){
			cout<<" Bremstrrulung energy not calculated properly";
			exit(1);
		}

		phi = randomObject->Uniform(0,2*PI);


		//----------   DETECTOR WORK BEGIN ---------
		// Detectors need angles in degrees. So do conversion here.
		double thetaDegrees = rad2degree(theta);
		double phiDegrees = rad2degree(phi);

		rosenbluthCrossSection *= rosenbluth_formula( theta, incidentElectronEnergy);
		if ( electronDetector->detect(thetaDegrees, phiDegrees, e3, crossSection) ) {
			totalCrossSection += crossSection;
			totalRosenbluthCrossSection += rosenbluthCrossSection;
			events_needed--;

			// photonEnergyHistogram->Fill((photonEnergy)*1e-6, crossSection); 
			// missingEnergyHistogram->Fill((deltaE)*1e-6, crossSection); 
		}

	}


	// Save the histogram of photon energies.
	// TFile *file = new TFile(DESTINATION_FILE,"UPDATE");
	// photonEnergyHistogram->GetYaxis()->SetTitle("cross-section(ub)");
	// photonEnergyHistogram->GetXaxis()->SetTitle("Photon Energy(MeV)");
	// photonEnergyHistogram->Write("E_g");
	// missingEnergyHistogram->GetYaxis()->SetTitle("cross-section(ub)");
	// missingEnergyHistogram->GetXaxis()->SetTitle("delta E (MeV)");
	// missingEnergyHistogram->Write("E_m");
	// file->Close();

   electronDetector->showHistogram();

	cout<<"Successfully generated "<<count<<" events."<<endl;
	cout<<"Total cross section was ( in ub ) "<<totalCrossSection<<endl;
	cout<<"Total Rosenbluth cross section was ( in ub ) "<<totalRosenbluthCrossSection<<endl;
	cout<<"Delta is : "<<(totalCrossSection/totalRosenbluthCrossSection)-1.0<<endl;
	cout<<"Histograms have been written in file: "<<DESTINATION_FILE<<endl;

}





