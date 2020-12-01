// See readme file on Github for references
#include "helper.C"
#include "detector.C"


// Implementing equation (6.1) in [1]
// returns mott cross section of electron scattered at 'theta' 
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

// See (6.12) in [1]
double dipole_form_factor(double Q2){
	double temp = 1 + Q2*LIGHT_SPEED*LIGHT_SPEED/(0.71e18);
	return TMath::Power(temp,-2);
}


// Implementing equation (6.10) in [1]
// returns rosenbluth cross section of electron scattered at 'theta' 
double rosenbluth_formula(double theta, double incidentElectronEnergy){
	
	double mott = mott_formula(theta,incidentElectronEnergy);
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* vectors = elastic4Vectors(incidentElectronEnergy, theta);
	double q2 = vectors[4].Dot(vectors[4]);
	double Q2 = -1*q2;
	double tau = Q2/(4*MASS_PROTON*MASS_PROTON*LIGHT_SPEED*LIGHT_SPEED);
	double dipole_electric = dipole_form_factor(Q2);
	double dipole_magnetic = 2.79*dipole_form_factor(Q2);
	double ans = mott*(
						(dipole_electric*dipole_electric + tau*dipole_magnetic*dipole_magnetic)/(1+tau) +
	 					2*tau*TMath::Power(dipole_magnetic*TMath::Tan(theta/2),2));
	return ABS(ans);

}



// This method for called for each event. Based on 'incidentElectronEnergy', 'scatteredElectronEnergy' and 'theta', this method
// produces a bremmstrahulg photon using the modifiedequivalent radiator approximation
// It returns the 4-vectors of incident and final particles after accounting for scattering and bremsstrahlung ( per event )
// Energy give to photon is stored in the variable 'lostEnergy' ( passed by reference )
// 'crossSection' is also modified/multiplied by the probability of this event to happen, also passed by reference
ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* 
 scatteringWithRadiators( double incidentElectronEnergy, double scatteredElectronEnergy, double theta, TRandom*randomObject, double & crossSection, double & lostEnergy ){
 	
 	bool incomingRadiator = true;
 	// decide with 50-50% chance if photon is emitted from initial state radiation or the final state radiation
 	if ( randomObject->Uniform(-1.0,1.0) > 0.0 ) incomingRadiator = false;


 	// calculate the parameters needed for calculation of "equivalent length" of external radiators using equation 58 and 62 of [3]
 	// All variable names are defined as in [3]
	double k = SQRT(incidentElectronEnergy*incidentElectronEnergy - MASS_ELECTRON*MASS_ELECTRON);
	double k_dash = SQRT(scatteredElectronEnergy*scatteredElectronEnergy - MASS_ELECTRON*MASS_ELECTRON);

	double lambda_e_tilda = lambda_e(incidentElectronEnergy) + extra_lambda(k, k_dash, theta);
	double lambda_e_dash_tilda = lambda_e(scatteredElectronEnergy) + extra_lambda(k, k_dash, theta);


	if ( incomingRadiator )
		crossSection *= simulateIncomingRadiator(incidentElectronEnergy, randomObject, lambda_e_tilda, lostEnergy);


	// The 4-vectors involved in elastic scattering are calculated here because some parameters are needed for 'simulateOutgoingRadiator'
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* vectors = 
						elastic4Vectors(incidentElectronEnergy, theta);

	// do actual elastic scattering using rosenbluth formula. This step and the previous step are quite redundant. But 'elastic4Vectors'
	// does not give the crossSection. 
	crossSection *= rosenbluth_formula( theta, incidentElectronEnergy);


	if ( !incomingRadiator )
		crossSection *= simulateOutgoingRadiator(vectors, randomObject, lambda_e_dash_tilda, lostEnergy);


	// 'lostEnergy' and 'crossSection' have already been modified by reference.
	// We just need to output the calculated 4-vectors now.
	return vectors;

}

// entry point to the program
// call this root macro to start simulation
void simulateScattering(double incidentElectronEnergy, int count){

	TRandom *randomObject = new TRandom(0);
	double originalIncidentEnergy = incidentElectronEnergy; // taking backup of 'incidentElectronEnergy' because it will be change in each event
	double theta = -1.0;
	double phi = -1.0;
	double crossSection = -1.0;

	// Histogram to plot the energies of bremsstrahlung photons. See readme at Github.
	TH1D *photonEnergyHistogram = new TH1D("E_g","E_g",nbins_E_g,xmin_E_g,xmax_E_g);


	// Placing two detectors to detect electrons and protons.
	Detector *electronDetector = new Detector(rad2degree(MINIMUM_DETECTOR_ANGLE),rad2degree(MAXIMUM_DETECTOR_ANGLE),0,360,originalIncidentEnergy);
	Detector *protonDetector = new Detector(rad2degree(MINIMUM_DETECTOR_ANGLE),rad2degree(MAXIMUM_DETECTOR_ANGLE),0,360,originalIncidentEnergy,"proton");

	double totalCrossSection = 0 ;

	for( int events_needed = count ; events_needed > 0 ; ){


		// start with the normalized cross section per event and change it to microbarns. It will be 
		// modified by reference in equivalent radiators
		crossSection = (2*PI)*( COS(MINIMUM_DETECTOR_ANGLE) - COS(MAXIMUM_DETECTOR_ANGLE))/(count);
		crossSection *= GeV2ub;
		incidentElectronEnergy = originalIncidentEnergy;


		// generate random theta between the required range
		double cosMax = COS( MINIMUM_DETECTOR_ANGLE - 2*ONE_DEGREE);
		double cosMin = COS( MAXIMUM_DETECTOR_ANGLE + 2*ONE_DEGREE);
		theta = ACOS( randomObject->Uniform(cosMin,cosMax) );

		// Calculate the 'expected' energy of elastically scattered electron because we use it to calculate radiator width
		double scatteredElectronEnergy = calculateScatteredElectronEnergy(incidentElectronEnergy, theta);

		// define the energy of bremsstrahung photon. It will be changed by reference in radiator methods
		double lostEnergy = 0;

		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* vectors = 
		scatteringWithRadiators(incidentElectronEnergy, scatteredElectronEnergy, theta, randomObject, crossSection, lostEnergy);

		if ( lostEnergy < 0){
			cout<<" Bremstrrulung energy not calculated properly";
			exit(1);
		}


		phi = randomObject->Uniform(0,2*PI);

		double protonTheta = TMath::ATan(vectors[3].Py()/vectors[3].Pz());

		if ( protonTheta < 0 )
			protonTheta = protonTheta + PI; // Since ATan gives
		// output from -pi/2 to pi/2 I am shifting the angles b/w ( -90,0) to (90,180)






		//----------   DETECTOR WORK BEGIN ---------

		// TODO:  the angles for proton are NOT CORRECT. If I calculate correct proton theta, I will
		// have to place the proton detector at correct place. For now I am just passing the electron theta
		// to proton detectors. Proton phi is also ad-hoc.

		// Detectors need angles in degrees. So do conversion here.
		double thetaDegrees = rad2degree(theta);
		double phiDegrees = rad2degree(phi);
		double thetaProton = rad2degree(protonTheta);
		// phi for proton can be found from phi electron by rotating around 180 degrees
		double phiProton = (phiDegrees + 180);
		if ( phiProton >= 360) phiProton -= 360;  // I need to be more careful here. Phi might have sign problem. 
		// But since the scattering is phi symmmetric, I am leaving it for now.

		
		if ( electronDetector->detect(thetaDegrees, phiDegrees, vectors[2].E(), crossSection) ) {
			totalCrossSection += crossSection;
			events_needed--;
			photonEnergyHistogram->Fill((lostEnergy)*1e-6, crossSection); 
		}

		// using thetaDegrees here as a hack, otherwise proton detector won't work
		protonDetector->detect(thetaDegrees, phiDegrees, vectors[3].E(), crossSection);
		//---------------DETECTOR WORK END----------------------------



	}


	// Save the histogram of photon energies.
    TFile *file = new TFile(DESTINATION_FILE,"UPDATE");
    photonEnergyHistogram->GetYaxis()->SetTitle("cross-section(ub)");
    photonEnergyHistogram->GetXaxis()->SetTitle("Total Energy(MeV)");
	photonEnergyHistogram->Write("E_g");
    file->Close();

    // Save the histogram of electron and proton energies
    electronDetector->showHistogram();
    protonDetector->showHistogram();



	cout<<"Successfully generated "<<count<<" events."<<endl;
	
	cout<<"Total cross section was ( in ub ) "<<totalCrossSection<<endl;
	cout<<"Histograms have been written in file: "<<DESTINATION_FILE<<endl;

}





