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

	double scatteredElectronEnergy = calculateScatteredElectronEnergy2(incidentElectronEnergy, theta);
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
// Energy give to photon is stored in the variable 'photonEnergy' ( passed by reference )
// 'crossSection' is also modified/multiplied by the probability of this event to happen, also passed by reference
ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* 
 scatteringWithRadiators( double e1, double e3, double theta, TRandom*randomObject, double & crossSection, double & photonEnergy ){
 	
 	bool incomingRadiator = true;
 	double e3_elastic = calculateScatteredElectronEnergy2(e1,theta);
 	double deltaE = e3_elastic - e3;
 	// decide with 50-50% chance if photon is emitted from initial state radiation or the final state radiation
 	if ( randomObject->Uniform(-1.0,1.0) > 0.0 ) incomingRadiator = false;


 	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* vectors;

	if ( incomingRadiator ){
		//cout<<"->->->->->->->->->Incoming radiator: "<<endl;
		double lambda_e_tilda = lambda_incident(e1, e3_elastic, theta);
		//double lambda_e_tilda = lambda_e(incidentElectronEnergy);
		double e1_required = incidentEnergyFromScatteredEnergy(e3,theta);
		photonEnergy = e1 - e1_required;
		if( photonEnergy < 0 ){
			cout<<"WARNING: Just change lost energy from "<<photonEnergy<<" to ZERO."<<endl;
			photonEnergy = 0;
			e1_required = e1;
		}

		crossSection *= simulateRadiator(e1, randomObject, lambda_e_tilda, photonEnergy);
		crossSection *= rosenbluth_formula( theta, e1_required);
		vectors = elastic4Vectors(e1_required, theta);
		//cout<<"Deviation is: "<<scatteredElectronEnergy - vectors[2].E()<<endl;
	}


	if ( !incomingRadiator ){
		//cout<<"Outgoing radiator->->->->->->->: "<<endl;
		crossSection *= rosenbluth_formula( theta, e1);
		vectors = elastic4Vectors(e1, theta);
		double lambda_e_dash_tilda = lambda_scattered(e1, e3_elastic, theta);
		//double lambda_e_dash_tilda = lambda_e(scatteredElectronEnergy);
		photonEnergy = e3_elastic - e3;
		crossSection *= simulateRadiator(e3_elastic, randomObject, lambda_e_dash_tilda, photonEnergy);

	}

	vectors[2].SetE(e3);
	crossSection *= DELTAEMAX;


	return vectors;

}


// entry point to the program
// call this root macro to start simulation
void simulateScatteringRange(double incidentElectronEnergy, int count){

	TRandom *randomObject = new TRandom(0);
	double originalIncidentEnergy = incidentElectronEnergy; // taking backup of 'incidentElectronEnergy' because it will be change in each event
	double theta = -1.0;
	double phi = -1.0;
	double crossSection = -1.0;
	double rosenbluthCrossSection = -1.0;

	// Histogram to plot the energies of bremsstrahlung photons. See readme at Github.
	TH1D *photonEnergyHistogram = new TH1D("E_g","E_g",nbins_E_g,xmin_E_g,xmax_E_g);
	TH1D *missingEnergyHistogram = new TH1D("E_m","E_m",nbins_E_g,xmin_E_g,xmax_E_g);


	// Placing two detectors to detect electrons and protons.
	Detector *electronDetector = new Detector(rad2degree(MINIMUM_DETECTOR_ANGLE),rad2degree(MAXIMUM_DETECTOR_ANGLE),0,360,originalIncidentEnergy);
	Detector *protonDetector = new Detector(rad2degree(MINIMUM_DETECTOR_ANGLE),rad2degree(MAXIMUM_DETECTOR_ANGLE),0,360,originalIncidentEnergy,"proton");

	double totalCrossSection = 0 ;
	double totalRosenbluthCrossSection = 0 ;

	//cout.precision(17);
	//cout<<incidentElectronEnergy<<endl;

	for( int events_needed = count ; events_needed > 0 ; ){


		// start with the normalized cross section per event and change it to microbarns. It will be 
		// modified by reference in equivalent radiators
		crossSection = (2*PI)*( COS(MINIMUM_DETECTOR_ANGLE) - COS(MAXIMUM_DETECTOR_ANGLE))/(count);
		crossSection *= GeV2ub;
		rosenbluthCrossSection = crossSection;
		incidentElectronEnergy = originalIncidentEnergy;

		// generate random theta between the required range
		double cosMax = COS( MINIMUM_DETECTOR_ANGLE - 0.05*ONE_DEGREE);
		double cosMin = COS( MAXIMUM_DETECTOR_ANGLE + 0.05*ONE_DEGREE);
		theta = ACOS( randomObject->Uniform(cosMin,cosMax) );

		// Calculate the 'expected' energy of elastically scattered electron because we use it to calculate radiator width
		double e3_el = calculateScatteredElectronEnergy2(incidentElectronEnergy, theta);
		double photonEnergy = -1;
		//double e3 = randomObject->Uniform(e3_el - DELTAEMAX,e3_el);
		double e3 = randomObject->Uniform(E3MAX - DELTAEMAX,E3MAX);
		double deltaE = e3_el - e3;

		if( deltaE < 0 ){
			crossSection = 0;
			if ( electronDetector->detect(rad2degree(theta), 1.57, e3, crossSection) ) {
				events_needed--;
				photonEnergyHistogram->Fill(0.0, crossSection); 
				missingEnergyHistogram->Fill(0.0, crossSection); 
			}
			continue;
		}

		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* vectors = 
		scatteringWithRadiators(incidentElectronEnergy, e3, theta, randomObject, crossSection, photonEnergy);

		//We have the value of four momentum transfer now. We use it to account for virtual corrections to cross section.
		double q2 = calculateq2(incidentElectronEnergy, theta);
		crossSection = (1-hardCorrections(q2))*crossSection;


		if ( photonEnergy < 0){
			cout<<" xxxxxxxxxx- Bremstrrulung energy not calculated properly -xxxxxxxxx"<<endl;
			cout<<e3_el<<" :max"<<endl;
			cout<<e3<<" :chosen"<<endl;
			cout<<photonEnergy<<" :photon"<<endl;
			exit(1);
		}


	//	phi = randomObject->Uniform(0,2*PI);
		phi = 1.57;

	//	double protonTheta = TMath::ATan(vectors[3].Py()/vectors[3].Pz());

	//	if ( protonTheta < 0 )
	//		protonTheta = protonTheta + PI; // Since ATan gives
		// output from -pi/2 to pi/2 I am shifting the angles b/w ( -90,0) to (90,180)






		//----------   DETECTOR WORK BEGIN ---------

		// TODO:  the angles for proton are NOT CORRECT. If I calculate correct proton theta, I will
		// have to place the proton detector at correct place. For now I am just passing the electron theta
		// to proton detectors. Proton phi is also ad-hoc.

		// Detectors need angles in degrees. So do conversion here.
		double thetaDegrees = rad2degree(theta);
		double phiDegrees = rad2degree(phi);
//		double thetaProton = rad2degree(protonTheta);
		// phi for proton can be found from phi electron by rotating around 180 degrees
//		double phiProton = (phiDegrees + 180);
//		if ( phiProton >= 360) phiProton -= 360;  // I need to be more careful here. Phi might have sign problem. 
		// But since the scattering is phi symmmetric, I am leaving it for now.

		rosenbluthCrossSection *= rosenbluth_formula( theta, originalIncidentEnergy);
		if ( electronDetector->detect(thetaDegrees, phiDegrees, vectors[2].E(), crossSection) ) {
		//if ( electronDetector->detect(thetaDegrees, phiDegrees, vectors[2].E(), 1) ) {
			totalCrossSection += crossSection;
			totalRosenbluthCrossSection += rosenbluthCrossSection;
			events_needed--;
			//photonEnergyHistogram->Fill((photonEnergy)*1e-6, crossSection); 
			//missingEnergyHistogram->Fill((deltaE)*1e-6, crossSection); 
		}

		// using thetaDegrees here as a hack, otherwise proton detector won't work
		//protonDetector->detect(thetaDegrees, phiDegrees, vectors[3].E(), crossSection);
		//---------------DETECTOR WORK END----------------------------



	}





	// Save the histogram of photon energies.
    TFile *file = new TFile(DESTINATION_FILE,"UPDATE");
    photonEnergyHistogram->GetYaxis()->SetTitle("cross-section(ub)");
    photonEnergyHistogram->GetXaxis()->SetTitle("Photon Energy(MeV)");
	photonEnergyHistogram->Write("E_g");

	missingEnergyHistogram->GetYaxis()->SetTitle("cross-section(ub)");
    missingEnergyHistogram->GetXaxis()->SetTitle("delta E (MeV)");
	missingEnergyHistogram->Write("E_m");

    file->Close();

    // Save the histogram of electron and proton energies
    electronDetector->showHistogram();
    protonDetector->showHistogram();



	cout<<"Successfully generated "<<count<<" events."<<endl;
	
	cout<<"Total cross section was ( in ub ) "<<totalCrossSection<<endl;
	cout<<"Total Rosenbluth cross section was ( in ub ) "<<totalRosenbluthCrossSection<<endl;
	cout<<"Delta is : "<<(totalCrossSection/totalRosenbluthCrossSection)-1.0<<endl;
	cout<<"Histograms have been written in file: "<<DESTINATION_FILE<<endl;

}





