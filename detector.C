#ifndef DETECTOR_H
#define DETECTOR_H

#include "helper.C" 

class Detector {     
  private:           

    double thetaMin;
    double thetaMax;
    double phiMin;
    double phiMax;
    double minEnergy;
    double maxEnergy;
    double incidentEnergy;

    TH1D *totalEnergyHistogram;
    TH1D *thetaHistogram;

    TString title;


    //TODO: should we make species an enum
    TString species;  

  public:

  	Detector(){
  	}

    Detector(double thetaMin, double thetaMax, double phiMin, double phiMax , double incidentEnergy, TString species = "electron"){
      initialize(thetaMin,thetaMax,phiMin,phiMax, incidentEnergy, species);
    }

  	void initialize(double thetaMin, double thetaMax, double phiMin, double phiMax, double incidentEnergy, TString species){
  		this->thetaMin = thetaMin;
  		this->thetaMax = thetaMax;
      this->phiMax = phiMax;
      this->phiMin = phiMin;
      this->incidentEnergy = incidentEnergy;
      this->species = species;
      if( species == "proton"){
        totalEnergyHistogram = new TH1D("E_p","E_p",nbins_E_p,xmin_E_p,xmax_E_p);
        thetaHistogram = new TH1D("theta_p","theta_p",nbins_theta_p,xmin_theta_p,xmax_theta_p);

      } else{
        totalEnergyHistogram = new TH1D("E_l","E_l",nbins_E_l,xmin_E_l,xmax_E_l);
        thetaHistogram = new TH1D("theta_l","theta_l",nbins_theta_p,xmin_theta_p,xmax_theta_p);
      
      }
      
  	}


  	bool detect(double theta, double phi, double energy, double cross_section){

  		if ( theta <= thetaMax && theta >= thetaMin && phi <= phiMax && phi >= phiMin && energy >= DETECTOR_THRESHOLD ){

        totalEnergyHistogram->Fill((energy)*1e-6, cross_section);  
        thetaHistogram->Fill(theta*PI/180,cross_section);
        return true;

  		}else{
        return false;
      } 

  	}


    void showHistogram(){
      
      decorateHistogram();
      
      TFile *file = new TFile(DESTINATION_FILE,"UPDATE");
      if( species == "proton"){

        totalEnergyHistogram->Write("E_p");  
        thetaHistogram->Write("theta_p");
      } 
      else{
        totalEnergyHistogram->Write("E_l");  
        thetaHistogram->Write("theta_l");
      } 

      file->Close();
      
    }

    void decorateHistogram(){

      totalEnergyHistogram->GetYaxis()->SetTitle("cross-section(ub)");
      totalEnergyHistogram->GetXaxis()->SetTitle("Total Energy(MeV)");

      thetaHistogram->GetYaxis()->SetTitle("cross-section(ub)");
      thetaHistogram->GetXaxis()->SetTitle("Theta(radians)");
      if( species == "electron"){
        totalEnergyHistogram->SetNameTitle("E_l","E_l");
        thetaHistogram->SetNameTitle("theta_l","theta_l");
      }
        
      else{
        totalEnergyHistogram->SetNameTitle("E_p","E_p");
        thetaHistogram->SetNameTitle("theta_p","theta_p");
      }


    }

};


#endif
