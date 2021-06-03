/*
*  This is ROOT macro that converts the .dat file produced from ESEPP into proper histogram.
*  Somehow, the ROOT files directly produced by ESEPP were not working with TBrowser().
*
*/

void ntuple(){
   ifstream in;

   // the name of the dat file from ESEPP goes here
   in.open("../esepp/output/e17_e-.dat");


   // total cross section from ESEPP should be entered here
   double total_cross_section = 3.9354e-05 ;
   double total_events = 100000000;
   double cpe = total_cross_section/total_events;



   Float_t E_l,theta_l,phi_l,E_p,theta_p,phi_p,E_g,theta_g,phi_g;
   Int_t nlines = 0;
   // the histogram for the ESEPP file is stored in the following location
   TFile *f = new TFile("../debug2/e17.root","RECREATE");

   /*
   TH1D *E_g_histogram = new TH1D("E_g_histogram","E_g distribution",500,0,2500);
   TH1D *E_l_histogram = new TH1D("E_l_histogram","E_l distribution",500,0,2500);
   TH1D *E_p_histogram = new TH1D("E_p_histogram","E_p distribution",500,0,2500);
   */
   
   //TH1D *E_g_histogram = new TH1D("E_g","E_g distribution",500,0,2500);
   TH1D *E_l_histogram = new TH1D("E_l","E_l distribution",500,1400,1900);
   //TH1D *E_p_histogram = new TH1D("E_p","E_p distribution",500,0,2500);

   TH1D *Theta_p_histogram = new TH1D("Theta_p_histogram","Theta_p distribution",120,0,1.6);
   TH1D *Theta_l_histogram = new TH1D("Theta_l_histogram","Theta_l distribution",120,0,1.6);

   TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","E_l");

   while (1) {

      for ( long i = 0; i <= 10; ++i){
         if ( nlines == (i*total_events)/10){
            cout<<"Parsed "<<i*10<<" percent events ... "<<endl;
            continue;
         }
      }


      in >> E_l >> theta_l >> phi_l >> E_p >> theta_p >> phi_p >> E_g >> theta_g >> phi_g;
      if (!in.good()) break;
      if (nlines < 5) cout<<E_l<<":"<<theta_l<<":"<<phi_l<<":"<<E_p<<":"<<theta_p<<":"<<phi_p<<":"<<E_g<<":"<<theta_g<<":"<<phi_g<<endl;
      //E_g_histogram->Fill(E_g,cpe);
      //E_g_histogram->Scale(cpe);
      E_l_histogram->Fill(E_l,cpe);
      //E_l_histogram->Scale(cpe);
      //E_p_histogram->Fill(E_p,cpe);
      //E_p_histogram->Scale(cpe);
      //Theta_p_histogram->Fill(theta_p,cpe);
      //Theta_p_histogram->Scale(cpe);
      //Theta_l_histogram->Fill(theta_l,cpe);

      ntuple->Fill(E_l);
      nlines++;
   }
   cout<<"Found lines: "<<nlines<<endl;

   in.close();

   f->Write();

}