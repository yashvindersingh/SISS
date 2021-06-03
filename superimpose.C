/*
*  Use this macro to compare two histograms
*/
#include "helper.C"
#include "constants.C"

void superimpose(){

    TString eseppFileName = "../debug2/e17.root";
    TString ermFileName = "../debug2/e16.root";
    TString title = "Title of the graph";

    TFile *file;
	TFile *esepp_file = new TFile(eseppFileName,"READ");
    TFile *erm_file = new TFile(ermFileName,"READ");

    TH1D *eseppProjection = (TH1D*)esepp_file->Get("E_l;1");
    TH1D *ermProjection = (TH1D*)erm_file->Get("E_l;1");
    TH1D *hnew = (TH1D*)myProjection->Clone("hnew");


    TCanvas *superimposeCanvas = new TCanvas("superimpose","superimpose",800,800);
    superimposeCanvas->cd();
    gPad->SetLogy();
    eseppProjection->SetLineColor(kGreen);
    eseppProjection->SetNameTitle("Comparison","Comparison");
    eseppProjection->GetYaxis()->SetRangeUser(10e-11,1);
    eseppProjection->Draw();
    ermProjection->SetLineColor(kRed);
    ermProjection->Draw("SAME");

    hnew->Divide(eseppProjection);
    hnew->SetLineColor(kBlue);
    hnew->Draw("SAME");

}