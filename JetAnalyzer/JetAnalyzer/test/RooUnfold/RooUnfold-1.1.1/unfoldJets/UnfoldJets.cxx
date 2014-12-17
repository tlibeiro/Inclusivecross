#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include <sstream>
#include <string>
#include "TRandom.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TArrayD.h"
#include "TMatrixD.h"
#include "TChain.h"
#include "TTree.h"
#include "TF1.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"
#include "UnfoldJets.h"
#include "matchOneToOne.hh"
#include "/uscms/home/tlibeiro/inclusivexsec/CMSSW_5_3_12_patch2/src/JetAnalyzer/JetAnalyzer/test/AnalyzerProcedures.h"
#endif
using namespace std;
using namespace UnfoldJetNS;
void groomHist(TH1D *h, int kCyan);
void normalizeTH2D(TH2D* h);
//==============================================================================
// Main program when run stand-alone
//==============================================================================
void Unfold();
int main (int argc, char** argv) {
  Unfold();
  return 0;
}
//==============================================================================
// Unfolding
//==============================================================================
void Unfold() {
#ifdef __CINT__
	gSystem->Load("libRooUnfold");
#endif

	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptTitle(kFALSE);
	//gROOT->SetBatch(kTRUE);

  const string indir("/uscms/home/tlibeiro/inclusivexsec/CMSSW_5_3_12_patch2/src/JetAnalyzer/JetAnalyzer/test");
  response= new RooUnfoldResponse ("response","Test");
	cout << "==================================== TRAIN ====================================" << endl;
  TChain* chain = new TChain("UnfoldJetTree");
  chain->Add((indir+"/AnalyzerHistos_allpt[0-9]_Unfold.root").c_str());
//	TFile* infile = new TFile("../../AnalyzerHistos_allpt6_Unfold.root","read");
	TTree* unfTree = chain;
	///Read Ntuple =============================
	const unsigned maxJets(15);
	int njets(0), ngen(0);
	float jtpt [maxJets];  float gjtpt [maxJets]; double mcWeight;
	float jteta[maxJets];  float gjteta[maxJets];
	float jtphi[maxJets];  float gjtphi[maxJets];
	float jtm  [maxJets];  float gjty  [maxJets];
	unfTree->SetBranchAddress("jtpt" , jtpt );  unfTree->SetBranchAddress("gjtpt" , gjtpt ); 
	unfTree->SetBranchAddress("jtphi", jtphi);  unfTree->SetBranchAddress("gjtphi", gjtphi);
	unfTree->SetBranchAddress("jteta", jteta);  unfTree->SetBranchAddress("gjteta", gjteta);
	unfTree->SetBranchAddress("jtm"  , jtm);    unfTree->SetBranchAddress("gjty"  , gjty);
	unfTree->SetBranchAddress("njets", &njets); unfTree->SetBranchAddress("ngen"  , &ngen);
  unfTree->SetBranchAddress("mcweight",&mcWeight);
	///End Read Ntuple ==========================
  RooUnfoldResponse* responseEta[nybins-1];
  ///Initialize Histos and Response class======
	for(unsigned i(0);i<nybins-1;++i) {
		sprintf(name,"traintrue_eta%d",i);
		hTrainTrue[i]= new TH1D (name,name,nx[i],&x[i][0]);
		hTrainTrue[i]->SetLineColor(kBlue);
		sprintf(name,"trainmeas_eta%d",i);
		hTrain= new TH1D (name,name,nx[i],&x[i][0]);
		hTrain->SetLineColor(kRed);
		sprintf(name,"resmat_eta%d",i);
		hResmat= new TH2D (name,name,nx[i],&x[i][0],nx[i],&x[i][0]);
		sprintf(name,"trainfake_eta%d",i);
		hTrainFake= new TH1D (name,name,nx[i],&x[i][0]);
		hTrainFake->SetLineColor(93);
		sprintf(name,"response_eta%d",i);
		responseEta[i] = new RooUnfoldResponse (name,name);
		responseEta[i]->Setup(hTrain,hTrainTrue[i]);
	}

	unsigned nentries(unfTree->GetEntries());
	//nentries=1000;
	cout<<"Reading Unfold tree: Entries "<<nentries<<endl;
	vector<UnfoldJetNS::Jet> jets, genjets;
	vector<int> matchedGen, matchedJet;
	const double dRMatch(0.5);//Gen-Jet match 
	DeltaRDistance d; 
	////Run on the training Tree===================
	for(unsigned entry(0);entry<nentries;++entry)
	{
		jets.clear(); genjets.clear();
		unfTree->GetEntry(entry);
		for(unsigned njet(0);njet<njets;++njet)
			jets.push_back(UnfoldJetNS::Jet(jtpt[njet], jteta[njet],
						jtphi[njet], jtm[njet]));
		for(unsigned ngjet(0);ngjet<ngen;++ngjet)
			genjets.push_back(UnfoldJetNS::Jet(gjtpt[ngjet], gjteta[ngjet],
						gjtphi[ngjet],0));
		matchOneToOne(jets,genjets,DeltaRDistance() ,&matchedGen,1000.0);  
		matchOneToOne(genjets,jets,DeltaRDistance() ,&matchedJet,1000.0);  
		if(false)
			cout<<"Entry "<<entry
				<< " == "<<" Gen "<<ngen<<" Jet "<<njets 
				<< endl;

		for(unsigned ngjet(0);ngjet<ngen;++ngjet) ///Matched Jets to Gen 
		{
			int match = matchedJet[ngjet];
			for(unsigned i(0);i<nybins-1;++i) 
				if(abs(genjets[ngjet].eta)>ybins[i] &&
						abs(genjets[ngjet].eta)<=ybins[i+1]) {
					hTrainTrue[i]->Fill(genjets[ngjet].pt,mcWeight);
					if(match>-1) 
					{
						const double dr = d(genjets[ngjet],jets[match]);
						if(dr<=dRMatch) {
							//	hTrain  ->Fill(jets[match].pt,mcWeight);
							//	hResmat ->Fill(jets[match].pt,genjets[ngjet].pt,mcWeight);
							responseEta[i]->Fill(jets[match].pt,genjets[ngjet].pt,mcWeight);
						}
						else { 
							responseEta[i]->Miss(genjets[ngjet].pt,mcWeight);
							responseEta[i]->Fake(jets[match].pt   ,mcWeight);
							//	hTrainFake->Fill(jets[match].pt ,mcWeight);
						}
					}
					else {
						responseEta[i]->Miss(genjets[ngjet].pt,mcWeight);
					}
				}
		}
		if(!(entry%200000))
			cout<<"Processing Event "<<(float)entry<<endl;
	}//entry loop

//	response->Mresponse().Print();
//	TMatrixD mat = response->Mresponse();
//	mat.Draw("colz");


	cout << "==================================== UNFOLD ===================================" << endl;
	TFile* fDat = new TFile((indir+"/AnalyzerHistosData.root").c_str(),"read"); 
	for(unsigned i(0);i<nybins-1;++i) {
		ostringstream oss;
		oss<<"hptUnf_eta_"<<ybins[i]<<"_"<<ybins[i+1];
		TH1F*  hDat = (TH1F*)fDat->Get(oss.str().c_str());
		TH1D*  hGen = (TH1D*)hTrainTrue[i];
		RooUnfoldBayes   unfold1  (responseEta[i], hDat, 3);    // OR
		RooUnfoldSvd     unfold2  (responseEta[i], hDat, 4);   // OR
		RooUnfoldBinByBin unfold3 (responseEta[i], hDat);

		TH1D* hUnf1= (TH1D*) unfold1.Hreco();
		TH1D* hUnf2= (TH1D*) unfold2.Hreco();
		TH1D* hUnf3= (TH1D*) unfold3.Hreco();

		unfold1.PrintTable(cout,hGen);
		unfold2.PrintTable(cout,hGen);
		cout<<"Drawing histos eta low edge "<<ybins[i]<<endl;
		TCanvas* c1 = new TCanvas();
		hUnf1->SetLineColor(kRed);
		hUnf1->GetYaxis()->SetRangeUser(1e-3,1e8);
		hUnf1->Draw();
		gPad->SetLogx(1);gPad->SetLogy(1);
		hDat->SetLineColor(kGreen);
		hDat->Draw("same");
		hGen->SetLineColor(kBlack);
		hGen->Draw("same");
		//		hTrain->SetLineColor(kBlue);
		//		hTrain->Draw("same");
		//	hUnf2->SetLineColor(kMagenta);
		//  hUnf2->Draw("same");
    TLegend* leg = new TLegend(0.711567,0.653719,0.930511,0.965385);
    leg->AddEntry(hDat,"Data","pel");
    leg->AddEntry(hGen,"MC Gen","pel");
    leg->AddEntry(hUnf1,"Unfolded Spectrum","pel");
    leg->Draw("same");

		TCanvas* c2 = new TCanvas();
    TMatrixD mat = responseEta[i]->Mresponse();
    mat.Draw("colz");
    	}

	//cout << "==================================== DRAW ===================================" << endl;
	//
	//	TH1D *Mont1= (TH1D*)hDat->Clone("Mont");
	//	TH1D *Mont2= (TH1D*)hUnf1->Clone("Mont");
	////	TH1D *Mont3= (TH1D*)hUnf2->Clone("Mont");
	//	TH1D *Mont4= (TH1D*)hUnf3->Clone("Mont");
	//
	//	leg1 = new TLegend(0.711567,0.653719,0.930511,0.965385);
	//	leg1->SetBorderSize(0);
	//	leg1->SetFillColor(kWhite);
	//	leg1->SetTextSize(0.03);
	//	leg1->AddEntry(hGen, "Truth","pel");
	//	leg1->AddEntry(hDat, "Measured","pel");
	//	leg1->AddEntry(hUnf1, "Unfold Bayes","pel");
	////	leg1->AddEntry(hUnf2, "Unfold SVD","pel");
	//	leg1->AddEntry(hUnf3, "Unfold Bin-by-Bin","pel");    
	//
	//	leg2 = new TLegend(0.711567,0.653719,0.930511,0.965385);
	//	leg2->SetBorderSize(0);
	//	leg2->SetFillColor(kWhite);
	//	leg2->SetTextSize(0.03);
	//	leg2->AddEntry(Mont1, "Measured/Truth","p");
	//	leg2->AddEntry(Mont2, "Unfold Bayes/Truth","p");
	////	leg2->AddEntry(Mont3, "Unfold SVD/Truth","p");
	//	leg2->AddEntry(Mont4, "Unfold Bin-by-Bin/Truth","p");
	//
	//	Mont1->Divide(hDat,hGen,1.0,1.0,"b");
	//	Mont2->Divide(hUnf1,hGen,1.0,1.0,"b");
	////	Mont3->Divide(hUnf2,hGen,1.0,1.0,"b");
	//	Mont4->Divide(hUnf3,hGen,1.0,1.0,"b");
	//
	//  TCanvas* cUnfold = new TCanvas("cUnfold","cUnfold",800,600);
	//  hUnf1->Draw("eX0C");
	////  hUnf2->Draw("same");
	//  hUnf3->Draw("same eX0C");
	//  leg1->Draw("same");
	//
	//  TCanvas* cRatio = new TCanvas("cRatio","cRatio",800,600);
	//  Mont1->Draw("peX0C");
	//  Mont2->Draw("peX0C same");
	////  Mont3->Draw("peX0C same");
	//  Mont4->Draw("peX0C same");
	//  leg2->Draw("same");
	//
	//  TCanvas *c2 = new TCanvas("c2","",600,600);
	//	char name[100];
	//	normalizeTH2D(hResmat);
	//	sprintf(name,"Reco %s","xtitle");
	//	hResmat->SetXTitle(name);
	//	sprintf(name,"Gen %s","y title");
	//	hResmat->GetYaxis()->SetTitleOffset(1.5);
	//	hResmat->SetYTitle(name);
	//	gPad->SetLogx(1); gPad->SetLogy(1);
	//	hResmat->Draw("colz");
	//


	///Write Histos and Response object to file
	TFile* outfile = new TFile("UnfoldJets.root","recreate");
  for(unsigned i(0);i<nybins-1;++i)
	{
		sprintf(name,"response_eta%i",i);
		outfile->WriteTObject(responseEta[i],name);
		hTrainTrue[i]->Write();
	}
	const TNamed* objs[] = { //unfold, 
		hTrain, hTrainFake ,// hTrue, hMeas, hFake,
		// hReco, hRes,
		hResmat//, hCorr, hMeasCorr, hUnfErr, hToyErr, ntChi2,
		// hParmChi2, hParmErr, hParmRes, hParmRms 
	};
	const unsigned numobj(sizeof(objs)/sizeof(objs[0])); 
	for (unsigned i(0);i<numobj;++i) {
		if (objs[i]) outfile->WriteTObject (objs[i], objs[i]->GetName());
	}
	outfile->Close();
}

void groomHist(TH1D *h, int kCyan ){
	h->SetMarkerColor(kCyan);
	h->SetLineColor(kCyan);
	h->SetMarkerStyle(22);//24
	h->SetMarkerSize(1.0);
	h->SetFillColor(kCyan);
	h->SetFillStyle(0);
	h->GetXaxis()->SetLabelSize(0.05);
	h->GetXaxis()->SetTitleSize(0.06);
	h->GetYaxis()->SetTitleSize(0.05);
	h->GetYaxis()->SetTitleOffset(0.9);
	h->GetYaxis()->SetLabelSize(0.05);
	h->SetTitleFont(42, "XYZ");
	h->SetLabelFont(42, "XYZ");
}

void normalizeTH2D(TH2D *h)
{
	int xbin(h->GetNbinsX()), ybin(h->GetNbinsY());
	for (int i(1); i <= ybin; i++){
		double sum(0.);
		for (int j(1); j <= xbin; j++){
			sum += h->GetBinContent(j,i);
		}
		for (int j(1); j <= xbin; j++){
			if (sum > 0){
				double val = h->GetBinContent(j,i)/sum;
				val=val*100;
				val=ceil(val);
				double nearest=val/100;
				h->SetBinContent(j,i,nearest);
			}
		}
	}
}

//		///Checking for fakes, jets not matched to any gens
//		for(unsigned njet(0) ;njet<njet ;++njet) ///Matched Check for jets not matched to gens
//		{ 
//			bool jetmatched(false);
//			for(unsigned ngjet(0);ngjet<ngen;++ngjet) ///Matched Check for jets not matched to gens
//			{  
//				int match = matchedJet[ngjet];
//				if(match>-1) 
//					if(njet==match)
//						jetmatched=true; 
//			}
//			if(!jetmatched)  {
//				response->Fake(jets[njet].pt);
//        hTrainFake->Fill(jets[njet].pt);
//			}
//		}


