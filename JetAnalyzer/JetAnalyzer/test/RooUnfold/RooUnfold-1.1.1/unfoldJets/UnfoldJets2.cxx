#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include "TRandom.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TLatex.h"
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
#include "../src/RooUnfoldResponse.h"
#include "../src/RooUnfoldBayes.h"
#include "../src/RooUnfoldBinByBin.h"
#include "../src/RooUnfoldSvd.h"
#include "../src/RooUnfoldTUnfold.h"
//#include "UnfoldJets.h"
#include "/uscms/home/tlibeiro/inclusivexsec/CMSSW_5_3_12_patch2/src/JetAnalyzer/JetAnalyzer/test/AnalyzerProcedures.h"
#endif
using namespace std;
void groomHist(TH1D *h, int kCyan);
void normalizeTH2D(TH2D* h);
void ScaleHistByBW(TH1* h);
TH1D* makeNLOSpectrum(char* infile,char* title);
void makeLOSpectrum(TH1F** h);
class DeltaRDistance;
void generateEvents(unsigned events, TF1* fres, TH1D* hnlo,
TH1D* hTrue, TH1D* hReco,TH2D* hMatrix,RooUnfoldResponse* resp=0);
//==============================================================================
// Main program when run stand-alone
//==============================================================================
void UnfoldJets2();
int main (int argc, char** argv) {
  UnfoldJets2();
  return 0;
}
//==============================================================================
// Unfolding
//==============================================================================
void UnfoldJets2() {
#ifdef __CINT__
	gSystem->Load("../libRooUnfold");
#endif
	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptTitle(kFALSE);
	//gROOT->SetBatch(kTRUE);
// Data Members
  TLegend            *lTrain, *lTest, *lErrors;
  TH1                *hTrain, *hTrainFake, *hTrue, *hMeas, *hReco, *hFake, *hRes, *hPulls;
  TH1D*              hTrainTrue[10];
  TH1D*              hTrainReco[10];
  TH1                *hUnfErr, *hToyErr, *hParmChi2, *hParmErr, *hParmRes, *hParmRms;
  TH2D               *hResmat, *hCorr, *hMeasCorr;
  TH2D               *hMat[10];
  RooUnfoldResponse* response;
  RooUnfold*         unfold;



  const string indir("/uscms/home/tlibeiro/inclusivexsec/CMSSW_5_3_12_patch2/src/JetAnalyzer/JetAnalyzer/test");
  const char* nloDir("/uscms_data/d3/tlibeiro/inclusivexsec/Generators/generateIncJets/inputXsecs");
  const char* fitDir("./");
   sprintf(name,"%s/ResolutionFits.root",fitDir);
   TFile* resFile = new TFile(name,"read");
   TF1*  fres[6];//resfits
	const unsigned eventsGenerate(1e5);

  response= new RooUnfoldResponse ("response","Test");
	cout << "==================================== TRAIN ====================================" << endl;
  RooUnfoldResponse* responseEta[nybins-1];
  TH1F* hlo[6];//lo pythia spectrum
  TH1D* hnlo[6];//nlo spectrum
  TH1D* hnlot[6][2];//nlo spectrum for comparison
  TH1D* hTrainTruet[6][2];// for comparison
  TH1D* hTrainRecot[6][2];// for comparison
  TH2D* hMatt[6][2];// for comparison
////get pythia spectrum 
  makeLOSpectrum(hlo);

  //////get nlo spectrum and fits=============
  ///Initialize Histos and Response class======
  for(unsigned i(0);i<6;++i) {
    char infile[500]; char title[500];
    sprintf(infile,"%s/nnpdf2.1y%i.cppInput",nloDir,i);
    sprintf(title ,"nloy%i",i);
    hnlo[i] = makeNLOSpectrum(infile,title);
    //comparing unfolding with different spectra
    sprintf(infile,"%s/ct10y%i.cppInput",nloDir,i);
    sprintf(title ,"nlot%dy%i",0,i);
    hnlot[i][0] = makeNLOSpectrum(infile,title);
    //
    sprintf(infile,"%s/abm11y%i.cppInput",nloDir,i);
    sprintf(title ,"abm11%dy%i",0,i);
    hnlot[i][1] = makeNLOSpectrum(infile,title);
    //
    sprintf(name,"fresfit_eta%d",i);
    fres[i] = (TF1*)resFile->Get(name);
		//
		sprintf(name,"traintrue_eta%d",i);
		hTrainTrue[i]= new TH1D (name,name,n1x[i],&x1[i][0]);
		hTrainTrue[i]->SetLineColor(kBlue);
		//
		sprintf(name,"trainmeas_eta%d",i);
		hTrainReco[i]= new TH1D (name,name,n1x[i],&x1[i][0]);
		hTrainReco[i]->SetLineColor(kRed);
		sprintf(name,"resmat_eta%d",i);
		hMat[i]= new TH2D (name,name,n1x[i],&x1[i][0],n1x[i],&x1[i][0]);
		//
    sprintf(name,"response_eta%d",i);
		responseEta[i] = new RooUnfoldResponse (name,name);
		responseEta[i]->Setup(hTrainReco[i],hTrainTrue[i]);
    /////for comparision
		sprintf(name,"resmatt%d_eta%d",0,i);
		hMatt[i][0]= new TH2D (name,name,n1x[i],&x1[i][0],n1x[i],&x1[i][0]);
		//
		sprintf(name,"traintruet%d_eta%d",0,i);
		hTrainTruet[i][0]= new TH1D (name,name,n1x[i],&x1[i][0]);
		hTrainTruet[i][0]->SetLineColor(kBlue);
		//
		sprintf(name,"trainmeast%d_eta%d",0,i);
		hTrainRecot[i][0]= new TH1D (name,name,n1x[i],&x1[i][0]);
		hTrainRecot[i][0]->SetLineColor(kRed);
    //
		sprintf(name,"resmatt%d_eta%d",1,i);
		hMatt[i][1]= new TH2D (name,name,n1x[i],&x1[i][0],n1x[i],&x1[i][0]);
		//
		sprintf(name,"traintruet%d_eta%d",1,i);
		hTrainTruet[i][1]= new TH1D (name,name,n1x[i],&x1[i][0]);
		hTrainTruet[i][1]->SetLineColor(kBlue);
		//
		sprintf(name,"trainmeast%d_eta%d",1,i);
		hTrainRecot[i][1]= new TH1D (name,name,n1x[i],&x1[i][0]);
		hTrainRecot[i][1]->SetLineColor(kRed);
  }
  ///generate events
	for(unsigned neta(0);neta<6;++neta) {
   generateEvents(eventsGenerate,fres[neta],hnlo[neta],hTrainTrue[neta],hTrainReco[neta],
		hMat[neta],responseEta[neta]);
   generateEvents(eventsGenerate,fres[neta],hnlot[neta][0],hTrainTruet[neta][0],
		hTrainRecot[neta][0],hMatt[neta][0]);
   generateEvents(eventsGenerate,fres[neta],hnlot[neta][1],hTrainTruet[neta][1],
		hTrainRecot[neta][1],hMatt[neta][1]);
	}

//	response->Mresponse().Print();
//	TMatrixD mat = response->Mresponse();
//	mat.Draw("colz");


	cout << "==================================== UNFOLD ===================================" << endl;
	TFile* fDat = new TFile((indir+"/AnalyzerIncHistosData276TeV.root").c_str(),"read"); 
	TFile* fMC  = new TFile((indir+"/AnalyzerIncHistosMC.root").c_str(),"read"); 
	TH1D* hclo[6];// unfolding closure
	TH1D* htom[6];// unfolded true / mesaured spectrum
	TH1D* hclot[6][2];// unfolding closure for comparision
	for(unsigned i(0);i<1;++i) {

		TH1D*  hDat = hTrainReco[i];
		TH1D*  hGen = (TH1D*)hTrainTrue[i];
//    ScaleHistByBW(hDat);
		RooUnfoldBayes   unfold1  (responseEta[i], hDat, 3);    // OR
//		RooUnfoldSvd     unfold2  (responseEta[i], hDat, 4);   // OR
//		RooUnfoldBinByBin unfold3 (responseEta[i], hDat);

		TH1D* hUnf1= (TH1D*) unfold1.Hreco();
//		TH1D* hUnf2= (TH1D*) unfold2.Hreco();
//		TH1D* hUnf3= (TH1D*) unfold3.Hreco();
		//for comparision 
//    RooUnfoldBinByBin   unfoldt1  (responseEta[i], hTrainRecot[i][0]);   
//		TH1D* hUnft1= (TH1D*) unfoldt1.Hreco();
//    RooUnfoldBinByBin   unfoldt2  (responseEta[i], hTrainRecot[i][1]);   
//		TH1D* hUnft2= (TH1D*) unfoldt2.Hreco();
   // 
//    RooUnfoldSvd   unfoldt1  (responseEta[i], hTrainRecot[i][0], 4);   
//		TH1D* hUnft1= (TH1D*) unfoldt1.Hreco();
//    RooUnfoldSvd   unfoldt2  (responseEta[i], hTrainRecot[i][1], 4);   
//		TH1D* hUnft2= (TH1D*) unfoldt2.Hreco();
//   // 
    RooUnfoldBayes   unfoldt1  (responseEta[i], hTrainRecot[i][0], 3);   
		TH1D* hUnft1= (TH1D*) unfoldt1.Hreco();
    RooUnfoldBayes   unfoldt2  (responseEta[i], hTrainRecot[i][1], 3);   
		TH1D* hUnft2= (TH1D*) unfoldt2.Hreco();

//		unfold1.PrintTable(cout,hGen);
//		unfold2.PrintTable(cout,hGen);
		cout<<"Drawing histos eta low edge "<<ybins[i]<<endl;
		TCanvas* c1 = new TCanvas();
		hUnf1->SetLineColor(kRed);
		hUnf1->SetMarkerStyle(21);
		hUnf1->SetMarkerSize(0.3);
		hUnf1->GetYaxis()->SetRangeUser(1e-8,1e8);
		hUnf1->Draw("hist");
		gPad->SetLogx(1);gPad->SetLogy(1);
		hDat->SetLineColor(kGreen);
		hDat->Draw("same");
		hGen->SetLineColor(kBlack);
		hGen->Draw("e same");
		//		hTrain->SetLineColor(kBlue);
		//		hTrain->Draw("same");
		//	hUnf2->SetLineColor(kMagenta);
		//  hUnf2->Draw("same");
    TLegend* leg = new TLegend(0.711567,0.653719,0.930511,0.965385);
    leg->AddEntry(hDat,"Data","pel");
    leg->AddEntry(hGen,"NLO ","pel");
    leg->AddEntry(hUnf1,"Unfolded Spectrum","pel");
    leg->Draw("same");
    //closure
	  TCanvas* ptclos = new TCanvas("ptclos","ptclos",600,600); //mc data closure
		ptclos->cd(); gPad->SetLogx(1);
    sprintf(name,"hclosure_eta%.2f",ybins[i+1]);
    hclo[i] = (TH1D*)hGen->Clone(name);
		hclo[i]->Reset(); 
    hclo[i]->Sumw2();
    hclo[i]->Divide(hUnf1,hGen);
    hclo[i]->GetYaxis()->SetRangeUser(0.9,1.1);
    hclo[i]->GetXaxis()->SetRangeUser(18,1000);
    hclo[i]->Draw();
	  TCanvas* ptclos1 = new TCanvas("ptclos1","ptclos1",600,600); //mc data closure
		ptclos1->cd(); gPad->SetLogx(1);
    sprintf(name,"hclosure1_eta%.2f",ybins[i+1]);
    hclot[i][0] = (TH1D*)hGen->Clone(name);
		hclot[i][0]->Reset(); 
    hclot[i][0]->Divide(hUnft1,hTrainTruet[i][0]);
    hclot[i][0]->GetYaxis()->SetRangeUser(0.9,1.1);
    hclot[i][0]->GetXaxis()->SetRangeUser(18,1000);
    hclot[i][0]->SetMarkerStyle(21);
    hclot[i][0]->Draw("");
    TCanvas* ptclos2 = new TCanvas("ptclos2","ptclos2",600,600); //mc data closure
		ptclos2->cd(); gPad->SetLogx(1);
    sprintf(name,"hclosure2_eta%.2f",ybins[i+1]);
    hclot[i][1] = (TH1D*)hGen->Clone(name);
		hclot[i][1]->Reset(); 
    hclot[i][1]->Divide(hUnft2,hTrainTruet[i][1]);
    hclot[i][1]->GetYaxis()->SetRangeUser(0.9,1.1);
    hclot[i][1]->GetXaxis()->SetRangeUser(18,1000);
    hclot[i][1]->SetMarkerStyle(21);
    hclot[i][1]->Draw("");
   ///testing pythia spectrum
	TCanvas* ct2 = new TCanvas();
		hUnft2->SetLineColor(kRed);
		hUnft2->SetMarkerStyle(21);
		hUnft2->SetMarkerSize(0.3);
		hUnft2->GetYaxis()->SetRangeUser(1e-8,1e8);
		hUnft2->Draw("hist");
		gPad->SetLogx(1);gPad->SetLogy(1);
		hTrainRecot[i][1]->SetLineColor(kGreen);
		hTrainRecot[i][1]->Draw("same");
		hTrainTruet[i][1]->SetLineColor(kBlack);
		hTrainTruet[i][1]->Draw("e same");
    TLegend* legt2 = new TLegend(0.711567,0.653719,0.930511,0.965385);
    legt2->AddEntry(hTrainRecot[i][1],"ABM11 Data","pel");
    legt2->AddEntry(hTrainTruet[i][1],"ABM11 Gen ","pel");
    legt2->AddEntry(hUnft2,"Unfolded Spectrum","pel");
    legt2->Draw("same");

//    cout<<fixed<<setprecision(3)<<scientific;
//		for(unsigned bin(0);bin<hGen->GetNbinsX();++bin)
//			cout<<" bin "<<bin
//			<<" gen "<<hGen->GetBinContent(bin)
//			<<" unfold "<<hUnf1->GetBinContent(bin)
//			<<" gen err "<<hGen->GetBinError(bin)
//			<<" unfold err  "<<hUnf1->GetBinError(bin)
//      <<" clos "<<hclo[i]->GetBinContent(bin)
//			<<" clos err "<<hclo[i]->GetBinError(bin)
//      <<endl;

    sprintf(name,"%.2f<|y|<%.2f",ybins[i],ybins[i+1]);
    TLatex latex;latex.SetTextSize(0.04);latex.SetTextFont(42);
    latex.DrawTextNDC(.7,.8,name);
//   //true over measured
//	  TCanvas* ctom = new TCanvas("ctom","ctom",600,600); //mc data closure
//		ctom->cd(); gPad->SetLogx(1);
//    sprintf(name,"htrueovmeasured_eta%.2f",ybins[i+1]);
//    htom[i] = (TH1D*)hGen->Clone(name);
//		htom[i]->Reset(); 
//    htom[i]->Divide(hUnf1,hDat);
//    htom[i]->GetYaxis()->SetRangeUser(0,5);
//    htom[i]->Draw();
//    sprintf(name,"%.2f<|y|<%.2f",ybins[i],ybins[i+1]);
//    latex.DrawTextNDC(.7,.8,name);
//
//
		TCanvas* c2 = new TCanvas();
    TMatrixD mat = responseEta[i]->Mresponse();
    mat.Draw("colz");
    TCanvas* c3 = new TCanvas();
    c3->SetLogy();c3->SetLogx();
    hMat[i]->Draw("colz");
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
		hTrainReco[i]->Write();
    hMat[i]->Write();
	}
	const TNamed* objs[] = { //unfold, 
//		hTrain, hTrainFake ,// hTrue, hMeas, hFake,
		// hReco, hRes,
//		hResmat//, hCorr, hMeasCorr, hUnfErr, hToyErr, ntChi2,
		// hParmChi2, hParmErr, hParmRes, hParmRms 
	};
	const unsigned numobj(sizeof(objs)/sizeof(objs[0])); 
	for (unsigned i(0);i<numobj;++i) {
		if (objs[i]) outfile->WriteTObject (objs[i], objs[i]->GetName());
	}
	outfile->Close();
}

void groomHist(TH1D *h, int kCyan ) {
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

//Scale bin content by bin width 
void ScaleHistByBW(TH1* h) {
	const unsigned numBins(h->GetNbinsX());
	for(unsigned bin(1);bin<=numBins;++bin)
	{ 
		const double binWidth  = h->GetBinWidth(bin);
		const double newBinVal = h->GetBinContent(bin)/(binWidth);
		const double newBinErr = h->GetBinError(bin)/(binWidth);
		h->SetBinContent(bin,newBinVal);
		h->SetBinError  (bin,newBinErr);
		if(false)
			cout<<"-------------------------------"
				<<" entries "<<h->GetBinContent(bin)<<" +- "
				<<h->GetBinError(bin)
				<<endl;
	}
}


void generateEvents(unsigned events, TF1* fres, TH1D* hnlo,
TH1D* hTrue, TH1D* hReco,TH2D* hMatrix,RooUnfoldResponse* resp)
{
    hTrue->Sumw2(); hReco->Sumw2();
	  const unsigned nbins(hnlo->GetNbinsX());
   // cout<<nbins<<' '<<events<<endl;
		for(unsigned bin(1);bin<=nbins;++bin)
		{
			double xsec  = hnlo->GetBinContent(bin);
			double binlo = hnlo->GetBinLowEdge(bin);
			double binhi = hnlo->GetBinWidth(bin)+binlo;
			for(unsigned i(0);i<events;++i)
			{
			 double genpt = gRandom->Uniform(binlo,binhi);
       double resol = gRandom->Gaus(1.0,fres->Eval(genpt));
       double recpt = genpt*resol;
//       cout<<genpt<<" "<<recpt<<" "<<xsec<<endl;
			 hTrue->Fill(genpt,xsec/events);
			 hReco->Fill(recpt,xsec/events);
			 hMatrix->Fill(recpt,genpt,xsec/events);
			 if(resp)
				 resp->Fill(recpt,genpt,xsec/events);
			}
		}
//		hTrue->Scale(1.0/events);
//		hReco->Scale(1.0/events);
		normalizeTH2D(hMatrix);
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

//==============================================================================
// Helper Functions And Classes
//==============================================================================
void colorIt(TH1D *hMidA, int kCyan ){
	hMidA->SetMarkerColor(kCyan);
	hMidA->SetLineColor(kCyan);
	hMidA->SetMarkerStyle(22);//24
	hMidA->SetMarkerSize(1.0);
	hMidA->SetFillColor(kCyan);
	hMidA->SetFillStyle(0);
	hMidA->GetXaxis()->SetLabelSize(0.05);
	hMidA->GetXaxis()->SetTitleSize(0.06);
	hMidA->GetYaxis()->SetTitleSize(0.05);
	hMidA->GetYaxis()->SetTitleOffset(0.9);
	hMidA->GetYaxis()->SetLabelSize(0.05);
	hMidA->SetTitleFont(42, "XYZ");
	hMidA->SetLabelFont(42, "XYZ");
};

class DeltaRDistance
{
  public:
    template<class JetType>
    double operator()(const JetType& jet1, const JetType& jet2) const
    {
      const double deltaEta = jet1.eta - jet2.eta;
      double deltaPhi = jet1.phi - jet2.phi;
      if (deltaPhi < -M_PI)
        deltaPhi += 2.0*M_PI;
      if (deltaPhi > M_PI)
        deltaPhi -= 2.0*M_PI;
      const double distance = sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
      return (distance);
    }
};

	TH1D* makeNLOSpectrum(char* infile,char* title) {
		ifstream rstream;
		rstream.open(infile);
		vector<double> ptMinV, xsecV; double xmax(0);
		if(rstream.is_open())
			while(!rstream.eof()) {
				int nBin(0);
				double ptMin(0), ptMax(0), xsec(0), knlo(0);
				rstream>>nBin>>ptMin>>ptMax>>xsec>>knlo;
				//cout<<" "<<nBin<<" "<<ptMin<<" "<<ptMax<<" "<<xsec<<" "<<knlo<<endl;
				if(ptMin && ptMax && xsec) {
					ptMinV.push_back(ptMin);
					xsecV .push_back(xsec); xmax = ptMax;
				}
			}
		//cout<<"num bins "<<ptMinV.size()<<endl;
		const unsigned nBins = ptMinV.size();
		double x[nBins+1];
		for(unsigned i(0);i<nBins;++i)
			x[i]  = ptMinV[i];
		x[nBins]=xmax;
		TH1D* h = new TH1D(title,title,nBins,x);
		//fill histogram
		for(unsigned i(0);i<nBins;++i) {
			h->SetBinContent(i+1,xsecV[i]);
			h->SetBinError(i+1,0);
		}
		return h;
	};

void makeLOSpectrum(TH1F** h) {
	char* indir = "/uscms/home/tlibeiro/inclusivexsec/CMSSW_5_3_12_patch2/src/JetAnalyzer/JetAnalyzer/test";
	sprintf(name,"%s/AnalyzerIncHistosMC.root",indir);
	TFile* f1 = new TFile(name,"read");
	for(unsigned i(0);i<nybins-1;++i)
	{
		sprintf(name,"hgenUnf_eta%i",i);
		h[i] = (TH1F*)f1->Get(name);
		cout<<h[i]->GetName()<<' '<<h[i]<<endl;
		ScaleHistByBW(h[i]);
   h[i]->Scale(1/lumi);
	}
};
