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
#include "TROOT.h"
#include "TF1.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "../src/RooUnfoldResponse.h"
#include "../src/RooUnfoldBayes.h"
#include "../src/RooUnfoldBinByBin.h"
#include "../src/RooUnfoldSvd.h"
#include "../src/RooUnfoldTUnfold.h"
//#include "UnfoldJets.h"
#include "../../../AnalyzerProcedures.h"
#endif
using namespace std;
void applyNP(TH1* h,TF1* f);
void normalizeTH2D(TH2D* h);
void normalizeTH2DX(TH2D* h);
void ScaleHistByBW(TH1* h);
TH1D* makeNLOSpectrum(char* infile,char* title);
void makeLOSpectrum(TH1F** h);
class DeltaRDistance;
void drawSpec(TCanvas* c,TH1* hreco, TH1* hgen, TH1* hunf, TLegend* leg,bool saveplots,char* label=0);
void drawClosure(TCanvas* c,TH1* hunf, TH1* hgen, TH1* hclo,unsigned eta,bool saveplots,char* label);
void drawClosure(TCanvas* c,unsigned eta,bool saveplots,
								 TH1* hunf1, TH1* hgen1, TH1* hclo1, char* label1,
                 TH1* hunf2, TH1* hgen2, TH1* hclo2, char* label2);
void drawClosure2(TCanvas* c,unsigned eta, bool saveplots,
								 TH1* hunf1, TH1* hgen1, TH1* hclo1, char* label1,
                 TH1* hunf2, TH1* hgen2, TH1* hclo2, char* label2);
void drawSpecCl(TCanvas* c,TH1* hdt, TH1* hmc,TH1* hclo,unsigned eta);
void drawMatrix(TH2* h,unsigned eta);
void spFillHisto(TH1* hin, TH1* hout, bool verbose=false);

void generateEvents(unsigned events, TF1* fres, TH1* hnlo,
TH1D* hTrue, TH1D* hReco,TH2D* hMatrix,double dataRes,
double lumi,RooUnfoldResponse* resp=0);
//==============================================================================
// Main program when run stand-alone
//==============================================================================
void UnfoldJets3();
int main (int argc, char** argv) {
  gSystem->Load("../libRooUnfold");
  UnfoldJets3();
  return 0;
}
//==============================================================================
// Unfolding
//==============================================================================
void UnfoldJets3() {
	gSystem->Load("../libRooUnfold");
	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptTitle(kFALSE);
// Data Members
  TH1D*              hTrainTrue[10];
  TH1D*              hTrainReco[10];
  TH2D               *hMat[10];
  RooUnfoldResponse* response;
  RooUnfold*         unfold;

  bool saveplots(true);
  if(saveplots) 	gROOT->SetBatch(kTRUE);
  gROOT->SetBatch(kFALSE);
  const string indir("/uscms/home/tlibeiro/inclusivexsec/CMSSW_5_3_12_patch2/src/JetAnalyzer/JetAnalyzer/test");
  const char* nloDir("/uscms_data/d3/tlibeiro/inclusivexsec/Generators/generateIncJets/inputXsecs");
  const char* fitDir("./");
   sprintf(name,"%s/ResolutionFits.root",fitDir);
   TFile* resFile = new TFile(name,"read");
   TF1*  fres[6];//resfits
   double cData[6]= {1.03,1.03,1.06,1.05,1.11,1.10};//jet res data/mc
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
  //for data
  RooUnfoldResponse* responseEtaDt[nybins-1]; //response eta with data res
  TH1D* hTrainTrueDt[6];// for comparison
  TH1D* hTrainRecoDt[6];// for comparison
  TH2D* hMatDt[6];// for comparison
////get pythia spectrum 
  makeLOSpectrum(hlo);

  ///Initialize Histos and Response class//////////////////////////////////////
  for(unsigned i(0);i<6;++i) {
    char infile[500]; char title[500];
//    sprintf(infile,"%s/nnpdf2.1y%i.cppInput",nloDir,i);
//    sprintf(title ,"nloy%i",i);
    sprintf(infile,"%s/ct10y%i.cppInput",nloDir,i);
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
		/////training histos
//    unsigned nbins(findBinsFromArray(18,maxPt[i],x1[i],n1x[i])+5);
    unsigned ngen(n1x[i]), nrec(nx[i]);
    const double *xgen = x1[i]; //x bins gen
    const double *xrec = x[i];// x bins reco
		sprintf(name,"traintrue_eta%d",i);
		hTrainTrue[i]= new TH1D (name,name,ngen,xgen);
		hTrainTrue[i]->SetLineColor(kBlue);
//    cout<<" histo bins "<<hTrainTrue[i]->GetNbinsX()<<" "<<hTrainTrue[i]->GetBinLowEdge(1)<<" "
//		<<hTrainTrue[i]->GetBinLowEdge(hTrainTrue[i]->GetNbinsX())+hTrainTrue[i]->GetBinWidth(hTrainTrue[i]->GetNbinsX())<<endl;
		//
		sprintf(name,"trainmeas_eta%d",i);
		hTrainReco[i]= new TH1D (name,name,nrec,xrec);
		hTrainReco[i]->SetLineColor(kRed);
		sprintf(name,"resmat_eta%d",i);
		hMat[i]= new TH2D (name,name,nrec,xrec,ngen,xgen);
		//
    sprintf(name,"response_eta%d",i);
		responseEta[i] = new RooUnfoldResponse (name,name);
		responseEta[i]->Setup(hTrainReco[i],hTrainTrue[i]);
    /////for comparision
		sprintf(name,"resmatt%d_eta%d",0,i);
		hMatt[i][0]= new TH2D (name,name,nrec,xrec,ngen,xgen);
		//
		sprintf(name,"traintruet%d_eta%d",0,i);
		hTrainTruet[i][0]= new TH1D (name,name,ngen,xgen);
		hTrainTruet[i][0]->SetLineColor(kBlue);
		//
		sprintf(name,"trainmeast%d_eta%d",0,i);
		hTrainRecot[i][0]= new TH1D (name,name,nrec,xrec);
		hTrainRecot[i][0]->SetLineColor(kRed);
    //
		sprintf(name,"resmatt%d_eta%d",1,i);
		hMatt[i][1]= new TH2D (name,name,nrec,xrec,ngen,xgen);
		//
		sprintf(name,"traintruet%d_eta%d",1,i);
		hTrainTruet[i][1]= new TH1D (name,name,ngen,xgen);
		hTrainTruet[i][1]->SetLineColor(kBlue);
		//
		sprintf(name,"trainmeast%d_eta%d",1,i);
		hTrainRecot[i][1]= new TH1D (name,name,nrec,xrec);
		hTrainRecot[i][1]->SetLineColor(kRed);
		//Data   -------------------------
		/////
		sprintf(name,"resmatDt%d_eta%d",0,i);
		hMatDt[i]= new TH2D (name,name,nrec,xrec,ngen,xgen);
		//
		sprintf(name,"traintrueDt%d_eta%d",0,i);
		hTrainTrueDt[i]= new TH1D (name,name,ngen,xgen);
		hTrainTrueDt[i]->SetLineColor(kBlue);
		//
		sprintf(name,"trainmeasDt%d_eta%d",0,i);
		hTrainRecoDt[i]= new TH1D (name,name,nrec,xrec);
		hTrainRecoDt[i]->SetLineColor(kRed);
		//
		sprintf(name,"responseDt_eta%d",i);
		responseEtaDt[i] = new RooUnfoldResponse (name,name);
		responseEtaDt[i]->Setup(hTrainRecoDt[i],hTrainTrueDt[i]);
	}
	///generate toy events and smear/////////////////////////////////////////////
	for(unsigned neta(0);neta<6;++neta) {
		generateEvents(eventsGenerate,fres[neta],hnlo[neta],hTrainTrue[neta],hTrainReco[neta],
				hMat[neta],0,lumi,responseEta[neta]);
		generateEvents(eventsGenerate,fres[neta],hnlot[neta][0],hTrainTruet[neta][0],
				hTrainRecot[neta][0],hMatt[neta][0],0,lumi);
		generateEvents(eventsGenerate,fres[neta],hnlot[neta][1],hTrainTruet[neta][1],
				hTrainRecot[neta][1],hMatt[neta][1],0,lumi);
	}
	///////Unfold/////////////////////////////////////////////////////////////////
	cout << "============== UNFOLD ===============" << endl;
	TFile* fDat = new TFile((indir+"/AnalyzerIncHistosData276TeV_8tevcorr.root").c_str(),"read"); 
	TFile* fMC  = new TFile((indir+"/AnalyzerIncHistosMC_8tevcorr.root").c_str(),"read"); 
	TCanvas* cptSpec[6];
	TH1D* hclo[6];     TCanvas* cclo[6];// unfolding closure
	TH1D* htom[6][2];  TCanvas* ctom[6];// unfolded true / mesaured spectrum
	TH1D* hclot[6][2]; TCanvas* cclot[6][2];// unfolding closure for comparision
	TLegend* leg[6];   TCanvas* cMat[6];///to draw matrix
	TCanvas* cptSpecCl[6];
	TCanvas* cClosCom[6][2]; //compare closures from different spectrum, unfolding routines
	TH1D* hclo2[6][2]; //for comparing unfolding routines
	//unfolded histograms //0-bayes, 1-svd, 2- bin by bin
	TH1D* hUnf[3];//unfolded
	TH1D* hUnft[3];//test pdfs

//	/// Check Unfolding/////////////////////////////-------------------------------
//	for(unsigned i(0);i<6;++i) {
//
//		TH1D*  hDat = hTrainReco[i];
//		TH1D*  hGen = (TH1D*)hTrainTrue[i];
//
//		RooUnfoldBayes   unfold1  (responseEta[i], hDat, 3);    // OR
//		RooUnfoldSvd     unfold2  (responseEta[i], hDat, 4);   // OR
//		RooUnfoldBinByBin unfold3 (responseEta[i], hDat);
//
//		hUnf[0]= (TH1D*) unfold1.Hreco();
//		//		hUnf[1]= (TH1D*) unfold2.Hreco();
//		hUnf[2]= (TH1D*) unfold3.Hreco();
//
//		//Unfolding with  different pdfs
//		RooUnfoldBinByBin   unfoldt1_bb  (responseEta[i], hTrainRecot[i][0]);   
//		TH1D* hUnft1_bb = (TH1D*) unfoldt1_bb.Hreco();
//		RooUnfoldBinByBin   unfoldt2_bb  (responseEta[i], hTrainRecot[i][1]);   
//		TH1D* hUnft2_bb = (TH1D*) unfoldt2_bb.Hreco();
//		// 
//		//    RooUnfoldSvd   unfoldt1  (responseEta[i], hTrainRecot[i][0], 4);   
//		//		TH1D* hUnft1= (TH1D*) unfoldt1.Hreco();
//		//    RooUnfoldSvd   unfoldt2  (responseEta[i], hTrainRecot[i][1], 4);   
//		//		TH1D* hUnft2= (TH1D*) unfoldt2.Hreco();
//		////   // 
//		RooUnfoldBayes   unfoldt1  (responseEta[i], hTrainRecot[i][0], 3);   
//		TH1D* hUnft1= (TH1D*) unfoldt1.Hreco();
//		RooUnfoldBayes   unfoldt2  (responseEta[i], hTrainRecot[i][1], 3);   
//		TH1D* hUnft2= (TH1D*) unfoldt2.Hreco();
//
//		//		unfold1.PrintTable(cout,hGen);
//		//		unfold2.PrintTable(cout,hGen);
//
//		char label[100];
//		char label1[100], label2[100];
//		sprintf(label,"%s","NNPDF");
//		drawClosure(cclo[i],hUnf[0],hGen,hclo[i],i,saveplots,label); 
//		sprintf(label1,"%s","ct10");
//		drawClosure(cclot[i][0],hUnft1,hTrainTruet[i][0],hclot[i][0],i,saveplots,label1); 
//		sprintf(label2,"%s","abm11");
//		drawClosure(cclot[i][1],hUnft2,hTrainTruet[i][1],hclot[i][0],i,saveplots,label2); 
//		//draw on the same canvas
//		drawClosure(cClosCom[i][0],i,saveplots,hUnft1,hTrainTruet[i][0],hclot[i][0],label1,
//				hUnft2,hTrainTruet[i][1],hclot[i][0],label2);
//
//		drawClosure2(cClosCom[i][1],i,saveplots,hUnf[0],hDat,hclo2[i][0],"Bayes",
//				hUnf[2],hDat,hclo2[i][1],"BinByBin");
//
//
//		sprintf(name,"c_%s",hMat[i]->GetName());
//		cMat[i] = new TCanvas(name,name,600,600);
//		cMat[i]->SetLogy();cMat[i]->SetLogx();
//		drawMatrix(hMat[i],i);
//		hMat[i]->Draw("colz");
//		sprintf(name,"%.1f<|y|<%.1f",ybins[i],ybins[i+1]);
//		TLatex latex;latex.SetTextSize(0.04);latex.SetTextFont(42);
//		latex.DrawTextNDC(.3,.8,name);
//
//
//		cout<<"Drawing histos eta low edge "<<ybins[i]<<endl;
//		drawSpec(cptSpec[i],hDat,hGen,hUnf[0],leg[i],saveplots,"nnpdf");
//		//save plots  
//		if(saveplots) {
//			cout<<cMat[i]->GetName()<<endl;
//			sprintf(name,"plots/%s.eps",cMat[i]->GetName());
//			cMat[i]->SaveAs(name);
//		}
//	}

	///generate toy events and smear with data res /////////////////////////////////////////
	for(unsigned neta(0);neta<6;++neta) {
		generateEvents(eventsGenerate,fres[neta],hnlo[neta],hTrainTrueDt[neta],hTrainRecoDt[neta],hMatDt[neta],
		0,lumi,responseEtaDt[neta]);
	}
	gROOT->SetBatch(kFALSE);
	///// Unfold Data ------------------------------------------------------------
	TFile* npfile= new TFile("npCorrections.root","read");
	TFile* jteff = new TFile("JetEfficiency.root","read");
  unsigned kiter[6] = {1,2,2,3,4,5};
  unsigned kmark[6] = {24,20,25,21,26,22};
  unsigned kcolr[6] = {kBlack,kTeal+2,kBlack,kBlue,kMagenta,kRed};
  unsigned kstyl[6] = {1,1,1,2,9,10};
	TH1D* hclo_tBay[7][6];//index : [testnum][neta]
	TH1D* hUnf_tBay[7][6];//index : [testnum][neta]
  TCanvas* cUnf_tBay[6];
  TCanvas* cclo_tBay[6];
  TH1D* hDat[6];
  TH1D* hUnfData[6];
	for(unsigned i(0);i<6;++i) {
		ostringstream oss;
		oss<<"hptUnf_eta_"<<ybins[i]<<"_"<<ybins[i+1];
		TH1D*  hDat_fullrange  = (TH1D*)fDat->Get(oss.str().c_str());
		TH1D*  hGen  = (TH1D*)hnlo[i];
		TH1D*  hclo1 = (TH1D*)hTrainTrueDt[i]->Clone("hclo1"); hclo1->Reset();
		TH1D*  hclo2 = (TH1D*)hTrainTrueDt[i]->Clone("hclo2"); hclo2->Reset();
		TH1D*  hclo3 = (TH1D*)hTrainTrueDt[i]->Clone("hclo3"); hclo3->Reset();
		///
		///set the histograms with the corect range
    sprintf(name,"hDat_%d",i);
		hDat[i] = new TH1D(name,name,n2x[i],&x2[i][0]);
		spFillHisto(hDat_fullrange,hDat[i]);
		///correct for jet efficiency
		sprintf(name,"heff_eta%d_fixed",i);
		TH1D* eff = (TH1D*)jteff->Get(name);
		TH1D* heff = (TH1D*)hDat[i]->Clone("hjeteff");
		spFillHisto(eff,heff);
		hDat[i]->Divide(heff);
		hDat[i]->Scale(1/0.99);

   //   RooUnfoldResponse response (hTrainRecoDt[i],hTrainTrueDt[i],hMatDt[i]);
    RooUnfoldBayes    unfold (responseEtaDt[i], hDat[i], 4);
    hUnf[0]= (TH1D*) unfold.Hreco();
//		RooUnfoldBayes    unfold1 (responseEtaDt[i], hDat[i], 3); 
//		RooUnfoldSvd      unfold2 (responseEtaDt[i], hDat[i], 3);  
//		RooUnfoldBinByBin unfold3 (responseEtaDt[i], hDat[i]);
cout<<"hdat bins "<<hDat[i]->GetNbinsX()<<" hunf bins "<<hUnf[0]->GetNbinsX()<<endl; 

//		hUnf[0]= (TH1D*) unfold1.Hreco();
		//	  hUnf[1]= (TH1D*) unfold2.Hreco();
//		hUnf[2]= (TH1D*) unfold3.Hreco();

		cout<<"Unfolding Data, eta low edge "<<ybins[i]<<endl;
		TH1D* htmp = (TH1D*)hTrainTrueDt[i]->Clone("htmp"); htmp->Reset();
		spFillHisto(hGen,htmp);
		//NPCorrs
		sprintf(name,"fratio_eta%d",i);
		TF1* f = (TF1*)npfile->Get(name);
		applyNP(htmp,f);
		applyNP(hGen,f);

		ScaleHistByBW(hUnf[0]);
		hUnf[0]->Scale(1.0/lumi);
		hclo1->Divide(hUnf[0],htmp);


		//////draw histogram
		sprintf(name,"c_hUnf0_eta%d",i);
		cptSpecCl[i] = new TCanvas(name,name,800,600);
		cptSpecCl[i]->SetLogx(); cptSpecCl[i]->SetLogy();
		drawSpecCl(cptSpecCl[i],hUnf[0],hGen,hclo1,i);
		//   hUnf[0]->Draw("");
		//   hnlo[i]->Draw("hist same");
		//hDat->Draw("same");
		//print latex 
		sprintf(name,"%.2f<|y|<%.2f",ybins[i],ybins[i+1]);
		TLatex latex;latex.SetTextSize(0.04);latex.SetTextFont(42);
		latex.DrawTextNDC(.5,.7,name);
		sprintf(name,"plots/hUnfBay_eta%d.eps",i);
		if(saveplots)
			cptSpecCl[i]->SaveAs(name);
		sprintf(name,"hUnf_eta%d",i);
		hUnfData[i] = (TH1D*)hUnf[0]->Clone(name);



		//        ////test bayesian unfolding ------------------------------------
		//		TLegend* leg_clo = new TLegend(0.6,0.7,0.8,0.9);
		//		leg_clo->SetFillColor(0);  leg_clo->SetTextFont(42);
		//		leg_clo->SetBorderSize(0); leg_clo->SetTextSize(0.025);
		//		TLegend* leg_unf = (TLegend*)leg_clo->Clone("leg_unf");
		//		//
		//		RooUnfoldBayes    unf (responseEtaDt[i], hDat[i], 3);
		//		sprintf(name,"hclo_tBay%d_eta%d",0,i);
		//		hclo_tBay[1][i] = (TH1D*)hTrainTrueDt[i]->Clone(name); 
		//		hclo_tBay[1][i]->Reset();
		//		hclo_tBay[1][i]->Sumw2();
		//		hUnf_tBay[1][i] = (TH1D*) unf.Hreco();
		//		hUnf_tBay[1][i]->Sumw2();
		//		htmp->Sumw2();
		//		ScaleHistByBW(hUnf_tBay[1][i]);
		//		hUnf_tBay[1][i]->Scale(1.0/lumi);
		//		gStyle->SetErrorX(0);
		//		ScaleHistByBW(hDat[i]);
		//		hDat[i]->Scale(1./lumi);
		//		sprintf(name,"cclo_tBay_eta%d",i);
		//		cclo_tBay[i] = new TCanvas(name,name,600,600);
		//		hclo_tBay[6][i] = (TH1D*)hTrainTrueDt[i]->Clone("hclo_noUnf");
		//		hclo_tBay[6][i]->Reset();
		//		spFillHisto(hTrainRecoDt[i],hclo_tBay[6][i]);
		//		ScaleHistByBW(hclo_tBay[6][i]);
		//		hclo_tBay[6][i]->Scale(1./lumi);
		//		hclo_tBay[6][i]->Divide(hUnf_tBay[1][i]);
		//		hclo_tBay[6][i]->Draw("pe");
		//		hclo_tBay[6][i]->SetMarkerStyle(34);
		//		hclo_tBay[6][i]->SetMarkerColor(kBlack);
		//		hclo_tBay[6][i]->GetXaxis()->SetRangeUser(74,maxPt[i]-1);
		//		hclo_tBay[6][i]->GetYaxis()->SetRangeUser(0,4);
		//		hclo_tBay[6][i]->GetYaxis()->SetTitle("Reco/Unfolded");
		//		leg_clo->AddEntry(hclo_tBay[6][i],"","pel");
		//		leg_clo->Draw("same");
		//
		//		TH1D* hunftmp[6];
		//		for(unsigned niter(0);niter<6;++niter) {
		//			sprintf(name,"hDattmp_%d",i);
		//			TH1D* hdattmp = (TH1D*)hDat[i]->Clone(name);
		//			spFillHisto(hDat_fullrange,hdattmp);
		//			hdattmp->Divide(heff);
		//			RooUnfoldBayes    unf (responseEtaDt[i], hdattmp, kiter[niter]);
		//			hUnf_tBay[niter][i] = (TH1D*) unf.Hreco();
		//			hUnf_tBay[niter][i]->Sumw2();
		//			ScaleHistByBW(hUnf_tBay[niter][i]);
		//			hUnf_tBay[niter][i]->Scale(1.0/lumi);
		//		}
		//
		//		for(unsigned niter(3);niter<6;++niter) {
		//			hunftmp[niter] = (TH1D*)hUnf_tBay[niter][i]->Clone("tmp");
		//			hunftmp[niter]->Reset();
		//			hunftmp[niter]->Divide(hUnf_tBay[niter][i],hUnf_tBay[niter-1][i]);
		//			hunftmp[niter]->SetMarkerStyle(kmark[niter]);
		//			hunftmp[niter]->SetMarkerColor(kcolr[niter]);
		//			hunftmp[niter]->SetLineStyle(kstyl[niter]);
		//			hunftmp[niter]->SetLineColor(kcolr[niter]);
		//			//      sprintf(name,"cUnf_tBay_eta%d_niter%d",i,niter);
		//			if(niter==3) {
		//				sprintf(name,"cUnf_tBay_eta%d",i);
		//				cUnf_tBay[i] = new TCanvas(name,name,600,600);
		//				hunftmp[niter]->GetXaxis()->SetRangeUser(74,maxPt[i]-1);
		//				hunftmp[niter]->GetYaxis()->SetRangeUser(0,2);
		//				hunftmp[niter]->GetYaxis()->SetTitle("Ratio of Unfolded Spectrum");
		//				hunftmp[niter]->Draw("pe");
		//			}
		//			else {
		//				hunftmp[niter]->Draw("pe same");
		//			}
		//			sprintf(name," niter %d/niter %d",kiter[niter],kiter[niter-1]);
		//			leg_unf->AddEntry(hunftmp[niter],name,"pel");
		//		}
		//		leg_unf->Draw("same");
		//		sprintf(name,"plots/iter_eta%d.eps",i);
		//    if(saveplots)
		//		cUnf_tBay[i]->SaveAs(name);
		//		////test bayesian unfolding ------------------------------------

	} ////eta loop
	//	//
	///Write Histos and Response object to file
	TFile* outfile = new TFile("UnfoldJets.root","recreate");
	for(unsigned i(0);i<nybins-1;++i)
	{
		sprintf(name,"response_eta%i",i);
		hTrainTrue[i]->Write();
		hTrainReco[i]->Write();
		hMat[i]->Write();
		cptSpecCl[i]->Write();
		hUnfData[i]->Write();
	}
	outfile->Close();
}


//==============================================================================
// Helper Functions And Classes
//==============================================================================
void drawSpec(TCanvas* c,TH1* hreco, TH1* hgen, TH1* hunf,TLegend* leg,bool saveplots,char* label) {
	sprintf(name,"c_%s",hreco->GetName());
	c = new TCanvas(name,name,600,600);
	c->cd();
	hunf->SetLineColor(kRed);
	hunf->SetMarkerStyle(21);
	hunf->SetMarkerSize(0.3);
	hunf->GetYaxis()->SetRangeUser(1e-8,1e8);
	hunf->Draw("hist");
	c->SetLogx(1);c->SetLogy(1);
	hreco->SetLineColor(kGreen);
	hreco->Draw("same");
	hgen->SetLineColor(kBlack);
	hgen->Draw("e same");

	leg = new TLegend(0.6,0.7,0.87,0.89);
	leg->SetFillColor(0);  leg->SetTextFont(42);
	leg->SetBorderSize(0); leg->SetTextSize(0.025);
	leg->AddEntry(hreco,"Reco","pel");
	leg->AddEntry(hgen ,"Gen ","pel");
	leg->AddEntry(hunf ,"Unfolded Spectrum","pel");
	leg->Draw("same");

	if(label) {
		TLatex latex;latex.SetTextSize(0.04);latex.SetTextFont(42);
		latex.DrawTextNDC(.7,.6,label);
	}
	if(saveplots) {
		sprintf(name,"plots/%s_%s.eps",c->GetName(),label);
		c->SaveAs(name);
	}
}

void drawClosure(TCanvas* c,TH1* hunf, TH1* hgen, TH1* hclo,unsigned eta,bool saveplots,char* label) {
	sprintf(name,"c_%s",hgen->GetName());
	cout<<name<<endl;
	c = new TCanvas(name,name,600,600);
	c->cd();
	c->SetLogx(1);
	sprintf(name,"hclos_%s",hgen->GetName());
	hclo = (TH1D*)hgen->Clone(name);
	hclo->Reset(); 
	hclo->Sumw2();
	hclo->Divide(hunf,hgen);
	hclo->GetYaxis()->SetTitle("Unfolded/Gen");
	hclo->GetXaxis()->SetTitle("Jet Pt Gev/c");
	hclo->GetYaxis()->SetRangeUser(0.9,1.1);
	hclo->GetXaxis()->SetRangeUser(18,maxPt[eta]-1);
	hclo->GetYaxis()->SetTitleOffset(1.5);
	hclo->SetLineColor(kBlack);
	hclo->SetMarkerColor(kBlack);
	hclo->SetMarkerStyle(20);
	hclo->SetMarkerSize(1.0);
	hclo->Draw();

	sprintf(name,"%.1f<|y|<%.1f",ybins[eta],ybins[eta+1]);
	TLatex latex;latex.SetTextSize(0.04);latex.SetTextFont(42);
	latex.DrawTextNDC(.6,.8,name);
	latex.DrawTextNDC(.6,.7,label);

	if(saveplots) {
		sprintf(name,"plots/%s_%s.eps",c->GetName(),label);
		c->SaveAs(name);
	}

}

void drawClosure(TCanvas* c,unsigned eta, bool saveplots,
		TH1* hunf1, TH1* hgen1, TH1* hclo1, char* label1,
		TH1* hunf2, TH1* hgen2, TH1* hclo2, char* label2)
{
	sprintf(name,"c_%s",hgen1->GetName());
	cout<<name<<endl;
	c = new TCanvas(name,name,600,600);
	c->cd();
	c->SetLogx(1);
	sprintf(name,"hclos_%s",hgen1->GetName());
	hclo1 = (TH1D*)hgen1->Clone(name);
	hclo1->Reset(); hclo1->Sumw2();
	hclo2 = (TH1D*)hgen2->Clone(name);
	hclo2->Reset(); hclo2->Sumw2();
	//make hclo and draw
	hclo1->Divide(hunf1,hgen1);
	hclo1->GetYaxis()->SetTitle("Unfolded/Gen");
	hclo1->GetXaxis()->SetTitle("Jet Pt Gev/c");
	hclo1->GetYaxis()->SetRangeUser(0.9,1.1);
	hclo1->GetXaxis()->SetRangeUser(18,maxPt[eta]-1);
	hclo1->GetYaxis()->SetTitleOffset(1.5);
	hclo1->SetLineColor(kBlack);
	hclo1->SetMarkerColor(kBlack);
	hclo1->SetMarkerStyle(20);
	hclo1->SetMarkerSize(1.0);
	hclo1->Draw();
	hclo2->Divide(hunf2,hgen2);
	hclo2->GetYaxis()->SetTitle("Unfolded/Gen");
	hclo2->GetXaxis()->SetTitle("Jet Pt Gev/c");
	hclo2->GetYaxis()->SetRangeUser(0.9,1.1);
	hclo2->GetXaxis()->SetRangeUser(18,maxPt[eta]-1);
	hclo2->GetYaxis()->SetTitleOffset(1.5);
	hclo2->SetLineColor(kBlue);
	hclo2->SetMarkerColor(kBlue);
	hclo2->SetMarkerStyle(33);
	hclo2->SetMarkerSize(1.0);
	hclo2->Draw("same");
	///find the slope
	TF1* f = new TF1("fline","pol1",18,maxPt[eta]-1);
	hclo1->Fit(f,"R same");
	hclo2->Draw("same");

	TLegend* leg = new TLegend(0.6,0.7,0.87,0.89);    
	leg->SetFillColor(0);  leg->SetTextFont(42);
	leg->SetBorderSize(0); leg->SetTextSize(0.035);
	leg->AddEntry(hclo1,label1,"pel");
	leg->AddEntry(hclo2,label2,"pel");
	leg->Draw("same");

	sprintf(name,"%.1f<|y|<%.1f",ybins[eta],ybins[eta+1]);
	TLatex latex;latex.SetTextSize(0.04);latex.SetTextFont(42);
	latex.DrawTextNDC(.4,.8,name);
	//    latex.DrawTextNDC(.6,.7,label1);
	if(saveplots) {
		sprintf(name,"plots/%s_%s_%s.eps",c->GetName(),label1,label2);
		c->SaveAs(name);
	}

}

void drawClosure2(TCanvas* c,unsigned eta, bool saveplots,
		TH1* hunf1, TH1* hgen1, TH1* hclo1, char* label1,
		TH1* hunf2, TH1* hgen2, TH1* hclo2, char* label2)
{
	sprintf(name,"c_%s_%s",hgen1->GetName(),"unfComp");
	cout<<name<<endl;
	c = new TCanvas(name,name,600,600);
	c->cd();
	c->SetLogx(1);
	sprintf(name,"hclos_%s",hgen1->GetName());
	hclo1 = (TH1D*)hgen1->Clone(name);
	hclo1->Reset(); hclo1->Sumw2();
	hclo2 = (TH1D*)hgen2->Clone(name);
	hclo2->Reset(); hclo2->Sumw2();
	//make hclo and draw
	hclo1->Divide(hunf1,hgen1);
	hclo1->GetYaxis()->SetTitle("Unfolded/Reco");
	hclo1->GetXaxis()->SetTitle("Jet Pt Gev/c");
	hclo1->GetYaxis()->SetRangeUser(0.3,1.1);
	hclo1->GetXaxis()->SetRangeUser(74,maxPt[eta]-1);
	hclo1->GetYaxis()->SetTitleOffset(1.5);
	hclo1->SetLineColor(kBlack);
	hclo1->SetMarkerColor(kBlack);
	hclo1->SetMarkerStyle(20);
	hclo1->SetMarkerSize(1.0);
	hclo1->Draw();
	hclo2->Divide(hunf2,hgen2);
	hclo2->GetYaxis()->SetTitle("Unfolded/Gen");
	hclo2->GetXaxis()->SetTitle("Jet Pt Gev/c");
	hclo2->GetYaxis()->SetRangeUser(0.3,1.1);
	hclo2->GetXaxis()->SetRangeUser(74,maxPt[eta]-1);
	hclo2->GetYaxis()->SetTitleOffset(1.5);
	hclo2->SetLineColor(kBlue);
	hclo2->SetMarkerColor(kBlue);
	hclo2->SetMarkerStyle(24);
	hclo2->SetMarkerSize(1.0);
	hclo2->Draw("same");

	TLegend* leg = new TLegend(0.6,0.7,0.87,0.89);    
	leg->SetFillColor(0);  leg->SetTextFont(42);
	leg->SetBorderSize(0); leg->SetTextSize(0.035);
	leg->AddEntry(hclo1,label1,"pel");
	leg->AddEntry(hclo2,label2,"pel");
	leg->Draw("same");

	sprintf(name,"%.1f<|y|<%.1f",ybins[eta],ybins[eta+1]);
	TLatex latex;latex.SetTextSize(0.04);latex.SetTextFont(42);
	latex.DrawTextNDC(.4,.8,name);
	//    latex.DrawTextNDC(.6,.7,label1);

	if(saveplots) {
		sprintf(name,"plots/%s_%s_%s.eps",c->GetName(),label1,label2);
		c->SaveAs(name);
	}


}



void drawMatrix(TH2* h,unsigned eta)
{
	h->GetXaxis()->SetRangeUser(74,maxPt[eta]-1);
	h->GetYaxis()->SetRangeUser(0,maxPt[eta]-1);
	h->GetYaxis()->SetTitleOffset(1.2);
	h->GetXaxis()->SetTitleOffset(1.2);
	h->GetYaxis()->SetTitle("Gen Pt Gev/c");
	h->GetXaxis()->SetTitle("Reco Pt Gev/c");

}



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
	sprintf(name,"%s/AnalyzerIncHistosMC_8tevcorr.root",indir);
	TFile* f1 = new TFile(name,"read");
	for(unsigned i(0);i<nybins-1;++i)
	{
		sprintf(name,"hgenUnf_eta%i",i);
		h[i] = (TH1F*)f1->Get(name);
		//	cout<<h[i]->GetName()<<' '<<h[i]<<endl;
		ScaleHistByBW(h[i]);
		h[i]->Scale(1/lumi);
	}
};

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

void normalizeTH2DX(TH2D *h)
{
	int xbin(h->GetNbinsX()), ybin(h->GetNbinsY());
	for (int i(1); i <= xbin; i++){
		double sum(0.);
		for (int j(1); j <= ybin; j++){
			sum += h->GetBinContent(i,j);
		}
		for (int j(1); j <= ybin; j++){
			if (sum > 0){
				double val = h->GetBinContent(i,j)/sum;
				val=val*100;
				val=ceil(val);
				double nearest=val/100;
				h->SetBinContent(i,j,nearest);
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


void generateEvents(unsigned events, TF1* fres, TH1* hnlo,
		TH1D* hTrue, TH1D* hReco,TH2D* hMatrix,double dataRes, 
		double lumi, RooUnfoldResponse* resp)
{
	hTrue->Sumw2(); hReco->Sumw2();
	const unsigned nbins(hnlo->GetNbinsX());
	// cout<<nbins<<' '<<events<<endl;
	for(unsigned bin(1);bin<=nbins;++bin)
	{
		double xsec  = hnlo->GetBinContent(bin)*hnlo->GetBinWidth(bin)*lumi;
		double binlo = hnlo->GetBinLowEdge(bin);
		double binhi = hnlo->GetBinWidth(bin)+binlo;
		for(unsigned i(0);i<events;++i)
		{
			double genpt = gRandom->Uniform(binlo,binhi);
			double resol = gRandom->Gaus(1.0,fres->Eval(genpt));
			if(dataRes) resol = resol; //for data resolution
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

void drawSpecCl(TCanvas* c,TH1* hdt, TH1* hmc,TH1* hclo, unsigned eta)
{
	hclo->GetXaxis()->SetTitle("Jet Pt Gev/c");
	hclo->GetXaxis()->SetTitleSize(0.1);
	hclo->GetXaxis()->SetLabelSize(0.1);
	hclo->GetYaxis()->SetLabelSize(0.1);
	hclo->SetMarkerStyle(20);
	hclo->GetXaxis()->SetRangeUser(74,maxPt[eta]-1);
	hclo->SetMarkerColor(kBlack);
	hclo->GetYaxis()->SetRangeUser(0.,2.0);
	hdt->GetXaxis()->SetTitle("");
	hdt->GetXaxis()->SetRangeUser(74,maxPt[eta]-1);
	hmc->GetXaxis()->SetRangeUser(74,maxPt[eta]-1);
	cout<<"\n eta "<<eta<<" maxpt "<<maxPt[eta]<<endl;
	hdt->GetYaxis()->SetRangeUser(1e-3,1e8);
	hdt->SetMarkerStyle(21);
	hdt->SetMarkerColor(kBlack);
	hdt->GetXaxis()->SetTitleSize(0.00001);
	hdt->GetXaxis()->SetLabelSize(0.00001);
	c->cd();
	sprintf(name,"tpad_%s",hdt->GetName());
	double ytop(0.25);
	TPad* tpad = new TPad(name,name,0,ytop, 1.0, 1.0);
	tpad->SetLogx(); tpad->SetLogy();
	tpad->Draw();
	tpad->cd();
	tpad->SetBottomMargin(0.1);
	hdt->Draw();
	hmc->Draw("hist same");

	c->cd();
	double ybot(0.3);
	sprintf(name,"bpad_%s",hclo->GetName());
	TPad* bpad = new TPad("bpad", "bpad", 0, 0.03, 1.0, 0.3);
	bpad->SetLogx();
	bpad->Draw();
	bpad->cd();
	bpad->SetTopMargin(0.1);
	bpad->SetBottomMargin(0.2);
	bpad->SetGrid(1,1);
	hclo->Draw();
	tpad->cd();
}

void spFillHisto(TH1* hin, TH1* hout,bool verbose)
{
	hout->Reset();
	const  unsigned numbins(hout->GetNbinsX());
	const  unsigned numbinsIn(hin->GetNbinsX());
	cout<<"numbins "<<numbins<<' '<<hout->GetName()<<endl;
	for(unsigned i(1);i<=numbins;++i)
		for(unsigned j(1);j<=numbinsIn;++j)
			if(((hin->GetBinLowEdge(j)+hin->GetBinWidth(j))==
						(hout->GetBinLowEdge(i)+hout->GetBinWidth(i))	)	
					&& hin->GetBinLowEdge(j)==hout->GetBinLowEdge(i) 
				)
			{
				hout->SetBinContent(i,hin->GetBinContent(j));
				hout->SetBinError  (i,hin->GetBinError(j));
				if(verbose)
					cout<<"Setting bin content "<<hout->GetBinLowEdge(i)<<" "<<hin->GetBinContent(j)
						<<" +- "<<hin->GetBinError(j)<<endl;
			}
			else {
				if(true && verbose)
					cout<<"No match bin hout "<<hout->GetBinLowEdge(i)<<"-"<<hout->GetBinLowEdge(i)+hout->GetBinWidth(i)
						<<"  hin "<<hin->GetBinLowEdge(j)<<"-"<<hin->GetBinLowEdge(j)+hin->GetBinWidth(j)
						<<endl;
			}
}

void applyNP(TH1* h,TF1* f)
{
	const  unsigned numbins(h->GetNbinsX());
	for(unsigned i(1);i<=numbins;++i)
	{
		double np = f->Eval(h->GetBinCenter(i));
		h->SetBinContent(i,h->GetBinContent(i)*np);
		//		cout<<"NPCorr "<<h->GetBinCenter(i)<<" "<<np<<endl;
	}
}

void testFunction(double xlo, double xhi)
{
	cout<<findBinsFromArray(xlo,xhi,x[0],nx[0])<<endl; 
}


