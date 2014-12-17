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
void addError(double* in, double* out, unsigned size);
void fillHistoFromArr(double* arr,TH1D* h);
void initArray(double* arr,unsigned size);
void applyNP(TH1* h,TF1* f);
void normalizeTH2D(TH2D* h);
void normalizeTH2DX(TH2D* h);
void ScaleHistByBW(TH1* h);
void histoFromErr(TH1* hin, TH1* hout);
void drawErrors(TCanvas* c, TH1* rec,TH1* unf, unsigned eta);
TH1D* makeNLOSpectrum(char* infile,char* title);
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
		double lumi,unsigned eta,RooUnfoldResponse* resp=0);
void generateEventsTest(unsigned events, TF1* fres, TH1* hnlo,
		TH1D* hTrue, TH1D* hReco,TH2D* hMatrix, 
		double lumi,unsigned eta,double scaleResolutionBy,RooUnfoldResponse* resp=0);
void generateEventsMC(unsigned events, TF1* fres, TH1* hnlo,
		TH1D* hTrue, TH1D* hReco,TH2D* hMatrix,double dataRes,
		double lumi,unsigned eta,RooUnfoldResponse* resp=0);
TH2D* CorrelationHist (const TMatrixD& cov, const char* name, const char* title,
                      const  double* bins, const unsigned nbins, unsigned i);

void getUnfoldUncertainty(TH1* hcn,TH1* hup,TH1* hdn,TH1* hunc);
void getParUncertainty(TH1D* hcn,TH1D** hvar,TH1D* hunc,unsigned size);
void getJesUncertainty(TH1* hcn, TH1* hup, TH1* hdn, double* uncert);
//==============================================================================
// Unfolding
//==============================================================================
void unfoldJesUnc() {
	bool saveplots(true);
//	if(saveplots) 	gROOT->SetBatch(kTRUE);

	//	const int nsrc = 19;
	//	const char* srcnames[nsrc] = {
	//		"AbsoluteStat", "AbsoluteScale","AbsoluteFlavMap", "AbsoluteMPFBias",
	//		"HighPtExtra","SinglePionECAL","SinglePionHCAL","FlavorQCD","Time",
	//		"RelativeJEREC1","RelativeJEREC2","RelativeJERHF","RelativePtBB",
	//		"RelativePtEC1","RelativePtEC2","RelativePtHF","RelativeFSR",
	//		"RelativeStatEC2","RelativeStatHF"
	//   //   ,"PileUpDataMC" ,"PileUpPtBB" ,"PileUpPtEC1"  ,"PileUpPtEC2"
	//   //   ,"PileUpPtHF"   ,"PileUpBias"
	//		} ;
	const int nsrc = 22;
	//winter14_v5
	const char* srcnames[nsrc] ={
		"AbsoluteStat","AbsoluteScale","AbsoluteFlavMap","AbsoluteMPFBias",
		"HighPtExtra","SinglePionECAL","SinglePionHCAL","FlavorQCD",
		"TimeEta","TimePt","RelativeJEREC1","RelativeJEREC2","RelativeJERHF",
		"RelativePtBB","RelativePtEC1","RelativePtEC2","RelativePtHF",
		"RelativeFSR","RelativeStatEC2","RelativeStatHF","PileUpMuZero",
		"TimeRunD"
	};

	const string indir("/uscms/home/tlibeiro/inclusivexsec/CMSSW_5_3_14/src/JetAnalyzer/JetAnalyzer/test");
	const char* nloDir("/uscms_data/d3/tlibeiro/inclusivexsec/Generators/generateIncJets/inputXsecs");
	const char* fitDir("./");
	sprintf(name,"%s/ResolutionFits.root",fitDir);
	TFile* resFile = new TFile(name,"read");
	TF1*  frescn[6];//resfits
	const unsigned eventsGenerate(1e6);

	cout << "==================================== TRAIN ====================================" << endl;
	//for data
	RooUnfoldResponse* responseEtacn[nybins-1]; 
	RooUnfoldResponse* responseEtaup[nybins-1]; 
	RooUnfoldResponse* responseEtadn[nybins-1]; 
	TH1D* hnlo[6];//nlo spectrum
	TH1D* hTrainTruecn[6];  TH1D* hTrainRecocn[6];  TH2D* hMatcn[6];
	///Initialize Histos and Response class//////////////////////////////////////
	for(unsigned i(0);i<6;++i) {
		char infile[500]; char title[500];
		sprintf(infile,"%s/ct10y%i.cppInput",nloDir,i);
		sprintf(title ,"nloy%i",i);
		hnlo[i] = makeNLOSpectrum(infile,title);
		//fres for smearing        
		sprintf(name,"fresfit_eta%d",i);
		frescn[i] = (TF1*)resFile->Get(name);
		/////training histos
		unsigned ngen(n1x[i]), nrec(nx[i]);
		const double *xgen = x1[i]; //x bins gen
		const double *xrec = x[i];// x bins reco
		//Data   -------------------------
		/////
		sprintf(name,"resmatcn%d_eta%d",0,i);
		hMatcn[i]= new TH2D (name,name,ngen,xgen,ngen,xgen);
		//
		sprintf(name,"traintruecn%d_eta%d",0,i);
		hTrainTruecn[i]= new TH1D (name,name,ngen,xgen);
		hTrainTruecn[i]->SetLineColor(kBlue);
		//
		sprintf(name,"trainmeascn%d_eta%d",0,i);
		hTrainRecocn[i]= new TH1D (name,name,nrec,xrec);
		hTrainRecocn[i]->SetLineColor(kRed);
		//
		sprintf(name,"responsecn_eta%d",i);
		responseEtacn[i] = new RooUnfoldResponse (name,name);
		responseEtacn[i]->Setup(hTrainRecocn[i],hTrainTruecn[i]);
	}
	///generate toy events and smear with data res /////////////////////////////////////////
	for(unsigned neta(0);neta<6;++neta) {
		generateEvents(eventsGenerate,frescn[neta],hnlo[neta],hTrainTruecn[neta],hTrainRecocn[neta],hMatcn[neta],
				0,lumi,neta,responseEtacn[neta]);
	}
	///// Unfold Data ------------------------------------------------------------
	cout << "================================ UNFOLD ================================" << endl;
	TCanvas* cjeserr[6];
	TH1D*  hjesunc[6][nsrc+2];//src+total+totalUncertainty (adding  all sources before varying pt)
	TFile* fjes = new TFile((indir+"/AnalyzerIncHistosMC_winter14_jecUnc.root").c_str(),"read"); 
	TH1D* hjesup[6];	TH1D* hjesdn[6];	TH1D* hjescn[6];
	TH1D* hjesup_allsrc[6];	TH1D* hjesdn_allsrc[6];	//calculated with all sources combined
	for(unsigned i(0);i<6;++i) {
		double unctotal[100]; ///uncertainties per bin
		initArray(unctotal,100);
		for(unsigned isrc(0);isrc<nsrc;++isrc)
		{
			//get the spectreum
			sprintf(name,"hptUnf_eta%d",i);
			hjescn[i]  = (TH1D*)fjes->Get(name);
			sprintf(name,"%s_dn_eta%d",srcnames[isrc],i);
			hjesdn[i]  = (TH1D*)fjes->Get(name);
			sprintf(name,"%s_up_eta%d",srcnames[isrc],i);
			hjesup[i]  = (TH1D*)fjes->Get(name);
			if(isrc==0)//checking if totalsrc = uncertainty with all sources
			{
				sprintf(name,"TotalUnc_up_eta%d",i);
				hjesup_allsrc[i] = (TH1D*)fjes->Get(name);
				sprintf(name,"TotalUnc_dn_eta%d",i);
				hjesdn_allsrc[i] = (TH1D*)fjes->Get(name);
			}
			///set the histograms with the corect range
			sprintf(name,"hcn_%s_%d",srcnames[isrc],i);
			TH1D* hcn = new TH1D(name,name,n2x[i],&x2[i][0]);
			spFillHisto( hjescn[i],hcn);
			sprintf(name,"hup_%s_%d",srcnames[isrc],i);
			TH1D* hup = new TH1D(name,name,n2x[i],&x2[i][0]);
			spFillHisto( hjesup[i],hup);
			sprintf(name,"hdn_%s_%d",srcnames[isrc],i);
			TH1D* hdn = new TH1D(name,name,n2x[i],&x2[i][0]);
			spFillHisto( hjesdn[i],hdn);
//checking if totalsrc = uncertainty with all sources
      TH1D *hup_allsrc, *hdn_allsrc;
			if(isrc==0)
			{
				sprintf(name,"hup_allsrc_%d",i);
				hup_allsrc = new TH1D(name,name,n2x[i],&x2[i][0]);
				spFillHisto( hjesup_allsrc[i],hup_allsrc);
				sprintf(name,"hdn_allsrc_%d",i);
				hdn_allsrc = new TH1D(name,name,n2x[i],&x2[i][0]);
				spFillHisto( hjesdn_allsrc[i],hdn_allsrc);
			}
			//unfold
			TH1D*  hUnfcn; TH1D*  hUnfup;TH1D*  hUnfdn;
			RooUnfoldBayes    unfoldcn (responseEtacn[i], hcn, 4);
			RooUnfoldBayes    unfoldup (responseEtacn[i], hup, 4);
			RooUnfoldBayes    unfolddn (responseEtacn[i], hdn, 4);
			hUnfcn= (TH1D*) unfoldcn.Hreco();
			hUnfup= (TH1D*) unfoldup.Hreco();
			hUnfdn= (TH1D*) unfolddn.Hreco();
			//checking if totalsrc = uncertainty with all sources
			RooUnfoldBayes    unfoldup_allsrc (responseEtacn[i], hup_allsrc, 4);
			RooUnfoldBayes    unfolddn_allsrc (responseEtacn[i], hdn_allsrc, 4);
			TH1D *hUnfup_allsrc, *hUnfdn_allsrc;
			if(isrc==0)
			{
				hUnfup_allsrc= (TH1D*) unfoldup_allsrc.Hreco();
				hUnfdn_allsrc= (TH1D*) unfolddn_allsrc.Hreco();
			}
			//caculate error
			double uncert[100]; ///uncertainties per bin
			initArray(uncert,100);
			getJesUncertainty(hUnfcn,hUnfup,hUnfdn,uncert);
			sprintf(name,"%s_up_eta%d",srcnames[isrc],i);
			hjesunc[i][isrc] = (TH1D*)hUnfup->Clone(name);
			fillHistoFromArr(uncert,hjesunc[i][isrc]);
			//total error
			addError(uncert,unctotal,100);
			if(isrc==0)
			{
				sprintf(name,"Total_up_eta%d",i);
				hjesunc[i][nsrc] = (TH1D*)hUnfup->Clone(name);
			}
			//checking if totalsrc = uncertainty with all sources
			double uncert_allsrc[100];   initArray(uncert_allsrc,100);
			if(isrc==0)
			{
				getJesUncertainty(hUnfcn,hUnfup_allsrc,hUnfdn_allsrc,uncert_allsrc);
				sprintf(name,"Total_allsrc_up_eta%d",i);
				hjesunc[i][nsrc+1] = (TH1D*)hUnfup->Clone(name);
				fillHistoFromArr(uncert_allsrc,hjesunc[i][nsrc+1]);
			}
		}///nsrc loop
		//total error
		fillHistoFromArr(unctotal,hjesunc[i][nsrc]);
	}//neta
	///Write Histos and Response object to file
	TFile* outfile = new TFile("JESUncertainties.root","recreate");
	for(unsigned i(0);i<6;++i)
		for(unsigned isrc(0);isrc<=nsrc+1;++isrc)
			hjesunc[i][isrc]->Write();
	outfile->Close();
}


//new definition of unfolding uncertainty 
//fitting parameters varied by 1-sigma,not 5% as done in the unfoldUnc()
void unfoldUncNew() {
	// Data Members
	TH1D*              hTrainTrue[10];
	TH1D*              hTrainReco[10];
	TH2D               *hMat[10];
	RooUnfoldResponse* response;
	RooUnfold*         unfold;

	//  bool saveplots(true);
	//  if(saveplots) 	gROOT->SetBatch(kTRUE);
	const string indir("/uscms/home/tlibeiro/inclusivexsec/CMSSW_5_3_14/src/JetAnalyzer/JetAnalyzer/test");
	const char* nloDir("/uscms_data/d3/tlibeiro/inclusivexsec/Generators/generateIncJets/inputXsecs");
	const char* fitDir("./");
	sprintf(name,"%s/ResolutionFits.root",fitDir);
	TFile* resFile = new TFile(name,"read");
	TF1*  frescn[6];//resfits mc
	TF1*  fresup[6];//resfits
	TF1*  fresdn[6];//resfits
	TF1*  fresdtcn[6];//resfits data
	TF1*  fresdtup[6];//resfits
	TF1*  fresdtdn[6];//resfits
	const unsigned eventsGenerate(1e6);

	cout << "==================================== TRAIN ====================================" << endl;
	//for mc fit
	RooUnfoldResponse* responseEtacn[nybins-1]; 
	RooUnfoldResponse* responseEtaup[nybins-1]; 
	RooUnfoldResponse* responseEtadn[nybins-1]; 
	//for data fit
	RooUnfoldResponse* responseEtadtcn[nybins-1]; 
	RooUnfoldResponse* responseEtadtup[nybins-1]; 
	RooUnfoldResponse* responseEtadtdn[nybins-1]; 
	TH1D* hnlo[6];//nlo spectrum
	TH1D* hTrainTruecn[6];  TH1D* hTrainRecocn[6];  TH2D* hMatcn[6];
	TH1D* hTrainTrueup[6];  TH1D* hTrainRecoup[6];  TH2D* hMatup[6];
	TH1D* hTrainTruedn[6];  TH1D* hTrainRecodn[6];  TH2D* hMatdn[6];
	///Initialize Histos and Response class//////////////////////////////////////
	for(unsigned i(0);i<6;++i) {
		char infile[500]; char title[500];
		sprintf(infile,"%s/ct10y%i.cppInput",nloDir,i);
		sprintf(title ,"nloy%i",i);
		hnlo[i] = makeNLOSpectrum(infile,title);
		//fres for smearing        
		sprintf(name,"fresfit_eta%d",i);
		frescn[i] = (TF1*)resFile->Get(name);
		sprintf(name,"fresfit_up_eta%d",i);
		fresup[i] = (TF1*)resFile->Get(name);
		sprintf(name,"fresfit_dn_eta%d",i);
		fresdn[i] = (TF1*)resFile->Get(name);
		//fres data for smearing        
		sprintf(name,"fresfitdt_eta%d",i);
		fresdtcn[i] = (TF1*)resFile->Get(name);
		sprintf(name,"fresfitdt_up_eta%d",i);
		fresdtup[i] = (TF1*)resFile->Get(name);
		sprintf(name,"fresfitdt_dn_eta%d",i);
		fresdtdn[i] = (TF1*)resFile->Get(name);
		/////training histos
		unsigned ngen(n1x[i]), nrec(nx[i]);
		const double *xgen = x1[i]; //x bins gen
		const double *xrec = x[i];// x bins reco
		//Data   -------------------------
		/////
		sprintf(name,"resmatcn%d_eta%d",0,i);
		hMatcn[i]= new TH2D (name,name,ngen,xgen,ngen,xgen);
		//
		sprintf(name,"traintruecn%d_eta%d",0,i);
		hTrainTruecn[i]= new TH1D (name,name,ngen,xgen);
		hTrainTruecn[i]->SetLineColor(kBlue);
		//
		sprintf(name,"trainmeascn%d_eta%d",0,i);
		hTrainRecocn[i]= new TH1D (name,name,nrec,xrec);
		hTrainRecocn[i]->SetLineColor(kRed);
		//
		sprintf(name,"responsecn_eta%d",i);
		responseEtacn[i] = new RooUnfoldResponse (name,name);
		responseEtacn[i]->Setup(hTrainRecocn[i],hTrainTruecn[i]);
		sprintf(name,"responsedtcn_eta%d",i);
		responseEtadtcn[i] = new RooUnfoldResponse (name,name);
		responseEtadtcn[i]->Setup(hTrainRecocn[i],hTrainTruecn[i]);
		//up variation----------------------------
		sprintf(name,"resmatup%d_eta%d",0,i);
		hMatup[i]= new TH2D (name,name,ngen,xgen,ngen,xgen);
		//
		sprintf(name,"traintrueup%d_eta%d",0,i);
		hTrainTrueup[i]= new TH1D (name,name,ngen,xgen);
		hTrainTrueup[i]->SetLineColor(kBlue);
		//
		sprintf(name,"trainmeasup%d_eta%d",0,i);
		hTrainRecoup[i]= new TH1D (name,name,nrec,xrec);
		hTrainRecoup[i]->SetLineColor(kRed);
		//
		sprintf(name,"responseup_eta%d",i);
		responseEtaup[i] = new RooUnfoldResponse (name,name);
		responseEtaup[i]->Setup(hTrainRecoup[i],hTrainTrueup[i]);
		sprintf(name,"responsedtup_eta%d",i);
		responseEtadtup[i] = new RooUnfoldResponse (name,name);
		responseEtadtup[i]->Setup(hTrainRecoup[i],hTrainTrueup[i]);
		//down variation----------------------------
		sprintf(name,"resmatdn%d_eta%d",0,i);
		hMatdn[i]= new TH2D (name,name,ngen,xgen,ngen,xgen);
		//
		sprintf(name,"traintruedn%d_eta%d",0,i);
		hTrainTruedn[i]= new TH1D (name,name,ngen,xgen);
		hTrainTruedn[i]->SetLineColor(kBlue);
		//
		sprintf(name,"trainmeasdn%d_eta%d",0,i);
		hTrainRecodn[i]= new TH1D (name,name,nrec,xrec);
		hTrainRecodn[i]->SetLineColor(kRed);
		//
		sprintf(name,"responsedn_eta%d",i);
		responseEtadn[i] = new RooUnfoldResponse (name,name);
		responseEtadn[i]->Setup(hTrainRecodn[i],hTrainTruedn[i]);
		sprintf(name,"responsedtdn_eta%d",i);
		responseEtadtdn[i] = new RooUnfoldResponse (name,name);
		responseEtadtdn[i]->Setup(hTrainRecodn[i],hTrainTruedn[i]);
	}
	///generate toy events and smear with data res /////////////////////////////////////////
	for(unsigned neta(0);neta<6;++neta) {
		generateEvents(eventsGenerate,frescn[neta],hnlo[neta],hTrainTruecn[neta],hTrainRecocn[neta],hMatcn[neta],
				0,lumi,neta,responseEtacn[neta]);
		generateEvents(eventsGenerate,fresup[neta],hnlo[neta],hTrainTrueup[neta],hTrainRecoup[neta],hMatup[neta],
				0,lumi,neta,responseEtaup[neta]);
		generateEvents(eventsGenerate,fresdn[neta],hnlo[neta],hTrainTruedn[neta],hTrainRecodn[neta],hMatdn[neta],
				0,lumi,neta,responseEtadn[neta]);
		/////datat resolution uncertainty
		generateEvents(eventsGenerate,fresdtcn[neta],hnlo[neta],hTrainTruecn[neta],hTrainRecocn[neta],hMatcn[neta],
				0,lumi,neta,responseEtadtcn[neta]);
		generateEvents(eventsGenerate,fresdtup[neta],hnlo[neta],hTrainTruedn[neta],hTrainRecodn[neta],hMatdn[neta],
				0,lumi,neta,responseEtadtup[neta]);
		generateEvents(eventsGenerate,fresdtdn[neta],hnlo[neta],hTrainTruedn[neta],hTrainRecodn[neta],hMatdn[neta],
				0,lumi,neta,responseEtadtdn[neta]);
	}
	////vary fit parameters 
	RooUnfoldResponse* responsePar[6][2];///neta, nparameter variation
	for(unsigned neta(0);neta<6;++neta) {
		TH1D* htmp = (TH1D*) hTrainRecodn[neta]->Clone("htmp");
		TH2D* htmpMat = (TH2D*) hMatdn[neta]->Clone("htmpMat");
		///+- 1 sigma variation
		double sg[3] = {
			frescn[neta]->GetParError(0),
			frescn[neta]->GetParError(1),
			frescn[neta]->GetParError(2)
		};
		///up variation
		sprintf(name,"response_eta%d_par%d",neta,0);
		responsePar[neta][0] = new RooUnfoldResponse (name,name);
		responsePar[neta][0]->Setup(hTrainRecodn[neta],hTrainTruedn[neta]);

		double p0(frescn[neta]->GetParameter(0)+sg[0]);
		double p1(frescn[neta]->GetParameter(1)+sg[1]);
		double p2(frescn[neta]->GetParameter(2)+sg[2]);
		//      cout<<"Parameters "<<frescn[neta]->GetName()<<"\n"
		//<<frescn[neta]->GetParameter(0)<<" +- "<<sg[0]<<"\n"
		//<<frescn[neta]->GetParameter(1)<<" +- "<<sg[1]<<"\n"
		//<<frescn[neta]->GetParameter(2)<<" +- "<<sg[2]<<"\n"
		//<<endl;
		frescn[neta]->SetParameter(0,p0);
		frescn[neta]->SetParameter(1,p1);
		frescn[neta]->SetParameter(2,p2);
		generateEvents(eventsGenerate,frescn[neta],hnlo[neta],htmp,htmp,htmpMat,
				0,lumi,neta,responsePar[neta][0]);
		///down variation
		sprintf(name,"response_eta%d_par%d",neta,1);
		responsePar[neta][1] = new RooUnfoldResponse (name,name);
		responsePar[neta][1]->Setup(hTrainRecodn[neta],hTrainTruedn[neta]);

		p0 = (frescn[neta]->GetParameter(0)-sg[0]);
		p1 = (frescn[neta]->GetParameter(1)-sg[1]);
		p2 = (frescn[neta]->GetParameter(2)-sg[2]);
		frescn[neta]->SetParameter(0,p0);
		frescn[neta]->SetParameter(1,p1);
		frescn[neta]->SetParameter(2,p2);
		generateEvents(eventsGenerate,frescn[neta],hnlo[neta],htmp,htmp,htmpMat,
				0,lumi,neta,responsePar[neta][1]);
	}
	///// Unfold Data ------------------------------------------------------------
	cout << "================================ UNFOLD ================================" << endl;
	TCanvas* cptSpecCl[6]; TCanvas* cUnfErr[6];
	TCanvas* cpar[6];
	//	TFile* fDat = new TFile((indir+"/AnalyzerIncHistosData276TeV_8tevcorr.root").c_str(),"read"); 
	TFile* fDat = new TFile((indir+"/AnalyzerIncHistosMC_winter14.root").c_str(),"read"); 
	TH1D* hDat[6];
	TH1D* hUnfDatacn[6];	TH1D* hUnfDataup[6];	TH1D* hUnfDatadn[6];
	TH1D* hUnfUnc[6];//unfolding uncertainty
	TH1D* hUnfUncdt[6];//unfolding uncertainty data smearing
	TH1D* hParUnc[6];//parameters uncertainty
	for(unsigned i(0);i<6;++i) {
		ostringstream oss;
		oss<<"hptUnf_eta_"<<ybins[i]<<"_"<<ybins[i+1];
		TH1D*  hDat_fullrange  = (TH1D*)fDat->Get(oss.str().c_str());
		TH1D*  hGen  = (TH1D*)hnlo[i];
		///set the histograms with the corect range
		sprintf(name,"hDat_%d",i);
		hDat[i] = new TH1D(name,name,nx[i],&x[i][0]);
		spFillHisto(hDat_fullrange,hDat[i]);

		TH1D*  hUnfcn; TH1D*  hUnfup;TH1D*  hUnfdn;
		RooUnfoldBayes    unfoldcn (responseEtacn[i], hDat[i], 4);
		RooUnfoldBayes    unfoldup (responseEtaup[i], hDat[i], 4);
		RooUnfoldBayes    unfolddn (responseEtadn[i], hDat[i], 4);
		hUnfcn= (TH1D*) unfoldcn.Hreco();
		hUnfup= (TH1D*) hUnfcn->Clone("hUnfup_copycn");
		hUnfdn= (TH1D*) hUnfcn->Clone("hUnfdn_copycn");
		//		hUnfup= (TH1D*) unfoldup.Hreco(); //disable +-10% smearng
		//		hUnfdn= (TH1D*) unfolddn.Hreco(); //disable +-10% smearng
		//////draw histogram----------------------------------------
		sprintf(name,"c_eta%d",i);
		cptSpecCl[i] = new TCanvas(name,name,800,600);
		sprintf(name,"hUnfUnc_eta%d",i);
		hUnfUnc[i] = (TH1D*)hUnfcn->Clone(name);
		getUnfoldUncertainty(hUnfcn,hUnfup,hUnfdn,hUnfUnc[i]);
		hUnfUnc[i]->GetXaxis()->SetRangeUser(74,maxPt[i]-1);
		hUnfUnc[i]->GetYaxis()->SetRangeUser(-0.4,0.6);
		hUnfUnc[i]->SetStats(kFALSE);
		hUnfUnc[i]->Draw("hist");

		////unfold data res uncertainty----------------------------
		TH1D*  hUnfdtcn; TH1D*  hUnfdtup;TH1D*  hUnfdtdn;
		RooUnfoldBayes    unfolddtcn (responseEtadtcn[i], hDat[i], 4);
		RooUnfoldBayes    unfolddtup (responseEtadtup[i], hDat[i], 4);
		RooUnfoldBayes    unfolddtdn (responseEtadtdn[i], hDat[i], 4);
		hUnfdtcn= (TH1D*) unfolddtcn.Hreco();
		hUnfdtup= (TH1D*) unfolddtup.Hreco();
		hUnfdtdn= (TH1D*) unfolddtdn.Hreco();

		//data unfold uncert
		sprintf(name,"hUnfUncdt_eta%d",i);
		hUnfUncdt[i] = (TH1D*)hUnfcn->Clone(name);
		getUnfoldUncertainty(hUnfdtcn,hUnfdtup,hUnfdtdn,hUnfUncdt[i]);
		hUnfUncdt[i]->GetXaxis()->SetRangeUser(74,maxPt[i]-1);
		hUnfUncdt[i]->GetYaxis()->SetRangeUser(-0.4,0.6);
		hUnfUncdt[i]->SetStats(kFALSE);
		hUnfUncdt[i]->Draw("same hist");


		//determine parameter uncertainty
		TH1D* hvr[2];
		RooUnfoldBayes    unfoldpar0 (responsePar[i][0], hDat[i], 4);
		hvr[0] = (TH1D*)unfoldpar0.Hreco();
		RooUnfoldBayes    unfoldpar1 (responsePar[i][1], hDat[i], 4);
		hvr[1] = (TH1D*)unfoldpar1.Hreco();
		sprintf(name,"hParUnc_eta%d",i);
		hParUnc[i] = (TH1D*)hUnfcn->Clone(name);
		getParUncertainty(hUnfcn,hvr,hParUnc[i],2);
		//////draw histogram
		sprintf(name,"cpar_eta%d",i);
		cpar[i] = new TCanvas(name,name,800,600);
		hParUnc[i]->GetXaxis()->SetRangeUser(74,maxPt[i]-1);
		hParUnc[i]->GetYaxis()->SetRangeUser(-0.4,0.6);
		hParUnc[i]->SetStats(kFALSE);
		hParUnc[i]->Draw("hist");
	}

	///Write Histos and Response object to file
	TFile* outfile = new TFile("UnfoldUncertainties.root","recreate");
	for(unsigned i(0);i<nybins-1;++i)
	{
		hUnfUnc[i]->Write();
		hUnfUncdt[i]->Write();
		hParUnc[i]->Write();
		//		hTrainRecoup[i]->Write();
		//		hTrainRecodn[i]->Write();
		//		hTrainRecocn[i]->Write();
		//		hUnfDatacn[i]->Write();
		//		hUnfDataup[i]->Write();
		//		hUnfDatadn[i]->Write();
	}
	outfile->Close();

}


void unfoldUnc() {
	// Data Members
	TH1D*              hTrainTrue[10];
	TH1D*              hTrainReco[10];
	TH2D               *hMat[10];
	RooUnfoldResponse* response;
	RooUnfold*         unfold;

	//  bool saveplots(true);
	//  if(saveplots) 	gROOT->SetBatch(kTRUE);
	const string indir("/uscms/home/tlibeiro/inclusivexsec/CMSSW_5_3_14/src/JetAnalyzer/JetAnalyzer/test");
	const char* nloDir("/uscms_data/d3/tlibeiro/inclusivexsec/Generators/generateIncJets/inputXsecs");
	const char* fitDir("./");
	sprintf(name,"%s/ResolutionFits.root",fitDir);
	TFile* resFile = new TFile(name,"read");
	TF1*  frescn[6];//resfits mc
	TF1*  fresup[6];//resfits
	TF1*  fresdn[6];//resfits
	TF1*  fresdtcn[6];//resfits data
	TF1*  fresdtup[6];//resfits
	TF1*  fresdtdn[6];//resfits
	const unsigned eventsGenerate(1e6);

	cout << "==================================== TRAIN ====================================" << endl;
	//for mc fit
	RooUnfoldResponse* responseEtacn[nybins-1]; 
	RooUnfoldResponse* responseEtaup[nybins-1]; 
	RooUnfoldResponse* responseEtadn[nybins-1]; 
	//for data fit
	RooUnfoldResponse* responseEtadtcn[nybins-1]; 
	RooUnfoldResponse* responseEtadtup[nybins-1]; 
	RooUnfoldResponse* responseEtadtdn[nybins-1]; 
	TH1D* hnlo[6];//nlo spectrum
	TH1D* hTrainTruecn[6];  TH1D* hTrainRecocn[6];  TH2D* hMatcn[6];
	TH1D* hTrainTrueup[6];  TH1D* hTrainRecoup[6];  TH2D* hMatup[6];
	TH1D* hTrainTruedn[6];  TH1D* hTrainRecodn[6];  TH2D* hMatdn[6];
	///Initialize Histos and Response class//////////////////////////////////////
	for(unsigned i(0);i<6;++i) {
		char infile[500]; char title[500];
		sprintf(infile,"%s/ct10y%i.cppInput",nloDir,i);
		sprintf(title ,"nloy%i",i);
		hnlo[i] = makeNLOSpectrum(infile,title);
		//fres for smearing        
		sprintf(name,"fresfit_eta%d",i);
		frescn[i] = (TF1*)resFile->Get(name);
		sprintf(name,"fresfit_up_eta%d",i);
		fresup[i] = (TF1*)resFile->Get(name);
		sprintf(name,"fresfit_dn_eta%d",i);
		fresdn[i] = (TF1*)resFile->Get(name);
		//fres data for smearing        
		sprintf(name,"fresfitdt_eta%d",i);
		fresdtcn[i] = (TF1*)resFile->Get(name);
		sprintf(name,"fresfitdt_up_eta%d",i);
		fresdtup[i] = (TF1*)resFile->Get(name);
		sprintf(name,"fresfitdt_dn_eta%d",i);
		fresdtdn[i] = (TF1*)resFile->Get(name);
		/////training histos
		unsigned ngen(n1x[i]), nrec(nx[i]);
		const double *xgen = x1[i]; //x bins gen
		const double *xrec = x[i];// x bins reco
		//Data   -------------------------
		/////
		sprintf(name,"resmatcn%d_eta%d",0,i);
		hMatcn[i]= new TH2D (name,name,ngen,xgen,ngen,xgen);
		//
		sprintf(name,"traintruecn%d_eta%d",0,i);
		hTrainTruecn[i]= new TH1D (name,name,ngen,xgen);
		hTrainTruecn[i]->SetLineColor(kBlue);
		//
		sprintf(name,"trainmeascn%d_eta%d",0,i);
		hTrainRecocn[i]= new TH1D (name,name,nrec,xrec);
		hTrainRecocn[i]->SetLineColor(kRed);
		//
		sprintf(name,"responsecn_eta%d",i);
		responseEtacn[i] = new RooUnfoldResponse (name,name);
		responseEtacn[i]->Setup(hTrainRecocn[i],hTrainTruecn[i]);
		sprintf(name,"responsedtcn_eta%d",i);
		responseEtadtcn[i] = new RooUnfoldResponse (name,name);
		responseEtadtcn[i]->Setup(hTrainRecocn[i],hTrainTruecn[i]);
		//up variation----------------------------
		sprintf(name,"resmatup%d_eta%d",0,i);
		hMatup[i]= new TH2D (name,name,ngen,xgen,ngen,xgen);
		//
		sprintf(name,"traintrueup%d_eta%d",0,i);
		hTrainTrueup[i]= new TH1D (name,name,ngen,xgen);
		hTrainTrueup[i]->SetLineColor(kBlue);
		//
		sprintf(name,"trainmeasup%d_eta%d",0,i);
		hTrainRecoup[i]= new TH1D (name,name,nrec,xrec);
		hTrainRecoup[i]->SetLineColor(kRed);
		//
		sprintf(name,"responseup_eta%d",i);
		responseEtaup[i] = new RooUnfoldResponse (name,name);
		responseEtaup[i]->Setup(hTrainRecoup[i],hTrainTrueup[i]);
		sprintf(name,"responsedtup_eta%d",i);
		responseEtadtup[i] = new RooUnfoldResponse (name,name);
		responseEtadtup[i]->Setup(hTrainRecoup[i],hTrainTrueup[i]);
		//down variation----------------------------
		sprintf(name,"resmatdn%d_eta%d",0,i);
		hMatdn[i]= new TH2D (name,name,ngen,xgen,ngen,xgen);
		//
		sprintf(name,"traintruedn%d_eta%d",0,i);
		hTrainTruedn[i]= new TH1D (name,name,ngen,xgen);
		hTrainTruedn[i]->SetLineColor(kBlue);
		//
		sprintf(name,"trainmeasdn%d_eta%d",0,i);
		hTrainRecodn[i]= new TH1D (name,name,nrec,xrec);
		hTrainRecodn[i]->SetLineColor(kRed);
		//
		sprintf(name,"responsedn_eta%d",i);
		responseEtadn[i] = new RooUnfoldResponse (name,name);
		responseEtadn[i]->Setup(hTrainRecodn[i],hTrainTruedn[i]);
		sprintf(name,"responsedtdn_eta%d",i);
		responseEtadtdn[i] = new RooUnfoldResponse (name,name);
		responseEtadtdn[i]->Setup(hTrainRecodn[i],hTrainTruedn[i]);
	}
	///generate toy events and smear with data res /////////////////////////////////////////
	for(unsigned neta(0);neta<6;++neta) {
		generateEvents(eventsGenerate,frescn[neta],hnlo[neta],hTrainTruecn[neta],hTrainRecocn[neta],hMatcn[neta],
				0,lumi,neta,responseEtacn[neta]);
		generateEvents(eventsGenerate,fresup[neta],hnlo[neta],hTrainTrueup[neta],hTrainRecoup[neta],hMatup[neta],
				0,lumi,neta,responseEtaup[neta]);
		generateEvents(eventsGenerate,fresdn[neta],hnlo[neta],hTrainTruedn[neta],hTrainRecodn[neta],hMatdn[neta],
				0,lumi,neta,responseEtadn[neta]);
		/////datat resolution uncertainty
		generateEvents(eventsGenerate,fresdtcn[neta],hnlo[neta],hTrainTruecn[neta],hTrainRecocn[neta],hMatcn[neta],
				0,lumi,neta,responseEtadtcn[neta]);
		generateEvents(eventsGenerate,fresdtup[neta],hnlo[neta],hTrainTruedn[neta],hTrainRecodn[neta],hMatdn[neta],
				0,lumi,neta,responseEtadtup[neta]);
		generateEvents(eventsGenerate,fresdtdn[neta],hnlo[neta],hTrainTruedn[neta],hTrainRecodn[neta],hMatdn[neta],
				0,lumi,neta,responseEtadtdn[neta]);
	}
	////vary fit parameters 
	RooUnfoldResponse* responsePar[6][16];///neta, nparameter variation
	double nvar0[8] = {0.05,0.00,0.05,0.05,0.05,0.00,0.00,0.00};
	double nvar1[8] = {0.05,0.05,0.00,0.05,0.00,0.00,0.05,0.00};
	double nvar2[8] = {0.05,0.05,0.05,0.00,0.00,0.05,0.00,0.00};
	for(unsigned neta(0);neta<6;++neta) {
		TH1D* htmp = (TH1D*) hTrainRecodn[neta]->Clone("htmp");
		TH2D* htmpMat = (TH2D*) hMatdn[neta]->Clone("htmpMat");
		for(unsigned par (0);par<8;++par ) {
			///up variation
			sprintf(name,"response_eta%d_par%d",neta,par);
			responsePar[neta][par] = new RooUnfoldResponse (name,name);
			responsePar[neta][par]->Setup(hTrainRecodn[neta],hTrainTruedn[neta]);

			double p0(frescn[neta]->GetParameter(0)*(1+nvar0[par]));
			double p1(frescn[neta]->GetParameter(1)*(1+nvar1[par]));
			double p2(frescn[neta]->GetParameter(2)*(1+nvar2[par]));
			frescn[neta]->SetParameter(0,p0);
			frescn[neta]->SetParameter(1,p1);
			frescn[neta]->SetParameter(2,p2);
			generateEvents(eventsGenerate,frescn[neta],hnlo[neta],htmp,htmp,htmpMat,
					0,lumi,neta,responsePar[neta][par]);
			///down variation
			sprintf(name,"response_eta%d_par%d",neta,par+8);
			responsePar[neta][par+8] = new RooUnfoldResponse (name,name);
			responsePar[neta][par+8]->Setup(hTrainRecodn[neta],hTrainTruedn[neta]);

			p0 = (frescn[neta]->GetParameter(0)*(1-nvar0[par]));
			p1 = (frescn[neta]->GetParameter(1)*(1-nvar1[par]));
			p2 = (frescn[neta]->GetParameter(2)*(1-nvar2[par]));
			frescn[neta]->SetParameter(0,p0);
			frescn[neta]->SetParameter(1,p1);
			frescn[neta]->SetParameter(2,p2);
			generateEvents(eventsGenerate,frescn[neta],hnlo[neta],htmp,htmp,htmpMat,
					0,lumi,neta,responsePar[neta][par+8]);
		}
	}
	///// Unfold Data ------------------------------------------------------------
	cout << "================================ UNFOLD ================================" << endl;
	TCanvas* cptSpecCl[6]; TCanvas* cUnfErr[6];
	TCanvas* cpar[6];
	//	TFile* fDat = new TFile((indir+"/AnalyzerIncHistosData276TeV_8tevcorr.root").c_str(),"read"); 
	TFile* fDat = new TFile((indir+"/AnalyzerIncHistosMC_winter14.root").c_str(),"read"); 
	TH1D* hDat[6];
	TH1D* hUnfDatacn[6];	TH1D* hUnfDataup[6];	TH1D* hUnfDatadn[6];
	TH1D* hUnfUnc[6];//unfolding uncertainty
	TH1D* hUnfUncdt[6];//unfolding uncertainty data smearing
	TH1D* hParUnc[6];//parameters uncertainty
	for(unsigned i(0);i<6;++i) {
		ostringstream oss;
		oss<<"hptUnf_eta_"<<ybins[i]<<"_"<<ybins[i+1];
		TH1D*  hDat_fullrange  = (TH1D*)fDat->Get(oss.str().c_str());
		TH1D*  hGen  = (TH1D*)hnlo[i];
		///set the histograms with the corect range
		sprintf(name,"hDat_%d",i);
		hDat[i] = new TH1D(name,name,nx[i],&x[i][0]);
		spFillHisto(hDat_fullrange,hDat[i]);

		TH1D*  hUnfcn; TH1D*  hUnfup;TH1D*  hUnfdn;
		RooUnfoldBayes    unfoldcn (responseEtacn[i], hDat[i], 4);
		RooUnfoldBayes    unfoldup (responseEtaup[i], hDat[i], 4);
		RooUnfoldBayes    unfolddn (responseEtadn[i], hDat[i], 4);
		hUnfcn= (TH1D*) unfoldcn.Hreco();
		hUnfup= (TH1D*) hUnfcn->Clone("hUnfup_copycn");
		hUnfdn= (TH1D*) hUnfcn->Clone("hUnfdn_copycn");
		//		hUnfup= (TH1D*) unfoldup.Hreco(); //disable +-10% smearng
		//		hUnfdn= (TH1D*) unfolddn.Hreco(); //disable +-10% smearng
		//////draw histogram----------------------------------------
		sprintf(name,"c_eta%d",i);
		cptSpecCl[i] = new TCanvas(name,name,800,600);
		sprintf(name,"hUnfUnc_eta%d",i);
		hUnfUnc[i] = (TH1D*)hUnfcn->Clone(name);
		getUnfoldUncertainty(hUnfcn,hUnfup,hUnfdn,hUnfUnc[i]);
		hUnfUnc[i]->GetXaxis()->SetRangeUser(74,maxPt[i]-1);
		hUnfUnc[i]->GetYaxis()->SetRangeUser(-0.4,0.6);
		hUnfUnc[i]->SetStats(kFALSE);
		hUnfUnc[i]->Draw("hist");

		////unfold data res uncertainty----------------------------
		TH1D*  hUnfdtcn; TH1D*  hUnfdtup;TH1D*  hUnfdtdn;
		RooUnfoldBayes    unfolddtcn (responseEtadtcn[i], hDat[i], 4);
		RooUnfoldBayes    unfolddtup (responseEtadtup[i], hDat[i], 4);
		RooUnfoldBayes    unfolddtdn (responseEtadtdn[i], hDat[i], 4);
		hUnfdtcn= (TH1D*) unfolddtcn.Hreco();
		hUnfdtup= (TH1D*) unfolddtup.Hreco();
		hUnfdtdn= (TH1D*) unfolddtdn.Hreco();

		//data unfold uncert
		sprintf(name,"hUnfUncdt_eta%d",i);
		hUnfUncdt[i] = (TH1D*)hUnfcn->Clone(name);
		getUnfoldUncertainty(hUnfdtcn,hUnfdtup,hUnfdtdn,hUnfUncdt[i]);
		hUnfUncdt[i]->GetXaxis()->SetRangeUser(74,maxPt[i]-1);
		hUnfUncdt[i]->GetYaxis()->SetRangeUser(-0.4,0.6);
		hUnfUncdt[i]->SetStats(kFALSE);
		hUnfUncdt[i]->Draw("same hist");


		//determine parameter uncertainty
		TH1D* hvr[16];
		for(unsigned par(0);par<16;++par)
		{
			RooUnfoldBayes    unfoldpar (responsePar[i][par], hDat[i], 4);
			hvr[par] = (TH1D*)unfoldpar.Hreco();
		}
		sprintf(name,"hParUnc_eta%d",i);
		hParUnc[i] = (TH1D*)hUnfcn->Clone(name);
		getParUncertainty(hUnfcn,hvr,hParUnc[i],16);
		//////draw histogram
		sprintf(name,"cpar_eta%d",i);
		cpar[i] = new TCanvas(name,name,800,600);
		hParUnc[i]->GetXaxis()->SetRangeUser(74,maxPt[i]-1);
		hParUnc[i]->GetYaxis()->SetRangeUser(-0.4,0.6);
		hParUnc[i]->SetStats(kFALSE);
		hParUnc[i]->Draw("hist");
	}

	///Write Histos and Response object to file
	TFile* outfile = new TFile("UnfoldUncertainties.root","recreate");
	for(unsigned i(0);i<nybins-1;++i)
	{
		hUnfUnc[i]->Write();
		hUnfUncdt[i]->Write();
		hParUnc[i]->Write();
		//		hTrainRecoup[i]->Write();
		//		hTrainRecodn[i]->Write();
		//		hTrainRecocn[i]->Write();
		//		hUnfDatacn[i]->Write();
		//		hUnfDataup[i]->Write();
		//		hUnfDatadn[i]->Write();
	}
	outfile->Close();

}


void unfoldData() {
	// Data Members
	TH1D*              hTrainTrue[10];
	TH1D*              hTrainReco[10];
	TH2D               *hMat[10];
	RooUnfoldResponse* response;
	RooUnfold*         unfold;

	bool saveplots(false);
	if(saveplots) 	gROOT->SetBatch(kTRUE);
	const string indir("/uscms/home/tlibeiro/inclusivexsec/CMSSW_5_3_14/src/JetAnalyzer/JetAnalyzer/test");
	const char* nloDir("/uscms_data/d3/tlibeiro/inclusivexsec/Generators/generateIncJets/inputXsecs");
	const char* fitDir("./");
	sprintf(name,"%s/ResolutionFits_GaussFit.root",fitDir);
	TFile* resFile = new TFile(name,"read");
	TF1*  fres[6];//resfits
	const unsigned eventsGenerate(1e6);

	response= new RooUnfoldResponse ("response","Test");
	cout << "==================================== TRAIN ====================================" << endl;

	//for data
	RooUnfoldResponse* responseEtaDt[nybins-1]; //response eta with data res
	TH1D* hnlo[6];//nlo spectrum
	TH1D* hTrainTrueDt[6];// for comparison
	TH1D* hTrainRecoDt[6];// for comparison
	TH2D* hMatDt[6];// for comparison
	///Initialize Histos and Response class//////////////////////////////////////
	for(unsigned i(0);i<6;++i) {
		char infile[500]; char title[500];
		//    sprintf(infile,"%s/nnpdf2.1y%i.cppInput",nloDir,i);
		//    sprintf(title ,"nloy%i",i);
		sprintf(infile,"%s/ct10y%i.cppInput",nloDir,i);
		sprintf(title ,"nloy%i",i);
		hnlo[i] = makeNLOSpectrum(infile,title);
		//fres for smearing        
		sprintf(name,"fresfitdt_eta%d",i);
		fres[i] = (TF1*)resFile->Get(name);
		/////training histos
		unsigned ngen(n1x[i]), nrec(nx[i]);
		const double *xgen = x1[i]; //x bins gen
		const double *xrec = x[i];// x bins reco
		//Data   -------------------------
		/////
		sprintf(name,"resmatDt%d_eta%d",0,i);
		hMatDt[i]= new TH2D (name,name,ngen,xgen,ngen,xgen);
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
	///generate toy events and smear with data res /////////////////////////////////////////
	for(unsigned neta(0);neta<6;++neta) {
		generateEvents(eventsGenerate,fres[neta],hnlo[neta],hTrainTrueDt[neta],hTrainRecoDt[neta],hMatDt[neta],
				0,lumi,neta,responseEtaDt[neta]);
	}
	///// Unfold Data ------------------------------------------------------------
	cout << "================================ UNFOLD ================================" << endl;
	TCanvas* cptSpecCl[6]; TCanvas* cUnfErr[6];
	TFile* fDat = new TFile((indir+"/AnalyzerIncHistosData276TeV_winter14.root").c_str(),"read"); 
	TFile* fMC = new TFile((indir+"/AnalyzerIncHistosMC_winter14.root").c_str(),"read"); 
	//	TFile* fDat = new TFile((indir+"/AnalyzerIncHistosData276TeV_ft53v10.root").c_str(),"read"); 
	TFile* npfile= new TFile("npCorrections.root","read");
	TFile* jteff = new TFile("JetEfficiency.root","read");
	TCanvas* cclo_tBay[6];
	TH1D* hDat[6];
	TH1D* hMC[6];
	TH1D* hUnfData[6];
	for(unsigned i(0);i<6;++i) {
		ostringstream oss;
		oss<<"hptUnf_eta_"<<ybins[i]<<"_"<<ybins[i+1];
		TH1D*  hDat_fullrange  = (TH1D*)fDat->Get(oss.str().c_str());
		TH1D*  hMC_fullrange  = (TH1D*)fMC->Get(oss.str().c_str());
		TH1D*  hGen  = (TH1D*)hnlo[i];
		//closure for three unfolding techniques
		TH1D*  hclo1 = (TH1D*)hTrainTrueDt[i]->Clone("hclo1"); hclo1->Reset();
		TH1D*  hclo2 = (TH1D*)hTrainTrueDt[i]->Clone("hclo2"); hclo2->Reset();
		TH1D*  hclo3 = (TH1D*)hTrainTrueDt[i]->Clone("hclo3"); hclo3->Reset();
		///set the histograms with the corect range
		sprintf(name,"hDat_%d",i);
		hDat[i] = new TH1D(name,name,n2x[i],&x2[i][0]);
		spFillHisto(hDat_fullrange,hDat[i]);
		sprintf(name,"hMC_%d",i);
		hMC[i] = new TH1D(name,name,n2x[i],&x2[i][0]);
		spFillHisto(hMC_fullrange,hMC[i]);
		///correct for jet efficiency
		sprintf(name,"heff_eta%d_fixed",i);
		TH1D* eff = (TH1D*)jteff->Get(name);
		TH1D* heff = (TH1D*)hDat[i]->Clone("hjeteff");
		spFillHisto(eff,heff);
		hDat[i]->Divide(heff);
		hDat[i]->Scale(1/0.99);//hlt efficiency

		TH1D*  hUnf[3];
		RooUnfoldBayes    unfold (responseEtaDt[i], hDat[i], 4);
		hUnf[0]= (TH1D*) unfold.Hreco();
		TH1D*  hUnfMC[3]; // for unfolding error and unfolding test
		RooUnfoldBayes    unfoldMC (responseEtaDt[i], hMC[i], 4);
		hUnfMC[0]= (TH1D*) unfoldMC.Hreco();

		//		RooUnfoldSvd      unfold2 (responseEtaDt[i], hDat[i], 3);  
		//		RooUnfoldBinByBin unfold3 (responseEtaDt[i], hDat[i]);
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
		///scale unfolded histos
		ScaleHistByBW(hUnf[0]);
		hUnf[0]->Scale(1.0/lumi);
		hclo1->Divide(hUnf[0],htmp);
		//errors before and after unfolding
		TH1D *hdaterr, *hunferr;
		TH1D* hDattmp = (TH1D*)hDat[i]->Clone("hDattmp");
		ScaleHistByBW(hDattmp);
		hDattmp->Scale(1/lumi);
		//////////make error histos
		sprintf(name,"%s_err",hDat[i]->GetName());
		hdaterr = (TH1D*)hDat[i]->Clone(name);
		sprintf(name,"%s_err",hUnf[0]->GetName());
		hunferr = (TH1D*)hUnf[0]->Clone(name);
		histoFromErr(hDattmp,hdaterr);
		histoFromErr(hUnf[0],hunferr);
		sprintf(name,"cUnfErr_eta%d",i);
		cUnfErr[i] = new TCanvas(name,name,600,600);
		drawErrors(cUnfErr[i],hdaterr,hunferr,i);
		sprintf(name,"plots/%s.eps",cUnfErr[i]->GetName());
		if(saveplots)
			cUnfErr[i]->SaveAs(name);
		//////draw histogram
		sprintf(name,"c_hUnf0_eta%d",i);
		cptSpecCl[i] = new TCanvas(name,name,800,600);
		cptSpecCl[i]->SetLogx(); cptSpecCl[i]->SetLogy();
		drawSpecCl(cptSpecCl[i],hUnf[0],hGen,hclo1,i);
		///////print latex 
		sprintf(name,"%.2f<|y|<%.2f",ybins[i],ybins[i+1]);
		TLatex latex;latex.SetTextSize(0.05);latex.SetTextFont(42);
		latex.DrawTextNDC(.5,.7,name);
		sprintf(name,"plots/hUnfBay_eta%d.eps",i);
		if(saveplots)
			cptSpecCl[i]->SaveAs(name);
		sprintf(name,"hUnf_eta%d",i);
		hUnfData[i] = (TH1D*)hUnf[0]->Clone(name);
	}

	///Write Histos and Response object to file
	TFile* outfile = new TFile("UnfoldJetsData_GaussFit.root","recreate");
	for(unsigned i(0);i<nybins-1;++i)
	{
		cptSpecCl[i]->Write();
		hUnfData[i]->Write();
		hnlo[i]->Write();
		hDat[i]->Write();
	}
	outfile->Close();

}

void unfoldtrainTest2() { //testing for unfolding factors, for regular use run unfoldtrainTest()
	// Data Members
	TH1D*              hTrainTrue[10];
	TH1D*              hTrainReco[10];
	TH2D               *hMat[10];
	RooUnfoldResponse* response;
	RooUnfold*         unfold;

	bool saveplots(true);
	if(saveplots) 	gROOT->SetBatch(kTRUE);
//	gROOT->SetBatch(kFALSE);
	const string indir("/uscms/home/tlibeiro/inclusivexsec/CMSSW_5_3_14/src/JetAnalyzer/JetAnalyzer/test");
	const char* nloDir("/uscms_data/d3/tlibeiro/inclusivexsec/Generators/generateIncJets/inputXsecs");
	const char* fitDir("./");
	sprintf(name,"%s/ResolutionFits.root",fitDir);
	TFile* resFile = new TFile(name,"read");
	TF1*  fres[6];//resfits
	const unsigned eventsGenerate(1e6);

	response= new RooUnfoldResponse ("response","Test");
	cout << "==================================== TRAIN ====================================" << endl;
	RooUnfoldResponse* responseEta[nybins-1];
	TH1D* hnlo[6];//nlo spectrum
	TH1D* hnlot[6][2];//nlo spectrum for comparison
	TH1D* hTrainTruet[6][3];// for comparison
	TH1D* hTrainRecot[6][3];// for comparison
	TH2D* hMatt[6][3];// for comparison
	///Initialize Histos and Response class//////////////////////////////////////
	for(unsigned i(0);i<6;++i) {
		char infile[500]; char title[500];
		//    sprintf(infile,"%s/nnpdf2.1y%i.cppInput",nloDir,i);
		//    sprintf(title ,"nloy%i",i);
		sprintf(infile,"%s/ct10y%i.cppInput",nloDir,i);
		sprintf(title ,"nloy%i",i);
		hnlo[i] = makeNLOSpectrum(infile,title);
		//comparing unfolding with different spectra
		sprintf(infile,"%s/nnpdf2.1y%i.cppInput",nloDir,i);
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
		const double *xrec = x1[i];// x bins reco
		sprintf(name,"traintrue_testUnf_eta%d",i);
		hTrainTrue[i]= new TH1D (name,name,ngen,xgen);
		hTrainTrue[i]->SetLineColor(kBlue);
		//    cout<<" histo bins "<<hTrainTrue[i]->GetNbinsX()<<" "<<hTrainTrue[i]->GetBinLowEdge(1)<<" "
		//		<<hTrainTrue[i]->GetBinLowEdge(hTrainTrue[i]->GetNbinsX())+hTrainTrue[i]->GetBinWidth(hTrainTrue[i]->GetNbinsX())<<endl;
		//
		sprintf(name,"trainmeas_testUnf_eta%d",i);
		hTrainReco[i]= new TH1D (name,name,nrec,xrec);
		hTrainReco[i]->SetLineColor(kRed);
		sprintf(name,"resmat_testUnf_eta%d",i);
		hMat[i]= new TH2D (name,name,ngen,xgen,ngen,xgen);
		//
		sprintf(name,"response_testUnf_eta%d",i);
		responseEta[i] = new RooUnfoldResponse (name,name);
		responseEta[i]->Setup(hTrainReco[i],hTrainTrue[i]);
		/////for comparision
		sprintf(name,"resmatt%d_testUnf_eta%d",0,i);
		hMatt[i][0]= new TH2D (name,name,ngen,xgen,ngen,xgen);
		//
		sprintf(name,"traintruet%d_testUnf_eta%d",0,i);
		hTrainTruet[i][0]= new TH1D (name,name,ngen,xgen);
		hTrainTruet[i][0]->SetLineColor(kBlue);
		//
		sprintf(name,"trainmeast%d_testUnf_eta%d",0,i);
		hTrainRecot[i][0]= new TH1D (name,name,nrec,xrec);
		hTrainRecot[i][0]->SetLineColor(kRed);
		//
		sprintf(name,"resmatt%d_testUnf_eta%d",1,i);
		hMatt[i][1]= new TH2D (name,name,ngen,xgen,ngen,xgen);
		//
		sprintf(name,"traintruet%d_eta%d",1,i);
		hTrainTruet[i][1]= new TH1D (name,name,ngen,xgen);
		hTrainTruet[i][1]->SetLineColor(kBlue);
		//
		sprintf(name,"trainmeast%d_testUnf_eta%d",1,i);
		hTrainRecot[i][1]= new TH1D (name,name,nrec,xrec);
		hTrainRecot[i][1]->SetLineColor(kRed);
		//
		sprintf(name,"traintruet%d_testUnf_eta%d",2,i);
		hTrainTruet[i][2]= new TH1D (name,name,ngen,xgen);
		hTrainTruet[i][2]->SetLineColor(kBlue);
		//
		sprintf(name,"trainmeast%d_testUnf_eta%d",2,i);
		hTrainRecot[i][2]= new TH1D (name,name,nrec,xrec);
		hTrainRecot[i][2]->SetLineColor(kRed);
		//
		sprintf(name,"resmatt%d_testUnf_eta%d",2,i);
		hMatt[i][2]= new TH2D (name,name,ngen,xgen,ngen,xgen);

	}
	///generate toy events and smear/////////////////////////////////////////////
	for(unsigned neta(0);neta<6;++neta) {
		generateEventsTest(eventsGenerate,fres[neta],hnlo[neta],hTrainTrue[neta],hTrainReco[neta],
				hMat[neta],lumi,neta,0,responseEta[neta]);
		generateEventsTest(eventsGenerate,fres[neta],hnlot[neta][0],hTrainTruet[neta][0],
				hTrainRecot[neta][0],hMatt[neta][0],lumi,neta,0);
		generateEventsTest(eventsGenerate,fres[neta],hnlot[neta][1],hTrainTruet[neta][1],
				hTrainRecot[neta][1],hMatt[neta][1],lumi,neta,0);
	}
	///////Unfold/////////////////////////////////////////////////////////////////
	cout << "============== UNFOLD ===============" << endl;
	TCanvas* cptSpec[6];
	TH1D* hclo[6];     TCanvas* cclo[6];// unfolding closure
	TH1D* htom[6][2];  TCanvas* ctom[6];// unfolded true / mesaured spectrum
	TH1D* hclot[6][2]; TCanvas* cclot[6][2];// unfolding closure for comparision
	TLegend* leg[6];   TCanvas* cMat[6];///to draw matrix
	TH1D* hcloMC[6];   TCanvas* cMC[6];
	TCanvas* cptSpecCl[6];
	TCanvas* cClosCom[6][2]; //compare closures from different spectrum, unfolding routines
	TH1D* hclo2[6][2]; //for comparing unfolding routines
	//unfolded histograms //0-bayes, 1-svd, 2- bin by bin
	TH1D* hUnf[4];//unfolded
	TH1D* hUnft[3];//test pdfs

	TH1D* hrecMC[6]; TH1D* hgenMC[6]; ///for unfolding check
	TFile* fMC = new TFile((indir+"/AnalyzerIncHistosMC_winter14.root").c_str(),"read"); 
	TH2D* hcov[6]; // covariance matrix
	TCanvas* ccov[6];
	/// Check Unfolding/////////////////////////////-------------------------------
	for(unsigned i(0);i<6;++i) {

		///set the genMC and recMC
		sprintf(name,"hrecMC_eta%d",i);
		hrecMC[i] = (TH1D*)fMC->Get(name);
		sprintf(name,"hgenMC_eta%d",i);
		hgenMC[i] = (TH1D*)fMC->Get(name);

		TH1D*  hDat = hTrainReco[i];
		TH1D*  hGen = (TH1D*)hTrainTrue[i];

		RooUnfoldBayes   unfold1  (responseEta[i], hDat, 4);    // OR
		RooUnfoldSvd     unfold2  (responseEta[i], hDat, 4);   // OR
		RooUnfoldBinByBin unfold3 (responseEta[i], hDat);

		hUnf[0]= (TH1D*) unfold1.Hreco();
		//		hUnf[1]= (TH1D*) unfold2.Hreco();
		hUnf[2]= (TH1D*) unfold3.Hreco();
		///test Pyhtia spectrum ======
		hrecMC[i]->Reset();
    sprintf(name,"hrecPy_eta%d",i);
    TH1D* hrecPy = (TH1D*)hTrainReco[i]->Clone(name);
    hrecPy->Reset();
    sprintf(name,"hgenPy_eta%d",i);
    TH1D* hgenPy = (TH1D*)hTrainTrue[i]->Clone(name);
    spFillHisto(hgenMC[i],hgenPy);
		generateEventsMC(eventsGenerate,fres[i],hgenPy,hTrainTruet[i][2],
				hrecPy,hMatt[i][2],0,lumi,i);
		RooUnfoldBayes unfoldMC (responseEta[i], hrecPy,4);
		hUnf[3] = (TH1D*) unfoldMC.Hreco();

		//Unfolding with  different pdfs
		RooUnfoldBinByBin   unfoldt1_bb  (responseEta[i], hTrainRecot[i][0]);   
		TH1D* hUnft1_bb = (TH1D*) unfoldt1_bb.Hreco();
		RooUnfoldBinByBin   unfoldt2_bb  (responseEta[i], hTrainRecot[i][1]);   
		TH1D* hUnft2_bb = (TH1D*) unfoldt2_bb.Hreco();
		////svd 
		//    RooUnfoldSvd   unfoldt1  (responseEta[i], hTrainRecot[i][0], 4);   
		//		TH1D* hUnft1= (TH1D*) unfoldt1.Hreco();
		//    RooUnfoldSvd   unfoldt2  (responseEta[i], hTrainRecot[i][1], 4);   
		//		TH1D* hUnft2= (TH1D*) unfoldt2.Hreco();
		//// bin by bin 
		RooUnfoldBayes   unfoldt1  (responseEta[i], hTrainRecot[i][0], 3);   
		TH1D* hUnft1= (TH1D*) unfoldt1.Hreco();
		RooUnfoldBayes   unfoldt2  (responseEta[i], hTrainRecot[i][1], 3);   
		TH1D* hUnft2= (TH1D*) unfoldt2.Hreco();

		///unfolding with different pdfs
		char label[100];
		char label1[100], label2[100];
		sprintf(label1,"%s","nnpdf2.1");
		sprintf(label2,"%s","abm11");
		sprintf(label,"%s","CT10");
		drawClosure(cclo[i],hUnf[0],hGen,hclo[i],i,saveplots,label); 
		//		drawClosure(cclot[i][0],hUnft1,hTrainTruet[i][0],hclot[i][0],i,saveplots,label1); 
		//		drawClosure(cclot[i][1],hUnft2,hTrainTruet[i][1],hclot[i][0],i,saveplots,label2); 
		drawClosure(cMC[i],hUnf[3],hTrainTruet[i][2],hcloMC[i],i,saveplots,"Pythia"); //pythia test unfolding
		//draw on the same canvas
		drawClosure(cClosCom[i][0],i,saveplots,hUnft1,hTrainTruet[i][0],hclot[i][0],label1,
				hUnft2,hTrainTruet[i][1],hclot[i][0],label2);
		///unfolding with different methods
		TH1D* htmp = (TH1D*)hUnf[0]->Clone("htmp");
		spFillHisto(hTrainReco[i],htmp);
		drawClosure2(cClosCom[i][1],i,saveplots,hUnf[0],htmp,hclo2[i][0],"Bayes",
				hUnf[2],htmp,hclo2[i][1],"BinByBin");

		sprintf(name,"c_%s",hMat[i]->GetName());
		cMat[i] = new TCanvas(name,name,600,600);
		cMat[i]->SetLogy();cMat[i]->SetLogx();
		drawMatrix(hMat[i],i); 
		hMat[i]->Draw("colz");

		sprintf(name,"%.1f < |y| < %.1f",ybins[i],ybins[i+1]);
		TLatex latex;latex.SetTextSize(0.04);latex.SetTextFont(42);
		latex.DrawTextNDC(.3,.8,name);
		//draw covariance matrix-----------------
		sprintf(name,"ccov_testUnf_%d",i);
		ccov[i] = new TCanvas(name,name,600,600);
		sprintf(name,"hcov_%d",i);
		//error = 2
		const unsigned nrec(hTrainReco[i]->GetNbinsX());
		const double *xrec = hTrainReco[i]->GetXaxis()->GetXbins()->GetArray();// x bins reco
		hcov[i]= CorrelationHist (unfold1.Ereco((RooUnfold::ErrorTreatment)2),
				name, "Unfolded correlation matrix",
				xrec,nrec,i
				);
		ccov[i]->SetLogx(1); ccov[i]->SetLogy(1);
		hcov[i]->Draw("colz");
		sprintf(name,"%.1f < |y| < %.1f",ybins[i],ybins[i+1]);
		latex.DrawTextNDC(.3,.8,name);
		//end draw covariance matrix-----------------

		//unfolded, gen and reco spectrum
		cout<<"Drawing histos eta low edge "<<ybins[i]<<endl;
		drawSpec(cptSpec[i],hDat,hGen,hUnf[0],leg[i],saveplots,"");
		//save plots  
		if(saveplots) {
			cout<<cMat[i]->GetName()<<endl;
			sprintf(name,"plots/%s.eps",cMat[i]->GetName());
			cMat[i]->SaveAs(name);
			sprintf(name,"plots/%s.eps",ccov[i]->GetName());
			cout<<name<<" "<<ccov[i]<<endl;
			ccov[i]->SaveAs(name);
		}
	}

	//////
	///Write Histos and Response object to file
	TFile* outfile = new TFile("UnfoldJetsTrainingTest2.root","recreate");
	for(unsigned i(0);i<nybins-1;++i)
	{
		hTrainTrue[i]->Write();
		hTrainReco[i]->Write();
		hMat[i]->Write();
		hcov[i]->Write();
	}
	outfile->Close();
}



void unfoldtrainTest() {
	// Data Members
	TH1D*              hTrainTrue[10];
	TH1D*              hTrainReco[10];
	TH2D               *hMat[10];
	RooUnfoldResponse* response;
	RooUnfold*         unfold;

	bool saveplots(true);
	if(saveplots) 	gROOT->SetBatch(kTRUE);
	//	gROOT->SetBatch(kFALSE);
	const string indir("/uscms/home/tlibeiro/inclusivexsec/CMSSW_5_3_14/src/JetAnalyzer/JetAnalyzer/test");
	const char* nloDir("/uscms_data/d3/tlibeiro/inclusivexsec/Generators/generateIncJets/inputXsecs");
	const char* fitDir("./");
	sprintf(name,"%s/ResolutionFits.root",fitDir);
	TFile* resFile = new TFile(name,"read");
	TF1*  fres[6];//resfits
	const unsigned eventsGenerate(1e6);

	response= new RooUnfoldResponse ("response","Test");
	cout << "==================================== TRAIN ====================================" << endl;
	RooUnfoldResponse* responseEta[nybins-1];
	TH1D* hnlo[6];//nlo spectrum
	TH1D* hnlot[6][2];//nlo spectrum for comparison
	TH1D* hTrainTruet[6][3];// for comparison
	TH1D* hTrainRecot[6][3];// for comparison
	TH2D* hMatt[6][3];// for comparison
	///Initialize Histos and Response class//////////////////////////////////////
	for(unsigned i(0);i<6;++i) {
		char infile[500]; char title[500];
		//    sprintf(infile,"%s/nnpdf2.1y%i.cppInput",nloDir,i);
		//    sprintf(title ,"nloy%i",i);
		sprintf(infile,"%s/ct10y%i.cppInput",nloDir,i);
		sprintf(title ,"nloy%i",i);
		hnlo[i] = makeNLOSpectrum(infile,title);
		//comparing unfolding with different spectra
		sprintf(infile,"%s/nnpdf2.1y%i.cppInput",nloDir,i);
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
		unsigned ngen(n1x[i]), nrec(n1x[i]);
		const double *xgen = x1[i]; //x bins gen
		const double *xrec = x1[i];// x bins reco
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
		hMat[i]= new TH2D (name,name,ngen,xgen,ngen,xgen);
		//
		sprintf(name,"response_eta%d",i);
		responseEta[i] = new RooUnfoldResponse (name,name);
		responseEta[i]->Setup(hTrainReco[i],hTrainTrue[i]);
		/////for comparision
		sprintf(name,"resmatt%d_eta%d",0,i);
		hMatt[i][0]= new TH2D (name,name,ngen,xgen,ngen,xgen);
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
		hMatt[i][1]= new TH2D (name,name,ngen,xgen,ngen,xgen);
		//
		sprintf(name,"traintruet%d_eta%d",1,i);
		hTrainTruet[i][1]= new TH1D (name,name,ngen,xgen);
		hTrainTruet[i][1]->SetLineColor(kBlue);
		//
		sprintf(name,"trainmeast%d_eta%d",1,i);
		hTrainRecot[i][1]= new TH1D (name,name,nrec,xrec);
		hTrainRecot[i][1]->SetLineColor(kRed);
		//
		sprintf(name,"traintruet%d_eta%d",2,i);
		hTrainTruet[i][2]= new TH1D (name,name,ngen,xgen);
		hTrainTruet[i][2]->SetLineColor(kBlue);
		//
		sprintf(name,"trainmeast%d_eta%d",2,i);
		hTrainRecot[i][2]= new TH1D (name,name,nrec,xrec);
		hTrainRecot[i][2]->SetLineColor(kRed);
		//
		sprintf(name,"resmatt%d_eta%d",2,i);
		hMatt[i][2]= new TH2D (name,name,ngen,xgen,ngen,xgen);

	}
	///generate toy events and smear/////////////////////////////////////////////
	for(unsigned neta(0);neta<6;++neta) {
		generateEvents(eventsGenerate,fres[neta],hnlo[neta],hTrainTrue[neta],hTrainReco[neta],
				hMat[neta],0,lumi,neta,responseEta[neta]);
		generateEvents(eventsGenerate,fres[neta],hnlot[neta][0],hTrainTruet[neta][0],
				hTrainRecot[neta][0],hMatt[neta][0],0,lumi,neta);
		generateEvents(eventsGenerate,fres[neta],hnlot[neta][1],hTrainTruet[neta][1],
				hTrainRecot[neta][1],hMatt[neta][1],0,lumi,neta);
	}
	///////Unfold/////////////////////////////////////////////////////////////////
	cout << "============== UNFOLD ===============" << endl;
	TCanvas* cptSpec[6];
	TH1D* hclo[6];     TCanvas* cclo[6];// unfolding closure
	TH1D* htom[6][2];  TCanvas* ctom[6];// unfolded true / mesaured spectrum
	TH1D* hclot[6][2]; TCanvas* cclot[6][2];// unfolding closure for comparision
	TLegend* leg[6];   TCanvas* cMat[6];///to draw matrix
	TH1D* hcloMC[6];   TCanvas* cMC[6];
	TCanvas* cptSpecCl[6];
	TCanvas* cClosCom[6][2]; //compare closures from different spectrum, unfolding routines
	TH1D* hclo2[6][2]; //for comparing unfolding routines
	//unfolded histograms //0-bayes, 1-svd, 2- bin by bin
	TH1D* hUnf[4];//unfolded
	TH1D* hUnft[3];//test pdfs

	TH1D* hrecMC[6]; TH1D* hgenMC[6]; ///for unfolding check
	TFile* fMC = new TFile((indir+"/AnalyzerIncHistosMC_winter14.root").c_str(),"read"); 
	TH2D* hcov[6]; // covariance matrix
	TCanvas* ccov[6];
	/// Check Unfolding/////////////////////////////-------------------------------
	for(unsigned i(0);i<6;++i) {

		///set the genMC and recMC
		sprintf(name,"hrecMC_eta%d",i);
		hrecMC[i] = (TH1D*)fMC->Get(name);
		sprintf(name,"hgenMC_eta%d",i);
		hgenMC[i] = (TH1D*)fMC->Get(name);

		TH1D*  hDat = hTrainReco[i];
		TH1D*  hGen = (TH1D*)hTrainTrue[i];

		RooUnfoldBayes   unfold1  (responseEta[i], hDat, 3);    // OR
		RooUnfoldSvd     unfold2  (responseEta[i], hDat, 4);   // OR
		RooUnfoldBinByBin unfold3 (responseEta[i], hDat);

		hUnf[0]= (TH1D*) unfold1.Hreco();
		//		hUnf[1]= (TH1D*) unfold2.Hreco();
		hUnf[2]= (TH1D*) unfold3.Hreco();
		///test Pyhtia spectrum ======
		hrecMC[i]->Reset();
    sprintf(name,"hrecPy_eta%d",i);
    TH1D* hrecPy = (TH1D*)hTrainReco[i]->Clone(name);
    hrecPy->Reset();
    sprintf(name,"hgenPy_eta%d",i);
    TH1D* hgenPy = (TH1D*)hTrainTrue[i]->Clone(name);
    spFillHisto(hgenMC[i],hgenPy);
		generateEventsMC(eventsGenerate,fres[i],hgenPy,hTrainTruet[i][2],
				hrecPy,hMatt[i][2],0,lumi,i);
		RooUnfoldBayes unfoldMC (responseEta[i], hrecPy,4);
		hUnf[3] = (TH1D*) unfoldMC.Hreco();

		//Unfolding with  different pdfs
		RooUnfoldBinByBin   unfoldt1_bb  (responseEta[i], hTrainRecot[i][0]);   
		TH1D* hUnft1_bb = (TH1D*) unfoldt1_bb.Hreco();
		RooUnfoldBinByBin   unfoldt2_bb  (responseEta[i], hTrainRecot[i][1]);   
		TH1D* hUnft2_bb = (TH1D*) unfoldt2_bb.Hreco();
		////svd 
		//    RooUnfoldSvd   unfoldt1  (responseEta[i], hTrainRecot[i][0], 4);   
		//		TH1D* hUnft1= (TH1D*) unfoldt1.Hreco();
		//    RooUnfoldSvd   unfoldt2  (responseEta[i], hTrainRecot[i][1], 4);   
		//		TH1D* hUnft2= (TH1D*) unfoldt2.Hreco();
		//// bin by bin 
		RooUnfoldBayes   unfoldt1  (responseEta[i], hTrainRecot[i][0], 3);   
		TH1D* hUnft1= (TH1D*) unfoldt1.Hreco();
		RooUnfoldBayes   unfoldt2  (responseEta[i], hTrainRecot[i][1], 3);   
		TH1D* hUnft2= (TH1D*) unfoldt2.Hreco();

		///unfolding with different pdfs
		char label[100];
		char label1[100], label2[100];
		sprintf(label1,"%s","nnpdf2.1");
		sprintf(label2,"%s","abm11");
		sprintf(label,"%s","CT10");
		drawClosure(cclo[i],hUnf[0],hGen,hclo[i],i,saveplots,label); 
		//		drawClosure(cclot[i][0],hUnft1,hTrainTruet[i][0],hclot[i][0],i,saveplots,label1); 
		//		drawClosure(cclot[i][1],hUnft2,hTrainTruet[i][1],hclot[i][0],i,saveplots,label2); 
		drawClosure(cMC[i],hUnf[3],hTrainTruet[i][2],hcloMC[i],i,saveplots,"Pythia"); //pythia test unfolding
		//draw on the same canvas
		drawClosure(cClosCom[i][0],i,saveplots,hUnft1,hTrainTruet[i][0],hclot[i][0],label1,
				hUnft2,hTrainTruet[i][1],hclot[i][0],label2);
		///unfolding with different methods
		TH1D* htmp = (TH1D*)hUnf[0]->Clone("htmp");
		spFillHisto(hTrainReco[i],htmp);
		drawClosure2(cClosCom[i][1],i,saveplots,hUnf[0],htmp,hclo2[i][0],"Bayes",
				hUnf[2],htmp,hclo2[i][1],"BinByBin");

		sprintf(name,"c_%s",hMat[i]->GetName());
		cMat[i] = new TCanvas(name,name,600,600);
		cMat[i]->SetLogy();cMat[i]->SetLogx();
		drawMatrix(hMat[i],i); 
		hMat[i]->Draw("colz");

		sprintf(name,"%.1f < |y| < %.1f",ybins[i],ybins[i+1]);
		TLatex latex;latex.SetTextSize(0.04);latex.SetTextFont(42);
		latex.DrawTextNDC(.3,.8,name);
		//draw covariance matrix-----------------
		sprintf(name,"ccov_%d",i);
		ccov[i] = new TCanvas(name,name,600,600);
		sprintf(name,"hcov_%d",i);
		//error = 2
		const unsigned nrec(hTrainReco[i]->GetNbinsX());
		const double *xrec = hTrainReco[i]->GetXaxis()->GetXbins()->GetArray();// x bins reco
		hcov[i]= CorrelationHist (unfold1.Ereco((RooUnfold::ErrorTreatment)2),
				name, "Unfolded correlation matrix",
				xrec,nrec,i
				);
		ccov[i]->SetLogx(1); ccov[i]->SetLogy(1);
		hcov[i]->Draw("colz");
		sprintf(name,"%.1f < |y| < %.1f",ybins[i],ybins[i+1]);
		latex.DrawTextNDC(.3,.8,name);
		//end draw covariance matrix-----------------

		//unfolded, gen and reco spectrum
		cout<<"Drawing histos eta low edge "<<ybins[i]<<endl;
		drawSpec(cptSpec[i],hDat,hGen,hUnf[0],leg[i],saveplots,"");
		//save plots  
		if(saveplots) {
			cout<<cMat[i]->GetName()<<endl;
			sprintf(name,"plots/%s.eps",cMat[i]->GetName());
			cMat[i]->SaveAs(name);
			sprintf(name,"plots/%s.eps",ccov[i]->GetName());
			cout<<name<<" "<<ccov[i]<<endl;
			ccov[i]->SaveAs(name);
		}
	}

	//	cout<<"Test unfolding "<<endl;
	//	////test bayesian unfolding ------------------------------------
	//	TFile* fMC = new TFile(("AnalyzerIncHistosMC.root"),"read"); 
	//	TH1D* hMC[6];
	//	unsigned kiter[6] = {1,2,2,3,4,5};
	//	unsigned kmark[6] = {24,20,25,21,26,22};
	//	unsigned kcolr[6] = {kBlack,kTeal+2,kBlack,kBlue,kMagenta,kRed};
	//	unsigned kstyl[6] = {1,1,1,2,9,10};
	//	TH1D* hclo_tBay[7][6];//index : [testnum][neta]
	//	TH1D* hUnf_tBay[7][6];//index : [testnum][neta]
	//	TCanvas* cUnf_tBay[6];
	//
	//	for(unsigned i(0);i<6;++i) {
	//		gStyle->SetErrorX(0);
	//		TH1D* hunftmp[6];
	//		ostringstream oss;
	//		oss<<"hptUnf_eta_"<<ybins[i]<<"_"<<ybins[i+1];
	//		TH1D*  hDat_fullrange  = (TH1D*)fMC->Get(oss.str().data());
	//		TLegend* leg_unf = new TLegend(0.4,0.6,0.6,0.8);
	//		leg_unf->SetFillColor(0);  leg_unf->SetTextFont(42);
	//		leg_unf->SetBorderSize(0); leg_unf->SetTextSize(0.035);
	//
	//		for(unsigned niter(0);niter<6;++niter) {
	//			sprintf(name,"hDattmp_%d",i);
	//			TH1D* hdattmp = (TH1D*)hTrainTrue[i]->Clone(name);
	//			spFillHisto(hDat_fullrange,hdattmp);
	//			RooUnfoldBayes    unf (responseEta[i], hDat_fullrange, kiter[niter]);
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
	//			hunftmp[niter]->SetStats(kFALSE);	hunftmp[niter]->SetTitle("");
	//			//      sprintf(name,"cUnf_tBay_eta%d_niter%d",i,niter);
	//			if(niter==3) {
	//				sprintf(name,"cUnf_tBay_eta%d",i);
	//				cUnf_tBay[i] = new TCanvas(name,name,600,600);
	//				hunftmp[niter]->GetXaxis()->SetRangeUser(74,maxPt[i]-1);
	//				hunftmp[niter]->GetYaxis()->SetRangeUser(0.5,1.5);
	//				hunftmp[niter]->GetYaxis()->SetTitle("Ratio of Unfolded Spectrum");
	//				hunftmp[niter]->GetYaxis()->SetTitleOffset(1.4);
	//				hunftmp[niter]->Draw("pe");
	//			}
	//			else {
	//				hunftmp[niter]->Draw("pe same");
	//			}
	//			sprintf(name," niter %d/niter %d",kiter[niter],kiter[niter-1]);
	//			leg_unf->AddEntry(hunftmp[niter],name,"pel");
	//		}
	//		leg_unf->Draw("same");
	//		TF1* fline = new TF1("fline","[0]",74,maxPt[i]-1);
	//		fline->SetParameter(0,1);
	//		fline->SetLineColor(kBlack);
	//		fline->Draw("same");
	//		sprintf(name,"%.2f<|y|<%.2f",ybins[i],ybins[i+1]);
	//		TLatex latex;latex.SetTextSize(0.035);latex.SetTextFont(42);
	//		latex.DrawTextNDC(.2,.7,name);
	//
	//
	//		sprintf(name,"plots/iter_eta%d.eps",i);
	//		if(saveplots)
	//			cUnf_tBay[i]->SaveAs(name);
	//	}
	//	////test bayesian unfolding ------------------------------------



	//////
	///Write Histos and Response object to file
	TFile* outfile = new TFile("UnfoldJetsTrainingTest.root","recreate");
	for(unsigned i(0);i<nybins-1;++i)
	{
		hTrainTrue[i]->Write();
		hTrainReco[i]->Write();
		hMat[i]->Write();
		hcov[i]->Write();
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
	hunf->SetMarkerColor(kRed);
	hunf->SetMarkerStyle(20);
	hunf->SetMarkerSize(0.5);
	c->SetLogx(1);c->SetLogy(1);

	hreco->SetLineColor(kSpring+4);
	hreco->SetMarkerColor(kSpring+4);
	hreco->SetMarkerStyle(20);
	hreco->SetMarkerSize(0.5);

	hgen->SetLineColor(kBlack);
	hgen->SetMarkerColor(kBlack);
	hgen->SetMarkerStyle(20);
	hgen->SetMarkerSize(0.1);
	hgen->GetXaxis()->SetMoreLogLabels(kTRUE);
	hgen->GetXaxis()->SetNoExponent(kTRUE);
	hgen->SetTitle(""); hgen->SetStats(kFALSE);
	hgen->GetXaxis()->SetTitle("Jet P_{T} (GeV)");

	hgen->GetYaxis()->SetRangeUser(1e2,1e8);
	hunf->GetYaxis()->SetRangeUser(1e2,1e8);
	hreco->GetYaxis()->SetRangeUser(1e2,1e8);
	hgen->GetXaxis()->SetRangeUser(30,300);
	hunf->GetXaxis()->SetRangeUser(30,300);
	hreco->GetXaxis()->SetRangeUser(30,300);


	hgen->Draw("hist");
	hreco->Draw("p same");
	hunf->Draw("p same");

	leg = new TLegend(0.6,0.7,0.87,0.89);
	leg->SetFillColor(0);  leg->SetTextFont(42);
	leg->SetBorderSize(0); leg->SetTextSize(0.025);
	leg->AddEntry(hreco,"Reco Spectrum","lpe");
	leg->AddEntry(hgen ,"Gen Spectrum","l");
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
	hclo->GetXaxis()->SetTitle("Jet P_{T} (GeV)");
	hclo->GetXaxis()->SetMoreLogLabels(kTRUE);
	hclo->GetXaxis()->SetNoExponent(kTRUE);
	hclo->GetYaxis()->SetRangeUser(0.9,1.1);
	hclo->GetXaxis()->SetRangeUser(74,maxPt[eta]-1);
	hclo->GetYaxis()->SetTitleOffset(1.5);
	hclo->SetLineColor(kBlack);
	hclo->SetMarkerColor(kBlack);
	hclo->SetMarkerStyle(20);
	hclo->SetMarkerSize(1.0);
	hclo->SetTitle("");hclo->SetStats(kFALSE);
	hclo->Draw();

	sprintf(name,"%.1f<|y|<%.1f",ybins[eta],ybins[eta+1]);
	TLatex latex;latex.SetTextSize(0.04);latex.SetTextFont(42);
	latex.DrawTextNDC(.6,.8,name);
	latex.DrawTextNDC(.6,.7,label);

	if(label=="Pythia")
	{
		///find the slope
		TF1* f = new TF1("fline","[0]",74,maxPt[eta]-1);
		f->SetParameter(0,1);
		//hclo1->Fit(f,"R same");
		f->Draw("same");
		hclo->Draw("same");
	}


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
	hclo1->SetTitle(""); hclo1->SetStats(kFALSE);
	hclo2->SetTitle(""); hclo2->SetStats(kFALSE);
	//make hclo and draw
	hclo1->Divide(hunf1,hgen1);
	hclo1->GetYaxis()->SetTitle("Unfolded/Gen");
	hclo1->GetXaxis()->SetTitle("Jet P_{T} (GeV)");
	hclo1->GetYaxis()->SetRangeUser(0.9,1.1);
	hclo1->GetXaxis()->SetRangeUser(74,maxPt[eta]-1);
	hclo1->GetYaxis()->SetTitleOffset(1.5);
	hclo1->GetXaxis()->SetMoreLogLabels(kTRUE);
	hclo1->GetXaxis()->SetNoExponent(kTRUE);
	hclo1->SetLineColor(kBlack);
	hclo1->SetMarkerColor(kBlack);
	hclo1->SetMarkerStyle(20);
	hclo1->SetMarkerSize(1.0);
	hclo1->Draw();
	hclo2->Divide(hunf2,hgen2);
	hclo2->GetYaxis()->SetTitle("Unfolded/Gen");
	hclo2->GetXaxis()->SetTitle("Jet P_{T} (GeV)");
	hclo2->GetYaxis()->SetRangeUser(0.9,1.1);
	hclo2->GetXaxis()->SetRangeUser(74,maxPt[eta]-1);
	hclo2->GetYaxis()->SetTitleOffset(1.5);
	hclo2->SetLineColor(kBlue);
	hclo2->SetMarkerColor(kBlue);
	hclo2->SetMarkerStyle(33);
	hclo2->SetMarkerSize(1.0);
	hclo2->Draw("same");
	///find the slope
	TF1* f = new TF1("fline","[0]",74,maxPt[eta]-1);
	f->SetParameter(0,1);
	//hclo1->Fit(f,"R same");
	f->Draw("same");
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
	sprintf(name,"c_%s_%s",hunf1->GetName(),"unfComp");
	cout<<name<<endl;
	c = new TCanvas(name,name,600,600);
	c->cd();
	c->SetLogx(1);
	sprintf(name,"hclos_%s",hgen1->GetName());
	hclo1 = (TH1D*)hgen1->Clone(name);
	hclo1->Reset(); hclo1->Sumw2();
	hclo2 = (TH1D*)hgen2->Clone(name);
	hclo2->Reset(); hclo2->Sumw2();
	hclo1->SetTitle(""); hclo1->SetStats(kFALSE);
	hclo2->SetTitle(""); hclo2->SetStats(kFALSE);
	//make hclo and draw
	hclo1->Divide(hunf1,hgen1);
	hclo1->GetYaxis()->SetTitle("Unfolded/Reco");
	hclo1->GetXaxis()->SetTitle("Jet P_{T} (GeV)");
	hclo1->GetYaxis()->SetRangeUser(0.3,1.1);
	hclo1->GetXaxis()->SetRangeUser(74,maxPt[eta]-1);
	hclo1->GetYaxis()->SetTitleOffset(1.5);
	hclo1->GetXaxis()->SetMoreLogLabels(kTRUE);
	hclo1->GetXaxis()->SetNoExponent(kTRUE);
	hclo1->SetLineColor(kBlack);
	hclo1->SetMarkerColor(kBlack);
	hclo1->SetMarkerStyle(20);
	hclo1->SetMarkerSize(1.0);
	hclo1->Draw();
	hclo2->Divide(hunf2,hgen2);
	hclo2->GetYaxis()->SetTitle("Unfolded/Gen");
	hclo2->GetXaxis()->SetTitle("Jet P_{T} (GeV)");
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
	h->SetTitle("");  h->SetStats(kFALSE);
	h->GetXaxis()->SetRangeUser(74,maxPt[eta]-1);

	h->GetXaxis()->SetMoreLogLabels(kTRUE);
	h->GetXaxis()->SetNoExponent(kTRUE);
	h->GetYaxis()->SetRangeUser(0,maxPt[eta]-1);
	h->GetYaxis()->SetTitleOffset(1.2);
	h->GetXaxis()->SetTitleOffset(1.2);
	h->GetYaxis()->SetTitle("Gen P_{T} (GeV)");
	h->GetXaxis()->SetTitle("Reco P_{T} (GeV)");

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
		double lumi,unsigned eta,RooUnfoldResponse* resp)
{
	hTrue->Sumw2(); hReco->Sumw2();
	const unsigned nbins(hnlo->GetNbinsX());
	// cout<<nbins<<' '<<events<<endl;
	//apply NLO
	TFile* npfile= new TFile("npCorrections.root","read");
	sprintf(name,"fratio_eta%d",eta);
	TF1* f = (TF1*)npfile->Get(name);
	TH1D* nlotmp = (TH1D*)hnlo->Clone("nptmp");
	applyNP(nlotmp,f);
	//finished apply np corr
	for(unsigned bin(1);bin<=nbins;++bin)
	{
		double xsec  = nlotmp->GetBinContent(bin)*nlotmp->GetBinWidth(bin)*lumi;
		double binlo = nlotmp->GetBinLowEdge(bin);
		double binhi = nlotmp->GetBinWidth(bin)+binlo;
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

	delete nlotmp;
	delete npfile;
	delete f;
}

void generateEventsTest(unsigned events, TF1* fres, TH1* hnlo,
		TH1D* hTrue, TH1D* hReco,TH2D* hMatrix, 
		double lumi,unsigned eta,double scaleResolutionBy,RooUnfoldResponse* resp)
{
	hTrue->Sumw2(); hReco->Sumw2();
	const unsigned nbins(hnlo->GetNbinsX());
	// cout<<nbins<<' '<<events<<endl;
	//apply NLO
	TFile* npfile= new TFile("npCorrections.root","read");
	sprintf(name,"fratio_eta%d",eta);
	TF1* f = (TF1*)npfile->Get(name);
	TH1D* nlotmp = (TH1D*)hnlo->Clone("nptmp");
	applyNP(nlotmp,f);
	//finished apply np corr
	for(unsigned bin(1);bin<=nbins;++bin)
	{
		double xsec  = nlotmp->GetBinContent(bin)*nlotmp->GetBinWidth(bin)*lumi;
		double binlo = nlotmp->GetBinLowEdge(bin);
		double binhi = nlotmp->GetBinWidth(bin)+binlo;
		for(unsigned i(0);i<events;++i)
		{
			double genpt = gRandom->Uniform(binlo,binhi);
			double resol = gRandom->Gaus(1.0,fres->Eval(genpt));
			double recpt = genpt*(resol+scaleResolutionBy);
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

	delete nlotmp;
	delete npfile;
	delete f;
}

void generateEventsMC(unsigned events, TF1* fres, TH1* hnlo,
		TH1D* hTrue, TH1D* hReco,TH2D* hMatrix,double dataRes, 
		double lumi,unsigned eta,RooUnfoldResponse* resp)
{
	hTrue->Sumw2(); hReco->Sumw2();
	const unsigned nbins(hnlo->GetNbinsX());
	// cout<<nbins<<' '<<events<<endl;
	//apply NLO
	TFile* npfile= new TFile("npCorrections.root","read");
	sprintf(name,"fratio_eta%d",eta);
	TF1* f = (TF1*)npfile->Get(name);
	TH1D* nlotmp = (TH1D*)hnlo->Clone("nptmp");
	//	applyNP(nlotmp,f);
	//finished apply np corr
	for(unsigned bin(1);bin<=nbins;++bin)
	{
		double xsec  = nlotmp->GetBinContent(bin);
		//		double xsec  = nlotmp->GetBinContent(bin)*nlotmp->GetBinWidth(bin);
		double binlo = nlotmp->GetBinLowEdge(bin);
		double binhi = nlotmp->GetBinWidth(bin)+binlo;
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

	delete nlotmp;
	delete npfile;
	delete f;
}

void drawSpecCl(TCanvas* c,TH1* hdt, TH1* hmc,TH1* hclo, unsigned eta)
{
	hdt->SetTitle("");    hdt->SetStats(kFALSE);
	hmc->SetTitle("");    hmc->SetStats(kFALSE);
	hclo->SetTitle("");   hclo->SetStats(kFALSE);
	hclo->GetXaxis()->SetTitle("Jet P_{T} (GeV)");
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
	hclo->GetXaxis()->SetNoExponent(kTRUE);
	hclo->GetXaxis()->SetMoreLogLabels(kTRUE);
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

void histoFromErr(TH1* hin, TH1* hout)
{
	hout->Reset();
	const  unsigned numbins(hout->GetNbinsX());
	for(unsigned i(1);i<=numbins;++i) {
		double bin(hin->GetBinContent(i));
		double bct(hin->GetBinCenter(i));
		double err(hin->GetBinError(i));
		if(bin)
			hout->SetBinContent(i,err*1.0/bin);
		else 
			hout->SetBinContent(i,0);
	}
}

void drawErrors(TCanvas* c, TH1* rec,TH1* unf,unsigned eta)
{
	c->SetLogx(); c->SetLogy();
	c->cd();
	rec->SetTitle("");   rec->SetStats(kFALSE);
	unf->SetTitle("");   unf->SetStats(kFALSE);
	unf->GetXaxis()->SetTitle("Jet P_{T} (GeV)");
	unf->GetYaxis()->SetTitle("Relative Error");
	unf->GetXaxis()->SetRangeUser(74,maxPt[eta]-1);
	//	rec->GetYaxis()->SetRangeUser(0.,2.0);
	rec->SetLineColor(kBlue);
	unf->SetLineColor(kRed);
	unf->Draw(""); 
	rec->Draw("same");

	TLegend* leg = new TLegend(0.4,0.5,0.6,0.7);
	leg->SetFillColor(0);  leg->SetTextFont(42);
	leg->SetBorderSize(0); leg->SetTextSize(0.03);
	leg->AddEntry(rec,"Reco","l");
	leg->AddEntry(unf,"Unfolded","l");
	leg->Draw("same");

	sprintf(name,"%.2f<|y|<%.2f",ybins[eta],ybins[eta+1]);
	TLatex latex;latex.SetTextSize(0.04);latex.SetTextFont(42);
	latex.DrawTextNDC(.4,.8,name);
}


void initArray(double* arr,unsigned size)
{
	assert(arr);
	for(unsigned i(0);i<size;++i)
		arr[i]=0;
}

void getParUncertainty(TH1D* hcn,TH1D** hvr,TH1D* hunc,unsigned size)
{
	hunc->Reset();
	unsigned bins (hcn->GetNbinsX());
	double dx[100];
	initArray(dx,100);
	for(unsigned par(0);par<size;++par) {
		assert(bins==hvr[par]->GetNbinsX());
		for(unsigned bin(1);bin<=bins;++bin)
		{
			double xcn(hcn->GetBinContent(bin));
			double xpr(hvr[par]->GetBinContent(bin));
			dx[bin] = max(fabs(xpr-xcn),dx[bin]);
		}//nparameters
	}
	cout<<fixed<<setprecision(3);
	for(unsigned bin(1);bin<=bins;++bin)
	{
		double xcn(hcn->GetBinContent(bin));
		if(xcn)
			cout<<hcn->GetBinCenter(bin)<<" x central "<<xcn<<" dx "<<dx[bin]<<" "<<dx[bin]/xcn<<endl;
		if(xcn)      hunc->SetBinContent(bin,dx[bin]/xcn);
		else         hunc->SetBinContent(bin,0);
	} //fill the histo
}



void getJesUncertainty(TH1* hcn, TH1* hup, TH1* hdn, double* uncert)
{
	cout<<setprecision(3)<<fixed;
	const unsigned bins = hcn->GetNbinsX();
	assert(hup->GetNbinsX()==bins);
	assert(hdn->GetNbinsX()==bins);
	for(unsigned bin(1);bin<=bins;++bin)
	{
		double xcn = hcn->GetBinContent(bin);
		double xup = hup->GetBinContent(bin);
		double xdn = hdn->GetBinContent(bin);
		double maxerr = max(fabs(xup-xcn),fabs(xdn-xcn));
		if(xcn)
			uncert[bin-1]  = maxerr/xcn; 
		else 
			uncert[bin-1] = 0;
		if(false && xcn)
			cout<<" "<<hcn->GetBinCenter(bin)
				<<" "<<uncert[bin-1]<<" "<<" maxerr "<<maxerr/xcn
				<<endl; 
	}//bins loop
}

void getUnfoldUncertainty(TH1* hcn,TH1* hup,TH1* hdn,TH1* hunc)
{
	hunc->Reset();
	unsigned bins (hcn->GetNbinsX());
	assert(bins==hup->GetNbinsX());
	assert(bins==hdn->GetNbinsX());
	for(unsigned bin(1);bin<=bins;++bin)
	{
		double xcn(hcn->GetBinContent(bin));
		double xup(hup->GetBinContent(bin));
		double xdn(hdn->GetBinContent(bin));
		//      double dx = max(fabs(xcn-xup),fabs(xcn-xdn));
		double dx = max(fabs(xdn-xcn),fabs(xup-xcn));
		if(xcn && true)
			cout<<hcn->GetBinCenter(bin)<<" x central "<<xcn<<" dx "<<dx<<" "<<dx/xcn<<endl;
		if(xcn)      hunc->SetBinContent(bin,dx/xcn);
		else         hunc->SetBinContent(bin,0);
	}
}

void fillHistoFromArr(double* arr,TH1D* h)
{
	h->Reset();
	unsigned bins(h->GetNbinsX());
	for(unsigned bin(1);bin<=bins;++bin)
	{
		if(arr[bin-1])
			h->SetBinContent(bin,arr[bin-1]);
		else 
			h->SetBinContent(bin,0);
	}
}


void addError(double* in, double* out,unsigned size)
{
	cout<<setprecision(3)<<fixed;
	for(unsigned i(0);i<size;++i)
	{
		if(true && in[i])
			cout<<out[i]<<"+-"<<in[i]
				<<"="<<sqrt(out[i]*out[i]+in[i]*in[i])<<endl;
		out[i]=sqrt(out[i]*out[i]+in[i]*in[i]);
	}
}

void applyNP(TH1* h,TF1* f)
{
	const  unsigned numbins(h->GetNbinsX());
	for(unsigned i(1);i<=numbins;++i)
	{
		double np = f->Eval(h->GetBinCenter(i));
		double x(h->GetBinContent(i));
		h->SetBinContent(i,h->GetBinContent(i)*np);
		//			cout<<h->GetBinCenter(i)<<" "<<x<<" x "<<np<<" = "<< h->GetBinContent(i)<<endl;
	}
}

void testFunction(double xlo, double xhi)
{
	cout<<findBinsFromArray(xlo,xhi,x[0],nx[0])<<endl; 
}

TH2D* CorrelationHist (const TMatrixD& cov, const char* name, const char* title,
		const  double* bins, const unsigned nbins,unsigned i)
{
	Int_t nb= cov.GetNrows();
	TH2D* h= new TH2D (name, title, nbins,bins, nbins, bins);
	h->SetAxisRange (-1.0, 1.0, "Z");
	for(int i=0; i < nb; i++)
		for(int j=0; j < nb; j++) {
			Double_t Viijj= cov(i,i)*cov(j,j);
			if (Viijj>0.0) h->SetBinContent (i+1, j+1, cov(i,j)/sqrt(Viijj));
		}
	h->SetTitle("");
	h->SetStats(kFALSE);
	h->GetXaxis()->SetMoreLogLabels(kTRUE);
	h->GetYaxis()->SetMoreLogLabels(kTRUE);
	h->GetXaxis()->SetNoExponent(kTRUE);
	h->GetYaxis()->SetNoExponent(kTRUE);
	h->GetYaxis()->SetRangeUser(74,maxPt[i]-1);
	h->GetXaxis()->SetRangeUser(74,maxPt[i]-1);
	h->GetXaxis()->SetTitle("Jet P_{T} (GeV)");
	h->GetYaxis()->SetTitle("Jet P_{T} (GeV)");
	return h;
}

