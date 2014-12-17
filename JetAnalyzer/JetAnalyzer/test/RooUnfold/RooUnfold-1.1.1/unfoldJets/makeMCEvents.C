#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "TRandom.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TArrayD.h"
#include "TMatrixD.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "matchOneToOne.hh"
#include "../../../AnalyzerProcedures.h"
/////Jet Class 
namespace makeMCEvents {
	class Jet {
		public : 
			double pt;       double mass;     
			double eta;      double y;
			double phi;      
			//constructors
			Jet(double p, double e, double ph, double m) 
			{
				pt = p;   	mass = m; 
				eta = e;  	y=-1;
				phi = ph; 
			};

			Jet()
				: pt(1e-8),eta(1e-8),phi(1e-8),mass(1e-8)	
			{ 
				y=-1;
			};
			~Jet(){ };
			Jet(const TLorentzVector& vec ){
				if(vec.Pt()>0)
				{
					pt = vec.Pt();eta = vec.Eta();phi = vec.Phi();mass = vec.M(); 
				}
				else 
				{
					pt = 1e-8;eta = 1e-8;phi = 1e-8;mass = 1e-8; 
				}
			};
			inline TLorentzVector  getLorentzVector(){
				TLorentzVector vec;
				vec.SetPtEtaPhiM(pt,eta,phi,mass);
				return (vec);
			};
			inline void setPtEtaPhiY(double p, double e,double ph, double rap ){
				pt=p;  eta=e; phi=ph; y=rap;
				TLorentzVector tmp; 
				tmp.SetPtEtaPhiM(pt,eta,phi,10);
				double pz = tmp.Pz();
				double eg = pz*((pow(2.71828,2*y)+1)/(pow(2.71828,2*y)-1));
				tmp.SetPtEtaPhiE(pt,eta,phi,eg);
				mass=tmp.M();
				if(true)
					std::cout
						<<"Calculated mass "<<mass
						<<"\nCalculated px|py|pz|eg "<<tmp.Px()<<'|'<<tmp.Py()<<'|'<<pz<<'|'<<eg
						<<"\ngiven|calculated y  "<<y<<"|"<<tmp.Rapidity()
						<<std::endl;
			};
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
//
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
		for(unsigned i(0);i<nBins;++i)
			h->SetBinContent(i+1,xsecV[i]);
		return h;
	};
}
using namespace std;
using namespace makeMCEvents;
void Unfold();
int main (int argc, char** argv) {
	Unfold();
	return 0;
}
void Unfold() {

	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptTitle(kFALSE);
	//gROOT->SetBatch(kTRUE);

  const char* nloDir("/uscms_data/d3/tlibeiro/inclusivexsec/Generators/generateIncJets/inputXsecs");
  const char* fitDir("./");
   sprintf(name,"%s/ResolutionFits.root",fitDir);
   TFile* resFile = new TFile(name,"read");

	const unsigned eventsGenerate(1e5);
	TH1F* hflat = new TH1F("hflat","hflat",1020,0,1020); hflat->Sumw2();
  TH1D* hnlo[6];//nlo spectrum
  TH1D* hsmr[6];//smeared spectrum
  TF1*  fres[6];//resfits

//	cout << "Generating Random Numbers " << endl;
//	for(unsigned i(0);i<eventsGenerate;++i)
//	{
//		hflat->Fill(gRandom->Uniform(20,1000));
//	//	hflat->Fill(gRandom->Gaus(500,50));
//		//cout<<gRandom->Uniform(20,1000)<<endl;
//	}
      ////get nlo spectrum and fits
  for(unsigned neta(0);neta<6;++neta) {
    char infile[500]; char title[500];
    sprintf(infile,"%s/fnl4332y%i.tab.xsec.cppInput",nloDir,neta);
    sprintf(title ,"nloy%i",neta);
    hnlo[neta] = makeNLOSpectrum(infile,title);
    sprintf(name,"fresfit_eta%d",neta);
    fres[neta] = (TF1*)resFile->Get(name);
  }
  ///generate events
	for(unsigned neta(0);neta<6;++neta) {
		char title[500];
		sprintf(title,"hsmr_eta%d",neta);
		hsmr[neta] = (TH1D*)hnlo[0]->Clone(title);
		hsmr[neta]->Reset();
		hsmr[neta]->Sumw2();
		const unsigned nbins(hsmr[neta]->GetNbinsX());
		for(unsigned bin(1);bin<=nbins;++bin)
		{
			double xsec  = hnlo[neta]->GetBinContent(bin);
			double binlo = hnlo[neta]->GetBinLowEdge(bin);
			double binhi = hnlo[neta]->GetBinWidth(bin)+binlo;
			for(unsigned i(0);i<eventsGenerate;++i)
			{
			 double genpt = gRandom->Uniform(binlo,binhi);
       double resol = gRandom->Gaus(1.0,fres[neta]->Eval(genpt));
       double recpt = genpt*resol;
			 hsmr[neta]->Fill(recpt,xsec);
			}
		}
    hsmr[neta]->Scale(1.0/eventsGenerate);
	}

//	TCanvas* cflat = new TCanvas("cflat","cflat",600,600);
//	hflat->Draw();
	TCanvas* csp = new TCanvas("cnlo","cnlo",800,600);
	csp->SetLogy(1);
  hnlo[0]->SetMarkerStyle(21);
	hnlo[0]->Draw("P");
//	TCanvas* csmr = new TCanvas("csmr","csmr",800,600);
//	csmr->SetLogy(1);
	hsmr[0]->Draw("hist E same");
}
