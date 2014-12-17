#include <iostream>
#include <math.h> 
#include "TNamed.h"
#include "TString.h"
#include "RooUnfold.h"
#include "TVectorDfwd.h"
#include "TMatrixDfwd.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TNtuple.h"
#include "RooUnfoldResponse.h"
#define lne 2.71828182846
  
namespace UnfoldJetNS {
// Data Members
  TCanvas*           canvas;
  TLegend            *lTrain, *lTest, *lErrors;
  TLegend            *leg1, *leg2;
  TH1                *hTrain, *hTrainFake, *hTrue, *hMeas, *hReco, *hFake, *hRes, *hPulls;
  TH1D*              hTrainTrue[10];
  TH1D*              hTrainReco[10];
  TH1                *hUnfErr, *hToyErr, *hParmChi2, *hParmErr, *hParmRes, *hParmRms;
  TH2D               *hResmat, *hCorr, *hMeasCorr;
  TH2D               *hMat[10];
  RooUnfoldResponse* response;
  RooUnfold*         unfold;
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
}

/////Jet Class 
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

}
