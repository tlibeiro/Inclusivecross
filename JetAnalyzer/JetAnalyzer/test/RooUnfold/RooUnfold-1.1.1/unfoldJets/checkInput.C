#include <iostream>
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
#include "TLorentzVector.h"
#include "matchOneToOne.hh"
#include "../../../AnalyzerProcedures.h"
/////Jet Class 
namespace checkInput {
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
}
using namespace std;
using namespace checkInput;
void Unfold();
int main (int argc, char** argv) {
	Unfold();
	return 0;
}
void Unfold() {

	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptTitle(kFALSE);
	//gROOT->SetBatch(kTRUE);

	cout << "==================================== TRAIN ====================================" << endl;
	TFile* infile = new TFile("../../../AnalyzerHistos_allpt_Unfold.root","read");
	TTree* unfTree = (TTree*)infile->Get("UnfoldJetTree");
	TChain* chain = new TChain("UnfoldJetTree");
  chain->Add("../../../AnalyzerHistos_allpt[0-7]_Unfold.root");
  unfTree = (TTree*)chain;
	///Read Ntuple =============================
	const unsigned maxJets(15);
	int njets(0), ngen(0);
	float jtpt [maxJets];  float gjtpt [maxJets];
	float jteta[maxJets];  float gjteta[maxJets];
	float jtphi[maxJets];  float gjtphi[maxJets];
	float jtm  [maxJets];  float gjty  [maxJets];
  double mcWeight(-1);
	unfTree->SetBranchAddress("jtpt" , jtpt );  unfTree->SetBranchAddress("gjtpt" , gjtpt ); 
	unfTree->SetBranchAddress("jtphi", jtphi);  unfTree->SetBranchAddress("gjtphi", gjtphi);
	unfTree->SetBranchAddress("jteta", jteta);  unfTree->SetBranchAddress("gjteta", gjteta);
	unfTree->SetBranchAddress("jtm"  , jtm);    unfTree->SetBranchAddress("gjty"  , gjty);
	unfTree->SetBranchAddress("njets", &njets); unfTree->SetBranchAddress("ngen"  , &ngen);
  unfTree->SetBranchAddress("mcweight", &mcWeight);
	///End Read Ntuple ==========================
	///Initialize Histos ========================

	unsigned nentries(unfTree->GetEntries());
	//nentries=1000;
	cout<<"Reading Unfold tree: Entries "<<nentries<<endl;
	vector<checkInput::Jet> jets, genjets;
	vector<int> matchedGen, matchedJet;
	const double dRMatch(0.5);//Gen-Jet match 
	DeltaRDistance d; 
  TH2D* hpteta1 = new TH2D ("hpteta1", "hpteta1",n1x[0],&x1[0][0],100,-3,3);
  TH2D* hpteta2 = new TH2D ("hpteta2", "hpteta2",n1x[0],&x1[0][0],100,-3,3);
  TH1D* hptgen = new TH1D ("hptgen", "hptgen",n1x[0],&x1[0][0]);
  TH1D* hptrec = new TH1D ("hptrec", "hptrec",n1x[0],&x1[0][0]);
	////Run on the training Tree===================
	for(unsigned entry(0);entry<nentries;++entry)
	{
		jets.clear(); genjets.clear();
		unfTree->GetEntry(entry);
		if(false)
			cout<<"Entry "<<entry
				<< " == "<<" Gen "<<ngen<<" Jet "<<njets 
				<< endl;
//		for(unsigned njet(0);njet<njets;++njet)
//			if(abs(jteta[njet])<2)
//				jets.push_back(checkInput::Jet(jtpt[njet], jteta[njet],
//							jtphi[njet], jtm[njet]));
//		for(unsigned ngjet(0);ngjet<ngen;++ngjet)
//			if(abs(gjteta[ngen])<2)
//				genjets.push_back(checkInput::Jet(gjtpt[ngjet], gjteta[ngjet],
//							gjtphi[ngjet],0));
//		matchOneToOne(jets,genjets,DeltaRDistance() ,&matchedGen,1000.0);  
//		matchOneToOne(genjets,jets,DeltaRDistance() ,&matchedJet,1000.0);  

		for(unsigned ngjet(0);ngjet<ngen;++ngjet) ///Matched Jets to Gen 
     {
			 hpteta1->Fill(gjtpt[ngjet],gjteta[ngjet]);
			 hptgen->Fill(gjtpt[ngjet],mcWeight);
		}
		for(unsigned njet(0);njet<njets;++njet) ///Matched Jets to Gen 
			{
      hpteta2->Fill(jtpt[njet],jteta[njet]);
			hptrec->Fill(jtpt[njet],mcWeight);
			}
		if(!(entry%200000))
			cout<<"Processing Event "<<(float)entry<<endl;
	}//entry loop
 
TCanvas* cgen = new TCanvas();
TCanvas* crec = new TCanvas();
TCanvas* cspec = new TCanvas();
 hpteta1->GetXaxis()->SetRangeUser(10,1000);
 hpteta2->GetXaxis()->SetRangeUser(10,1000);
// cgen->cd(); hpteta1->Draw("colz");
// crec->cd(); hpteta2->Draw("colz");
 cspec->cd(); 
	cspec->SetLogx(1);	cspec->SetLogy(1);
	hptgen->Draw();hptrec->Draw("same");
}
