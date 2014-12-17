#ifndef NTUPLEJETFILLER_HH
#define NTUPLEJETFILLER_HH
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <cassert>
#include <memory>
#include <algorithm>
#include <vector>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TProfile.h"
#include <TMath.h>
#include "TRandom.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

using namespace std;
using namespace edm;
class NtupleJetFiller {
	public: 
		static const unsigned maxJets=20;
		NtupleJetFiller(const std::string& name,
				TTree* t,const std::string& JEC_GlobalTag,const bool runMC);

		edm::Handle<vector<reco::PFJet> > jetCollection;
		const std::string name;
    const bool runOnMC;
		TTree* tree;
		void fillBranch(const edm::Event& iEvent,const edm::InputTag& inputTag,	
										const edm::InputTag& srcPrimaryVertex, const bool runOnMC);

			float jetpt [maxJets];      float jetpt_uncorr[maxJets];
			float jetm  [maxJets];   
			float jeteta[maxJets];     
			float jetphi[maxJets];    
			float jety  [maxJets];
      //
			void SetBranches();
			void SetBranch(float* x,const char* name,unsigned maxJetsOrPhotons);
			void SetBranch(int* x  ,const char* name,unsigned maxJetsOrPhotons);
    void SetBranchSingle( float* x, std::string name);
    void SetBranchSingle(unsigned  int* x, std::string name);
			unsigned nPV;
			FactorizedJetCorrector* jec;
      unsigned numGoodJets;
			const TLorentzVector getCorrectedJet(const reco::PFJet& jet);
			double getJEC(double curJetEta,double curJetPt,double curJetE);
};
//definitions
NtupleJetFiller::NtupleJetFiller(const std::string& n,
		TTree* t,const std::string& JEC_GlobalTag,const bool runMC)
:name(n),runOnMC(runMC)
{
	assert(t);
	tree = t;
  numGoodJets=0;
	for(unsigned i(0);i<maxJets;++i)
	{
		jetpt[i] =-1;     jetpt_uncorr[i]=-1; 
		jetm[i]  =-1;     
		jeteta[i]=-1;
		jetphi[i]=-1;
		jety[i]  =-1;
	}
	cout<<"arrays set "<<jetpt[0]<<endl;
	SetBranches();
	// ---- setting up the jec on-the-fly from text files...
	const std::string& fDir = JEC_GlobalTag;
	std::vector< JetCorrectorParameters > jecPars;
	std::vector< std::string > jecStr;
	const bool applyJECToGroomedJets(true);

	if(applyJECToGroomedJets) {
		jecStr.push_back( fDir + "_L1FastJet_AK7PFchs.txt" );
		jecStr.push_back( fDir + "_L2Relative_AK7PFchs.txt" );
		jecStr.push_back( fDir + "_L3Absolute_AK7PFchs.txt" );
		if (!runOnMC)
			jecStr.push_back( fDir + "_L2L3Residual_AK7PFchs.txt" );

		for (unsigned int i = 0; i < jecStr.size(); ++i){
			JetCorrectorParameters* ijec = new JetCorrectorParameters( jecStr[i] );
			jecPars.push_back( *ijec );
		}
		jec = new FactorizedJetCorrector(jecPars);
	}
	//end corrections
}
//
void NtupleJetFiller::SetBranches()
{
	//set branches
//	SetBranchSingle(&numGoodJets    ,"nref");
//	SetBranchSingle(&nPV            ,"nPV");
	SetBranch(jetpt    ,"jtpt" ,maxJets);
	SetBranch(jeteta   ,"jteta",maxJets);
	SetBranch(jetphi   ,"jtphi",maxJets);
	SetBranch(jetm     ,"jtm"  ,maxJets);
	SetBranch(jety     ,"jty"  ,maxJets);
	SetBranch(jetpt_uncorr,"rawpt",maxJets);
}
//
void NtupleJetFiller::fillBranch(const edm::Event& iEvent,const edm::InputTag& inputTag, 
		const edm::InputTag& srcPrimaryVertex, const bool runOnMC) 
{
	cout<<"filling brnca"<<endl;
	cout<<"Input tag "<<inputTag<<endl;
//	Handle<vector<reco::PFJet> > jets;
//	iEvent.getByLabel("ak7PFJets",jets);	
//	edm::Handle <edm::View<reco::Vertex> > recVtxs;
//	iEvent.getByLabel("offlinePrimaryVertices",recVtxs);
	nPV = 0.;
	double nPVval = 0;
	for(unsigned int ind=0;ind<recVtxs->size();ind++){
		if (!((*recVtxs)[ind].isFake()) && ((*recVtxs)[ind].ndof()>=4)
				&& (fabs((*recVtxs)[ind].z())<=24.0) &&
				((*recVtxs)[ind].position().Rho()<=2.0) ) {
			nPVval += 1;
		}
	}
	nPV = nPVval;

		const unsigned numJets(jetCollection->size());
		for(unsigned i(0);i<maxJets && i<numJets;++i)
		{ cout<<"running over "<<i <<endl;
			double pt = jetCollection->at(i).pt();
	    if(pt>20) { 
			// save ungroomed jet
			jetpt_uncorr[numGoodJets] = pt ;
			const TLorentzVector& jet_corr = getCorrectedJet(jetCollection->at(i));
			jetm[numGoodJets]   = jet_corr.M();
			jetpt[numGoodJets]  = jet_corr.Pt();
			jeteta[numGoodJets] = jet_corr.Eta();
			jetphi[numGoodJets] = jet_corr.Phi();
			jety[numGoodJets]   = jet_corr.Rapidity();
			++numGoodJets;
			}
		}//jet loop
	//cout<<"filling brnca"<<endl;
}
/////////Helper functions //////////////
void NtupleJetFiller::SetBranch(float* x, const char* name,unsigned maxJetsOrPhotons)
{
	std::ostringstream oss;
	oss<<name<<"["<<maxJetsOrPhotons<<"]"<<"/F";
	tree->Branch( name, x, oss.str().c_str() );
}
//
void NtupleJetFiller::SetBranch(int* x, const char* name,unsigned maxJetsOrPhotons)
{
	std::ostringstream oss;
	oss<<name<<"["<<maxJetsOrPhotons<<"]"<<"/I";
	tree->Branch( name, x, oss.str().c_str() );
}
//
void NtupleJetFiller::SetBranchSingle( float* x, std::string name)
{
	tree->Branch( name.c_str(), x, ( name+"/F").c_str() );
}

void NtupleJetFiller::SetBranchSingle( unsigned int* x, std::string name)
{
	tree->Branch( name.c_str(), x, ( name+"/I").c_str() );
}
//
const TLorentzVector NtupleJetFiller::getCorrectedJet(const reco::PFJet& jet) {
	double jecVal = 1.0;
	jecVal = getJEC( jet.eta(), jet.pt(), jet.energy());
	return(TLorentzVector (jet.px() * jecVal,
				jet.py() * jecVal,
				jet.pz() * jecVal,
				jet.energy() * jecVal));
}
//
double NtupleJetFiller::getJEC(double curJetEta,
		double curJetPt,
		double curJetE){
	jec->setJetEta( curJetEta );
	jec->setJetPt ( curJetPt );
	jec->setJetE  ( curJetE );
	jec->setNPV   ( nPV );
	double corr = jec->getCorrection();
	return corr;
}

#endif 
