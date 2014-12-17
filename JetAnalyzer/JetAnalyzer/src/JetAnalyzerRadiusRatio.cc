// system include files
#include <memory>
#include <vector>
#include <iostream>
#include <algorithm>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
 #include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "matchOneToOne.hh"
#include "JetAnalyzer.h"

//#include "NtupleJetFiller.h"
//
// class declaration
//
using namespace reco;
using namespace std;
using namespace edm;
class JetAnalyzerRadiusRatio : public edm::EDAnalyzer {
   public:
      explicit JetAnalyzerRadiusRatio(const edm::ParameterSet&);
      ~JetAnalyzerRadiusRatio();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
			virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

			// ----------member data ---------------------------
			const edm::InputTag srcak7Jets, srcPrimaryVertex; 
			const edm::InputTag srcGenJets;
			const bool runOnMC;
      const bool saveEventInfo;
      const string corrLabel;
			const edm::ParameterSet iC;
			TTree* outTree;
			int run  ; 
			int event;
			int lumi ;
			int bunch;
			//			auto_ptr<NtupleJetFiller> jetfiller;
			string JEC_GlobalTag;
			static const unsigned maxJets=40;
			static const unsigned maxGenJets=40;

			edm::Handle<vector<reco::PFJet > > jetCollection;
			edm::Handle<vector<reco::GenJet> > genCollection;
			const std::string name;

			float jetpt [maxJets];   float jetpt_uncorr[maxJets];
			float jetm  [maxJets];   float jetphi[maxJets];
			float jeteta[maxJets];   float jety  [maxJets];  
      float jecl1   [maxJets];
      float jecl2l3 [maxJets];
      float jecres  [maxJets];
			float chargedHardSum [maxJets];  float neutralHardN [maxJets];
			float chargedHardN   [maxJets];  float emEnergy     [maxJets]; 
			float neutralHardSum [maxJets];  float muSum        [maxJets];
			float photonSum      [maxJets];  float eSum         [maxJets];
			float chemEnergy     [maxJets];  float nuemEnergy   [maxJets];
			int   trackN         [maxJets];  int   constituentN [maxJets];
			float matchGenPt     [maxJets];	 float matchGendr   [maxJets];
			float genpt          [maxGenJets];	 float genphi       [maxGenJets];
			float geneta         [maxGenJets];	 float genm         [maxGenJets];
      float qscale;

			const edm::InputTag metLabel;
			edm::InputTag l1TriggersTags;
			edm::InputTag hlTriggersTags;
			L1GtUtils   l1GtUtils;
			vector<std::string>   l1TriggerList ;
			vector<std::string>   hlTriggerList ;
			HLTConfigProvider hltConfigProv;
			double analyzeTrigger(const edm::Event& iEvent, 
					const edm::EventSetup& iSetup, 
					const std::string& triggerName);
			edm::Handle<edm::TriggerResults> hlTriggersObj;
			edm::Handle<trigger::TriggerEvent> summary;
			//
			unsigned nPV;
			float rhoVal, met,sumET;
			int l1prescl[maxJets]; int l1trg[maxJets];
			int hltprescl[maxJets];int hlttrg[maxJets];
			int hltl1prescl[maxJets];
			float hltPt[maxJets];
			FactorizedJetCorrector* jec;
			unsigned numGoodJets, numGoodGenJets;
			const TLorentzVector getCorrectedJet(const reco::PFJet& jet);
			void  getSubCorr(const reco::PFJet& jet, double* subc);
			double getJEC(double curJetEta,double curJetPt,double curJetE,double curJetArea);
      
      void Init();
			static const  bool useL1EventSetup = true;
			static const  bool useL1GtTriggerMenuLite = false;

};

JetAnalyzerRadiusRatio::JetAnalyzerRadiusRatio(const edm::ParameterSet& iConfig)
	:		srcak7Jets(iConfig.getParameter<edm::InputTag>("srcak7Jets")),
	srcGenJets(iConfig.getParameter<edm::InputTag>("srcGenJets")),
	srcPrimaryVertex(iConfig.getParameter<edm::InputTag>("srcPrimaryVertex")),
	JEC_GlobalTag(iConfig.getParameter<string>("JEC_GlobalTag")),
	runOnMC(iConfig.getParameter<bool>("runOnMC")),
	saveEventInfo(iConfig.getParameter<bool>("saveEventInfo")),
  corrLabel(iConfig.getParameter<string>("corrLabel")),
	metLabel(iConfig.getParameter<edm::InputTag>("metLabel")),
	iC(iConfig),outTree(0)
{
	l1TriggersTags = iConfig.getParameter<edm::InputTag>("l1TriggersTags");
	l1TriggerList  = iConfig.getParameter<vector<std::string>   >("l1TriggerList");
	hlTriggersTags = iConfig.getParameter<edm::InputTag>("hlTriggersTags");
	hlTriggerList  = iConfig.getParameter<vector<std::string>   >("hlTriggerList");
	/////////
	////set branches 
	ostringstream oss; oss<<maxJets; string njets=oss.str();
	edm::Service<TFileService> fs;
	outTree = fs->make<TTree>("JetTree","JetTree");
	if(saveEventInfo)
	{
		outTree->Branch("event_runNo",  &run,   "event_runNo/I");
		outTree->Branch("event_evtNo",  &event, "event_evtNo/I");
		outTree->Branch("event_lumi",   &lumi,  "event_lumi/I");
		outTree->Branch("event_bunch",  &bunch, "event_bunch/I");
		outTree->Branch("nPV",  &nPV,   "nPV/I");
		outTree->Branch("rho",  &rhoVal,"rhoVal/F");
		outTree->Branch("met",  &met,   "met/F");
		outTree->Branch("sumET",&sumET, "sumET/F");
	}
	else
	{
		outTree->Branch("njets", &numGoodJets, "njets/I");
		outTree->Branch("jtpt" ,jetpt       ,"jtpt[njets]/F");
		outTree->Branch("jteta",jeteta      ,"jteta[njets]/F");
		outTree->Branch("jtphi",jetphi      ,"jtphi[njets]/F");
		outTree->Branch("jtm"  ,jetm        ,"jtm[njets]/F");
		outTree->Branch("jty"  ,jety        ,"jty[njets]/F");
		outTree->Branch("rawpt",jetpt_uncorr,"rawpt[njets]/F");
		outTree->Branch("jecl1"  ,jecl1     ,"jecl1[njets]/F");
		outTree->Branch("jecl2l3",jecl2l3   ,"jecl2l3[njets]/F");
		outTree->Branch("jecres" ,jecres    ,"jecres[njets]/F");
		outTree->Branch("chargedHardSum",chargedHardSum,"chargedHardSum[njets]/F");
		outTree->Branch("chargedHardN"  ,chargedHardN  ,"chargedHardN[njets]/F");
		outTree->Branch("neutralHardSum",neutralHardSum,"neutralHardSum[njets]/F");
		outTree->Branch("neutralHardN"  ,neutralHardN  ,"neutralHardN[njets]/F");
		outTree->Branch("emEnergy"      ,emEnergy      ,"emEnergy[njets]/F");
		outTree->Branch("chemEnergy"    ,chemEnergy    ,"chemEnergy[njets]/F");
		outTree->Branch("nuemEnergy"    ,nuemEnergy    ,"nuemEnergy[njets]/F");
		outTree->Branch("muSum"         ,muSum         ,"muSum[njets]/F");
		outTree->Branch("eSum"          ,eSum          ,"eSum[njets]/F");
		outTree->Branch("photonSum"     ,photonSum     ,"photonSum[njets]/F");
		outTree->Branch("trackN"        ,trackN        ,"trackN[njets]/I");
		outTree->Branch("constituentN"  ,constituentN  ,"constituentN[njets]/I");
		if(runOnMC) {
			outTree->Branch("matchgenpt",matchGenPt,"matchgenpt[njets]/F");
			outTree->Branch("matchgendr",matchGendr,"matchgendr[njets]/F");
			outTree->Branch("ngen",   &numGoodGenJets, "ngen/I");
			outTree->Branch("qscale", &qscale, "qscale/F");
			oss.str("");oss<<maxGenJets; string ngen=oss.str();
			outTree->Branch("genpt" ,genpt       ,"genpt[ngen]/F");
			outTree->Branch("geneta",geneta      ,"geneta[ngen]/F");
			outTree->Branch("genphi",genphi      ,"genphi[ngen]/F");
			outTree->Branch("genm"  ,genm        ,"genm[ngen]/F");
		}//if running on mc , save gen info
	}
	if(saveEventInfo) {
		const unsigned nl1trgs(l1TriggerList.size()), 
					nhlttrgs(hlTriggerList.size());
		for(unsigned i(0);i<nl1trgs;++i) {
			outTree->Branch((l1TriggerList[i]+"_Prescl").c_str(),&l1prescl[i],
					(l1TriggerList[i]+"_Prescl/I").c_str());
			outTree->Branch((l1TriggerList[i]).c_str(),&l1trg[i],
					(l1TriggerList[i]+"/I").c_str());
		}
		for(unsigned i(0);i<nhlttrgs;++i) {
			outTree->Branch((hlTriggerList[i]+"_Prescl").c_str(),&hltprescl[i],
					(hlTriggerList[i]+"_Prescl/I").c_str());
			outTree->Branch((hlTriggerList[i]).c_str(),&hlttrg[i],
					(hlTriggerList[i]+"/I").c_str());
			outTree->Branch((hlTriggerList[i]+"L1_Prescl").c_str(),&hltl1prescl[i],
					(hlTriggerList[i]+"L1_Prescl/I").c_str());
			outTree->Branch((hlTriggerList[i]+"_ObjPt").c_str(),&hltPt[i],
					(hlTriggerList[i]+"_ObjPt/F").c_str());
		}
	}
	// ---- setting up the jec on-the-fly from text files...
	const std::string& fDir = JEC_GlobalTag;
	std::vector< JetCorrectorParameters > jecPars;
	std::vector< std::string > jecStr;
//	jecStr.push_back( fDir + "_L1FastJet_"+corrLabel+".txt" );
	jecStr.push_back( fDir + "_L2Relative_"+corrLabel+".txt" );
	jecStr.push_back( fDir + "_L3Absolute_"+corrLabel+".txt" );
	if (!runOnMC)
		jecStr.push_back( fDir + "_L2L3Residual_"+corrLabel+".txt" );

	for (unsigned int i = 0; i < jecStr.size(); ++i){
		cout<<jecStr[i]<<endl;
		JetCorrectorParameters* ijec = new JetCorrectorParameters( jecStr[i] );
		jecPars.push_back( *ijec );
	}
	jec = new FactorizedJetCorrector(jecPars);
	//end corrections
}

JetAnalyzerRadiusRatio::~JetAnalyzerRadiusRatio()
{
}

	void
JetAnalyzerRadiusRatio::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	const unsigned nl1trgs(l1TriggerList.size()),
				nhlttrgs(hlTriggerList.size());
	Init();
	// 	// level 1 triggers
	edm::Handle<L1GlobalTriggerReadoutRecord> l1TriggersObj;
	iEvent.getByLabel(l1TriggersTags,l1TriggersObj);
	edm::ESHandle<L1GtTriggerMenu> l1TriggerMenu;
	iSetup.get<L1GtTriggerMenuRcd>().get(l1TriggerMenu);
	DecisionWord const algoWord = l1TriggersObj->decisionWord();
	l1GtUtils.getL1GtRunCache(iEvent, iSetup, useL1EventSetup, useL1GtTriggerMenuLite);
	for (size_t l = 0; l< nl1trgs; ++l) {
		AlgorithmMap const algoMap = l1TriggerMenu->gtAlgorithmMap();
		AlgorithmMap::const_iterator it = algoMap.find(l1TriggerList.at(l));
		if (it != algoMap.end()) {
			int bitNumber = (it->second).algoBitNumber();
			bool l1trgval = algoWord.at(bitNumber);
			l1trg[l]=(int) l1trgval;
			if(l1trgval) {
				int errorCode = 0;
				int preScale = l1GtUtils.prescaleFactor(iEvent,it->second.algoName(), errorCode);
				l1prescl[l]=preScale;
				if(false) {
					cout
						<<l1TriggerList[l]<<' '<<l1trgval
						<<" prescale "<<preScale
						<<" errorcode "<<errorCode
						<<endl;
				}
			}
		}
	}///l1 loop
	//
	// high level triggers
	iEvent.getByLabel(hlTriggersTags, hlTriggersObj);
	edm::TriggerNames triggerNames = iEvent.triggerNames(*hlTriggersObj);
	iEvent.getByLabel("hltTriggerSummaryAOD",summary);
	for (size_t l = 0; l<hlTriggerList.size(); ++l) 
		for (size_t t = 0; t<hlTriggersObj->size(); ++t) 
			if (triggerNames.triggerName(t).find(hlTriggerList.at(l))!=string::npos) 
			{
				bool hlttrgval = hlTriggersObj->accept(t);
				hlttrg[l]=(int)hlttrgval;
				if (hlttrgval) {
					int preScale = hltConfigProv.prescaleValue(iEvent, iSetup,triggerNames.triggerName(t));
					hltprescl[l] = preScale ;
					pair<int,int> pscls =
						hltConfigProv.prescaleValues(iEvent,iSetup,triggerNames.triggerName(t)) ;
					hltl1prescl[l] = pscls.first;
					hltPt[l] = analyzeTrigger(iEvent,iSetup,triggerNames.triggerName(t));
				}
				if(false) {
					cout
						<<"-----------------------\n"
						<<triggerNames.triggerName(t)<<' '<<hlttrg[l]
						<<' '<<hlTriggerList.at(l)
						<<" prescale "<<hltl1prescl[l]
						<<" l "<<l<<" t "<<t
						<<endl;
				}
			}
	// write event information: run, event, bunch crossing, ....
	run   = iEvent.id().run();
	event = iEvent.id().event();
	lumi  = iEvent.luminosityBlock();
	bunch = iEvent.bunchCrossing();
	iEvent.getByLabel(srcak7Jets,jetCollection);	
	if(runOnMC)
		iEvent.getByLabel(srcGenJets,genCollection);	
	////// fro rho 
	edm::Handle<double> rho;
	const edm::InputTag eventrho("kt6PFJets", "rho");
	iEvent.getByLabel(eventrho,rho);
	rhoVal = *rho;
	/////for met 
	double mpfMET(-1),mpfSumET(-1);
	edm::Handle<edm::View<reco::MET> > pfmet;
	iEvent.getByLabel(metLabel, pfmet);
	if (pfmet->size() != 0) 
	{
		met      = (*pfmet)[0].et();
		sumET    = (*pfmet)[0].sumEt();
	}
	//	////triggers//////////
	///////////////get pile up infor
	edm::Handle <edm::View<reco::Vertex> > recVtxs;
	iEvent.getByLabel(srcPrimaryVertex,recVtxs);	
	nPV = 0.;
	double nPVval = 0;
	for(unsigned int ind=0;ind<recVtxs->size();ind++){
		if (!((*recVtxs)[ind].isFake()) && ((*recVtxs)[ind].ndof()>=4)
				&& (fabs((*recVtxs)[ind].z())<=24.0) ) 
			nPVval += 1;
	}
	nPV = nPVval;
	///////////////get pile up infor
	numGoodJets=0;
	vector<JA::Jet> recjets, genjets;
	const unsigned numJets(jetCollection->size());
	for(unsigned i(0);i<maxJets && i<numJets;++i)
	{
		reco::PFJet uncorrjet(jetCollection->at(i));
		const TLorentzVector& jet_corr = getCorrectedJet(uncorrjet);
		double pt = jet_corr.Pt();
		// get sub corrections
		double subc[3] = {0.,0.,0.};
		getSubCorr(uncorrjet,subc);
		cout<<fixed<<setprecision(2);
    if(false)   
		cout<<" jet pt/eta "<<uncorrjet.pt()<<"|"<<uncorrjet.eta()
      //  <<" "<<jet_corr.Pt()/uncorrjet.pt()
						<<" l1|l2l3|res "
						<<subc[0]<<"|"<<subc[1]<<"|"<<subc[2]
						<<endl;
		if(pt>5) { 
			// save ungroomed jet
			jetpt_uncorr[numGoodJets] = uncorrjet.pt() ;
			jetm[numGoodJets]   = jet_corr.M();
			jetpt[numGoodJets]  = jet_corr.Pt();
			jeteta[numGoodJets] = jet_corr.Eta();
			jetphi[numGoodJets] = jet_corr.Phi();
			jety[numGoodJets]   = jet_corr.Rapidity();
			jecl1[numGoodJets]  = subc[0];
			jecl2l3[numGoodJets]= subc[1];
			jecres[numGoodJets] = subc[2];
			chargedHardSum [numGoodJets]   = uncorrjet.chargedHadronEnergyFraction();
			chargedHardN   [numGoodJets]   = uncorrjet.chargedHadronMultiplicity();
			//			neutralHardSum [numGoodJets]   = uncorrjet.neutralHadronEnergyFraction();
			neutralHardSum [numGoodJets]   = (uncorrjet.neutralHadronEnergy() + 
					uncorrjet.HFHadronEnergy())/uncorrjet.energy();
			neutralHardN   [numGoodJets]   = uncorrjet.neutralHadronMultiplicity();
			emEnergy       [numGoodJets]   = uncorrjet.chargedEmEnergyFraction()+
				uncorrjet.neutralEmEnergyFraction();
			chemEnergy     [numGoodJets]   = uncorrjet.chargedEmEnergyFraction();
			nuemEnergy     [numGoodJets]   = uncorrjet.neutralEmEnergyFraction();
			muSum          [numGoodJets]   = uncorrjet.muonEnergyFraction();
			eSum           [numGoodJets]   = uncorrjet.electronEnergyFraction();
			photonSum      [numGoodJets]   = uncorrjet.photonEnergyFraction();
			constituentN   [numGoodJets]   = uncorrjet.chargedMultiplicity()+
				uncorrjet.neutralMultiplicity();
			//			trackN         [numGoodJets]   = uncorrjet.getTrackRefs().size();
			//---- get the vector of tracks -----
			trackN[numGoodJets] = 0;
			reco::TrackRefVector vTrks(uncorrjet.getTrackRefs());
			for(reco::TrackRefVector::const_iterator i_trk = vTrks.begin(); i_trk != vTrks.end(); i_trk++) 
				if((*i_trk)->pt()) ++trackN[numGoodJets];

			++numGoodJets;
			recjets.push_back(JA::Jet(jet_corr.Pt(),jet_corr.Eta(),
						jet_corr.Phi(),jet_corr.M()));
			if(false)
				cout<<"jet "<<i<<" "
					<<jet_corr.Pt()
					<<endl;
		}
	}//jet loop
	if(runOnMC) {
		edm::Handle<GenEventInfoProduct>    genEventScale;
		if (iEvent.getByLabel("generator", genEventScale))
			qscale = genEventScale->qScale();

		const unsigned numGenJets(genCollection->size());
		numGoodGenJets=0;
		for(unsigned i(0);i<maxGenJets && i<numGenJets;++i)
		{
			reco::GenJet jet(genCollection->at(i));
			double pt = jet.pt();
			if(pt>5) { 
				// save ungroomed jet
				genm[numGoodGenJets]   = jet.mass();
				genpt[numGoodGenJets]  = jet.pt();
				geneta[numGoodGenJets] = jet.eta();
				genphi[numGoodGenJets] = jet.phi();
				++numGoodGenJets;
				genjets.push_back(JA::Jet(jet.pt(),jet.eta(),
							jet.phi(),jet.mass()));
			}
		}//genjet loop
		//match recojets to gen jets
		JA::DeltaRDistance d; 
		vector<int> matchVec;
		matchOneToOne(recjets,genjets,JA::DeltaRDistance() ,&matchVec,1000.0);  
		for(unsigned i(0);i<numGoodJets;++i)
			if(matchVec[i]>-1)
			{
				matchGenPt[i] = genjets[matchVec[i]].pt ;
				matchGendr[i] = d(genjets[matchVec[i]],recjets[i]) ;
			}
	}///if running on mc write gen information
	//Fill
	outTree->Fill();
}

// ------------ helper functions  ------------
const TLorentzVector JetAnalyzerRadiusRatio::getCorrectedJet(const reco::PFJet& jet) {
	double jecVal = 1.0;
	jecVal = getJEC( jet.eta(), jet.pt(), jet.energy(),jet.jetArea());
	//  cout<<"max distance "<<jet.maxDistance()<<" "<<jet.pt()<<endl;
	return(TLorentzVector (jet.px() * jecVal,
				jet.py() * jecVal,
				jet.pz() * jecVal,
				jet.energy() * jecVal));
}
//
void JetAnalyzerRadiusRatio::getSubCorr(const reco::PFJet& jet, double* subc) {
	//initialize array
	subc[0]=0;  subc[1]=0;  subc[2]=0;

	double jecVal = 1.0;
	jecVal = getJEC( jet.eta(), jet.pt(), jet.energy(),jet.jetArea());
	jec->setJetEta( jet.eta());
	jec->setJetPt ( jet.pt() );
	jec->setJetE  ( jet.energy() );
	jec->setNPV   ( nPV );
	jec->setJetA  ( jet.jetArea() );
	jec->setRho   ( rhoVal );
	vector<float> v = jec->getSubCorrections();
	subc[0] = v[0];
	subc[1] = v[2]/v[0];
	subc[2] = ((!runOnMC) ? v[3]/v[2] : 1.);
	assert(jecVal == v[v.size()-1]);
   if(false)   
		cout<<"========================================\n"
<<" jet pt/eta "<<jet.pt()<<"|"<<jet.eta()<<"\n"
            <<"v0|v1|v2 "<<v[0]<<"|"<<v[1]<<"|"<<v[2]<<"\n"
						<<" l1|l2l3|res "
						<<subc[0]<<"|"<<subc[1]<<"|"<<subc[2]
						<<endl;

}

double JetAnalyzerRadiusRatio::getJEC(double curJetEta,
		double curJetPt,
		double curJetE, double curJetArea){
	jec->setJetEta( curJetEta );
	jec->setJetPt ( curJetPt );
	jec->setJetE  ( curJetE );
	jec->setNPV   ( nPV );
	jec->setJetA  ( curJetArea );
	jec->setRho   ( rhoVal );
	double corr = jec->getCorrection();
	return corr;
}
// ------------ method called once each job just before starting event loop  ------------
	void 
JetAnalyzerRadiusRatio::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
JetAnalyzerRadiusRatio::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
	void 
JetAnalyzerRadiusRatio::beginRun(edm::Run const& run, edm::EventSetup const& setup)
{
	l1GtUtils.getL1GtRunCache(run, setup, useL1EventSetup, useL1GtTriggerMenuLite);
	// high level triggers
	bool isChanged = true;
	hltConfigProv.init(run, setup, hlTriggersTags.process(), isChanged);
}

// ------------ method called when ending the processing of a run  ------------
	void 
JetAnalyzerRadiusRatio::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
	void 
JetAnalyzerRadiusRatio::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
	void 
JetAnalyzerRadiusRatio::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

double JetAnalyzerRadiusRatio::analyzeTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string& triggerName) {

	using namespace std;
	using namespace edm;
	using namespace reco;
	using namespace trigger;
	bool printOutput(false);
	if(printOutput)
		cout << endl;

	iEvent.getByLabel(hlTriggersTags, hlTriggersObj);
	iEvent.getByLabel("hltTriggerSummaryAOD",summary);
	const unsigned int n(hltConfigProv.size());
	const unsigned int triggerIndex(hltConfigProv.triggerIndex(triggerName));
	assert(triggerIndex==iEvent.triggerNames(*hlTriggersObj).triggerIndex(triggerName));

	// abort on invalid trigger name
	if (triggerIndex>=n) {
		if(printOutput)
			cout << "HLTEventAnalyzerAOD::analyzeTrigger: path "
				<< triggerName << " - not found!" << endl;
		return -1.;
	}

	const std::pair<int,int> prescales(hltConfigProv.prescaleValues(iEvent,iSetup,triggerName));
	if(printOutput)
		cout << "HLTEventAnalyzerAOD::analyzeTrigger: path "
			<< triggerName << " [" << triggerIndex << "] "
			<< "prescales L1T,HLT: " << prescales.first << "," << prescales.second
			<< endl;
	//	const std::pair<std::vector<std::pair<std::string,int> >,int> prescalesInDetail(hltConfigProv.prescaleValuesInDetail(iEvent,iSetup,triggerName));
	//	std::ostringstream message;
	//	for (unsigned int i=0; i<prescalesInDetail.first.size(); ++i) {
	//		message << " " << i << ":" << prescalesInDetail.first[i].first << "/" << prescalesInDetail.first[i].second;
	//	}
	//	cout << "HLTEventAnalyzerAOD::analyzeTrigger: path "
	//		<< triggerName << " [" << triggerIndex << "] "
	//		<< endl
	//		<< "prescales L1T: " << prescalesInDetail.first.size() <<  message.str()
	//		<< endl
	//		<< "prescale  HLT: " << prescalesInDetail.second
	//		<< endl;

	// modules on this trigger path
	const unsigned int m(hltConfigProv.size(triggerIndex));
	const vector<string>& moduleLabels(hltConfigProv.moduleLabels(triggerIndex));

	// Results from TriggerResults product
	if(printOutput)
		cout << " Trigger path status:"
			<< " WasRun=" << hlTriggersObj->wasrun(triggerIndex)
			<< " Accept=" << hlTriggersObj->accept(triggerIndex)
			<< " Error =" << hlTriggersObj->error(triggerIndex)
			<< endl;
	const unsigned int moduleIndex(hlTriggersObj->index(triggerIndex));
	if(printOutput)
		cout << " Last active module - label/type: "
			<< moduleLabels[moduleIndex] << "/" << hltConfigProv.moduleType(moduleLabels[moduleIndex])
			<< " [" << moduleIndex << " out of 0-" << (m-1) << " on this path]"
			<< endl;
	assert (moduleIndex<m);

	vector<double> hltobjpt;
	// Results from TriggerEvent product - Attention: must look only for
	// modules actually run in this path for this event!
	for (unsigned int j=0; j<=moduleIndex; ++j) {
		const string& moduleLabel(moduleLabels[j]);
		const string  moduleType(hltConfigProv.moduleType(moduleLabel));
		// check whether the module is packed up in TriggerEvent product
		const unsigned int filterIndex(summary->filterIndex(InputTag(moduleLabel,"","HLT")));
		if (filterIndex<summary->sizeFilters()) {
			if(printOutput)
				cout << " 'L3' filter in slot " << j << " - label/type " << moduleLabel << "/" << moduleType << endl;
			const Vids& VIDS (summary->filterIds(filterIndex));
			const Keys& KEYS(summary->filterKeys(filterIndex));
			const size_type nI(VIDS.size());
			const size_type nK(KEYS.size());
			assert(nI==nK);
			const size_type n(max(nI,nK));
			if(printOutput)
				cout << "   " << n  << " accepted 'L3' objects found: " << endl;
			const TriggerObjectCollection& TOC(summary->getObjects());
			for (size_type i=0; i!=n; ++i) {
				const TriggerObject& TO(TOC[KEYS[i]]);
				if(printOutput)
					cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
						<< TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass()
						<<"\n module Type "<<moduleType
						<< endl;
				if(moduleType.find("HLT1CaloJet")==0) {
					hltobjpt.push_back(TO.pt());
				}//
			}
		}///filters
	}///modules

	//  for(unsigned i(0);i<hltobjpt.size();++i) 	
	//		cout<<"\n HLT Objects "<<hltobjpt[i]<<endl;
	//sort the vector----------------------------
	sort(hltobjpt.begin(),hltobjpt.end());

	return (*(hltobjpt.end()-1));
}


void JetAnalyzerRadiusRatio::Init()
{
	numGoodJets=0; met=0; sumET=0;numGoodGenJets=0;
	qscale=-1;
	for(unsigned i(0);i<maxJets;++i)
	{
		jetpt[i] =-1;   jetpt_uncorr[i]=-1; 
		jetm[i]  =-1;  	jetphi[i]=-1;
		jeteta[i]=-1;  	jety[i]  =-1;
		chargedHardSum [i]=-1;  emEnergy  [i]=-1;
		chargedHardN   [i]=-1;  muSum     [i]=-1;
		neutralHardSum [i]=-1;  photonSum [i]=-1;
		neutralHardN   [i]=-1;  eSum      [i]=-1;
		chemEnergy     [i]=-1;  nuemEnergy [i]=-1;
		trackN         [i]=-1;  constituentN[i]=-1;
		matchGenPt     [i]=-1;  matchGendr     [i]=-1;
		jecl1          [i]=-1;
		jecl2l3        [i]=-1;
		jecres         [i]=-1;
		///trgs
		l1prescl[i]=-10; l1trg[i]=-10;
		hltprescl[i]=-10;hlttrg[i]=-10;
		hltl1prescl[i]=-10; hltPt[i]=-10.0;
	}
	for(unsigned i(0);i<maxGenJets;++i)
	{
		genpt [i]=-1; geneta[i]=-1;
		genphi[i]=-1; genm  [i]=-1;
	}
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetAnalyzerRadiusRatio);
