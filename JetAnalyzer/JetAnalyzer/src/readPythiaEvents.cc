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
//#include "readPythiaEvents.h"

// class declaration
using namespace reco;
using namespace std;
using namespace edm;
class readPythiaEvents : public edm::EDAnalyzer {
   public:
      explicit readPythiaEvents(const edm::ParameterSet&);
      ~readPythiaEvents();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
			virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

			// ----------member data ---------------------------
			const edm::InputTag srcGenJets;
			TTree* outTree;
			void Init();
      static const unsigned maxGenJets=40;
			unsigned numGoodGenJets;
      float genpt    [maxGenJets]; float genphi [maxGenJets];
			float geneta   [maxGenJets]; float genm   [maxGenJets];
      float geny     [maxGenJets]; float nConstituent[maxGenJets];
      float qscale;
      float weight;
      float q_binLimit, qtmp;
			edm::Handle<vector<reco::GenJet> > genCollection;
      unsigned jetsgt0, jetsgt5, jetsgt10, events;

};

readPythiaEvents::readPythiaEvents(const edm::ParameterSet& iConfig)
	:	srcGenJets(iConfig.getParameter<edm::InputTag>("srcGenJets")),
	  qtmp(iConfig.getParameter<double>("q_binLimit")),
		outTree(0)
{
jetsgt0=jetsgt5=jetsgt10=events=0;
	////set branches 
	edm::Service<TFileService> fs;
	outTree = fs->make<TTree>("JetTree","JetTree");
	outTree->Branch("ngen",   &numGoodGenJets, "ngen/I");
	outTree->Branch("qscale", &qscale, "qscale/F");
	outTree->Branch("weight", &weight, "weight/F");
	outTree->Branch("q_binLimit", &q_binLimit, "q_binLimit/F");
	ostringstream oss;
	oss.str("");oss<<maxGenJets; string ngen=oss.str();
	outTree->Branch("genpt" ,genpt       ,("genpt[" +ngen+"]/F").c_str());
	outTree->Branch("geneta",geneta      ,("geneta["+ngen+"]/F").c_str());
	outTree->Branch("geny"  ,geny        ,("geny["  +ngen+"]/F").c_str());
	outTree->Branch("genphi",genphi      ,("genphi["+ngen+"]/F").c_str());
	outTree->Branch("genm"  ,genm        ,("genm["  +ngen+"]/F").c_str());
	outTree->Branch("nConstituent"  ,nConstituent        ,("nConstituent["  +ngen+"]/F").c_str());
}

readPythiaEvents::~readPythiaEvents()
{
cout<<" All jets | jets pt > 5 | jets pt > 10 "<<endl;
cout<<jetsgt0<<" "<<jetsgt5<<" "<<jetsgt10<<" "<<events<<endl;
}

	void
readPythiaEvents::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	// write event information: run, event, bunch crossing, ....
		iEvent.getByLabel(srcGenJets,genCollection);	
jetsgt0+=genCollection->size();
    for(unsigned i(0); i<genCollection->size();++i)
		{
      if( genCollection->at(i).pt()>5) jetsgt5+=1;
      if( genCollection->at(i).pt()>10) jetsgt10+=1;
		}
events+=1;




	Init();
	const unsigned numJets(genCollection->size());
  edm::Handle<GenEventInfoProduct>    genEventScale;
  if (iEvent.getByLabel("generator", genEventScale))
		{
      qscale = genEventScale->qScale();
      weight = genEventScale->weight();
		}

  q_binLimit = qtmp;

		const unsigned numGenJets(genCollection->size());
 //   cout<<"Num Gen Jets "<<genCollection->size()<<endl;
		numGoodGenJets=0;
		for(unsigned i(0);i<maxGenJets && i<numGenJets;++i)
		{
			reco::GenJet jet(genCollection->at(i));
			double pt = jet.pt();
       if(pt>5) {
				// save ungroomed jet
				nConstituent[numGoodGenJets]   = jet.getGenConstituents().size();
				genm[numGoodGenJets]   = jet.mass();
				genpt[numGoodGenJets]  = jet.pt();
				geneta[numGoodGenJets] = jet.eta();
				genphi[numGoodGenJets] = jet.phi();
        TLorentzVector v; 
				v.SetPtEtaPhiM(jet.pt(),jet.eta(),jet.phi(),jet.mass());
				geny[numGoodGenJets] = v.Rapidity();
				++numGoodGenJets;
			}
		}//genjet loop
	//Fill
	outTree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
	void 
readPythiaEvents::beginJob()
{
}
// ------------ method called once each job just after ending the event loop  ------------
	void 
readPythiaEvents::endJob() 
{
}
// ------------ method called when starting to processes a run  ------------
	void 
readPythiaEvents::beginRun(edm::Run const& run, edm::EventSetup const& setup)
{
}
// ------------ method called when ending the processing of a run  ------------
	void 
readPythiaEvents::endRun(edm::Run const&, edm::EventSetup const&)
{
}
// ------------ method called when starting to processes a luminosity block  ------------
	void 
readPythiaEvents::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
// ------------ method called when ending the processing of a luminosity block  ------------
	void 
readPythiaEvents::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
//Init function
void readPythiaEvents::Init()
{
	numGoodGenJets=0;
	qscale=-1;
	weight=-1;
	q_binLimit=-1;
	for(unsigned i(0);i<maxGenJets;++i)
	{
		genpt [i]=-1; geneta[i]=-1;
		genphi[i]=-1; genm  [i]=-1;
		geny[i]=-1;   nConstituent[i]=-1;
	}
}

//define this as a plug-in
DEFINE_FWK_MODULE(readPythiaEvents);
