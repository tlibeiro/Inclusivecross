#ifndef runJetAnalyzert_h
#define runJetAnalyzert_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <cassert>
#include "JetTree8TeV.h"
#include "hltanalysisHltTree.h"
#include "ak5PFJetAnalyzert.h"
#include "anaMETmetTree.h"
#include "skimanalysisHltTree.h"
#include "JetTreeAK5.h"
#include "JetTreeAK7.h"

using namespace std;
class runJetAnalyzer {
public:
	runJetAnalyzer(const vector<string>& infiles, const string& ofile,
      bool runmc,bool ak5 ,TTree *tree=0); 
	void Loop();

  private:
	///functions and variables 
   float*   jtpt ;    int  *   nPV;   float*   genpt ;
   float*   jteta;    float*   rho;   float*   geneta;
   float*   jtphi;    float*   met;   float*   genphi;
   float*   jtm  ;    float*   sumET; float*   genm  ;
   float*   jty  ;     int *   nref;
   float*   jecl1;
   float*   jecl2l3;
   float*   jecres;
   float*   rawpt;
   float*   chargedHardSum; float*   emEnergy ;
   float*   chargedHardN  ; float*   photonSum;
   float*   neutralHardSum; float*   muSum    ;
   float*   neutralHardN  ; float*   eSum     ;
	 int*     trackN;         int*     constituentN;          
	 int* run; int* event; int* lumiBlock;
   int*     IchargedHardN;  float* nuemEnergy, *chemEnergy;
   float*   qscale;
   float    pthat, maxpthat;
   ////trigger specific variables 
  static const unsigned numMaxTrgs = 10;
  int* hltPass[numMaxTrgs]; int* hltPsc[numMaxTrgs];
  int* l1Pass[numMaxTrgs];  int* l1Psc[numMaxTrgs];
  int* hltl1Psc[numMaxTrgs]; //l1 prescale for the associated hlt
  float* hltobj[numMaxTrgs]; //hlt trigger objects
  string* trgnames;
  int*    HLTJetPt;
  int*    HLTJetPtTh;
  unsigned numtrgs;
  ///for unfolding ntuple
  static const unsigned maxGenJets =14;  
  static const unsigned maxRecoJets=12;
  float jetptUnf[maxRecoJets];   float gjetptUnf[maxGenJets]; 
  float jetetaUnf[maxRecoJets];  float gjetetaUnf[maxGenJets];
  float jetphiUnf[maxRecoJets];  float gjetphiUnf[maxGenJets];
  float jetmUnf[maxRecoJets];    float gjetyUnf[maxGenJets];  
  int ngenUnf, njetUnf;

	TFile* f_out;
	TFile* f_unf;
  TTree* treeUnf;
	string outfile;
	hltanalysisHltTree* hltAnalyzer;
	ak5PFJetAnalyzert* ak5Analyzer;
  skimanalysisHltTree* skimAnalyzer;
	JetTreeAK7* treeak7;
	JetTreeAK7* treeak5;
  anaMETmetTree* metAnalyzer;
	bool runningOnIncxsec;
  const bool runOnMC;
  double mcWeight;
  void Init();
  void setMCWeight();
  bool jetId(unsigned jetnum);
  bool jetIdTagPrb(unsigned jetnum);
  void fillArrayFromVec(const vector<float>& vec,float* arr,
												const unsigned arrSize);
  ////declare Stability Histos==============================
  TH1* stabHists[500]; unsigned numStabHists;
  void makeStabilityHists();
  const vector<string> inputfiles;
  bool  runak5;
};
#endif

#ifdef runJetAnalyzer_cxx
runJetAnalyzer::runJetAnalyzer(const vector<string>& infiles, const string& ofile,
		bool runmc, bool ak5, TTree *tree) 
: outfile(ofile), runOnMC(runmc), runak5(ak5),mcWeight(-1), inputfiles(infiles)
{
	TChain * ak5chain = new TChain("ak5PFJetAnalyzer/t","");
	///for hlt 
	TChain * hltchain = new TChain("hltanalysis/HltTree","");
	///for treeak7 and ak5
	TChain * treeak7chain = new TChain("jetanalyzer/JetTree","");
	TChain * treeak5chain = new TChain("jetanalyzerak5/JetTree","");
  //for met from HI trees
	TChain * metchain = new TChain("anaMET/metTree","");
	//for Filters
  TChain * skimchain = new TChain("skimanalysis/HltTree"); 
	const unsigned numinput(infiles.size());
	for(unsigned i(0);i<numinput;++i)
	{
		string infile(infiles[i]);
		if(infile.find("pnfs")!=string::npos)
			infile = "dcache:"+infile;
		cout<<"Running on HI \n";
//		cout<<"Files Added "<<ak5chain->Add(infile.c_str())<<endl;
//		cout<<"Files Added HLT "<<hltchain->Add(infile.c_str())<<endl;
//		cout<<"Files Added met "<<metchain->Add(infile.c_str())<<endl;
//		cout<<"Files Added skim "<<skimchain->Add(infile.c_str())<<endl;
		cout<<"Files Added IncXsec ak7 "<<treeak7chain->Add(infile.c_str())<<endl;
		cout<<"Files Added IncXsec ak5 "<<treeak5chain->Add(infile.c_str())<<endl;
		cout<<infile<<endl;
	}
///jet tree
  ak5Analyzer = new ak5PFJetAnalyzert();
	ak5Analyzer->Init(ak5chain);
	//hlt
	hltAnalyzer = new hltanalysisHltTree();
	hltAnalyzer->Init(hltchain);
  ///for met info
  metAnalyzer = new anaMETmetTree();
	metAnalyzer->Init(metchain);
  ///for skimanalysis, HI filters
  skimAnalyzer = new skimanalysisHltTree();
  skimAnalyzer->Init(skimchain);
	/// for IncXsec
	treeak7 = new JetTreeAK7();
	treeak5 = new JetTreeAK7();
	treeak5->Init(treeak5chain);
	if(runak5)
		treeak7 = treeak5;
	else
		treeak7->Init(treeak7chain);

	///determine which dataset 
	const unsigned ak5entries(ak5Analyzer->fChain->GetEntries()),
				hltentries(hltAnalyzer->fChain->GetEntries()),
				metentries(metAnalyzer->fChain->GetEntries()),
				skmentries(skimAnalyzer->fChain->GetEntries()),
				treeentries(treeak7->fChain->GetEntries());
	cout<<"Running on Entries "<<ak5entries<<endl;
	//		cout<<"Hlt tree "<<hltentries<<endl;
	//		cout<<"met tree "<<metentries<<endl;
	//		cout<<"skim tree "<<skmentries<<endl;
	cout<<"IncXsec Tree "<<treeentries<<endl;
	if(ak5entries && hltentries)
	{
		assert(!treeentries);
		runningOnIncxsec=false;
	}
	else if (treeentries)
	{
		assert(!ak5entries && !hltentries);
		runningOnIncxsec=true;
	}
	//Initialize Unfold ntuple
	ostringstream mjets; mjets<<maxRecoJets;
	treeUnf = new TTree("UnfoldJetTree","UnfoldJetTree");
	treeUnf->Branch("njets",&njetUnf    , "njetUnf/I");
	treeUnf->Branch("jtpt" ,jetptUnf       ,("jtptUnf[" +mjets.str()+"]/F").c_str());
	treeUnf->Branch("jteta",jetetaUnf      ,("jtetaUnf["+mjets.str()+"]/F").c_str());
	treeUnf->Branch("jtphi",jetphiUnf      ,("jtphiUnf["+mjets.str()+"]/F").c_str());
	treeUnf->Branch("jtm"  ,jetmUnf        ,("jtmUnf["+mjets.str()+"]/F").c_str());
	mjets.str(""); mjets<<maxGenJets;
	treeUnf->Branch("ngen" ,&ngenUnf       , "ngenUnf/I");
	treeUnf->Branch("gjtpt" ,gjetptUnf       ,("gjtptUnf[" +mjets.str()+"]/F").c_str());
	treeUnf->Branch("gjteta",gjetetaUnf      ,("gjtetaUnf["+mjets.str()+"]/F").c_str());
	treeUnf->Branch("gjtphi",gjetphiUnf      ,("gjtphiUnf["+mjets.str()+"]/F").c_str());
	treeUnf->Branch("gjty"  ,gjetyUnf        ,("gjtyUnf["+mjets.str()+"]/F").c_str());
	treeUnf->Branch("mcweight" ,&mcWeight    , "mcWeight/D");
	///get mc weight for the sample, needed for unfolding 
	pthat= -1; maxpthat=-1;
	setMCWeight();
}

#endif
