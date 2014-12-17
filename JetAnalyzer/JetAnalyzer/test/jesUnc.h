#ifndef jesAnalyzert_h
#define jesAnalyzert_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <cassert>
#include "JetTreeAK7.h"
#include "AnalyzerProcedures.h"

using namespace std;
class jesAnalyzer {
public:
	jesAnalyzer(const vector<string>& infiles, const string& ofile,
      bool runmc, TTree *tree=0); 
	void Loop();

  private:
	///functions and variables 
   float*   jtpt ;    int  *   nPV;   float*   genpt ;
   float*   jteta;    float*   rho;   float*   geneta;
   float*   jtphi;    float*   met;   float*   genphi;
   float*   jtm  ;    float*   sumET; float*   genm  ;
   float*   jty  ;     int *   nref;
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
 
	TFile* f_out;
	string outfile;
	JetTreeAK7* treeak7;
	bool runningOnIncxsec;
  const bool runOnMC;
  double mcWeight;
  void Init();
  void setMCWeight();
  bool jetId(unsigned jetnum);
  void fillArrayFromVec(const vector<float>& vec,float* arr,
												const unsigned arrSize);
  const vector<string> inputfiles;
public:
  void makejecUncHists();
};

jesAnalyzer::jesAnalyzer(const vector<string>& infiles, const string& ofile,
		bool runmc, TTree *tree) 
: outfile(ofile), runOnMC(runmc), mcWeight(-1), inputfiles(infiles)
{
	///for treeak7 and ak5
	TChain * treeak7chain = new TChain("jetanalyzer/JetTree","");
	const unsigned numinput(infiles.size());
	for(unsigned i(0);i<numinput;++i)
	{
		string infile(infiles[i]);
		if(infile.find("pnfs")!=string::npos)
			infile = "dcache:"+infile;
		cout<<"Files Added IncXsec ak7 "<<treeak7chain->Add(infile.c_str())<<endl;
		cout<<infile<<endl;
	}
	/// for IncXsec
	treeak7 = new JetTreeAK7();
	treeak7->Init(treeak7chain);
	///determine which dataset 
	const unsigned treeentries(treeak7->fChain->GetEntries());
	cout<<"IncXsec Tree "<<treeentries<<endl;
	///get mc weight for the sample, needed for unfolding 
	pthat= -1; maxpthat=-1;
	setMCWeight();
}

#endif
