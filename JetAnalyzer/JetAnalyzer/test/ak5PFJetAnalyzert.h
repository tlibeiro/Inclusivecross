#ifndef ak5PFJetAnalyzert_h
#define ak5PFJetAnalyzert_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <iostream>
#include <cassert>
#include "hltanalysisHltTree.h"
#include "JetTree8TeV.h"
// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

using namespace std;
class ak5PFJetAnalyzert {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           evt;
   Float_t         b;
   Int_t           nref;
   Float_t         rawpt[30];   //[nref]
   Float_t         jtpt[30];   //[nref]
   Float_t         jteta[30];   //[nref]
   Float_t         jty[30];   //[nref]
   Float_t         jtphi[30];   //[nref]
   Float_t         jtpu[30];   //[nref]
   Float_t         jtm[30];   //[nref]
   Float_t         discr_fr01[30];   //[nref]
   Float_t         trackMax[30];   //[nref]
   Float_t         trackSum[30];   //[nref]
   Int_t           trackN[30];   //[nref]
   Float_t         trackHardSum[30];   //[nref]
   Int_t           trackHardN[30];   //[nref]
   Float_t         chargedMax[30];   //[nref]
   Float_t         chargedSum[30];   //[nref]
   Int_t           chargedN[30];   //[nref]
   Float_t         chargedHardSum[30];   //[nref]
   Int_t           chargedHardN[30];   //[nref]
   Float_t         photonMax[30];   //[nref]
   Float_t         photonSum[30];   //[nref]
   Int_t           photonN[30];   //[nref]
   Float_t         photonHardSum[30];   //[nref]
   Int_t           photonHardN[30];   //[nref]
   Float_t         neutralMax[30];   //[nref]
   Float_t         neutralSum[30];   //[nref]
   Int_t           neutralN[30];   //[nref]
   Float_t         hcalSum[30];   //[nref]
   Float_t         ecalSum[30];   //[nref]
   Float_t         eMax[30];   //[nref]
   Float_t         eSum[30];   //[nref]
   Int_t           eN[30];   //[nref]
   Float_t         muMax[30];   //[nref]
   Float_t         muSum[30];   //[nref]
   Int_t           muN[30];   //[nref]
   Float_t         fHPD[30];   //[nref]
   Float_t         fRBX[30];   //[nref]
   Int_t           n90[30];   //[nref]
   Float_t         fSubDet1[30];   //[nref]
   Float_t         fSubDet2[30];   //[nref]
   Float_t         fSubDet3[30];   //[nref]
   Float_t         fSubDet4[30];   //[nref]
   Float_t         restrictedEMF[30];   //[nref]
   Int_t           nHCAL[30];   //[nref]
   Int_t           nECAL[30];   //[nref]
   Float_t         apprHPD[30];   //[nref]
   Float_t         apprRBX[30];   //[nref]
   Int_t           n2RPC[30];   //[nref]
   Int_t           n3RPC[30];   //[nref]
   Int_t           nRPC[30];   //[nref]
   Float_t         fEB[30];   //[nref]
   Float_t         fEE[30];   //[nref]
   Float_t         fHB[30];   //[nref]
   Float_t         fHE[30];   //[nref]
   Float_t         fHO[30];   //[nref]
   Float_t         fLong[30];   //[nref]
   Float_t         fShort[30];   //[nref]
   Float_t         fLS[30];   //[nref]
   Float_t         fHFOOT[30];   //[nref]
   Float_t         matchedPt[30];   //[nref]
   Float_t         matchedRawPt[30];   //[nref]
   Float_t         matchedR[30];   //[nref]
   Float_t         refpt[30];   //[nref]
   Float_t         refeta[30];   //[nref]
   Float_t         refy[30];   //[nref]
   Float_t         refphi[30];   //[nref]
   Float_t         refdphijt[30];   //[nref]
   Float_t         refdrjt[30];   //[nref]
   Float_t         refparton_pt[30];   //[nref]
   Int_t           refparton_flavor[30];   //[nref]
   Int_t           refparton_flavorForB[30];   //[nref]
   Float_t         genChargedSum[30];   //[nref]
   Float_t         genHardSum[30];   //[nref]
   Float_t         signalChargedSum[30];   //[nref]
   Float_t         signalHardSum[30];   //[nref]
   Int_t           subid[30];   //[nref]
   Int_t           ngen;
   Int_t           genmatchindex[38];   //[ngen]
   Float_t         genpt[38];   //[ngen]
   Float_t         geneta[38];   //[ngen]
   Float_t         geny[38];   //[ngen]
   Float_t         genphi[38];   //[ngen]
   Float_t         gendphijt[38];   //[ngen]
   Float_t         gendrjt[38];   //[ngen]
   Int_t           gensubid[38];   //[ngen]


   // List of branches
   TBranch        *b_evt;   //!
   TBranch        *b_b;   //!
   TBranch        *b_nref;   //!
   TBranch        *b_rawpt;   //!
   TBranch        *b_jtpt;   //!
   TBranch        *b_jteta;   //!
   TBranch        *b_jty;   //!
   TBranch        *b_jtphi;   //!
   TBranch        *b_jtpu;   //!
   TBranch        *b_jtm;   //!
   TBranch        *b_discr_fr01;   //!
   TBranch        *b_trackMax;   //!
   TBranch        *b_trackSum;   //!
   TBranch        *b_trackN;   //!
   TBranch        *b_trackHardSum;   //!
   TBranch        *b_trackHardN;   //!
   TBranch        *b_chargedMax;   //!
   TBranch        *b_chargedSum;   //!
   TBranch        *b_chargedN;   //!
   TBranch        *b_chargedHardSum;   //!
   TBranch        *b_chargedHardN;   //!
   TBranch        *b_photonMax;   //!
   TBranch        *b_photonSum;   //!
   TBranch        *b_photonN;   //!
   TBranch        *b_photonHardSum;   //!
   TBranch        *b_photonHardN;   //!
   TBranch        *b_neutralMax;   //!
   TBranch        *b_neutralSum;   //!
   TBranch        *b_neutralN;   //!
   TBranch        *b_hcalSum;   //!
   TBranch        *b_ecalSum;   //!
   TBranch        *b_eMax;   //!
   TBranch        *b_eSum;   //!
   TBranch        *b_eN;   //!
   TBranch        *b_muMax;   //!
   TBranch        *b_muSum;   //!
   TBranch        *b_muN;   //!
   TBranch        *b_fHPD;   //!
   TBranch        *b_fRBX;   //!
   TBranch        *b_n90;   //!
   TBranch        *b_fSubDet1;   //!
   TBranch        *b_fSubDet2;   //!
   TBranch        *b_fSubDet3;   //!
   TBranch        *b_fSubDet4;   //!
   TBranch        *b_restrictedEMF;   //!
   TBranch        *b_nHCAL;   //!
   TBranch        *b_nECAL;   //!
   TBranch        *b_apprHPD;   //!
   TBranch        *b_apprRBX;   //!
   TBranch        *b_n2RPC;   //!
   TBranch        *b_n3RPC;   //!
   TBranch        *b_nRPC;   //!
   TBranch        *b_fEB;   //!
   TBranch        *b_fEE;   //!
   TBranch        *b_fHB;   //!
   TBranch        *b_fHE;   //!
   TBranch        *b_fHO;   //!
   TBranch        *b_fLong;   //!
   TBranch        *b_fShort;   //!
   TBranch        *b_fLS;   //!
   TBranch        *b_fHFOOT;   //!
   TBranch        *b_matchedPt;   //!
   TBranch        *b_matchedRawPt;   //!
   TBranch        *b_matchedR;   //!
   TBranch        *b_refpt;   //!
   TBranch        *b_refeta;   //!
   TBranch        *b_refy;   //!
   TBranch        *b_refphi;   //!
   TBranch        *b_refdphijt;   //!
   TBranch        *b_refdrjt;   //!
   TBranch        *b_refparton_pt;   //!
   TBranch        *b_refparton_flavor;   //!
   TBranch        *b_refparton_flavorForB;   //!
   TBranch        *b_genChargedSum;   //!
   TBranch        *b_genHardSum;   //!
   TBranch        *b_signalChargedSum;   //!
   TBranch        *b_signalHardSum;   //!
   TBranch        *b_subid;   //!
   TBranch        *b_ngen;   //!
   TBranch        *b_genmatchindex;   //!
   TBranch        *b_genpt;   //!
   TBranch        *b_geneta;   //!
   TBranch        *b_geny;   //!
   TBranch        *b_genphi;   //!
   TBranch        *b_gendphijt;   //!
   TBranch        *b_gendrjt;   //!
   TBranch        *b_gensubid;   //!

//   ak5PFJetAnalyzert(const vector<string>& infiles, const string& outfile, 
//										 TTree *tree=0);
   ak5PFJetAnalyzert(TTree *tree=0);
   virtual ~ak5PFJetAnalyzert();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
virtual void     Loop();
   TFile* f_out;
   string outfile;
   TTree* hlt_tree;
   hltanalysisHltTree* hltAnalyzer;
   JetTree8TeV* tree8tev;
   bool runningOn8Tev;
};
#endif

#ifdef ak5PFJetAnalyzert_cxx
//ak5PFJetAnalyzert::ak5PFJetAnalyzert(const vector<string>& infiles, const string& ofile,
//																		TTree *tree) 
ak5PFJetAnalyzert::ak5PFJetAnalyzert(TTree *tree) 
: fChain(0)
{
// TChain * chain = new TChain("ak5PFJetAnalyzer/t","");
//  ///for hlt 
//  TChain * hltchain = new TChain("hltanalysis/HltTree","");
//  ///for tree8tev 
//  TChain * tree8chain = new TChain("jetanalyzer/JetTree","");
//  //
//	const unsigned numinput(infiles.size());
//	for(unsigned i(0);i<numinput;++i)
//	{
//		string infile(infiles[i]);
//		if(infile.find("pnfs")!=string::npos)
//			infile = "dcache:"+infile;
//		cout<<"Running on HI \n";
//		cout<<"Files Added "<<chain->Add(infile.c_str())<<endl;
//		cout<<"Files Added HLT "<<hltchain->Add(infile.c_str())<<endl;
//		cout<<"Files Addad 8 TeV "<<tree8chain->Add(infile.c_str())<<endl;
//		cout<<infile<<endl;
//	}
//
//	tree = chain;
//  Init(tree);
//  //hlt
//  hltAnalyzer = new hltanalysisHltTree(infiles);
//  hltAnalyzer->Init(hltchain);
//	/// for 8 tev 
//	tree8tev = new JetTree8TeV();
//	tree8tev->Init(tree8chain);
//	///determine which dataset 
//	const unsigned ak5entries(fChain->GetEntries()),
//				hltentries(hltAnalyzer->fChain->GetEntries()),
//				treeentries(tree8tev->fChain->GetEntries());
////	cout<<"Running on Entries "<<ak5entries<<endl;
////	cout<<"Hlt tree "<<hltentries<<endl;
////	cout<<"8 TeV Tree "<<treeentries<<endl;
//	if(ak5entries && hltentries)
//	{
//		assert(!treeentries);
//		runningOn8Tev=false;
//	}
//	else if (treeentries)
//	{
//		assert(!ak5entries && !hltentries);
//		runningOn8Tev=true;
//	}
}

ak5PFJetAnalyzert::~ak5PFJetAnalyzert()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

Int_t ak5PFJetAnalyzert::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
Long64_t ak5PFJetAnalyzert::LoadTree(Long64_t entry)
{
	// Set the environment to read one entry
	if (!fChain) return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if (centry < 0) return centry;
	if (fChain->GetTreeNumber() != fCurrent) {
		fCurrent = fChain->GetTreeNumber();
	}
	return centry;
}

void ak5PFJetAnalyzert::Init(TTree *tree)
{
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).

	// Set branch addresses and branch pointers
	if (!tree) return;
	fChain = tree;
	fCurrent = -1;
	fChain->SetMakeClass(1);

	fChain->SetBranchAddress("evt", &evt, &b_evt);
	fChain->SetBranchAddress("b", &b, &b_b);
	fChain->SetBranchAddress("nref", &nref, &b_nref);
	fChain->SetBranchAddress("rawpt", rawpt, &b_rawpt);
	fChain->SetBranchAddress("jtpt", jtpt, &b_jtpt);
	fChain->SetBranchAddress("jteta", jteta, &b_jteta);
	fChain->SetBranchAddress("jty", jty, &b_jty);
	fChain->SetBranchAddress("jtphi", jtphi, &b_jtphi);
	fChain->SetBranchAddress("jtpu", jtpu, &b_jtpu);
	fChain->SetBranchAddress("jtm", jtm, &b_jtm);
	fChain->SetBranchAddress("discr_fr01", discr_fr01, &b_discr_fr01);
	fChain->SetBranchAddress("trackMax", trackMax, &b_trackMax);
	fChain->SetBranchAddress("trackSum", trackSum, &b_trackSum);
	fChain->SetBranchAddress("trackN", trackN, &b_trackN);
	fChain->SetBranchAddress("trackHardSum", trackHardSum, &b_trackHardSum);
	fChain->SetBranchAddress("trackHardN", trackHardN, &b_trackHardN);
	fChain->SetBranchAddress("chargedMax", chargedMax, &b_chargedMax);
	fChain->SetBranchAddress("chargedSum", chargedSum, &b_chargedSum);
	fChain->SetBranchAddress("chargedN", chargedN, &b_chargedN);
	fChain->SetBranchAddress("chargedHardSum", chargedHardSum, &b_chargedHardSum);
	fChain->SetBranchAddress("chargedHardN", chargedHardN, &b_chargedHardN);
	fChain->SetBranchAddress("photonMax", photonMax, &b_photonMax);
	fChain->SetBranchAddress("photonSum", photonSum, &b_photonSum);
	fChain->SetBranchAddress("photonN", photonN, &b_photonN);
	fChain->SetBranchAddress("photonHardSum", photonHardSum, &b_photonHardSum);
	fChain->SetBranchAddress("photonHardN", photonHardN, &b_photonHardN);
	fChain->SetBranchAddress("neutralMax", neutralMax, &b_neutralMax);
	fChain->SetBranchAddress("neutralSum", neutralSum, &b_neutralSum);
	fChain->SetBranchAddress("neutralN", neutralN, &b_neutralN);
	fChain->SetBranchAddress("hcalSum", hcalSum, &b_hcalSum);
	fChain->SetBranchAddress("ecalSum", ecalSum, &b_ecalSum);
	fChain->SetBranchAddress("eMax", eMax, &b_eMax);
	fChain->SetBranchAddress("eSum", eSum, &b_eSum);
	fChain->SetBranchAddress("eN", eN, &b_eN);
	fChain->SetBranchAddress("muMax", muMax, &b_muMax);
	fChain->SetBranchAddress("muSum", muSum, &b_muSum);
	fChain->SetBranchAddress("muN", muN, &b_muN);
	fChain->SetBranchAddress("fHPD", fHPD, &b_fHPD);
	fChain->SetBranchAddress("fRBX", fRBX, &b_fRBX);
	fChain->SetBranchAddress("n90", n90, &b_n90);
	fChain->SetBranchAddress("fSubDet1", fSubDet1, &b_fSubDet1);
	fChain->SetBranchAddress("fSubDet2", fSubDet2, &b_fSubDet2);
	fChain->SetBranchAddress("fSubDet3", fSubDet3, &b_fSubDet3);
	fChain->SetBranchAddress("fSubDet4", fSubDet4, &b_fSubDet4);
	fChain->SetBranchAddress("restrictedEMF", restrictedEMF, &b_restrictedEMF);
	fChain->SetBranchAddress("nHCAL", nHCAL, &b_nHCAL);
	fChain->SetBranchAddress("nECAL", nECAL, &b_nECAL);
	fChain->SetBranchAddress("apprHPD", apprHPD, &b_apprHPD);
	fChain->SetBranchAddress("apprRBX", apprRBX, &b_apprRBX);
	fChain->SetBranchAddress("n2RPC", n2RPC, &b_n2RPC);
	fChain->SetBranchAddress("n3RPC", n3RPC, &b_n3RPC);
	fChain->SetBranchAddress("nRPC", nRPC, &b_nRPC);
	fChain->SetBranchAddress("fEB", fEB, &b_fEB);
	fChain->SetBranchAddress("fEE", fEE, &b_fEE);
	fChain->SetBranchAddress("fHB", fHB, &b_fHB);
	fChain->SetBranchAddress("fHE", fHE, &b_fHE);
	fChain->SetBranchAddress("fHO", fHO, &b_fHO);
	fChain->SetBranchAddress("fLong", fLong, &b_fLong);
	fChain->SetBranchAddress("fShort", fShort, &b_fShort);
	fChain->SetBranchAddress("fLS", fLS, &b_fLS);
	fChain->SetBranchAddress("fHFOOT", fHFOOT, &b_fHFOOT);
	fChain->SetBranchAddress("matchedPt", matchedPt, &b_matchedPt);
	fChain->SetBranchAddress("matchedRawPt", matchedRawPt, &b_matchedRawPt);
	fChain->SetBranchAddress("matchedR", matchedR, &b_matchedR);
	fChain->SetBranchAddress("refpt", refpt, &b_refpt);
	fChain->SetBranchAddress("refeta", refeta, &b_refeta);
	fChain->SetBranchAddress("refy", refy, &b_refy);
	fChain->SetBranchAddress("refphi", refphi, &b_refphi);
	fChain->SetBranchAddress("refdphijt", refdphijt, &b_refdphijt);
	fChain->SetBranchAddress("refdrjt", refdrjt, &b_refdrjt);
	fChain->SetBranchAddress("refparton_pt", refparton_pt, &b_refparton_pt);
	fChain->SetBranchAddress("refparton_flavor", refparton_flavor, &b_refparton_flavor);
	fChain->SetBranchAddress("genChargedSum", genChargedSum, &b_genChargedSum);
	fChain->SetBranchAddress("genHardSum", genHardSum, &b_genHardSum);
	fChain->SetBranchAddress("ngen", &ngen, &b_ngen);
	fChain->SetBranchAddress("genmatchindex", genmatchindex, &b_genmatchindex);
	fChain->SetBranchAddress("genpt", genpt, &b_genpt);
	fChain->SetBranchAddress("geneta", geneta, &b_geneta);
	fChain->SetBranchAddress("geny", geny, &b_geny);
	fChain->SetBranchAddress("genphi", genphi, &b_genphi);

}


#endif // #ifdef ak5PFJetAnalyzert_cxx
