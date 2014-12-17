//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr  4 11:13:53 2014 by ROOT version 5.32/00
// from TTree JetTree/JetTree
// found on file: /uscms/home/tlibeiro/eos/Run2013APP276/ntuples_incxsec_12_1_VjD.root
//////////////////////////////////////////////////////////

#ifndef JetTreeAK7_h
#define JetTreeAK7_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class JetTreeAK7 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           event_runNo;
   Int_t           event_evtNo;
   Int_t           event_lumi;
   Int_t           nPV;
   Float_t         rho;
   Float_t         met;
   Float_t         sumET;
   Int_t           nref;
   Float_t         jtpt[40];
   Float_t         jteta[40];
   Float_t         jtphi[40];
   Float_t         jtm[40];
   Float_t         jty[40];
   Float_t         jecl1[40];
   Float_t         jecl2l3[40];
   Float_t         jecres[40];
   Float_t         rawpt[40];
   Float_t         chargedHardSum[40];
   Float_t         chargedHardN[40];
   Float_t         neutralHardSum[40];
   Float_t         neutralHardN[40];
   Float_t         emEnergy[40];
   Float_t         chemEnergy[40];
   Float_t         nuemEnergy[40];
   Float_t         muSum[40];
   Float_t         eSum[40];
   Float_t         photonSum[40];
   Int_t           trackN[40];
   Int_t           constituentN[40];
   Int_t           ngen;
   Float_t         genpt[40];
   Float_t         geneta[40];
   Float_t         genphi[40];
   Float_t         genm[40];
   Float_t         matchgenpt[40];
   Float_t         matchgendr[40];
   Float_t         qscale;
   Int_t           L1_SingleJet16_BptxAND_Prescl;
   Int_t           L1_SingleJet16_BptxAND;
   Int_t           L1_SingleJet36_Prescl;
   Int_t           L1_SingleJet36;
   Int_t           HLT_PAJet20_NoJetID_v1_Prescl;
   Int_t           HLT_PAJet20_NoJetID_v1;
   Int_t           HLT_PAJet20_NoJetID_v1L1_Prescl;
   Float_t         HLT_PAJet20_NoJetID_v1_ObjPt;
   Int_t           HLT_PAJet40_NoJetID_v1_Prescl;
   Int_t           HLT_PAJet40_NoJetID_v1;
   Int_t           HLT_PAJet40_NoJetID_v1L1_Prescl;
   Float_t         HLT_PAJet40_NoJetID_v1_ObjPt;
   Int_t           HLT_PAJet60_NoJetID_v1_Prescl;
   Int_t           HLT_PAJet60_NoJetID_v1;
   Int_t           HLT_PAJet60_NoJetID_v1L1_Prescl;
   Float_t         HLT_PAJet60_NoJetID_v1_ObjPt;
   Int_t           HLT_PAJet80_NoJetID_v1_Prescl;
   Int_t           HLT_PAJet80_NoJetID_v1;
   Int_t           HLT_PAJet80_NoJetID_v1L1_Prescl;
   Float_t         HLT_PAJet80_NoJetID_v1_ObjPt;
   Int_t           HLT_PAJet100_NoJetID_v1_Prescl;
   Int_t           HLT_PAJet100_NoJetID_v1;
   Int_t           HLT_PAJet100_NoJetID_v1L1_Prescl;
   Float_t         HLT_PAJet100_NoJetID_v1_ObjPt;
   Int_t           HLT_PAJet120_NoJetID_v1_Prescl;
   Int_t           HLT_PAJet120_NoJetID_v1;
   Int_t           HLT_PAJet120_NoJetID_v1L1_Prescl;
   Float_t         HLT_PAJet120_NoJetID_v1_ObjPt;

   // List of branches
   TBranch        *b_event_runNo;   //!
   TBranch        *b_event_evtNo;   //!
   TBranch        *b_event_lumi;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_rhoVal;   //!
   TBranch        *b_met;   //!
   TBranch        *b_sumET;   //!
   TBranch        *b_nref;   //!
   TBranch        *b_jtpt;   //!
   TBranch        *b_jteta;   //!
   TBranch        *b_jtphi;   //!
   TBranch        *b_jtm;   //!
   TBranch        *b_jty;   //!
   TBranch        *b_jecl1;   //!
   TBranch        *b_jecl2l3;   //!
   TBranch        *b_jecres;   //!
   TBranch        *b_rawpt;   //!
   TBranch        *b_chargedHardSum;   //!
   TBranch        *b_chargedHardN;   //!
   TBranch        *b_neutralHardSum;   //!
   TBranch        *b_neutralHardN;   //!
   TBranch        *b_emEnergy;   //!
   TBranch        *b_chemEnergy;   //!
   TBranch        *b_nuemEnergy;   //!
   TBranch        *b_muSum;   //!
   TBranch        *b_eSum;   //!
   TBranch        *b_photonSum;   //!
   TBranch        *b_trackN;   //!
   TBranch        *b_constituentN;   //!
   TBranch        *b_ngen;   //!
   TBranch        *b_genpt;   //!
   TBranch        *b_geneta;   //!
   TBranch        *b_genphi;   //!
   TBranch        *b_genm;   //!
   TBranch        *b_matchgenpt;
   TBranch        *b_matchgendr;
   TBranch        *b_qscale;
   TBranch        *b_L1_SingleJet16_BptxAND_Prescl;   //!
   TBranch        *b_L1_SingleJet16_BptxAND;   //!
   TBranch        *b_L1_SingleJet36_Prescl;   //!
   TBranch        *b_L1_SingleJet36;   //!
   TBranch        *b_HLT_PAJet20_NoJetID_v1_Prescl;   //!
   TBranch        *b_HLT_PAJet20_NoJetID_v1;   //!
   TBranch        *b_HLT_PAJet20_NoJetID_v1L1_Prescl;   //!
   TBranch        *b_HLT_PAJet20_NoJetID_v1_ObjPt;   //!
   TBranch        *b_HLT_PAJet40_NoJetID_v1_Prescl;   //!
   TBranch        *b_HLT_PAJet40_NoJetID_v1;   //!
   TBranch        *b_HLT_PAJet40_NoJetID_v1L1_Prescl;   //!
   TBranch        *b_HLT_PAJet40_NoJetID_v1_ObjPt;   //!
   TBranch        *b_HLT_PAJet60_NoJetID_v1_Prescl;   //!
   TBranch        *b_HLT_PAJet60_NoJetID_v1;   //!
   TBranch        *b_HLT_PAJet60_NoJetID_v1L1_Prescl;   //!
   TBranch        *b_HLT_PAJet60_NoJetID_v1_ObjPt;   //!
   TBranch        *b_HLT_PAJet80_NoJetID_v1_Prescl;   //!
   TBranch        *b_HLT_PAJet80_NoJetID_v1;   //!
   TBranch        *b_HLT_PAJet80_NoJetID_v1L1_Prescl;   //!
   TBranch        *b_HLT_PAJet80_NoJetID_v1_ObjPt;   //!
   TBranch        *b_HLT_PAJet100_NoJetID_v1_Prescl;   //!
   TBranch        *b_HLT_PAJet100_NoJetID_v1;   //!
   TBranch        *b_HLT_PAJet100_NoJetID_v1L1_Prescl;   //!
   TBranch        *b_HLT_PAJet100_NoJetID_v1_ObjPt;   //!
   TBranch        *b_HLT_PAJet120_NoJetID_v1_Prescl;   //!
   TBranch        *b_HLT_PAJet120_NoJetID_v1;   //!
   TBranch        *b_HLT_PAJet120_NoJetID_v1L1_Prescl;   //!
   TBranch        *b_HLT_PAJet120_NoJetID_v1_ObjPt;   //!

   JetTreeAK7(TTree *tree=0);
   virtual ~JetTreeAK7();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef JetTreeAK7_cxx
JetTreeAK7::JetTreeAK7(TTree *tree) : fChain(0) 
{
//// if parameter tree is not specified (or zero), connect the file
//// used to generate this class and read the Tree.
//   if (tree == 0) {
//      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/uscms/home/tlibeiro/eos/Run2013APP276/ntuples_incxsec_12_1_VjD.root");
//      if (!f || !f->IsOpen()) {
//         f = new TFile("/uscms/home/tlibeiro/eos/Run2013APP276/ntuples_incxsec_12_1_VjD.root");
//      }
//      TDirectory * dir = (TDirectory*)f->Get("/uscms/home/tlibeiro/eos/Run2013APP276/ntuples_incxsec_12_1_VjD.root:/jetanalyzer");
//      dir->GetObject("JetTree",tree);
//
//   }
//   Init(tree);
}

JetTreeAK7::~JetTreeAK7()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t JetTreeAK7::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t JetTreeAK7::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void JetTreeAK7::Init(TTree *tree)
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

   fChain->SetBranchAddress("event_runNo", &event_runNo, &b_event_runNo);
   fChain->SetBranchAddress("event_evtNo", &event_evtNo, &b_event_evtNo);
   fChain->SetBranchAddress("event_lumi", &event_lumi, &b_event_lumi);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("rho", &rho, &b_rhoVal);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("sumET", &sumET, &b_sumET);
   fChain->SetBranchAddress("nref", &nref, &b_nref);
   fChain->SetBranchAddress("jtpt", jtpt, &b_jtpt);
   fChain->SetBranchAddress("jteta", jteta, &b_jteta);
   fChain->SetBranchAddress("jtphi", jtphi, &b_jtphi);
   fChain->SetBranchAddress("jtm", jtm, &b_jtm);
   fChain->SetBranchAddress("jty", jty, &b_jty);
   fChain->SetBranchAddress("jecl1",   jecl1, &b_jecl1);
   fChain->SetBranchAddress("jecl2l3", jecl2l3, &b_jecl2l3);
   fChain->SetBranchAddress("jecres",  jecres, &b_jecres);
   fChain->SetBranchAddress("rawpt", rawpt, &b_rawpt);
   fChain->SetBranchAddress("chargedHardSum", chargedHardSum, &b_chargedHardSum);
   fChain->SetBranchAddress("chargedHardN", chargedHardN, &b_chargedHardN);
   fChain->SetBranchAddress("neutralHardSum", neutralHardSum, &b_neutralHardSum);
   fChain->SetBranchAddress("neutralHardN", neutralHardN, &b_neutralHardN);
   fChain->SetBranchAddress("emEnergy", emEnergy, &b_emEnergy);
   fChain->SetBranchAddress("chemEnergy", chemEnergy, &b_chemEnergy);
   fChain->SetBranchAddress("nuemEnergy", nuemEnergy, &b_nuemEnergy);
   fChain->SetBranchAddress("muSum", muSum, &b_muSum);
   fChain->SetBranchAddress("eSum", eSum, &b_eSum);
   fChain->SetBranchAddress("photonSum", photonSum, &b_photonSum);
   fChain->SetBranchAddress("trackN", trackN, &b_trackN);
   fChain->SetBranchAddress("constituentN", constituentN, &b_constituentN);
   fChain->SetBranchAddress("ngen", &ngen, &b_ngen);
   fChain->SetBranchAddress("genpt", genpt, &b_genpt);
   fChain->SetBranchAddress("geneta", geneta, &b_geneta);
   fChain->SetBranchAddress("genphi", genphi, &b_genphi);
   fChain->SetBranchAddress("genm", genm, &b_genm);
   fChain->SetBranchAddress("matchgenpt", matchgenpt, &b_matchgenpt);
   fChain->SetBranchAddress("matchgendr", matchgendr, &b_matchgendr);
   fChain->SetBranchAddress("qscale", &qscale, &b_qscale);

   fChain->SetBranchAddress("L1_SingleJet16_BptxAND_Prescl", &L1_SingleJet16_BptxAND_Prescl, &b_L1_SingleJet16_BptxAND_Prescl);
   fChain->SetBranchAddress("L1_SingleJet16_BptxAND", &L1_SingleJet16_BptxAND, &b_L1_SingleJet16_BptxAND);
   fChain->SetBranchAddress("L1_SingleJet36_Prescl", &L1_SingleJet36_Prescl, &b_L1_SingleJet36_Prescl);
   fChain->SetBranchAddress("L1_SingleJet36", &L1_SingleJet36, &b_L1_SingleJet36);
   fChain->SetBranchAddress("HLT_PAJet20_NoJetID_v1_Prescl", &HLT_PAJet20_NoJetID_v1_Prescl, &b_HLT_PAJet20_NoJetID_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAJet20_NoJetID_v1", &HLT_PAJet20_NoJetID_v1, &b_HLT_PAJet20_NoJetID_v1);
   fChain->SetBranchAddress("HLT_PAJet20_NoJetID_v1L1_Prescl", &HLT_PAJet20_NoJetID_v1L1_Prescl, &b_HLT_PAJet20_NoJetID_v1L1_Prescl);
   fChain->SetBranchAddress("HLT_PAJet20_NoJetID_v1_ObjPt", &HLT_PAJet20_NoJetID_v1_ObjPt, &b_HLT_PAJet20_NoJetID_v1_ObjPt);
   fChain->SetBranchAddress("HLT_PAJet40_NoJetID_v1_Prescl", &HLT_PAJet40_NoJetID_v1_Prescl, &b_HLT_PAJet40_NoJetID_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAJet40_NoJetID_v1", &HLT_PAJet40_NoJetID_v1, &b_HLT_PAJet40_NoJetID_v1);
   fChain->SetBranchAddress("HLT_PAJet40_NoJetID_v1L1_Prescl", &HLT_PAJet40_NoJetID_v1L1_Prescl, &b_HLT_PAJet40_NoJetID_v1L1_Prescl);
   fChain->SetBranchAddress("HLT_PAJet40_NoJetID_v1_ObjPt", &HLT_PAJet40_NoJetID_v1_ObjPt, &b_HLT_PAJet40_NoJetID_v1_ObjPt);
   fChain->SetBranchAddress("HLT_PAJet60_NoJetID_v1_Prescl", &HLT_PAJet60_NoJetID_v1_Prescl, &b_HLT_PAJet60_NoJetID_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAJet60_NoJetID_v1", &HLT_PAJet60_NoJetID_v1, &b_HLT_PAJet60_NoJetID_v1);
   fChain->SetBranchAddress("HLT_PAJet60_NoJetID_v1L1_Prescl", &HLT_PAJet60_NoJetID_v1L1_Prescl, &b_HLT_PAJet60_NoJetID_v1L1_Prescl);
   fChain->SetBranchAddress("HLT_PAJet60_NoJetID_v1_ObjPt", &HLT_PAJet60_NoJetID_v1_ObjPt, &b_HLT_PAJet60_NoJetID_v1_ObjPt);
   fChain->SetBranchAddress("HLT_PAJet80_NoJetID_v1_Prescl", &HLT_PAJet80_NoJetID_v1_Prescl, &b_HLT_PAJet80_NoJetID_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAJet80_NoJetID_v1", &HLT_PAJet80_NoJetID_v1, &b_HLT_PAJet80_NoJetID_v1);
   fChain->SetBranchAddress("HLT_PAJet80_NoJetID_v1L1_Prescl", &HLT_PAJet80_NoJetID_v1L1_Prescl, &b_HLT_PAJet80_NoJetID_v1L1_Prescl);
   fChain->SetBranchAddress("HLT_PAJet80_NoJetID_v1_ObjPt", &HLT_PAJet80_NoJetID_v1_ObjPt, &b_HLT_PAJet80_NoJetID_v1_ObjPt);
   fChain->SetBranchAddress("HLT_PAJet100_NoJetID_v1_Prescl", &HLT_PAJet100_NoJetID_v1_Prescl, &b_HLT_PAJet100_NoJetID_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAJet100_NoJetID_v1", &HLT_PAJet100_NoJetID_v1, &b_HLT_PAJet100_NoJetID_v1);
   fChain->SetBranchAddress("HLT_PAJet100_NoJetID_v1L1_Prescl", &HLT_PAJet100_NoJetID_v1L1_Prescl, &b_HLT_PAJet100_NoJetID_v1L1_Prescl);
   fChain->SetBranchAddress("HLT_PAJet100_NoJetID_v1_ObjPt", &HLT_PAJet100_NoJetID_v1_ObjPt, &b_HLT_PAJet100_NoJetID_v1_ObjPt);
   fChain->SetBranchAddress("HLT_PAJet120_NoJetID_v1_Prescl", &HLT_PAJet120_NoJetID_v1_Prescl, &b_HLT_PAJet120_NoJetID_v1_Prescl);
   fChain->SetBranchAddress("HLT_PAJet120_NoJetID_v1", &HLT_PAJet120_NoJetID_v1, &b_HLT_PAJet120_NoJetID_v1);
   fChain->SetBranchAddress("HLT_PAJet120_NoJetID_v1L1_Prescl", &HLT_PAJet120_NoJetID_v1L1_Prescl, &b_HLT_PAJet120_NoJetID_v1L1_Prescl);
   fChain->SetBranchAddress("HLT_PAJet120_NoJetID_v1_ObjPt", &HLT_PAJet120_NoJetID_v1_ObjPt, &b_HLT_PAJet120_NoJetID_v1_ObjPt);
   Notify();
}

Bool_t JetTreeAK7::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void JetTreeAK7::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t JetTreeAK7::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef JetTreeAK7_cxx
