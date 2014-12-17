#ifndef JetTree8TeV_h
#define JetTree8TeV_h
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class JetTree8TeV {
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
   Float_t         jtpt[30];
   Float_t         jteta[30];
   Float_t         jtphi[30];
   Float_t         jtm[30];
   Float_t         jty[30];
   Float_t         rawpt[30];
   Float_t         chargedHardSum[30];
   Float_t         chargedHardN[30];
   Float_t         neutralHardSum[30];
   Float_t         neutralHardN[30];
   Float_t         emEnergy[30];
   Float_t         muSum[30];
   Float_t         eSum[30];
   Float_t         photonSum[30];
   Int_t           L1_SingleJet16_Prescl;
   Int_t          L1_SingleJet16;
   Int_t           L1_SingleJet36_Prescl;
   Int_t          L1_SingleJet36;
   Int_t           L1_SingleJet68_Prescl;
   Int_t          L1_SingleJet68;
   Int_t           L1_SingleJet92_Prescl;
   Int_t          L1_SingleJet92;
   Int_t           L1_SingleJet128_Prescl;
   Int_t          L1_SingleJet128;
   Int_t           HLT_PFJet40_Prescl;
   Int_t          HLT_PFJet40;
   Int_t           HLT_PFJet80_Prescl;
   Int_t          HLT_PFJet80;
   Int_t           HLT_PFJet140_Prescl;
   Int_t          HLT_PFJet140;
   Int_t           HLT_PFJet200_Prescl;
   Int_t          HLT_PFJet200;
   Int_t           HLT_PFJet260_Prescl;
   Int_t          HLT_PFJet260;
   Int_t           HLT_PFJet320_Prescl;
   Int_t          HLT_PFJet320;
   Int_t           HLT_PFJet400_Prescl;
   Int_t          HLT_PFJet400;
   Int_t           HLT_PFJet40L1_Prescl;
   Int_t           HLT_PFJet80L1_Prescl;
   Int_t           HLT_PFJet140L1_Prescl;
   Int_t           HLT_PFJet200L1_Prescl;
   Int_t           HLT_PFJet260L1_Prescl;
   Int_t           HLT_PFJet320L1_Prescl;
   Int_t           HLT_PFJet400L1_Prescl;
   Float_t          HLT_PFJet40_ObjPt;
   Float_t          HLT_PFJet80_ObjPt;
   Float_t          HLT_PFJet140_ObjPt;
   Float_t          HLT_PFJet200_ObjPt;
   Float_t          HLT_PFJet260_ObjPt;
   Float_t          HLT_PFJet320_ObjPt;
   Float_t          HLT_PFJet400_ObjPt;

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
   TBranch        *b_rawpt;   //!
   TBranch        *b_chargedHardSum;   //!
   TBranch        *b_chargedHardN;   //!
   TBranch        *b_neutralHardSum;   //!
   TBranch        *b_neutralHardN;   //!
   TBranch        *b_emEnergy;   //!
   TBranch        *b_muSum;   //!
   TBranch        *b_eSum;   //!
   TBranch        *b_photonSum;   //!
   TBranch        *b_L1_SingleJet16_Prescl;   //!
   TBranch        *b_L1_SingleJet16;   //!
   TBranch        *b_L1_SingleJet36_Prescl;   //!
   TBranch        *b_L1_SingleJet36;   //!
   TBranch        *b_L1_SingleJet68_Prescl;   //!
   TBranch        *b_L1_SingleJet68;   //!
   TBranch        *b_L1_SingleJet92_Prescl;   //!
   TBranch        *b_L1_SingleJet92;   //!
   TBranch        *b_L1_SingleJet128_Prescl;   //!
   TBranch        *b_L1_SingleJet128;   //!
   TBranch        *b_HLT_PFJet40_Prescl;   //!
   TBranch        *b_HLT_PFJet40;   //!
   TBranch        *b_HLT_PFJet80_Prescl;   //!
   TBranch        *b_HLT_PFJet80;   //!
   TBranch        *b_HLT_PFJet140_Prescl;   //!
   TBranch        *b_HLT_PFJet140;   //!
   TBranch        *b_HLT_PFJet200_Prescl;   //!
   TBranch        *b_HLT_PFJet200;   //!
   TBranch        *b_HLT_PFJet260_Prescl;   //!
   TBranch        *b_HLT_PFJet260;   //!
   TBranch        *b_HLT_PFJet320_Prescl;   //!
   TBranch        *b_HLT_PFJet320;   //!
   TBranch        *b_HLT_PFJet400_Prescl;   //!
   TBranch        *b_HLT_PFJet400;   //!
   TBranch        *b_HLT_PFJet40L1_Prescl;   //!
   TBranch        *b_HLT_PFJet80L1_Prescl;   //!
   TBranch        *b_HLT_PFJet140L1_Prescl;   //!
   TBranch        *b_HLT_PFJet200L1_Prescl;   //!
   TBranch        *b_HLT_PFJet260L1_Prescl;   //!
   TBranch        *b_HLT_PFJet320L1_Prescl;   //!
   TBranch        *b_HLT_PFJet400L1_Prescl;   //!
   TBranch        *b_HLT_PFJet40_ObjPt;
   TBranch        *b_HLT_PFJet80_ObjPt;
   TBranch        *b_HLT_PFJet140_ObjPt;
   TBranch        *b_HLT_PFJet200_ObjPt;
   TBranch        *b_HLT_PFJet260_ObjPt;
   TBranch        *b_HLT_PFJet320_ObjPt;
   TBranch        *b_HLT_PFJet400_ObjPt;


   JetTree8TeV(TTree *tree=0);
   virtual ~JetTree8TeV();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
};

#endif

#ifdef JetTree8TeV_cxx
JetTree8TeV::JetTree8TeV(TTree *tree) : fChain(0) 
{
//   if (tree == 0) {
//      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ntuples_incxsec.root");
//      if (!f || !f->IsOpen()) {
//         f = new TFile("ntuples_incxsec.root");
//      }
//      TDirectory * dir = (TDirectory*)f->Get("ntuples_incxsec.root:/jetanalyzer");
//      dir->GetObject("JetTree",tree);
//
//   }
//   Init(tree);
}

JetTree8TeV::~JetTree8TeV()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t JetTree8TeV::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t JetTree8TeV::LoadTree(Long64_t entry)
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

void JetTree8TeV::Init(TTree *tree)
{
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
   fChain->SetBranchAddress("rawpt", rawpt, &b_rawpt);
   fChain->SetBranchAddress("chargedHardSum", chargedHardSum, &b_chargedHardSum);
   fChain->SetBranchAddress("chargedHardN", chargedHardN, &b_chargedHardN);
   fChain->SetBranchAddress("neutralHardSum", neutralHardSum, &b_neutralHardSum);
   fChain->SetBranchAddress("neutralHardN", neutralHardN, &b_neutralHardN);
   fChain->SetBranchAddress("emEnergy", emEnergy, &b_emEnergy);
   fChain->SetBranchAddress("muSum", muSum, &b_muSum);
   fChain->SetBranchAddress("eSum", eSum, &b_eSum);
   fChain->SetBranchAddress("photonSum", photonSum, &b_photonSum);
   fChain->SetBranchAddress("L1_SingleJet16_Prescl", &L1_SingleJet16_Prescl, &b_L1_SingleJet16_Prescl);
   fChain->SetBranchAddress("L1_SingleJet16", &L1_SingleJet16, &b_L1_SingleJet16);
   fChain->SetBranchAddress("L1_SingleJet36_Prescl", &L1_SingleJet36_Prescl, &b_L1_SingleJet36_Prescl);
   fChain->SetBranchAddress("L1_SingleJet36", &L1_SingleJet36, &b_L1_SingleJet36);
   fChain->SetBranchAddress("L1_SingleJet68_Prescl", &L1_SingleJet68_Prescl, &b_L1_SingleJet68_Prescl);
   fChain->SetBranchAddress("L1_SingleJet68", &L1_SingleJet68, &b_L1_SingleJet68);
   fChain->SetBranchAddress("L1_SingleJet92_Prescl", &L1_SingleJet92_Prescl, &b_L1_SingleJet92_Prescl);
   fChain->SetBranchAddress("L1_SingleJet92", &L1_SingleJet92, &b_L1_SingleJet92);
   fChain->SetBranchAddress("L1_SingleJet128_Prescl", &L1_SingleJet128_Prescl, &b_L1_SingleJet128_Prescl);
   fChain->SetBranchAddress("L1_SingleJet128", &L1_SingleJet128, &b_L1_SingleJet128);
   fChain->SetBranchAddress("HLT_PFJet40_Prescl", &HLT_PFJet40_Prescl, &b_HLT_PFJet40_Prescl);
   fChain->SetBranchAddress("HLT_PFJet40", &HLT_PFJet40, &b_HLT_PFJet40);
   fChain->SetBranchAddress("HLT_PFJet80_Prescl", &HLT_PFJet80_Prescl, &b_HLT_PFJet80_Prescl);
   fChain->SetBranchAddress("HLT_PFJet80", &HLT_PFJet80, &b_HLT_PFJet80);
   fChain->SetBranchAddress("HLT_PFJet140_Prescl", &HLT_PFJet140_Prescl, &b_HLT_PFJet140_Prescl);
   fChain->SetBranchAddress("HLT_PFJet140", &HLT_PFJet140, &b_HLT_PFJet140);
   fChain->SetBranchAddress("HLT_PFJet200_Prescl", &HLT_PFJet200_Prescl, &b_HLT_PFJet200_Prescl);
   fChain->SetBranchAddress("HLT_PFJet200", &HLT_PFJet200, &b_HLT_PFJet200);
   fChain->SetBranchAddress("HLT_PFJet260_Prescl", &HLT_PFJet260_Prescl, &b_HLT_PFJet260_Prescl);
   fChain->SetBranchAddress("HLT_PFJet260", &HLT_PFJet260, &b_HLT_PFJet260);
   fChain->SetBranchAddress("HLT_PFJet320_Prescl", &HLT_PFJet320_Prescl, &b_HLT_PFJet320_Prescl);
   fChain->SetBranchAddress("HLT_PFJet320", &HLT_PFJet320, &b_HLT_PFJet320);
   fChain->SetBranchAddress("HLT_PFJet400_Prescl", &HLT_PFJet400_Prescl, &b_HLT_PFJet400_Prescl);
   fChain->SetBranchAddress("HLT_PFJet400", &HLT_PFJet400, &b_HLT_PFJet400);
   fChain->SetBranchAddress("HLT_PFJet40L1_Prescl", &HLT_PFJet40L1_Prescl, &b_HLT_PFJet40L1_Prescl);
   fChain->SetBranchAddress("HLT_PFJet80L1_Prescl", &HLT_PFJet80L1_Prescl, &b_HLT_PFJet80L1_Prescl);
   fChain->SetBranchAddress("HLT_PFJet140L1_Prescl", &HLT_PFJet140L1_Prescl, &b_HLT_PFJet140L1_Prescl);
   fChain->SetBranchAddress("HLT_PFJet200L1_Prescl", &HLT_PFJet200L1_Prescl, &b_HLT_PFJet200L1_Prescl);
   fChain->SetBranchAddress("HLT_PFJet260L1_Prescl", &HLT_PFJet260L1_Prescl, &b_HLT_PFJet260L1_Prescl);
   fChain->SetBranchAddress("HLT_PFJet320L1_Prescl", &HLT_PFJet320L1_Prescl, &b_HLT_PFJet320L1_Prescl);
   fChain->SetBranchAddress("HLT_PFJet400L1_Prescl", &HLT_PFJet400L1_Prescl, &b_HLT_PFJet400L1_Prescl);
   fChain->SetBranchAddress("HLT_PFJet40_ObjPt", &HLT_PFJet40_ObjPt, &b_HLT_PFJet40_ObjPt);
   fChain->SetBranchAddress("HLT_PFJet80_ObjPt", &HLT_PFJet80_ObjPt, &b_HLT_PFJet80_ObjPt);
   fChain->SetBranchAddress("HLT_PFJet140_ObjPt", &HLT_PFJet140_ObjPt, &b_HLT_PFJet140_ObjPt);
   fChain->SetBranchAddress("HLT_PFJet200_ObjPt", &HLT_PFJet200_ObjPt, &b_HLT_PFJet200_ObjPt);
   fChain->SetBranchAddress("HLT_PFJet260_ObjPt", &HLT_PFJet260_ObjPt, &b_HLT_PFJet260_ObjPt);
   fChain->SetBranchAddress("HLT_PFJet320_ObjPt", &HLT_PFJet320_ObjPt, &b_HLT_PFJet320_ObjPt);
   fChain->SetBranchAddress("HLT_PFJet400_ObjPt", &HLT_PFJet400_ObjPt, &b_HLT_PFJet400_ObjPt);

}


#endif // #ifdef JetTree8TeV_cxx
