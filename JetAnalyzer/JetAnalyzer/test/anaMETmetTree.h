//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan 28 11:10:37 2014 by ROOT version 5.32/00
// from TTree metTree/v1
// found on file: /uscms/home/tlibeiro/3days/pp_minbiasSkim_forest_53x_2013-08-15-0155.root
//////////////////////////////////////////////////////////

#ifndef anaMETmetTree_h
#define anaMETmetTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class anaMETmetTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           nEv;
   Int_t           nLumi;
   Int_t           nBX;
   Int_t           nRun;
   Int_t           nMET;
   Float_t         MET[1];   //[nMET]
   Float_t         METPhi[1];   //[nMET]
   Float_t         SumEt[1];   //[nMET]

   // List of branches
   TBranch        *b_nEv;   //!
   TBranch        *b_nLumi;   //!
   TBranch        *b_nBX;   //!
   TBranch        *b_nRun;   //!
   TBranch        *b_nMET;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_METPhi;   //!
   TBranch        *b_SumEt;   //!

   anaMETmetTree(TTree *tree=0);
   virtual ~anaMETmetTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef anaMETmetTree_cxx
anaMETmetTree::anaMETmetTree(TTree *tree) : fChain(0) 
{
//// if parameter tree is not specified (or zero), connect the file
//// used to generate this class and read the Tree.
//   if (tree == 0) {
//      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/uscms/home/tlibeiro/3days/pp_minbiasSkim_forest_53x_2013-08-15-0155.root");
//      if (!f || !f->IsOpen()) {
//         f = new TFile("/uscms/home/tlibeiro/3days/pp_minbiasSkim_forest_53x_2013-08-15-0155.root");
//      }
//      TDirectory * dir = (TDirectory*)f->Get("/uscms/home/tlibeiro/3days/pp_minbiasSkim_forest_53x_2013-08-15-0155.root:/anaMET");
//      dir->GetObject("metTree",tree);
//
//   }
//   Init(tree);
}

anaMETmetTree::~anaMETmetTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t anaMETmetTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t anaMETmetTree::LoadTree(Long64_t entry)
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

void anaMETmetTree::Init(TTree *tree)
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

   fChain->SetBranchAddress("nEv", &nEv, &b_nEv);
   fChain->SetBranchAddress("nLumi", &nLumi, &b_nLumi);
   fChain->SetBranchAddress("nBX", &nBX, &b_nBX);
   fChain->SetBranchAddress("nRun", &nRun, &b_nRun);
   fChain->SetBranchAddress("nMET", &nMET, &b_nMET);
   fChain->SetBranchAddress("MET", MET, &b_MET);
   fChain->SetBranchAddress("METPhi", METPhi, &b_METPhi);
   fChain->SetBranchAddress("SumEt", SumEt, &b_SumEt);
   Notify();
}

Bool_t anaMETmetTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void anaMETmetTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t anaMETmetTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef anaMETmetTree_cxx
