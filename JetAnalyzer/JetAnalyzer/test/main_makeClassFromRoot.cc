#include "TROOT.h"
#include <string>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
using namespace std;
using std::scientific;

int main(int argc, char **argv)
//int main_makeClassFromRoot()
{

	TROOT root("astring","bstring");
	root.SetBatch(kTRUE);
	string indirDCache="/pnfs/cms/WAX/11/store/user/tlibeiro/PatTuples"; 
	string indir = "/uscms/home/tlibeiro/3days";
//	string outdir= "TreeClasses";
//  const string infilename(indir+"/pp_minbiasSkim_forest_53x_2013-08-15-0155.root");
  const string infilename(indir+"/pt15_pp2013_P01_prod22_v81_merged_forest_0.root");
  
  TFile*  infile = new TFile(infilename.c_str(),"READ");
	TList* keys = infile->GetListOfKeys();
  const unsigned numkeys(keys->GetEntries());
  cout<<"Num Keys found: "<<numkeys<<endl;
//  cout<<"Writing  classes to dir: "<<outdir<<endl;
  for(unsigned i(0);i<numkeys;++i)
		{
			TDirectory* dir =(TDirectoryFile*)keys->At(i);
			gDirectory->cd(dir->GetName());
      const string dirname(dir->GetName());
      TList* list = (TList*)gDirectory->GetListOfKeys();
      unsigned numtrees = list->GetSize();
      cout<<"Number of Trees "<<numtrees<<endl;
			for(unsigned j(0); j<numtrees; ++j) {
      const unsigned numtrees =list->GetSize() ;
      const string treename =list->At(j)->GetName() ;
      cout<<"Making Class for Tree: "<<dirname+"/"+treename<<endl;
      TTree* tree2 = (TTree*)gFile->Get((dirname+"/"+treename).c_str());
      tree2->MakeClass((dirname+treename).c_str());
		  }
      gROOT->cd("/uscms/home/tlibeiro/3days/pt15_pp2013_P01_prod22_v81_merged_forest_0.root:/");
	  }
  cout<<"Done\n";
}

