#include "TROOT.h"
#include <string>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "runJetAnalyzer.h"
using namespace std;
using std::scientific;

int main(int argc, char **argv)
{
	string outdir= ".";
	string indirDCache="/pnfs/cms/WAX/11/store/user/tlibeiro//HI_fromYue"; 
	string indirEOS="/eos/uscms/store/user/tlibeiro/"; 
//	string indirEOS="/uscmst1b_scratch/lpc1/3DayLifetime/tlibeiro/"; 
	string infile1(indirEOS+"/HI_fromYue/PP2013_HiForest_PromptReco_JSon_Jet40Jet60_ppTrack_forestv84.root"),outfile(outdir+"/AnalyzerHistosData.root") ;
  string infile2(indirEOS+"/HI_fromYue/PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_forestv82.root");
  string infileMB(indirEOS+"/HI_fromYue/pp_minbiasSkim_forest_53x_2013-08-15-0155.root"),outMB(outdir+"/AnalyzerHistosMinBias.root");
//-------MC Datasets
//	string inpt15 (indirEOS+"/HI_fromYue/pt15_pp2013_P01_prod22_v81_merged_forest_0.root"),  outpt15 (outdir+"/AnalyzerHistos_pt15.root") ;
//	string inpt30 (indirEOS+"/HI_fromYue/pt30_pp2013_P01_prod22_v81_merged_forest_0.root"),  outpt30 (outdir+"/AnalyzerHistos_pt30.root") ;
//	string inpt50 (indirEOS+"/HI_fromYue/pt50_pp2013_P01_prod22_v81_merged_forest_0.root"),  outpt50 (outdir+"/AnalyzerHistos_pt50.root") ;
//	string inpt80 (indirEOS+"/HI_fromYue/pt80_pp2013_P01_prod22_v81_merged_forest_0.root"),  outpt80 (outdir+"/AnalyzerHistos_pt80.root") ;
//	string inpt120(indirEOS+"/HI_fromYue/pt120_pp2013_P01_prod22_v81_merged_forest_0.root"), outpt120(outdir+"/AnalyzerHistos_pt120.root") ;
//	string inpt170(indirEOS+"/HI_fromYue/pt170_pp2013_P01_prod22_v81_merged_forest_0.root"), outpt170(outdir+"/AnalyzerHistos_pt170.root") ;
//	string inpt220(indirEOS+"/HI_fromYue/pt220_pp2013_P01_prod22_v81_merged_forest_0.root"), outpt220(outdir+"/AnalyzerHistos_pt220.root") ;
//	string inpt280(indirEOS+"/HI_fromYue/pt280_pp2013_P01_prod22_v81_merged_forest_0.root"), outpt280(outdir+"/AnalyzerHistos_pt280.root") ;
	string inpt15 (indirEOS+"/QCDPP276_Winter14_V5/ntuples_incxsec_pt15.root" ),  outpt15 (outdir+"/AnalyzerIncHistos_pt15.root") ;
	string inpt30 (indirEOS+"/QCDPP276_Winter14_V5/ntuples_incxsec_pt30.root" ),  outpt30 (outdir+"/AnalyzerIncHistos_pt30.root") ;
	string inpt50 (indirEOS+"/QCDPP276_Winter14_V5/ntuples_incxsec_pt50.root" ),  outpt50 (outdir+"/AnalyzerIncHistos_pt50.root") ;
	string inpt80 (indirEOS+"/QCDPP276_Winter14_V5/ntuples_incxsec_pt80.root" ),  outpt80 (outdir+"/AnalyzerIncHistos_pt80.root") ;
	string inpt120(indirEOS+"/QCDPP276_Winter14_V5/ntuples_incxsec_pt120.root"), outpt120(outdir+"/AnalyzerIncHistos_pt120.root") ;
	string inpt170(indirEOS+"/QCDPP276_Winter14_V5/ntuples_incxsec_pt170.root"), outpt170(outdir+"/AnalyzerIncHistos_pt170.root") ;
	string inpt220(indirEOS+"/QCDPP276_Winter14_V5/ntuples_incxsec_pt220.root"), outpt220(outdir+"/AnalyzerIncHistos_pt220.root") ;
	string inpt280(indirEOS+"/QCDPP276_Winter14_V5/ntuples_incxsec_pt280.root"), outpt280(outdir+"/AnalyzerIncHistos_pt280.root") ;
	string inpt370(indirEOS+"/QCDPP276_Winter14_V5/ntuples_incxsec_pt370.root"), outpt370(outdir+"/AnalyzerIncHistos_pt370.root") ;
	string inpt460(indirEOS+"/QCDPP276_Winter14_V5/ntuples_incxsec_pt460.root"), outpt460(outdir+"/AnalyzerIncHistos_pt460.root") ;
	string inpt540(indirEOS+"/QCDPP276_Winter14_V5/ntuples_incxsec_pt540.root"), outpt540(outdir+"/AnalyzerIncHistos_pt540.root") ;
////8 TeV
	string infile8t(indirEOS+"/JetHT2012C/ntuples_incxsec_*.root"),outfile8t(outdir+"/AnalyzerHistosData8TeV.root") ;
////2.76 TeV inclusive x sec ntuples
	string infile276t(indirEOS+"/Run2013APP276_8tevCorr_Winter14_V5/ntuples_incxsec.root"),outfile276t(outdir+"/AnalyzerIncHistosData276TeV_winter14.root") ;

	vector<string> infiles8t;
	infiles8t.push_back(infile8t);
///276 tev
  vector<string> infiles276t;
	infiles276t.push_back(infile276t);

	//input  files
	vector<string> infiles;
	infiles.push_back(infile1);
	
  vector<string> inputfilesAllMC;
  inputfilesAllMC.push_back(inpt15);
  inputfilesAllMC.push_back(inpt30);
  inputfilesAllMC.push_back(inpt50);
  inputfilesAllMC.push_back(inpt80);
  inputfilesAllMC.push_back(inpt120);
  inputfilesAllMC.push_back(inpt170);
  inputfilesAllMC.push_back(inpt220);
  inputfilesAllMC.push_back(inpt280);
  inputfilesAllMC.push_back(inpt370);
  inputfilesAllMC.push_back(inpt460);
  inputfilesAllMC.push_back(inpt540);

  TROOT root("astring","bstring");
	root.SetBatch(kTRUE);
	//cout << scientific;
	cout<<setprecision(4);

	bool runOnMC=false;
  bool runak5=true;
	cout<<"Data =========================================="<<endl;
//	runJetAnalyzer datahi(infiles,outfile,runOnMC); //HI Forest ntuples
//	datahi.Loop();
	string otfiledt(outfile276t.substr(0,outfile276t.find(".root")));
	if(runak5)
		otfiledt = otfiledt+"_ak5.root";
	else
		otfiledt = outfile276t;
	runJetAnalyzer dataic(infiles276t,otfiledt,runOnMC,runak5); //IncXsec ntuples
	dataic.Loop();
	//	cout<<"Data 8TeV=========================================="<<endl;
	//	runJetAnalyzer pu5pf8t(infiles8t,outfile8t,runOnMC);
	//	pu5pf8t.Loop();


//	runOnMC=true;
//	for(unsigned i(0);i<inputfilesAllMC.size();++i) {
//		vector<string> input; input.push_back(inputfilesAllMC[i]);
//    string ak5("");
//    if(runak5) ak5="_ak5";
//		ostringstream oss; oss<<"AnalyzerIncHistos_allpt"<<i<<"_winter14"<<ak5<<".root";
//		runJetAnalyzer allmc (input,oss.str(),runOnMC,runak5);
//		allmc.Loop();
//	}
}

