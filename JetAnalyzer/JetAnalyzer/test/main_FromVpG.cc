#include "VplusGTreeAnalyzer.h"
#include "TROOT.h"
#include <string>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
using namespace std;
using std::scientific;

int main(int argc, char **argv)
{
	string indir = "/uscms_data/d3/tlibeiro/VPlusGamma";
	string outdir= "/uscms_data/d3/tlibeiro/VPlusGamma/VplusGAnalyzerOutput";
	string indirDCache="/pnfs/cms/WAX/11/store/user/tlibeiro/PatTuples"; 
	string infile_sg_wa(indir+"/PatTuples_WA_merged_wa_lo/VplusGTree_*.root"), outfile_sg_wa(outdir+"/VplusGTreeAnalyzer_WA.root");
	string infile_bk_gj200(indirDCache+"/PatTuples_GJets200/VplusGTree_*.root"), outfile_bk_gj200(outdir+"/VplusGTreeAnalyzer_GJets200.root");
	string infile_bk_gj400(indirDCache+"/PatTuples_GJets400/VplusGTree_*.root"), outfile_bk_gj400(outdir+"/VplusGTreeAnalyzer_GJets400.root");
	string infile_data1(indirDCache+"/PatTuples_Run2012A/VplusGTree_*.root"), outfile_data(outdir+"/VplusGTreeAnalyzer_Data.root") ;
	string infile_data2(indirDCache+"/PatTuples_Run2012B/VplusGTree_*.root");
	string infile_data3(indirDCache+"/PatTuples_Run2012C/VplusGTree_*.root");
	string infile_data4(indirDCache+"/PatTuples_Run2012D/VplusGTree_*.root");
	//signal files
  vector<string> infiles_sg_wa;
  infiles_sg_wa.push_back(infile_sg_wa);
	//data files
	vector<string> infiles_data;
//  infiles_data.push_back(argv[1]);
//  outfile_data=outdir+'/'+argv[2];
  infiles_data.push_back(infile_data1);
	infiles_data.push_back(infile_data2);
	infiles_data.push_back(infile_data3);
	infiles_data.push_back(infile_data4);
	//bkg file 
	vector<string> infiles_bk_gj200;
	infiles_bk_gj200.push_back(infile_bk_gj200);
	vector<string> infiles_bk_gj400;
	infiles_bk_gj400.push_back(infile_bk_gj400);

	TROOT root("astring","bstring");
	root.SetBatch(kTRUE);
	cout << scientific;
	cout<<setprecision(4);
	cout<<"Signal ======== =========================================="<<endl;
	VplusGTreeAnalyzer sg_wa(infiles_sg_wa,outfile_sg_wa);
	sg_wa.Loop();
	cout<<"Gamma+Jets HT 200 =========================================="<<endl;
	VplusGTreeAnalyzer bkGJets200(infiles_bk_gj200,outfile_bk_gj200);
	bkGJets200.Loop();
	cout<<"Gamma+Jets HT 400 =========================================="<<endl;
	VplusGTreeAnalyzer bkGJets400(infiles_bk_gj400,outfile_bk_gj400);
	bkGJets400.Loop();
	cout<<"Data Combined =========================================="<<endl;
	VplusGTreeAnalyzer dataCombined(infiles_data,outfile_data);
	dataCombined.Loop();
//	//  vector<string> testvec; testvec.push_back("VplusGTree.root");
//	//  VplusGTreeAnalyzer test(testvec,"test.root");
//	//  test.Loop();
}

