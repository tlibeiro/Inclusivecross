#define jesUnc8tev_cxx
#include <math.h>
#include <iomanip>
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include <vector>
#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TF1.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>
#include "AnalyzerProcedures.h"

using namespace std;
TH1D* makeNLOSpectrum(char* infile,char* title);
void generateEvents(unsigned events, TF1* fres, TH1* hnlo,
		double lumi, TTree* tree,double unc);

int main(int argc,char** argv) {
	string outdir= ".";

	TROOT root("astring","bstring");
	root.SetBatch(kTRUE);
	cout<<setprecision(4);

	bool runOnMC=false;
	runOnMC=true;
//	int num = strtol(argv[1],NULL,10);
//	cout<<"RUnning on "<<num<<endl;
	//	for(unsigned i(num);i<num+1 && i<inputfilesAllMC.size();++i) {
	//		vector<string> input; input.push_back(inputfilesAllMC[i]);
	//		ostringstream oss; oss<<"AnalyzerIncHistos_allpt"<<i<<"_8tevcorr.root";
	//		jesAnalyzer allmc (input,oss.str(),runOnMC );
	//		allmc.makejecUncHists();
	//	}
	return 0;
}

void makejecUncHists()
{
	// Instantiate uncertainty sources Summer13_V4, or Fall12_V7
	const int nsrc = 22;
	const char* srcnames[nsrc] =
	{"Absolute", "HighPtExtra",  "SinglePionECAL", "SinglePionHCAL",
		"FlavorQCD", "Time",
		"RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
		"RelativePtBB","RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeFSR",
		"RelativeStatEC2", "RelativeStatHF",
		"FlavorZJet","FlavorPhotonJet","FlavorPureGluon","FlavorPureQuark","FlavorPureCharm","FlavorPureBottom"};

	vector<JetCorrectionUncertainty*> vsrc(nsrc);
	for (int isrc = 0; isrc < nsrc; ++isrc) {
		const char *uncsrc = srcnames[isrc];
		JetCorrectorParameters *p = new JetCorrectorParameters("Summer13_V4_DATA_UncertaintySources_AK7PF.txt", uncsrc);
		JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
		vsrc[isrc] = unc;
	}

	cout<<setprecision(6);
	cout << scientific;
	///declare hists ===========================================
	TH1D* hptUnf [nybins-1]; //pt spectrum, originally used for unfolding
	TH1D* hjesUnc [nsrc+1][nybins-1][2];//nsrc+1 (total uncertainty), 0-down,1-up
	for(unsigned i(0);i<nybins-1;++i)
	{
		sprintf(name,"hptUnf_eta%d",i);
		hptUnf[i] = new TH1D(name,name,n1x[i],&x1[i][0]);	
		hptUnf[i]->Sumw2();
		for(unsigned j(0);j<2;++j) {
			const char* ud = j==0?"dn":"up";
			for (int isrc = 0; isrc < nsrc; isrc++) {
				sprintf(name,"%s_%s_eta%d",srcnames[isrc],ud,i);
				hjesUnc[isrc][i][j] = new TH1D(name,name,n1x[i],&x1[i][0]);
				hjesUnc[isrc][i][j]->Sumw2();
			}
			sprintf(name,"TotalUnc_%s_eta%d",ud,i);
			hjesUnc[nsrc][i][j] = new TH1D(name,name,n1x[i],&x1[i][0]);
			hjesUnc[nsrc][i][j]->Sumw2();
		}
	}
	unsigned nentries(1e5);
	//  nentries=5e2;
	cout<<"Running on jet scale uncertainty plots "<<nentries<<endl;
	for(unsigned i(0);i<nybins-1;++i) ///Cycle over eta bins 
	{
		for(unsigned jentry(0);jentry<nentries;++jentry) {
			//////------Determine event weight ---------------------
			double wt = 1.0;
			// Calculate uncertainty per source 
			double sum2_up(0), sum2_dw(0);
			for (int isrc = 0; isrc < nsrc; ++isrc) {
				JetCorrectionUncertainty *unc = vsrc[isrc];
				unc->setJetPt(1);
				unc->setJetEta(1);
				double sup = unc->getUncertainty(true); // up variation
				unc->setJetPt(1);
				unc->setJetEta(1);
				double sdw = unc->getUncertainty(false); // down variation
				//take the largest variation
				double maxerr = max(fabs(sup),fabs(sdw));
				if(false)
					cout<<"jet pt:eta "<<""<<":"<<""
						<<" "<<srcnames[isrc]<<" "
						<<" "<<sup<<"/"<<sdw<<endl;
				double ptup((1.+maxerr)*1.);
				double ptdn((1.-maxerr)*1.);
				hjesUnc[isrc][i][1]->Fill(ptup,wt);
				hjesUnc[isrc][i][0]->Fill(ptdn,wt);
			} // for isrc
			///Filling Spectrum histogram-----------------------
			hptUnf[i]->Fill(1.,wt);		
		if(!(jentry%200000))
			cout<<"Processing Event "<<(float)jentry<<endl;
		}///eta loop

	}//nentry loop
	///OUTPUT==============================================================
//	string unfFile(outfile.substr(0,outfile.find(".root")));
//	TFile* fstab = new TFile((unfFile+"_jecUnc8tev.root").c_str(),"RECREATE");
//	fstab->cd();
//	//write histos 
//	for(unsigned i(0);i<nybins-1;++i) {
//		hptUnf[i]->Write();
//		hjesUnc[nsrc][i][0]->Write();
//		hjesUnc[nsrc][i][1]->Write();
//		for (int isrc = 0; isrc < nsrc; ++isrc) {
//			hjesUnc[isrc][i][0]->Write();
//			hjesUnc[isrc][i][1]->Write();
//		}
//	}
//	fstab->Write();
//	fstab->Close();
	//delete histos 
	for(unsigned i(0);i<nybins-1;++i) {
		hptUnf[i]->Delete();
		hjesUnc[nsrc][i][0]->Delete();
		hjesUnc[nsrc][i][1]->Delete();
		for (int isrc = 0; isrc < nsrc; ++isrc) {
			hjesUnc[isrc][i][0]->Delete();
			hjesUnc[isrc][i][1]->Delete();
		}
	}
	//delte JEC
	for (int isrc = 0; isrc < nsrc; ++isrc) 
		delete  vsrc[isrc];

	cout<<"Jet scale Uncertainty Plots Done\n";
}

void fillArrayFromVec(const vector<float>& vec,
		float* arr, const unsigned arrSize)
{
	const unsigned size(vec.size());
	assert(size<=arrSize);
	for(unsigned i(0);i<size&&i<arrSize;++i)
		if(i<size)
			arr[i]=vec[i];
		else
			arr[i]=-1;
}


void generateEvents(unsigned events, TF1* fres, TH1* hnlo,
		double lumi, TTree* tree,double unc)
{
  tree->Reset();
	const unsigned nbins(hnlo->GetNbinsX());
	// cout<<nbins<<' '<<events<<endl;
	for(unsigned bin(1);bin<=nbins;++bin)
	{
		double xsec  = hnlo->GetBinContent(bin)*hnlo->GetBinWidth(bin)*lumi;
		double binlo = hnlo->GetBinLowEdge(bin);
		double binhi = hnlo->GetBinWidth(bin)+binlo;
		for(unsigned i(0);i<events;++i)
		{
			double genpt = gRandom->Uniform(binlo,binhi);
			double resol = gRandom->Gaus(1.0,fres->Eval(genpt));
			double recpt = genpt*resol;
			if(unc)
				recpt = recpt*unc;
			//       cout<<genpt<<" "<<recpt<<" "<<xsec<<endl;
//			tree->Fill(genpt,recpt,xsec/events);
		}
	}
}

TH1D* makeNLOSpectrum(char* infile,char* title) {
	ifstream rstream;
	rstream.open(infile);
	vector<double> ptMinV, xsecV; double xmax(0);
	if(rstream.is_open())
		while(!rstream.eof()) {
			int nBin(0);
			double ptMin(0), ptMax(0), xsec(0), knlo(0);
			rstream>>nBin>>ptMin>>ptMax>>xsec>>knlo;
			//cout<<" "<<nBin<<" "<<ptMin<<" "<<ptMax<<" "<<xsec<<" "<<knlo<<endl;
			if(ptMin && ptMax && xsec) {
				ptMinV.push_back(ptMin);
				xsecV .push_back(xsec); xmax = ptMax;
			}
		}
	//cout<<"num bins "<<ptMinV.size()<<endl;
	const unsigned nBins = ptMinV.size();
	double x[nBins+1];
	for(unsigned i(0);i<nBins;++i)
		x[i]  = ptMinV[i];
	x[nBins]=xmax;
	TH1D* h = new TH1D(title,title,nBins,x);
	//fill histogram
	for(unsigned i(0);i<nBins;++i) {
		h->SetBinContent(i+1,xsecV[i]);
		h->SetBinError(i+1,0);
	}
	return h;
};


