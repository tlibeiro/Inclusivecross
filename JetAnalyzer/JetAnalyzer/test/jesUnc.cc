#define jesUnc_cxx
#include <math.h>
#include <iomanip>
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include <vector>
#include "TFile.h"
#include "TROOT.h"
#include "TH1D.h"
#include <string>
#include <iostream>
#include "jesUnc.h"

using namespace std;

int main(int argc,char** argv) {
	string outdir= ".";
//	string indirEOS="/eos/uscms/store/user/tlibeiro/"; 
	string indirEOS="~/3days/"; 
//-------MC Datasets
	string inpt15 (indirEOS+"/QCDPP276_8tevCorr/ntuples_incxsec_pt15.root" ),  outpt15 (outdir+"/AnalyzerIncHistos_pt15.root") ;
	string inpt30 (indirEOS+"/QCDPP276_8tevCorr/ntuples_incxsec_pt30.root" ),  outpt30 (outdir+"/AnalyzerIncHistos_pt30.root") ;
	string inpt50 (indirEOS+"/QCDPP276_8tevCorr/ntuples_incxsec_pt50.root" ),  outpt50 (outdir+"/AnalyzerIncHistos_pt50.root") ;
	string inpt80 (indirEOS+"/QCDPP276_8tevCorr/ntuples_incxsec_pt80.root" ),  outpt80 (outdir+"/AnalyzerIncHistos_pt80.root") ;
	string inpt120(indirEOS+"/QCDPP276_8tevCorr/ntuples_incxsec_pt120.root"), outpt120(outdir+"/AnalyzerIncHistos_pt120.root") ;
	string inpt170(indirEOS+"/QCDPP276_8tevCorr/ntuples_incxsec_pt170.root"), outpt170(outdir+"/AnalyzerIncHistos_pt170.root") ;
	string inpt220(indirEOS+"/QCDPP276_8tevCorr/ntuples_incxsec_pt220.root"), outpt220(outdir+"/AnalyzerIncHistos_pt220.root") ;
	string inpt280(indirEOS+"/QCDPP276_8tevCorr/ntuples_incxsec_pt280.root"), outpt280(outdir+"/AnalyzerIncHistos_pt280.root") ;
	string inpt370(indirEOS+"/QCDPP276_8tevCorr/ntuples_incxsec_pt370.root"), outpt370(outdir+"/AnalyzerIncHistos_pt370.root") ;
	string inpt460(indirEOS+"/QCDPP276_8tevCorr/ntuples_incxsec_pt460.root"), outpt460(outdir+"/AnalyzerIncHistos_pt460.root") ;
	string inpt540(indirEOS+"/QCDPP276_8tevCorr/ntuples_incxsec_pt540.root"), outpt540(outdir+"/AnalyzerIncHistos_pt540.root") ;
////2.76 TeV inclusive x sec ntuples
	string infile276t(indirEOS+"/Run2013APP276_8tevCorr/ntuples_incxsec.root"),outfile276t(outdir+"/AnalyzerIncHistosData276TeV_8tevcorr.root") ;
///276 tev
  vector<string> infiles276t;
	infiles276t.push_back(infile276t);
	
	vector<string> infilespt15 ;infilespt15 .push_back(inpt15 );
	vector<string> infilespt30 ;infilespt30 .push_back(inpt30 );
	vector<string> infilespt50 ;infilespt50 .push_back(inpt50 );
	vector<string> infilespt80 ;infilespt80 .push_back(inpt80 );
	vector<string> infilespt120;infilespt120.push_back(inpt120);
	vector<string> infilespt170;infilespt170.push_back(inpt170);
	vector<string> infilespt220;infilespt220.push_back(inpt220);
	vector<string> infilespt280;infilespt280.push_back(inpt280);
	vector<string> infilespt370;infilespt280.push_back(inpt370);
	vector<string> infilespt460;infilespt280.push_back(inpt460);
	vector<string> infilespt540;infilespt280.push_back(inpt540);
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
	cout<<setprecision(4);

	bool runOnMC=false;
//	cout<<"Data =========================================="<<endl;
//	jesAnalyzer dataic(infiles276t,outfile276t,runOnMC); //IncXsec ntuples
//	dataic.makejecUncHists();
	 
	runOnMC=true;
  int num = strtol(argv[1],NULL,10);
  cout<<"RUnning on "<<num<<endl;
	for(unsigned i(num);i<num+1 && i<inputfilesAllMC.size();++i) {
		vector<string> input; input.push_back(inputfilesAllMC[i]);
		ostringstream oss; oss<<"AnalyzerIncHistos_allpt"<<i<<"_winter14.root";
		jesAnalyzer allmc (input,oss.str(),runOnMC );
		allmc.makejecUncHists();
	}
		
return 0;
}

void jesAnalyzer::makejecUncHists()
{

	const int nsrc = 22;
	//winter14_v5
	const char* srcnames[nsrc] ={
		"AbsoluteStat","AbsoluteScale","AbsoluteFlavMap","AbsoluteMPFBias",
		"HighPtExtra","SinglePionECAL","SinglePionHCAL","FlavorQCD",
		"TimeEta","TimePt","RelativeJEREC1","RelativeJEREC2","RelativeJERHF",
		"RelativePtBB","RelativePtEC1","RelativePtEC2","RelativePtHF",
		"RelativeFSR","RelativeStatEC2","RelativeStatHF","PileUpMuZero",
		"TimeRunD"
	};


	////winter14_v4
	//	const char* srcnames[nsrc] = {
	//		"AbsoluteStat", "AbsoluteScale","AbsoluteFlavMap", "AbsoluteMPFBias",
	//		"HighPtExtra","SinglePionECAL","SinglePionHCAL","FlavorQCD","Time",
	//		"RelativeJEREC1","RelativeJEREC2","RelativeJERHF","RelativePtBB",
	//		"RelativePtEC1","RelativePtEC2","RelativePtHF","RelativeFSR",
	//		"RelativeStatEC2","RelativeStatHF",
	//		"FlavorZJet","FlavorPhotonJet","FlavorPureGluon","FlavorPureQuark",
	//		"FlavorPureCharm","FlavorPureBottom",
	//		,"PileUpDataMC"	,"PileUpPtBB"	,"PileUpPtEC1"	,"PileUpPtEC2"
	//			,"PileUpPtHF"		,"PileUpBias"
	//	};
	// Instantiate uncertainty sources Summer13_V4, or Fall12_V7
	//		const int nsrc = 22;
	//		const char* srcnames[nsrc] =
	//		{"Absolute", "HighPtExtra",  "SinglePionECAL", "SinglePionHCAL",
	//			"FlavorQCD", "Time",
	//			"RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
	//			"RelativePtBB","RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeFSR",
	//			"RelativeStatEC2", "RelativeStatHF",
	//			"FlavorZJet","FlavorPhotonJet","FlavorPureGluon","FlavorPureQuark","FlavorPureCharm","FlavorPureBottom"};

	vector<JetCorrectionUncertainty*> vsrc(nsrc);
	for (int isrc = 0; isrc < nsrc; ++isrc) {
		const char *uncsrc = srcnames[isrc];
		JetCorrectorParameters *p = new JetCorrectorParameters("Winter14_V5_DATA_UncertaintySources_AK7PFchs.txt", uncsrc);
		JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
		vsrc[isrc] = unc;
	}

	cout<<setprecision(6);
	cout << scientific;
	///declare hists ===========================================
	TH1D* hptUnf [nybins-1]; //pt spectrum, originally used for unfolding
	TH1D* hjesUnc [nsrc+1][nybins-1][2];//nsrc+1 (total uncertainty), 0-down,1-up
	TProfile* huncFactors[nsrc+1][nybins-1][3];//errup, errdn, maxerr
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
				//average jec uncertainties
				sprintf(name,"huncFac_%s_%s_eta%d",srcnames[isrc],ud,i);
				huncFactors[isrc][i][j] = new TProfile(name,name,100,10,1000,-5.,5.);
				huncFactors[isrc][i][j]->Sumw2();
				if(j==0)
				{
					sprintf(name,"huncFac_%s_max_eta%d",srcnames[isrc],i);
					huncFactors[isrc][i][2] = new TProfile(name,name,100,10,1000,-5.,5.);
					huncFactors[isrc][i][2]->Sumw2();
				}
				//average jec uncertainties eend
			}//isrc
			sprintf(name,"TotalUnc_%s_eta%d",ud,i);
			hjesUnc[nsrc][i][j] = new TH1D(name,name,n1x[i],&x1[i][0]);
			hjesUnc[nsrc][i][j]->Sumw2();
			//average jec uncertainties
			sprintf(name,"huncFac_total_%s_eta%d",ud,i);
			huncFactors[nsrc][i][j] = new TProfile(name,name,100,10,1000,-5.,5.);
			huncFactors[nsrc][i][j]->Sumw2();
			if(j==0)
			{
				sprintf(name,"huncFac_total_max_eta%d",i);
				huncFactors[nsrc][i][2] = new TProfile(name,name,100,10,1000,-5.,5.);
				huncFactors[nsrc][i][2]->Sumw2();
			}
		}//up/down
	}//eta loop
	unsigned nentries(treeak7->fChain->GetEntries());
	//  nentries=5e2;
	Init();	//Initiate local variables
	cout<<"Running on jet scale uncertainty plots "<<nentries<<endl;
	for(unsigned jentry(0);jentry<nentries;++jentry) {
		////load entry 
		treeak7->fChain->GetEntry(jentry);
		// Dont process events outside pt window-----------------
		if(runOnMC && (*qscale)>maxpthat) continue;
		//////------Determine prescale for the event------------
		int prescl[numtrgs];
		for(unsigned ntrg(0);ntrg<numtrgs;++ntrg)  
			prescl[ntrg] = (*(hltPsc[ntrg]))*(*(hltl1Psc[ntrg])); //psc with hltl1pass

		for(unsigned ntrg(1);ntrg<numtrgs;++ntrg) //hlt20 is excluded  
			if((*(hltPass[ntrg]))>0.5)     {      //hlt pass 
				for(unsigned njet(0);njet<(*nref);++njet ) {
					bool cut(true);
					double chf(0), nhf(0), emf(0),cemf(0),nemf(0),
								 muf(0), elf(0), phf(0);
					nhf  = neutralHardSum[njet];
					muf  = muSum[njet];
					elf  = eSum[njet];
					//do jet id and pt cut here-----------------------
					cut  =
						jtpt[njet]>=20            &&
						(*nPV) >= 1               &&
						(*met)/(*sumET) <= 0.3    &&              
						(elf < 0.9) && 
						(muf < 0.9) &&
						(nhf < 0.9) &&
						jetId(njet)               
						;
					if(cut) 
						for(unsigned i(0);i<nybins-1;++i) ///Cycle over eta bins 
							if(fabs(jty[njet])>=ybins[i] && fabs(jty[njet])<ybins[i+1]) {
								//////------Determine event weight ---------------------
								double wt = 1.0;
								if(runOnMC) wt = mcWeight;
								else        wt = prescl[ntrg];

								// Calculate uncertainty per source 
								double sum2_up(0), sum2_dw(0);
								for (int isrc = 0; isrc < nsrc; isrc++) {
									JetCorrectionUncertainty *unc = vsrc[isrc];
									unc->setJetPt(jtpt[njet]);
									unc->setJetEta(jteta[njet]);
									double sup = unc->getUncertainty(true); // up variation
									unc->setJetPt(jtpt[njet]);
									unc->setJetEta(jteta[njet]);
									double sdw = unc->getUncertainty(false); // down variation
									//take the largest variation
									double maxerr = max(fabs(sup),fabs(sdw));
									/////total unc not used in macros down the chain
									sum2_up += pow(max(sup,sdw),2);
									sum2_dw += pow(min(sup,sdw),2);
									//                  cout<<setprecision(3); cout<<fixed;
									if(false)
										cout<<"jet pt:eta "<<jtpt[njet]<<":"<<jteta[njet]
											<<" "<<srcnames[isrc]<<" "
											<<" "<<sup<<"/"<<sdw<<endl;
									/////total unc not used in macros down the chain
									double ptup((1.+maxerr)*jtpt[njet]);
									double ptdn((1.-maxerr)*jtpt[njet]);
									if (ptup<HLTJetPtTh[ntrg+1] && ptup>=HLTJetPtTh[ntrg])
										hjesUnc[isrc][i][1]->Fill(ptup,wt);
									if (ptdn<HLTJetPtTh[ntrg+1] && ptdn>=HLTJetPtTh[ntrg])
										hjesUnc[isrc][i][0]->Fill(ptdn,wt);
									/////average jec uncertainties
									huncFactors[isrc][i][0]->Fill(jtpt[njet],sdw,1);
									huncFactors[isrc][i][1]->Fill(jtpt[njet],sup,1);
									huncFactors[isrc][i][2]->Fill(jtpt[njet],maxerr,1);
								} // for isrc
								//total uncertainty
								double ptup((1.+sqrt(sum2_up))*jtpt[njet]);
								double ptdn((1.-sqrt(sum2_dw))*jtpt[njet]);
								if (ptup<HLTJetPtTh[ntrg+1] && ptup>=HLTJetPtTh[ntrg])
									hjesUnc[nsrc][i][1]->Fill(ptup,wt);
								if (ptdn<HLTJetPtTh[ntrg+1] && ptdn>=HLTJetPtTh[ntrg])
									hjesUnc[nsrc][i][0]->Fill(ptdn,wt);
								/////average jec uncertainties
								huncFactors[nsrc][i][0]->Fill(jtpt[njet],sqrt(sum2_dw),1);
								huncFactors[nsrc][i][1]->Fill(jtpt[njet],sqrt(sum2_up),1);
								huncFactors[nsrc][i][2]->Fill(jtpt[njet],0,1);

								///Filling Spectrum histogram-----------------------
								if(    jtpt[njet]<HLTJetPtTh[ntrg+1]  //hlt pt bin
										&& jtpt[njet]>=HLTJetPtTh[ntrg])  //hlt pt bin
								{
									hptUnf[i]->Fill(jtpt[njet],wt);		
								}//hlt pt bin
							}///eta loop
				}//for loop on jets

			}////trigger loop
		if(!(jentry%200000))
			cout<<"Processing Event "<<(float)jentry<<endl;
	}//nentry loop
	///OUTPUT==============================================================
	string unfFile(outfile.substr(0,outfile.find(".root")));
	TFile* fstab = new TFile((unfFile+"_jecUnc.root").c_str(),"RECREATE");
	fstab->cd();
	//write histos 
	for(unsigned i(0);i<nybins-1;++i) {
		hptUnf[i]->Write();
		hjesUnc[nsrc][i][0]->Write();
		hjesUnc[nsrc][i][1]->Write();
		for (int isrc = 0; isrc <= nsrc; ++isrc) {
			hjesUnc[isrc][i][0]->Write();
			hjesUnc[isrc][i][1]->Write();
			huncFactors[isrc][i][0]->Write();
			huncFactors[isrc][i][1]->Write();
			huncFactors[isrc][i][2]->Write();
		}
	}
	fstab->Write();
	fstab->Close();
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


void jesAnalyzer::Init() {
	treeak7->fChain->SetBranchStatus("*",1);
	//	treeak5->fChain->SetBranchStatus("*",1);
	HLTJetPt  = HLTJetPtN; //HLTJetPtN8Tev;
	HLTJetPtTh= HLTJetPtS; //HLTJetPtS8Tev;
	trgnames  = trignames; //trignames8Tev;
	numtrgs   = numtrigs;  //numtrigs8Tev;
	jtpt = treeak7->jtpt  ;     chargedHardSum= treeak7->chargedHardSum ;
	jteta= treeak7->jteta ;     chargedHardN  = treeak7->chargedHardN   ;
	jtphi= treeak7->jtphi ;     neutralHardSum= treeak7->neutralHardSum ;
	jtm  = treeak7->jtm   ;     neutralHardN  = treeak7->neutralHardN   ;
	jty  = treeak7->jty   ;     emEnergy      = treeak7->emEnergy       ;
	rawpt= treeak7->rawpt ;     photonSum     = treeak7->photonSum      ;
	nPV  = &treeak7->nPV  ;     muSum         = treeak7->muSum          ;
	rho  = &treeak7->rho  ;     eSum          = treeak7->eSum           ;
	met  = &treeak7->met  ;	    trackN        = treeak7->trackN         ;
	sumET= &treeak7->sumET;     nuemEnergy    = treeak7->nuemEnergy     ;
	nref = &treeak7->nref ;     chemEnergy    = treeak7->chemEnergy     ;
	genpt = treeak7->genpt  ;   constituentN  = treeak7->constituentN   ;
	geneta= treeak7->geneta ;   
	genphi= treeak7->genphi ;
	genm  = treeak7->genm   ;  

	////trigger specific variables 
	run       = &treeak7->event_runNo;
	event     = &treeak7->event_evtNo;
	lumiBlock = &treeak7->event_lumi;
	qscale    = &treeak7->qscale    ;
	hltPass[0]= &treeak7->HLT_PAJet20_NoJetID_v1;      
	hltPass[1]= &treeak7->HLT_PAJet40_NoJetID_v1;      
	hltPass[2]= &treeak7->HLT_PAJet60_NoJetID_v1;
	hltPass[3]= &treeak7->HLT_PAJet80_NoJetID_v1;
	hltPass[4]= &treeak7->HLT_PAJet100_NoJetID_v1;
	hltPass[5]= &treeak7->HLT_PAJet120_NoJetID_v1;
	hltPsc[0]=  &treeak7->HLT_PAJet20_NoJetID_v1_Prescl;
	hltPsc[1]=  &treeak7->HLT_PAJet40_NoJetID_v1_Prescl;
	hltPsc[2]=  &treeak7->HLT_PAJet60_NoJetID_v1_Prescl;
	hltPsc[3]=  &treeak7->HLT_PAJet80_NoJetID_v1_Prescl;
	hltPsc[4]=  &treeak7->HLT_PAJet100_NoJetID_v1_Prescl;
	hltPsc[5]=  &treeak7->HLT_PAJet120_NoJetID_v1_Prescl;
	//
	hltl1Psc[0]=&treeak7->HLT_PAJet20_NoJetID_v1L1_Prescl; 
	hltl1Psc[1]=&treeak7->HLT_PAJet40_NoJetID_v1L1_Prescl; 
	hltl1Psc[2]=&treeak7->HLT_PAJet60_NoJetID_v1L1_Prescl; 
	hltl1Psc[3]=&treeak7->HLT_PAJet80_NoJetID_v1L1_Prescl; 
	hltl1Psc[4]=&treeak7->HLT_PAJet100_NoJetID_v1L1_Prescl;
	hltl1Psc[5]=&treeak7->HLT_PAJet120_NoJetID_v1L1_Prescl;
	//
	hltobj[0] =&treeak7->HLT_PAJet20_NoJetID_v1_ObjPt ;
	hltobj[1] =&treeak7->HLT_PAJet40_NoJetID_v1_ObjPt ;
	hltobj[2] =&treeak7->HLT_PAJet60_NoJetID_v1_ObjPt;
	hltobj[3] =&treeak7->HLT_PAJet80_NoJetID_v1_ObjPt;
	hltobj[4] =&treeak7->HLT_PAJet100_NoJetID_v1_ObjPt;
	hltobj[5] =&treeak7->HLT_PAJet120_NoJetID_v1_ObjPt;
}//end iniit

void jesAnalyzer::fillArrayFromVec(const vector<float>& vec,
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

bool jesAnalyzer::jetId(unsigned jetnum)
{
	double chf  = chargedHardSum[jetnum];
	unsigned chm= chargedHardN[jetnum];
	double nhf  = neutralHardSum[jetnum];
	double muf  = muSum[jetnum];
	double elf  = eSum[jetnum];
	double phf  = photonSum[jetnum];
	double cemf = chemEnergy[jetnum];
	double nemf = nuemEnergy[jetnum];
	double eta  = jteta[jetnum];
	int constituents  = constituentN[jetnum];
	bool id = constituents>1  && nhf<0.99 && phf<0.99 && ((fabs(eta)<=2.4 && 
				nhf<0.9 &&phf<0.9&& elf<0.99 &&chf>0 && chm>0 ) || fabs(eta)>2.4) ;
	//	bool id = constituents>1 && nhf<0.9 && ((fabs(eta)<=2.4 && 
	//				chf>0 && chm>0 && cemf<0.99) || fabs(eta)>2.4) ;
	//	bool id  = 
	//		muf < 0.9 && nhf < 0.9 &&  
	//		elf < 0.9 && phf < 0.9 &&
	//		nemf< 0.9 && 
	//		constituents>1 &&
	//		((fabs(eta)<=2.4  &&
	//			chargedHardN  [jetnum] > 0  &&
	//			chf > 0.0   &&  cemf< 0.99) ||
	//		 fabs(eta)>2.4);

	return	(id);
}

void jesAnalyzer::setMCWeight()
{
	double mcSampleWt[nMCSamples];
	for(unsigned i(0);i<nMCSamples;++i) {
		float xsec(0); 
		if(i<nMCSamples-1)
			xsec = mcSampleXsecmb[i]-mcSampleXsecmb[i+1];
		else 
			xsec = mcSampleXsecmb[i] - 0.;
		mcSampleWt[i]=lumi*(xsec)*1.0e9/mcSampleEvents[i];
		//cout<<"================================"<<endl;
		//cout<<"xsec "<<xsec<<" pt "<<mcSamplePt[i];
	}

	if(runOnMC){
		for(unsigned i(0);i<nMCSamples;++i)
		{
			sprintf(name,"pt%d.root",mcSamplePt[i]);
			string sam(name);
			assert(inputfiles.size());
			for(unsigned j(0);j<inputfiles.size();++j) 
				if(inputfiles[j].find(sam)!=string::npos) {
					cout<<"================================"<<endl;
					cout<<"Applying weight for sample "<<sam<<endl;
					if(mcWeight==-1) {
						mcWeight =  mcSampleWt[i];
						pthat = (float)mcSamplePt[i];
						maxpthat = i<nMCSamples-1?(float)mcSamplePt[i+1]:999999.;
					}
					else if (mcSampleWt[i]== mcWeight) 
						cout<<endl;
					else 
					{
						cout<<"\nMC weight is already set for this sample:"<<mcWeight
							<<"\nRunning on file:"<<inputfiles[j]
							<<"\nPlease make sure only one pt sample is run at a time.\n"
							<<endl;
						assert(0);
					}
				}
		}
		cout<<fixed;
		cout<<scientific;
		cout<<setprecision(3);
		cout<<"Sample weight "<<mcWeight<<endl;
		cout<<"pthat "<<pthat<<" maxpthat "<<maxpthat<<endl;
		cout<<"================================"<<endl;
	}
}


