#define runJetAnalyzer_cxx
#include "runJetAnalyzer.h"
#include "math.h"
#include "TRandom.h"
#include "TProfile.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdlib.h>
#include <sstream>      // std::ostringstream
#include "AnalyzerProcedures.h"
#include "stability.cc"

void runJetAnalyzer::Loop()
{
	unsigned nentries(treeak7->fChain->GetEntries());
	//Initiate local variables
	Init();
	//For Checking Data Stability
//	makeStabilityHists();
	cout<<setprecision(4);

	const int maxPtBins(50);
	//   	 //declare histos
	TH1D* hpt [numtrgs][nybins-1]; // for spectrum 
	TH1D* PFJet0[numtrgs][nybins-1]; // for hlt thresholds
	TH1D* PFJet1[numtrgs][nybins-1];
	TH1D* heta[numtrgs][nybins-1];  TH1D* hchf[numtrgs][nybins-1];
	TH1D* hphi[numtrgs][nybins-1];  TH1D* hnhf[numtrgs][nybins-1];
	TH1D* hjec[numtrgs][nybins-1];  TH1D* hmuf[numtrgs][nybins-1];
	TH1D* htrk[numtrgs][nybins-1];  TH1D* helf[numtrgs][nybins-1];
	TH1D* hmet[numtrgs][nybins-1];  TH1D* hphf[numtrgs][nybins-1];
  TH1D* hcef[numtrgs][nybins-1];  TH1D* hnef[numtrgs][nybins];
	TH1D* hres[nybins-1][maxPtBins]; // pt resolution
	TH1D* hgen[nybins-1]; 					 // gen spectrum
	TH1D* hgenUnf[nybins-1]; 					 // lo gen spectrum for unfolding
	TH1D* hptUnf[nybins-1]; // to unfold spectrum
	TH1D* hraw[nybins-1]; // raw  spectrum 
	TH1D* hcor[nybins-1]; // raw  spectrum 
	TH1D* htjec[nybins-1][maxPtBins];//test jec
	TH1D* htchf[nybins-1][maxPtBins];//test chf
	TH1D* htnhf[nybins-1][maxPtBins];//test nhf
	TH1D* hnPV[2];
  TH1D* hprb[numtrgs][nybins-1];//to calculate jet efficiency 
  TH1D* htag[numtrgs][nybins-1];//to calculate jet efficiency 
  TProfile* hjecl1  [nybins-1];
  TProfile* hjecl2l3[nybins-1];
  TProfile* hjecres [nybins-1];
	TH1D* hrecMC[nybins-1]; // to unfold spectrum
	TH1D* hgenMC[nybins-1]; // to unfold spectrum
	//instantialize histos-----------npv
	sprintf(name,"nPVequals1");
	hnPV[0] = new TH1D(name,name,2,0,2);
	sprintf(name,"nPVgreaterthan1");
	hnPV[1] = new TH1D(name,name,2,0,2);
	//instantialize histos-----------hlt efficiency plots
	for(unsigned ntrg(0);ntrg<numtrgs;++ntrg)
		for(unsigned i(0);i<nybins-1;++i)
		{
			//hlt histos
			sprintf(name,"PFJet0_HLT%i_eta%i",HLTJetPt[ntrg],i);
			PFJet0[ntrg][i] = new TH1D(name,name,nxhlt[i],&xhlt[i][0]);
			PFJet0[ntrg][i]->Sumw2();
			////
			sprintf(name,"PFJet1_HLT%i_eta%i",HLTJetPt[ntrg],i);
			PFJet1[ntrg][i] = new TH1D(name,name,nxhlt[i],&xhlt[i][0]);
			PFJet1[ntrg][i]->Sumw2();
		}
	///-------------Spectrum plot---------------------------
	for(unsigned ntrg(0);ntrg<numtrgs;++ntrg)
		for(unsigned i(0);i<nybins-1;++i)
		{
		sprintf(name,"hprb_HLT%i_eta%i",HLTJetPt[ntrg],i);
			hprb[ntrg][i] = new TH1D(name,name,nx[i],&x[i][0]);
			hprb[ntrg][i]->Sumw2();
		sprintf(name,"htag_HLT%i_eta%i",HLTJetPt[ntrg],i);
			htag[ntrg][i] = new TH1D(name,name,nx[i],&x[i][0]);
			htag[ntrg][i]->Sumw2();
      ///
			ostringstream oss; 
			oss<<"hpt_eta_"<<ybins[i]<<"_"<<ybins[i+1]
				<<"_"<<trgnames[ntrg];
			hpt[ntrg][i] = new TH1D(oss.str().c_str(),oss.str().c_str(),nx[i],&x[i][0]);
			hpt[ntrg][i]->Sumw2();
     
		   oss.str("");
			oss<<"hptUnf_eta_"<<ybins[i]<<"_"<<ybins[i+1];
			if(ntrg==0) {
				hptUnf[i] = new TH1D(oss.str().c_str(),oss.str().c_str(),n1x[i],&x1[i][0]);	
				hptUnf[i]->Sumw2();
				sprintf(name,"hrawpt_eta%d",i);
				hraw[i] = new TH1D(name,name,n1x[i],&x1[i][0]);	
				hraw[i]->Sumw2();
        sprintf(name,"hcorpt_eta%d",i);
				hcor[i] = new TH1D(name,name,n1x[i],&x1[i][0]);	
				hcor[i]->Sumw2();
        sprintf(name,"jecl1_eta%i",i);
        hjecl1 [i]= new TProfile(name,name,n1x[i],&x1[i][0],0,4);
        sprintf(name,"jecl2l3_eta%i",i);
        hjecl2l3 [i]= new TProfile(name,name,n1x[i],&x1[i][0],0,4);
        sprintf(name,"jecres_eta%i",i);
        hjecres [i]= new TProfile(name,name,n1x[i],&x1[i][0],0,4);
        sprintf(name,"hrecMC_eta%i",i);
			 hrecMC[i] = new TH1D(name,name,n1x[i],&x1[i][0]);	
			 hrecMC[i]->Sumw2();
       sprintf(name,"hgenMC_eta%i",i);
			 hgenMC[i] = new TH1D(name,name,n1x[i],&x1[i][0]);	
			 hgenMC[i]->Sumw2();
			}
			//histograms for jet id variables--------------------
			oss.str("");
			oss<<"hphi_eta_"<<ybins[i]<<"_"<<ybins[i+1]
				<<"_"<<trgnames[ntrg];
			hphi[ntrg][i] = new TH1D(oss.str().c_str(),oss.str().c_str(),30,-M_PI,M_PI);				
			oss.str("");
			oss<<"hchf_eta_"<<ybins[i]<<"_"<<ybins[i+1]
				<<"_"<<trgnames[ntrg];
			hchf[ntrg][i] = new TH1D(oss.str().c_str(),oss.str().c_str(),30,0,1);				
			oss.str("");
			oss<<"hnhf_eta_"<<ybins[i]<<"_"<<ybins[i+1]
				<<"_"<<trgnames[ntrg];
			hnhf[ntrg][i] = new TH1D(oss.str().c_str(),oss.str().c_str(),30,0,1);			
			oss.str("");
			oss<<"hnef_eta_"<<ybins[i]<<"_"<<ybins[i+1]
				<<"_"<<trgnames[ntrg];
			hnef[ntrg][i] = new TH1D(oss.str().c_str(),oss.str().c_str(),30,0,1);			
			oss.str("");
			oss<<"hcef_eta_"<<ybins[i]<<"_"<<ybins[i+1]
				<<"_"<<trgnames[ntrg];
			hcef[ntrg][i] = new TH1D(oss.str().c_str(),oss.str().c_str(),30,0,1);			
			oss.str("");
			oss<<"hphf_eta_"<<ybins[i]<<"_"<<ybins[i+1]
				<<"_"<<trgnames[ntrg];
			hphf[ntrg][i] = new TH1D(oss.str().c_str(),oss.str().c_str(),30,0,1);
			oss.str("");
			oss<<"heta_eta_"<<ybins[i]<<"_"<<ybins[i+1]
				<<"_"<<trgnames[ntrg];
			heta[ntrg][i] = new TH1D(oss.str().c_str(),oss.str().c_str(),30,-5,5);
			oss.str("");
			oss<<"hmuf_eta_"<<ybins[i]<<"_"<<ybins[i+1]
				<<"_"<<trgnames[ntrg];
			hmuf[ntrg][i] = new TH1D(oss.str().c_str(),oss.str().c_str(),30,0,1);	
			oss.str("");
			oss<<"helf_eta_"<<ybins[i]<<"_"<<ybins[i+1]
				<<"_"<<trgnames[ntrg];
			helf[ntrg][i] = new TH1D(oss.str().c_str(),oss.str().c_str(),30,0,1);;
			oss.str("");
			oss<<"hjec_eta_"<<ybins[i]<<"_"<<ybins[i+1]
				<<"_"<<trgnames[ntrg];
			hjec[ntrg][i] = new TH1D(oss.str().c_str(),oss.str().c_str(),80,0.5,1.5);	
			oss.str("");
			oss<<"htrk_eta_"<<ybins[i]<<"_"<<ybins[i+1]
				<<"_"<<trgnames[ntrg];
			htrk[ntrg][i] = new TH1D(oss.str().c_str(),oss.str().c_str(),30,0,30);	
			oss.str("");
			oss<<"hmet_eta_"<<ybins[i]<<"_"<<ybins[i+1]
				<<"_"<<trgnames[ntrg];
			hmet[ntrg][i] = new TH1D(oss.str().c_str(),oss.str().c_str(),30,0,5);	
			heta[ntrg][i]->Sumw2(); hmuf[ntrg][i]->Sumw2();
			hphi[ntrg][i]->Sumw2(); helf[ntrg][i]->Sumw2();
			hchf[ntrg][i]->Sumw2(); hphf[ntrg][i]->Sumw2();
			hnhf[ntrg][i]->Sumw2(); hjec[ntrg][i]->Sumw2();
			htrk[ntrg][i]->Sumw2(); hmet[ntrg][i]->Sumw2();
			hnef[ntrg][i]->Sumw2(); hcef[ntrg][i]->Sumw2();
		}
	//Initialize pt resolution histo
	for(unsigned i(0);i<nybins-1;++i)
		if(runOnMC) 
		{
			for(unsigned ptb(0);ptb<n1x[i];++ptb)
			{
				//Initialize resolution histos
				sprintf(name, "hres_eta%i_bin%i",i,ptb);
				hres[i][ptb] = new TH1D(name,name,400,0.2,1.8);
				hres[i][ptb]->Sumw2();
				sprintf(name, "htjec_eta%i_bin%i", i,ptb);
				htjec[i][ptb] = new TH1D(name,name,400,0,4);
				sprintf(name, "htchf_eta%i_bin%i", i,ptb);
				htchf[i][ptb] = new TH1D(name,name,400,0,4);
				sprintf(name, "htnhf_eta%i_bin%i", i,ptb);
				htnhf[i][ptb] = new TH1D(name,name,400,0,4);

			}
			//Initialize gen spectrum histo
			sprintf(name, "hgen_eta%i", i);
			//			hgen[i] =  new TH1D(name,name,n1x[0],&x1[0][0]);
			hgen[i] =  new TH1D(name,name,1000,0,1000);
			hgen[i]->Sumw2();
			sprintf(name, "hgenUnf_eta%i", i);
			hgenUnf[i] =  new TH1D(name,name,n1x[0],&x1[0][0]);
			hgenUnf[i]->Sumw2();
		}

//		nentries=1000;
	//////------Run on Entries -----------------------------
	cout<<"Running on entries "<<nentries<<endl;
	for(unsigned jentry(0);jentry<nentries;++jentry)
	{
		////Vectors for filling unfold ntuple
		vector<float> jptUnf;   vector<float> gjptUnf; 
		vector<float> jetaUnf;  vector<float> gjetaUnf;
		vector<float> jphiUnf;  vector<float> gjphiUnf;
		vector<float> jmUnf;    vector<float> gjyUnf;  	
		////load entry 
		treeak7->fChain->GetEntry(jentry);
		// Dont process events outside pt window-----------------
		if(runOnMC && (*qscale)>maxpthat) continue;
		//////------Determine prescale for the event------------
		int prescl[numtrgs];
		for(unsigned ntrg(0);ntrg<numtrgs;++ntrg)  
			prescl[ntrg] = (*(hltPsc[ntrg]))*(*(hltl1Psc[ntrg])); //psc with hltl1pass
		///-------Filling histograms----------------------------
		for(unsigned ntrg(0);ntrg<numtrgs;++ntrg)  
		{
			bool hlton  =(*(hltPass[ntrg]))>0.5 ;
			if  (trgnames[ntrg].find("Jet80")!=string::npos)
				hlton = (*(hltPass[ntrg]))>0.5 || 
					(*(hltPass[ntrg+1]))>0.5 ||  //hlt100 
					(*(hltPass[ntrg+2]))>0.5 ;   //hlt120
			if(hlton)//hlt pass 
			{	
				bool hltcut = false;
				if(ntrg!=numtrgs-1) {
					hltcut = (*(hltobj[ntrg]))>HLTJetPt[ntrg+1];
				}
				njetUnf=0;
				for(unsigned njet(0);njet<(*nref);++njet ) {
					bool cut(true);
					double chf(0), nhf(0), emf(0),cemf(0),nemf(0),
								 muf(0), elf(0), phf(0);
					nhf  = neutralHardSum[njet];
					phf  = photonSum[njet];
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
						//		((fabs(jteta[njet])<=2.4 && trackN[njet]>4) || 
						//				(fabs(jteta[njet])>2.4)) &&
						jetId(njet)               
						;
					if(cut) 
						for(unsigned i(0);i<nybins-1;++i) ///Cycle over eta bins 
							if(fabs(jty[njet])>=ybins[i] && fabs(jty[njet])<ybins[i+1]) {
								///--------HLT Thresholds--------------------------
								PFJet0[ntrg][i]->Fill(jtpt[njet]);
								if(hltcut) //passed trigger
									//	PFJet1[ntrg][i]->Fill(jtpt[njet],prescl[ntrg+1]);
									PFJet1[ntrg][i]->Fill(jtpt[njet]);
								///--------HLT Thresholds Done----------------------
								//////------Determine event weight ---------------------
								double wt = 1.0;
								if(runOnMC) wt = mcWeight;
								else        wt = prescl[ntrg];
								///look at raw and corrected spectrum ----------------
								if(trgnames[ntrg].find("Jet20")==string::npos)//exclude hlt20
								{
									hraw[i]->Fill(rawpt[njet],wt);
									hcor[i]->Fill(jtpt[njet],wt);
								}
								///Filling Spectrum histogram-----------------------
								bool ptcut = (fabs(jty[njet])< 2.0  && fabs(jty[njet])>= 1.5  &&jtpt[njet]> 272) ||
									(fabs(jty[njet])< 2.5  && fabs(jty[njet])>= 2.0  &&jtpt[njet]> 153) 
									;
								if(    jtpt[njet]<HLTJetPtTh[ntrg+1]  //hlt pt bin
										&& jtpt[njet]>=HLTJetPtTh[ntrg])  //hlt pt bin
									if(true)
									{

										hpt[ntrg][i] ->Fill(jtpt[njet],wt);
										if(trgnames[ntrg].find("Jet20")==string::npos)
											hptUnf[i]    ->Fill(jtpt[njet],wt);
										hjecl1[i]->Fill(jtpt[njet],jecl1[njet],1);
										hjecl2l3[i]->Fill(jtpt[njet],jecl2l3[njet],1);
										hjecres[i]->Fill(jtpt[njet],jecres[njet],1);
										heta[ntrg][i]->Fill((jteta[njet]),wt);
										hphi[ntrg][i]->Fill(jtphi[njet],wt);
										hchf[ntrg][i]->Fill(chargedHardSum[njet],wt);
										hcef[ntrg][i]->Fill(chemEnergy[njet],wt);
										hnef[ntrg][i]->Fill(nuemEnergy[njet],wt);
										hnhf[ntrg][i]->Fill(neutralHardSum[njet],wt);
										hmuf[ntrg][i]->Fill(muSum[njet],wt);
										helf[ntrg][i]->Fill(eSum[njet],wt);
										hphf[ntrg][i]->Fill(photonSum[njet],wt);
										hjec[ntrg][i]->Fill(jtpt[njet]/rawpt[njet],wt);
										htrk[ntrg][i]->Fill(trackN[njet],wt);
										if(njet==0)
											hmet[ntrg][i]->Fill((*met)/(*sumET),wt);
										cout<<fixed<<setprecision(3);
										if(false)
											cout
												<<" jet "<<njet<<" pt "<<jtpt[njet]<<" eta "<<jty[njet]
												<<" phi "<<jtphi[njet]<<" wt "<<wt<<'\n'
												<<" chf "<< chf <<" muf "<< muf<<'\n'
												<<" nhf "<< nhf <<" phf "<< phf<<'\n'
												<<" elf "<< elf <<" tracks "<<trackN[njet]<<'\n'
												<<" npv "<<(*nPV) <<" rho "<<(*rho)<<'\n'
												<<" met "<<(*met)/(*sumET) <<" jec "<<jtpt[njet]/rawpt[njet]<<'\n'
												<<" run "<<*run<<" event "<<*event<<" lumiBlock "<<*lumiBlock<<'\n'
												<<"total jets "<<(*nref)
												<<"==========================="
												<<endl;
										cout<<fixed<<setprecision(3)<<scientific;

									}// fill specctrum 
								///Filling Spectrum histogram-----------------------
							}///eta loop
				}//for loop on jets
				//-------- jet id efficiency --------------------------------------------
				if((!runOnMC) ){
					int irand1 = (gRandom->Uniform()>0.5)?0:1;
					int irand2 = (irand1+1)%2;
					bool cut = jetIdTagPrb(irand1) && 
						jetIdTagPrb(irand2) && 
						(*met)/(*sumET) < 0.3 &&
						jtpt[irand1]> 20 && 
						jtpt[irand2]> 20 && 
						(*nPV) >= 1   	;
					double dphi(deltaPhi(jtphi[irand1],jtphi[irand2]));
					//cut on dijet events
					if(cut && fabs(dphi)>2.7 && (*nref)==2) {
						for(unsigned i(0);i<nybins-1;++i) {
							if(fabs(jty[irand2])>=ybins[i] && fabs(jty[irand2])<ybins[i+1]) {
								if(jetId(irand1) ){
									if(    jtpt[irand2]<HLTJetPtTh[ntrg+1]  //hlt pt bin
											&& jtpt[irand2]>=HLTJetPtTh[ntrg]){  //hlt pt bin

										htag[ntrg][i]->Fill(jtpt[irand2],prescl[ntrg]);

										if(jetId(irand2))
											hprb[ntrg][i]->Fill(jtpt[irand2],prescl[ntrg]); 



										///testt test ------------------------------------
										double chf  = chargedHardSum[irand2];
										unsigned chm= chargedHardN[irand2];
										double nhf  = neutralHardSum[irand2];
										double muf  = muSum[irand2];
										double elf  = eSum[irand2];
										double phf  = photonSum[irand2];
										double cemf = chemEnergy[irand2];
										double nemf = nuemEnergy[irand2];
										double eta  = jteta[irand2];
										int constituents  = constituentN[irand2];

										if(jtpt[irand2]>=330 && jtpt[irand2]<362 && (!jetId(irand2)))
											cout
												<<" jet "<<irand2<<" pt "<<jtpt[irand2]<<" eta "<<jty[irand2]
												<<" phi "<<jtphi[irand2]<<" wt "<<prescl[ntrg]<<'\n'
												<<" chf "<< chf <<" muf "<< muf<<'\n'
												<<" nhf "<< nhf <<" phf "<< phf<<'\n'
												<<" elf "<< elf <<" tracks "<<trackN[irand2]<<'\n'
												<<" npv "<<(*nPV) <<" rho "<<(*rho)<<'\n'
												<<" met "<<(*met)/(*sumET) <<" jec "<<jtpt[irand2]/rawpt[irand2]<<'\n'
												<<" run "<<*run<<" event "<<*event<<" lumiBlock "<<*lumiBlock<<'\n'
												<<"total jets "<<(*nref)<<'\n'
												<<"jet1,2,3 pt "<<jtpt[0]<<','<<jtpt[1]<<','<<jtpt[2]<<'\n'
												<<"==========================="
												<<endl;
										///testt test ------------------------------------

									} //jtpt[irand2]<HLTJetPtTh[ntrg+1]
								}//jetId(irand1)
							}//fabs(jty[irand2])>=ybins[i]
						}//eta loop
					}//dphi>2.7
				}//if not running on MC
				//-------- jet id end efficiency ------------------------------------------
			}// if hlt has passed 
		}////trigger loop
		//=========Filled once per event==========================================
		//-----///for NPV efficiency----------------------------------------------
		if((*nPV)==1)  hnPV[0]->Fill(1);
		if((*nPV)>=1)  hnPV[1]->Fill(1);
		//		///------------------Resolution-------------------------------------
		if(runOnMC) {
			njetUnf=0;
			for(unsigned njet(0);njet<(*nref);++njet ) {
				bool cut(true);
				//do jet id and pt cut here-----------------------
				double chf(0), nhf(0), emf(0),cemf(0),nemf(0),
							 muf(0), elf(0), phf(0);
				nhf  = neutralHardSum[njet];
				muf  = muSum[njet];
				elf  = eSum[njet];
				phf  = photonSum[njet];

				cut  =
					jtpt[njet]>20             &&
					(*nPV) >= 1               &&
					(*met)/(*sumET) < 0.3     &&              
					jetId(njet) 							&& 
					(elf < 0.9) && 
					(muf < 0.9) && 
					(nhf < 0.9) &&
					treeak7->matchgendr[njet]<0.5 && 
					treeak7->matchgendr[njet]>=0 //dR to matched genjet
					;
				if(cut) 
					for(unsigned i(0);i<nybins-1;++i) ///Cycle over eta bins 
						if(fabs(jty[njet])>=ybins[i] && fabs(jty[njet])<ybins[i+1]) {
							for(unsigned ptb(0); ptb<n1x[i]; ++ptb)
								///Filling Resolution histogram-----------------------
								if(    treeak7->matchgenpt[njet]>x1[i][ptb]
										&& treeak7->matchgenpt[njet]<x1[i][ptb+1] )  //gen pt bin
								{
									hres[i][ptb]->Fill(jtpt[njet]/treeak7->matchgenpt[njet]);
									htjec[i][ptb]->Fill(jtpt[njet]/rawpt[njet]);
									htchf[i][ptb]->Fill(chargedHardSum[njet]);
									htnhf[i][ptb]->Fill(neutralHardSum[njet]);
								}//
							///Filling Resolution histogram-----------------------
						}///eta loop
				///Fill Jet distribution for Unfolding---------------------------
				//do jet id and pt cut here-----------------------
				cut =
					jtpt[njet]>20             &&
					(*nPV) >= 1               &&
					(*met)/(*sumET) < 0.3     &&
					jetId(njet)              
					;
				if(cut)
				{
					jptUnf .push_back(jtpt[njet]);		   
					jetaUnf.push_back(jteta[njet]);		   
					jphiUnf.push_back(jtphi[njet]);		   
					jmUnf  .push_back(jtm[njet]);		   
					++njetUnf;
				}
				///for unfolding check------------------------------------------
				for(unsigned i(0);i<nybins-1;++i) ///Cycle over eta bins 
					if(fabs(jty[njet])>=ybins[i] && fabs(jty[njet])<ybins[i+1]) 
						hrecMC[i]->Fill(jtpt[njet],mcWeight);


				///Done Fill Jet distribution for Unfolding----------------------
			}//for loop on jets
			// ----------- Unbiased GenJet Pt Distribution -----------------
			ngenUnf=0;
			for (unsigned ngen=0;ngen<treeak7->ngen;++ngen) {
				TLorentzVector gvec; 
				gvec.SetPtEtaPhiM(genpt[ngen],geneta[ngen],genphi[ngen],genm[ngen]);
				for (unsigned i(0); i<nybins-1; ++i)
					if (    fabs(gvec.Rapidity())>= ybins[i] 
							&&  fabs(gvec.Rapidity())<  ybins[i+1]) {
						hgen[i]->Fill(genpt[ngen],mcWeight);
						hgenUnf[i]->Fill(genpt[ngen],mcWeight);
						hgenMC[i]->Fill(genpt[ngen],mcWeight);
					}
				///Fill GenJet distribution for Unfolding---------------------------
				if(genpt[ngen]>18) {
					gjptUnf .push_back(genpt[ngen] );		   
					gjetaUnf.push_back(geneta[ngen]);		   
					gjphiUnf.push_back(genphi[ngen]);		   
					gjyUnf  .push_back(gvec.Rapidity());		   
					++ngenUnf;
				}
				///Done Fill GenJet distribution for Unfolding----------------------
			} // loop over gen jets
		}//runOnMC and not runningon8tev
		//Fill Unfold Tree -----------------------------------------------------
		//		if(runOnMC && runningOnIncxsec) 
		//			if(ngenUnf>maxGenJets ||
		//					njetUnf>maxRecoJets)
		//				cout<<"njets filled "<<jptUnf.size()
		//					<<" max "<<maxRecoJets
		//					<<" njetUnf "<<njetUnf
		//					<<" ngjets filled "<<gjptUnf.size()
		//					<<" max  "<<maxGenJets
		//					<<" ngenUnf "<<ngenUnf
		//					<<endl;

		fillArrayFromVec(gjyUnf  ,gjetyUnf  ,maxGenJets);
		fillArrayFromVec(gjphiUnf,gjetphiUnf,maxGenJets);
		fillArrayFromVec(gjetaUnf,gjetetaUnf,maxGenJets);
		fillArrayFromVec(gjptUnf ,gjetptUnf ,maxGenJets);
		fillArrayFromVec(jptUnf  ,jetptUnf  ,maxRecoJets);
		fillArrayFromVec(jphiUnf ,jetphiUnf ,maxRecoJets);
		fillArrayFromVec(jetaUnf ,jetetaUnf ,maxRecoJets);
		fillArrayFromVec(jmUnf   ,jetmUnf   ,maxRecoJets);
		if(ngenUnf || njetUnf)
			treeUnf->Fill();
		//checking prescales
		if(false)
			for(unsigned ntrg(0);ntrg<numtrgs;++ntrg)  
				if(*(hltPass[ntrg]))
					cout<<"---------------------\n"
						<<trgnames[ntrg]<<" prescale "
						//		<<" l1 "<<*(l1Psc[ntrg])
						<<" hlt "<<*(hltPsc[ntrg])
						<<" l1*hlt "
						//			<< (*(hltPsc[ntrg]))*(*(l1Psc[ntrg]))
						<<" jentry "<<jentry
						<<endl;

		cout<<fixed<<setprecision(3)<<scientific;
		if(!(jentry%200000))
			cout<<"Processing Event "<<(float)jentry<<endl;
	}//nentry loop

	///////Write histograms to file -----------------------------------
	f_out = new TFile(outfile.c_str(),"RECREATE");
	f_out->cd();
	hnPV[0]->Write();  
	hnPV[1]->Write();
	//write histos 
	for(unsigned i(0);i<nybins-1;++i) {
		for(unsigned ntrg(0);ntrg<numtrgs;++ntrg) {
			PFJet0[ntrg][i]->Write();
			PFJet1[ntrg][i]->Write();
			hpt[ntrg][i]->Write();
			hprb[ntrg][i]->Write();
			htag[ntrg][i]->Write();
			if(ntrg==0)
			{
				hptUnf[i]->Write();
				hjecl1[i]->Write();
				hjecl2l3[i]->Write();
				hjecres[i]->Write();
				hrecMC[i]->Write();
				hgenMC[i]->Write();
			}
			heta[ntrg][i]->Write();
			hphi[ntrg][i]->Write();
			hchf[ntrg][i]->Write();
			hnhf[ntrg][i]->Write();
			hnef[ntrg][i]->Write();
			hcef[ntrg][i]->Write();
			hmuf[ntrg][i]->Write();
			helf[ntrg][i]->Write();
			hphf[ntrg][i]->Write();
			hjec[ntrg][i]->Write();
			htrk[ntrg][i]->Write();
			hmet[ntrg][i]->Write();
		}
		if (runOnMC) {
			////jet resolution
			for(unsigned ptb(0); ptb<n1x[i]; ++ptb) {
				hres[i][ptb]->Write();
				htjec[i][ptb]->Write();
				htchf[i][ptb]->Write();
				htnhf[i][ptb]->Write();
			}
			//gen spectrum 
			hgen[i]->Write();
			hgenUnf[i]->Write();
		}
		//raw and corr spectrum
		hraw[i]->Write();
		hcor[i]->Write();
	}
	f_out->Write();
	f_out->Close();
	cout<<"Wrote "<<f_out->GetName()<<endl;
	//Write Unfold file ---------
	if(runOnMC && false) {
		string unfFile(outfile.substr(0,outfile.find(".root")));
		f_unf  = new TFile((unfFile+"_Unfold.root").c_str(),"RECREATE");
		f_unf->cd();
		treeUnf->Write();
		f_unf->Write();
		f_unf->Close();
	}
	///////Write jets for unfolding to file --------------------------
	//delete histos---------------------------------------------------
	hnPV[0]->Delete();  
	hnPV[1]->Delete();
	for(unsigned i(0);i<nybins-1;++i) {
		for(unsigned ntrg(0);ntrg<numtrgs;++ntrg) {
			hpt[ntrg][i]->Delete();
			hprb[ntrg][i]->Delete();
			htag[ntrg][i]->Delete();
			if(ntrg==0)
			{
				hptUnf[i]->Delete();
				hjecl1[i]->Delete();
				hjecl2l3[i]->Delete();
				hjecres[i]->Delete();
				hrecMC[i]->Delete();
				hgenMC[i]->Delete();
			}
			PFJet0[ntrg][i]->Delete();
			PFJet1[ntrg][i]->Delete();
			heta[ntrg][i]->Delete();
			hphi[ntrg][i]->Delete();
			hchf[ntrg][i]->Delete();
			hnhf[ntrg][i]->Delete();
			hcef[ntrg][i]->Delete();
			hnef[ntrg][i]->Delete();
			hmuf[ntrg][i]->Delete();
			helf[ntrg][i]->Delete();
			hphf[ntrg][i]->Delete();
			hjec[ntrg][i]->Delete();
			htrk[ntrg][i]->Delete();
			hmet[ntrg][i]->Delete();
		}
		if (runOnMC) {
			////jet resolution
			for(unsigned ptb(0); ptb<n1x[i]; ++ptb){
				hres[i][ptb]->Delete();
				htjec[i][ptb]->Delete();
				htchf[i][ptb]->Delete();
				htnhf[i][ptb]->Delete();
			}
			//gen spectrum 
			hgen[i]->Delete();
			hgenUnf[i]->Delete();
		}
		hraw[i]->Delete();
		hcor[i]->Delete();
	}
	cout<<"Clean up Done.\n";
}//END analyze

void runJetAnalyzer::Init() {
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
	jecl1= treeak7->jecl1 ;
	jecl2l3= treeak7->jecl2l3;
	jecres= treeak7->jecres  ;
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
	//Initialize Unfold ntuple variables
	ngenUnf=-1;
	njetUnf=-1;
	for(unsigned i(0);i<maxRecoJets;++i) 
	{
		jetptUnf[i] =-1; jetphiUnf[i] =-1;
		jetetaUnf[i]=-1; jetmUnf[i]=-1;
	}
	for(unsigned i(0);i<maxGenJets;++i) 
	{
		gjetptUnf[i] =-1; gjetphiUnf[i] =-1;
		gjetetaUnf[i]=-1; gjetyUnf[i]=-1;
	}
}//end iniit

void runJetAnalyzer::fillArrayFromVec(const vector<float>& vec,
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

bool runJetAnalyzer::jetId(unsigned jetnum)
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
	bool id = constituents>1  && phf<0.99 && nhf<0.99  && ((fabs(eta)<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && chf>0 && chm>0 ) || fabs(eta)>2.4) ;
	//	bool id = constituents>1   && ((fabs(eta)<=2.4
	//				 &&chf>0 && chm>0 && cemf<0.99 ) || fabs(eta)>2.4) ;
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

bool runJetAnalyzer::jetIdTagPrb(unsigned jetnum)
{
	double nhf  = neutralHardSum[jetnum];
	double muf  = muSum[jetnum];
	double elf  = eSum[jetnum];
	bool id  = 
		muf < 0.9 && nhf < 0.9 &&  
		elf < 0.9 ;

	return	(id);
}

void runJetAnalyzer::setMCWeight()
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
