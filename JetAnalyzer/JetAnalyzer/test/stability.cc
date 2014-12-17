#define stability_cxx
#include <math.h>
#include "TProfile.h"
#include <iomanip>
void runJetAnalyzer::makeStabilityHists()
{
	gROOT->SetBatch(kTRUE);
	////Get the run numbers in a vector/array==================
	TH1D* hrun = new TH1D("hrun","hrun",2000,211000.5,213000.5);
	
	if (runningOnIncxsec) 
		   treeak7->fChain->Draw("event_runNo>>hrun");
  else hltAnalyzer->fChain->Draw("Run>>hrun");

	vector<int> runNumbers, runEvents;
	cout<<setprecision(6);
  cout << scientific;
	for(unsigned i(0);i<2000;++i)
		if(hrun->GetBinContent(i+1)
     && 211821!=(int)hrun->GetBinCenter(i+1) 
		 && 211822!=(int)hrun->GetBinCenter(i+1)
			) 
		{
			runNumbers.push_back((int)hrun->GetBinCenter(i+1));
			runEvents .push_back((int)hrun->GetBinContent(i+1));
		}
	const unsigned numRuns(runNumbers.size());
	if(true)
		for(unsigned rn(0);rn<numRuns;++rn) {
			cout<<"Run "<<runNumbers[rn]<<" Index "<<rn+1
				<<" Entries "<<runEvents[rn]
				<<endl;
		}
	///declare hists ===========================================
	TProfile* Jet1pt [numtrgs][nybins-1]; TProfile* Jet2pt [numtrgs][nybins-1];
	TProfile* Jet1chf[numtrgs][nybins-1]; TProfile* Jet2chf[numtrgs][nybins-1];
	TProfile* Jet1nhf[numtrgs][nybins-1]; TProfile* Jet2nhf[numtrgs][nybins-1];
	TProfile* Jet1mas[numtrgs][nybins-1]; TProfile* Jet2mas[numtrgs][nybins-1];
	TProfile* Jet1phi[numtrgs][nybins-1]; TProfile* Jet2phi[numtrgs][nybins-1];
	TProfile* EvtRate[numtrgs][nybins-1];
	TH1D* hj1pt [numtrgs][nybins-1]; TH1D* hj2pt [numtrgs][nybins-1];
	TH1D* hj1phi[numtrgs][nybins-1]; TH1D* hj2phi[numtrgs][nybins-1];
	TH1D* hj1chf[numtrgs][nybins-1]; TH1D* hj2chf[numtrgs][nybins-1];
	TH1D* hj1nhf[numtrgs][nybins-1]; TH1D* hj2nhf[numtrgs][nybins-1];
  TH1D* hmet  [numtrgs][nybins-1]; TH1D* hrho  [numtrgs][nybins-1];

	numStabHists=0; 
	for(unsigned ntrg(0);ntrg<numtrgs;++ntrg)
		for(unsigned i(0);i<nybins-1;++i)
		{
			sprintf(name,"hj1pt_HLT%i_eta%i",HLTJetPt[ntrg],i);
			hj1pt  [ntrg][i]= new TH1D(name,name,40,0,600);
			sprintf(name,"hj1phi_HLT%i_eta%i",HLTJetPt[ntrg],i);
			hj1phi [ntrg][i]= new TH1D(name,name,40,-M_PI,M_PI);
			sprintf(name,"hj1chf_HLT%i_eta%i",HLTJetPt[ntrg],i);
			hj1chf [ntrg][i]= new TH1D(name,name,30,0,1);
			sprintf(name,"hj1nhf_HLT%i_eta%i",HLTJetPt[ntrg],i);
			hj1nhf [ntrg][i]= new TH1D(name,name,30,0,1);
			sprintf(name,"hj2pt_HLT%i_eta%i",HLTJetPt[ntrg],i);
			hj2pt  [ntrg][i]= new TH1D(name,name,40,0,600);
			sprintf(name,"hj2phi_HLT%i_eta%i",HLTJetPt[ntrg],i);
			hj2phi [ntrg][i]= new TH1D(name,name,40,-M_PI,M_PI);
			sprintf(name,"hj2chf_HLT%i_eta%i",HLTJetPt[ntrg],i);
			hj2chf [ntrg][i]= new TH1D(name,name,30,0,1);
			sprintf(name,"hj2nhf_HLT%i_eta%i",HLTJetPt[ntrg],i);
			hj2nhf [ntrg][i]= new TH1D(name,name,30,0,1);
			sprintf(name,"hmet_HLT%i_eta%i",HLTJetPt[ntrg],i);
			hmet [ntrg][i]= new TH1D(name,name,30,0,1);
			sprintf(name,"hrho_HLT%i_eta%i",HLTJetPt[ntrg],i);
			hrho [ntrg][i]= new TH1D(name,name,30,0,30);
			hj1pt[ntrg][i] ->Sumw2(); hj2pt [ntrg][i]->Sumw2();
			hj1phi[ntrg][i]->Sumw2(); hj2phi[ntrg][i]->Sumw2();
			hj1chf[ntrg][i]->Sumw2(); hj2chf[ntrg][i]->Sumw2();
			hj1nhf[ntrg][i]->Sumw2(); hj2nhf[ntrg][i]->Sumw2();
			hmet  [ntrg][i]->Sumw2(); hrho  [ntrg][i]->Sumw2();

			if(!runOnMC) {
				sprintf(name,"Jet1pt_HLT%i_eta%i",HLTJetPt[ntrg],i);
				Jet1pt  [ntrg][i]= new TProfile(name,name,numRuns,0.5,numRuns+0.5,50,500);
				sprintf(name,"Jet1chf_HLT%i_eta%i",HLTJetPt[ntrg],i);
				Jet1chf [ntrg][i]= new TProfile(name,name,numRuns,0.5,numRuns+0.5,0,1);
				sprintf(name,"Jet1nhf_HLT%i_eta%i",HLTJetPt[ntrg],i);
				Jet1nhf [ntrg][i]= new TProfile(name,name,numRuns,0.5,numRuns+0.5,0,1);
				sprintf(name,"Jet1mas_HLT%i_eta%i",HLTJetPt[ntrg],i);
				Jet1mas [ntrg][i]= new TProfile(name,name,numRuns,0.5,numRuns+0.5,0,50);
				sprintf(name,"Jet1phi_HLT%i_eta%i",HLTJetPt[ntrg],i);
				Jet1phi [ntrg][i]= new TProfile(name,name,numRuns,0.5,numRuns+0.5,-M_PI,M_PI);
				sprintf(name,"Jet2pt_HLT%i_eta%i",HLTJetPt[ntrg],i);
				Jet2pt  [ntrg][i]= new TProfile(name,name,numRuns,0.5,numRuns+0.5,50,500);
				sprintf(name,"Jet2chf_HLT%i_eta%i",HLTJetPt[ntrg],i);
				Jet2chf [ntrg][i]= new TProfile(name,name,numRuns,0.5,numRuns+0.5,0,1);
				sprintf(name,"Jet2nhf_HLT%i_eta%i",HLTJetPt[ntrg],i);
				Jet2nhf [ntrg][i]= new TProfile(name,name,numRuns,0.5,numRuns+0.5,0,1);
				sprintf(name,"Jet2mas_HLT%i_eta%i",HLTJetPt[ntrg],i);
				Jet2mas [ntrg][i]= new TProfile(name,name,numRuns,0.5,numRuns+0.5,0,50);
				sprintf(name,"Jet2phi_HLT%i_eta%i",HLTJetPt[ntrg],i);
				Jet2phi [ntrg][i]= new TProfile(name,name,numRuns,0.5,numRuns+0.5,-M_PI,M_PI);
				sprintf(name,"EvtRate_HLT%i_eta%i",HLTJetPt[ntrg],i);
				EvtRate [ntrg][i]= new TProfile(name,name,numRuns,0.5,numRuns+0.5,0,10);
			}
		}

   unsigned nentries(treeak7->fChain->GetEntries());

//		nentries = 1000;
	cout<<"Running on stability plots "<<nentries<<endl;
	for(unsigned jentry(0);jentry<nentries;++jentry) {
		////load entry 
			treeak7->fChain->GetEntry(jentry);
//      treeak5->fChain->GetEntry(jentry);
		//find the run index-------------------------------
		unsigned runIndex(0);
		for(unsigned rn(0);rn<numRuns && !runOnMC;++rn)
			if(runNumbers[rn]==treeak7->event_runNo)
				runIndex=rn+1;
		//		 assert(runIndex);
		//found the run index------------------------------
   // Dont process events outside pt window-----------------
   if(runOnMC && (*qscale)>maxpthat) continue;
		//////------Determine prescale for the event------------
		int prescl[numtrgs];
		for(unsigned ntrg(0);ntrg<numtrgs;++ntrg)  
				prescl[ntrg] = (*(hltPsc[ntrg]))*(*(hltl1Psc[ntrg])); //psc with hltl1pass

		for(unsigned ntrg(0);ntrg<numtrgs;++ntrg)  
			if((*(hltPass[ntrg]))>0.5)     {      //hlt pass 
				unsigned jetsThisEvent[nybins-1]; //for event rate
				for(unsigned i(0);i<nybins-1;++i) 
					jetsThisEvent[i]=0;

				for(unsigned njet(0);njet<(*nref);++njet ) {
					bool cut(true);
					double chf(0), nhf(0), emf(0),cemf(0),nemf(0),
								 muf(0), elf(0), phf(0);
					nhf  = neutralHardSum[njet];
					muf  = muSum[njet];
					elf  = eSum[njet];
					//do jet id and pt cut here-----------------------
					cut  =
						jtpt[njet]>20             &&
						(*nPV) >= 1               &&
						(*met)/(*sumET) < 0.3     &&              
						(elf < 0.9) && 
						(muf < 0.9) && 
						(nhf < 0.9) &&
						jetId(njet)               
						;

					if(cut) 
						for(unsigned i(0);i<nybins-1;++i) ///Cycle over eta bins 
							if(fabs(jty[njet])>=ybins[i] && fabs(jty[njet])<ybins[i+1]) {
								++jetsThisEvent[i];
						//////------Determine event weight ---------------------
								double wt = 1.0;
								if(runOnMC) wt = mcWeight;
								else        wt = prescl[ntrg];
								///Filling Spectrum histogram-----------------------
								if(    jtpt[njet]<HLTJetPtTh[ntrg+1]  //hlt pt bin
										&& jtpt[njet]>=HLTJetPtTh[ntrg])  //hlt pt bin
								{
									if(njet==0)///met/sumet and rho
									{		
										hmet[ntrg][i]->Fill((*met)/(*sumET),wt);
										hrho[ntrg][i]->Fill((*rho),wt);
									}
									if(njet==0) { ///leading jet
										if(!runOnMC) {
											Jet1pt [ntrg][i]->Fill(runIndex,jtpt[njet],1);    
											Jet1chf[ntrg][i]->Fill(runIndex,chargedHardSum[njet],1); 
											Jet1nhf[ntrg][i]->Fill(runIndex,neutralHardSum[njet],1); 
											Jet1mas[ntrg][i]->Fill(runIndex,jtm[njet],1); 
											Jet1phi[ntrg][i]->Fill(runIndex,jtphi[njet],1); }
											hj1pt [ntrg][i]->Fill(jtpt[njet],wt); 
											hj1phi[ntrg][i]->Fill(jtphi[njet],wt);
											hj1chf[ntrg][i]->Fill(chargedHardSum[njet],wt);
											hj1nhf[ntrg][i]->Fill(neutralHardSum[njet],wt);
									}
									else if(njet==1) {///second jet
										if(!runOnMC) {
											Jet2pt [ntrg][i]->Fill(runIndex,jtpt[njet],1);    
											Jet2chf[ntrg][i]->Fill(runIndex,chargedHardSum[njet],1); 
											Jet2nhf[ntrg][i]->Fill(runIndex,neutralHardSum[njet],1); 
											Jet2mas[ntrg][i]->Fill(runIndex,jtm[njet],1); 
											Jet2phi[ntrg][i]->Fill(runIndex,jtphi[njet],1); }
											hj2pt [ntrg][i]->Fill(jtpt[njet],wt); 
											hj2phi[ntrg][i]->Fill(jtphi[njet],wt);
											hj2chf[ntrg][i]->Fill(chargedHardSum[njet],wt);
											hj2nhf[ntrg][i]->Fill(neutralHardSum[njet],wt);
									}
								}//hlt pt bin
							}///eta loop
				}//for loop on jets

				if(!runOnMC)
					for(unsigned i(0);i<nybins-1;++i) ///Fill the Rate plot
						EvtRate[ntrg][i]->Fill(runIndex,jetsThisEvent[i],1);
			}////trigger loop
		if(!(jentry%200000))
			cout<<"Processing Event "<<(float)jentry<<endl;
	}//nentry loop

	///run Index histogram
	TH2F* hrIn = new TH2F("hrIn","hrIn",2000,211000.5,213000.5,10,0.5,10.5);
	for(unsigned i(0);i<numRuns && !runOnMC;++i)
		hrIn->Fill(runNumbers[i],i);

	///OUTPUT==============================================================
	string unfFile(outfile.substr(0,outfile.find(".root")));
	TFile* fstab = new TFile((unfFile+"_stabplots.root").c_str(),"RECREATE");
	fstab->cd();
	//write histos 
	hrIn->Write();
	for(unsigned i(0);i<nybins-1;++i) 
		for(unsigned ntrg(0);ntrg<numtrgs;++ntrg) {
			hj1pt [ntrg][i]->Write();
			hj1phi[ntrg][i]->Write();
			hj1chf[ntrg][i]->Write();
			hj1nhf[ntrg][i]->Write();
			hj2pt [ntrg][i]->Write();
			hj2phi[ntrg][i]->Write();
			hj2chf[ntrg][i]->Write();
			hj2nhf[ntrg][i]->Write();
			hmet[ntrg][i]->Write();
			hrho[ntrg][i]->Write();
			if(!runOnMC) {
				Jet1pt [ntrg][i]->Write();
				Jet1chf[ntrg][i]->Write();
				Jet1nhf[ntrg][i]->Write();
				Jet1mas[ntrg][i]->Write();
				Jet1phi[ntrg][i]->Write();
				Jet2pt [ntrg][i]->Write();
				Jet2chf[ntrg][i]->Write();
				Jet2nhf[ntrg][i]->Write();
				Jet2mas[ntrg][i]->Write(); 
				Jet2phi[ntrg][i]->Write(); 
				EvtRate[ntrg][i]->Write();
			}
		}
	fstab->Write();
	fstab->Close();
	//delete histos 
	hrIn->Delete();
	hrun->Delete();
	for(unsigned i(0);i<nybins-1;++i) 
		for(unsigned ntrg(0);ntrg<numtrgs;++ntrg) {
			hj1pt [ntrg][i]->Delete();
			hj1phi[ntrg][i]->Delete();
			hj1chf[ntrg][i]->Delete();
			hj1nhf[ntrg][i]->Delete();
			hj2pt [ntrg][i]->Delete();
			hj2phi[ntrg][i]->Delete();
			hj2chf[ntrg][i]->Delete();
			hj2nhf[ntrg][i]->Delete();
			hmet[ntrg][i]  ->Delete();
			hrho[ntrg][i]  ->Delete();

			if(!runOnMC) {
				Jet1pt [ntrg][i]->Delete();
				Jet1chf[ntrg][i]->Delete();
				Jet1nhf[ntrg][i]->Delete();
				Jet1mas[ntrg][i]->Delete();
				Jet1phi[ntrg][i]->Delete();
				Jet2pt [ntrg][i]->Delete();
				Jet2chf[ntrg][i]->Delete();
				Jet2nhf[ntrg][i]->Delete();
				Jet2mas[ntrg][i]->Delete();
				Jet2phi[ntrg][i]->Delete();
				EvtRate[ntrg][i]->Delete();
			}
		}
	cout<<"Stability Plots Done\n";
}
