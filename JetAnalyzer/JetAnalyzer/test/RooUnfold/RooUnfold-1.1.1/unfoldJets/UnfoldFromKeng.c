#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include <string>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"

#include "TArrayD.h"
#include "TMatrixD.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"
#endif

void colorIt(TH1D *hMidA, int kCyan ){
	hMidA->SetMarkerColor(kCyan);
	hMidA->SetLineColor(kCyan);
	hMidA->SetMarkerStyle(22);//24
	hMidA->SetMarkerSize(1.0);
	hMidA->SetFillColor(kCyan);
	hMidA->SetFillStyle(0);
	hMidA->GetXaxis()->SetLabelSize(0.05);
	hMidA->GetXaxis()->SetTitleSize(0.06);
	hMidA->GetYaxis()->SetTitleSize(0.05);
	hMidA->GetYaxis()->SetTitleOffset(0.9);
	hMidA->GetYaxis()->SetLabelSize(0.05);
	hMidA->SetTitleFont(42, "XYZ");
	hMidA->SetLabelFont(42, "XYZ");
}

void normalizeTH2D(TH2D *h)
{
	int xbin(h->GetNbinsX()), ybin(h->GetNbinsY());
	for (int i(1); i <= ybin; i++){
		double sum(0.);
		for (int j(1); j <= xbin; j++){
			sum += h->GetBinContent(j,i);
		}
		for (int j(1); j <= xbin; j++){
			if (sum > 0){
				double val = h->GetBinContent(j,i)/sum;
				val=val*100;
				val=ceil(val);
				double nearest=val/100;
				h->SetBinContent(j,i,nearest);

			}
		}
	}
}
//==============================================================================
// Unfolding
//==============================================================================

void Unfold(
		//Njets//
		//std::string varData="n_jet_30_unfold",std::string varMC="n_jet_30_Gen_VS_reco",std::string varGen="n_jet_30_gen",char *Xtitle="N_{jet}",char *Ytitle="d#sigma/dN_{jet}",
		//1st pT jet//
		//std::string varData="jetpt_1",std::string varMC="jet1pt_Gen_VS_reco",std::string varGen="Gjetpt_1",char *Xtitle="p_{T}^{1st} [GeV]",char *Ytitle="d#sigma/dp_{T}^{1st}",
		//2nd pT jet//
		//std::string varData="jetpt_2",std::string varMC="jet2pt_Gen_VS_reco",std::string varGen="Gjetpt_2",char *Xtitle="p_{T}^{2nd} [GeV]",char *Ytitle="d#sigma/dp_{T}^{2nd}",
		//1st Y jet//
		//std::string varData="jetY_1",std::string varMC="jet1Y_Gen_VS_reco",std::string varGen="GjetY_1",char *Xtitle="#eta^{1st}",char *Ytitle="d#sigma/d#eta^{1st}",
		//2nd Y jet//
		std::string varData="jetY_2",std::string varMC="jet2Y_Gen_VS_reco",std::string varGen="GjetY_2",char *Xtitle="#eta^{2nd}",char *Ytitle="d#sigma/d#eta^{2nd}",
		///////////////////////////////////////

		std::string file1="Unfolding_063013/DoubleMu2012_Unfolding_063013.root",
		std::string file2="Unfolding_063013/DYJetsToLL_Unfolding_063013.root",
		std::string file3="Unfolding_063013/TTjets_Unfolding_063013.root",
		std::string file4="Unfolding_063013/WWJetsTo2L2Nu_Unfolding_063013.root",
		std::string file5="Unfolding_063013/WZJetsTo2L2Q_Unfolding_063013.root",
		std::string file6="Unfolding_063013/WZJetsTo3LNu_Unfolding_063013.root",
		std::string file7="Unfolding_063013/ZZJetsTo2L2Nu_Unfolding_063013.root",
		std::string file8="Unfolding_063013/ZZJetsTo4L_Unfolding_063013.root",
		std::string file9="Unfolding_063013/ZZjetTo2L2Q_Unfolding_063013.root",
		int rebin = 1,
		int plot=0
		)
{
#ifdef __CINT__
	gSystem->Load("libRooUnfold");
#endif

	gStyle->SetOptStat(kFALSE);
	gStyle->SetOptTitle(kFALSE);
	//gROOT->SetBatch(kTRUE);

	cout << "==================================== TRAIN ====================================" << endl;

	TFile *f1   = TFile::Open(file1.data());
	TFile *f2   = TFile::Open(file2.data());
	TFile *f3   = TFile::Open(file3.data());
	TFile *f4   = TFile::Open(file4.data());
	TFile *f5   = TFile::Open(file5.data());
	TFile *f6   = TFile::Open(file6.data());
	TFile *f7   = TFile::Open(file7.data());
	TFile *f8   = TFile::Open(file8.data());
	TFile *f9   = TFile::Open(file9.data());

	TH1D *hDat_withBG     =(TH1D*)(f1->Get(varData.data()));

	TH1D *hRec      =(TH1D*)(f2->Get(varData.data()));
	TH1D *hGen      =(TH1D*)(f2->Get(varGen.data()));
	TH2D *hMat      =(TH2D*)(f2->Get(varMC.data()));

	TH1D *ttRec       =(TH1D*)(f3->Get(varData.data()));
	TH1D *wwjLNuRec   =(TH1D*)(f4->Get(varData.data()));
	TH1D *wzjLQRec    =(TH1D*)(f5->Get(varData.data()));
	TH1D *wzjLNuRec   =(TH1D*)(f6->Get(varData.data()));
	TH1D *zzjLNuRec   =(TH1D*)(f7->Get(varData.data()));
	TH1D *zzjLRec     =(TH1D*)(f8->Get(varData.data()));
	TH1D *zzjLQRec    =(TH1D*)(f9->Get(varData.data()));

	Double_t tt_scale= 234*19696/13847500.;
	Double_t wwjLNu_scale= 54.838*0.10608*19696/3866470.;
	Double_t wzjLQ_scale= 33.21*0.068258*19696/3207951.;
	Double_t wzjLNu_scale=33.21*0.032887 *19696/4035958.;
	Double_t zzjLNu_scale= 17.654*0.04039*19696/1909822.;
	Double_t zzjL_scale= 17.654*0.010196*19696/4795874.;
	Double_t zzjLQ_scale= 17.654*0.14118*19696/3873454.;

	ttRec->Scale(tt_scale);
	zzjLQRec->Scale(zzjLQ_scale);
	zzjLNuRec->Scale(zzjLNu_scale);
	zzjLRec->Scale(zzjL_scale);
	wzjLQRec->Scale(wzjLQ_scale);
	wzjLNuRec->Scale(wzjLNu_scale);
	wwjLNuRec->Scale(wwjLNu_scale);

	TH1D *mc_bg= (TH1D*)ttRec->Clone("mc_bg");
	//  mc_bg->Add(ttRec);
	mc_bg->Add(wwjLNuRec);
	mc_bg->Add(wzjLQRec);
	mc_bg->Add(wzjLNuRec);
	mc_bg->Add(zzjLNuRec);
	mc_bg->Add(zzjLQRec);
	mc_bg->Add(zzjLRec);

	hDat= (TH1D*) hDat_withBG->Clone();
	hDat->Add(mc_bg,-1.0);

	bin = hDat->GetNbinsX();
	//cout<<"n bin="<<bin<<endl;

	std::string str1 ("n_jet_30_unfold");
	std::string str2 ("jetpt_1");
	std::string str3 ("jetpt_2");
	std::string str4 ("jetY_1");
	std::string str5 ("jetY_2");
	if(str1.compare(varData)==0) plot=1;
	if(str2.compare(varData)==0) plot=2;
	if(str3.compare(varData)==0) plot=3;  
	if(str4.compare(varData)==0) plot=4;
	if(str5.compare(varData)==0) plot=5; 

	RooUnfoldResponse response (hRec, hGen, hMat); 

	cout << "==================================== UNFOLD ===================================" << endl;
	RooUnfoldBayes   unfold1 (&response, hDat, 3);    // OR
	RooUnfoldSvd     unfold2 (&response, hDat, 4);   // OR
	//RooUnfoldTUnfold unfold (&response, hDat);
	RooUnfoldBinByBin unfold3 (&response, hDat);

	TH1D* hUnf1= (TH1D*) unfold1.Hreco();
	TH1D* hUnf2= (TH1D*) unfold2.Hreco();
	TH1D* hUnf3= (TH1D*) unfold3.Hreco();

	//unfold.PrintTable (cout, hGen);

	cout << "==================================== DRAW ===================================" << endl;

	TH1D *Mont= (TH1D*)hDat->Clone("Mont");
	TH1D *Mont2= (TH1D*)hUnf1->Clone("Mont");
	TH1D *Mont3= (TH1D*)hUnf2->Clone("Mont");
	TH1D *Mont4= (TH1D*)hUnf3->Clone("Mont");

	leg1 = new TLegend(0.711567,0.653719,0.930511,0.965385);
	leg1->SetBorderSize(0);
	leg1->SetFillColor(kWhite);
	leg1->SetTextSize(0.05);
	leg1->AddEntry(hGen, "Truth","pel");
	leg1->AddEntry(hDat, "Measured","pel");
	leg1->AddEntry(hUnf1, "Unfold Bayes","pel");
	leg1->AddEntry(hUnf2, "Unfold SVD","pel");
	leg1->AddEntry(hUnf3, "Unfold Bin-by-Bin","pel");    

	leg2 = new TLegend(0.711567,0.653719,0.930511,0.965385);
	leg2->SetBorderSize(0);
	leg2->SetFillColor(kWhite);
	leg2->SetTextSize(0.05);
	leg2->AddEntry(Mont, "Measured/Truth","p");
	leg2->AddEntry(Mont2, "Unfold Bayes/Truth","p");
	leg2->AddEntry(Mont3, "Unfold SVD/Truth","p");
	leg2->AddEntry(Mont4, "Unfold Bin-by-Bin/Truth","p");

	Double_t Mc_scale= 3503.7/30383355.;
	Double_t data_scale= 1/19696.;
	hDat->Scale(data_scale);
	hUnf1->Scale(data_scale);
	hUnf2->Scale(data_scale);
	hUnf3->Scale(data_scale);
	hGen->Scale(Mc_scale);

	if(plot==2){
		double jetPt_Zinc1jet[12] = {30, 40, 52, 68, 88, 113, 144, 184, 234, 297, 377, 480};
		for(int i=1;i<=bin;i++){
			double ndata = hDat->GetBinContent(i)/(jetPt_Zinc1jet[i]-jetPt_Zinc1jet[i-1]);
			double ndataerr = hDat->GetBinError(i)/(jetPt_Zinc1jet[i]-jetPt_Zinc1jet[i-1]);
			hDat->SetBinContent(i,ndata);
			hDat->SetBinError(i,ndataerr);
			double nGen = hGen->GetBinContent(i)/(jetPt_Zinc1jet[i]-jetPt_Zinc1jet[i-1]);
			double nGenerr = hGen->GetBinError(i)/(jetPt_Zinc1jet[i]-jetPt_Zinc1jet[i-1]);
			hGen->SetBinContent(i,nGen);
			hGen->SetBinError(i,nGenerr);
			double nUnf1 = hUnf1->GetBinContent(i)/(jetPt_Zinc1jet[i]-jetPt_Zinc1jet[i-1]);
			double nUnf1err = hUnf1->GetBinError(i)/(jetPt_Zinc1jet[i]-jetPt_Zinc1jet[i-1]);
			hUnf1->SetBinContent(i,nUnf1);
			hUnf1->SetBinError(i,nUnf1err);
			double nUnf2 = hUnf2->GetBinContent(i)/(jetPt_Zinc1jet[i]-jetPt_Zinc1jet[i-1]);
			double nUnf2err = hUnf2->GetBinError(i)/(jetPt_Zinc1jet[i]-jetPt_Zinc1jet[i-1]);
			hUnf2->SetBinContent(i,nUnf2);
			hUnf2->SetBinError(i,nUnf2err);
			double nUnf3 = hUnf3->GetBinContent(i)/(jetPt_Zinc1jet[i]-jetPt_Zinc1jet[i-1]);
			double nUnf3err = hUnf3->GetBinError(i)/(jetPt_Zinc1jet[i]-jetPt_Zinc1jet[i-1]);
			hUnf3->SetBinContent(i,nUnf3);
			hUnf3->SetBinError(i,nUnf3err);
		}
	}
	if(plot==3){
		double jetPt_Zinc2jet[11] = {30, 40, 52, 68, 88, 113, 144, 184, 234, 297, 377};
		for(int i=1;i<=bin;i++){
			double ndata = hDat->GetBinContent(i)/(jetPt_Zinc2jet[i]-jetPt_Zinc2jet[i-1]);
			double ndataerr = hDat->GetBinError(i)/(jetPt_Zinc2jet[i]-jetPt_Zinc2jet[i-1]);
			hDat->SetBinContent(i,ndata);
			hDat->SetBinError(i,ndataerr);
			double nGen = hGen->GetBinContent(i)/(jetPt_Zinc2jet[i]-jetPt_Zinc2jet[i-1]);
			double nGenerr = hGen->GetBinError(i)/(jetPt_Zinc2jet[i]-jetPt_Zinc2jet[i-1]);
			hGen->SetBinContent(i,nGen);
			hGen->SetBinError(i,nGenerr);
			double nUnf1 = hUnf1->GetBinContent(i)/(jetPt_Zinc2jet[i]-jetPt_Zinc2jet[i-1]);
			double nUnf1err = hUnf1->GetBinError(i)/(jetPt_Zinc2jet[i]-jetPt_Zinc2jet[i-1]);
			hUnf1->SetBinContent(i,nUnf1);
			hUnf1->SetBinError(i,nUnf1err);
			double nUnf2 = hUnf2->GetBinContent(i)/(jetPt_Zinc2jet[i]-jetPt_Zinc2jet[i-1]);
			double nUnf2err = hUnf2->GetBinError(i)/(jetPt_Zinc2jet[i]-jetPt_Zinc2jet[i-1]);
			hUnf2->SetBinContent(i,nUnf2);
			hUnf2->SetBinError(i,nUnf2err);
			double nUnf3 = hUnf3->GetBinContent(i)/(jetPt_Zinc2jet[i]-jetPt_Zinc2jet[i-1]);
			double nUnf3err = hUnf3->GetBinError(i)/(jetPt_Zinc2jet[i]-jetPt_Zinc2jet[i-1]);
			hUnf3->SetBinContent(i,nUnf3);
			hUnf3->SetBinError(i,nUnf3err);
		}
	}
	if(plot==4 || plot==5){
		for(int i=1;i<=bin;i++){
			double ndata = hDat->GetBinContent(i)/0.2;
			double ndataerr = hDat->GetBinError(i)/0.2;
			hDat->SetBinContent(i,ndata);
			hDat->SetBinError(i,ndataerr);
			double nGen = hGen->GetBinContent(i)/0.2;
			double nGenerr = hGen->GetBinError(i)/0.2;
			hGen->SetBinContent(i,nGen);
			hGen->SetBinError(i,nGenerr);
			double nUnf1 = hUnf1->GetBinContent(i)/0.2;
			double nUnf1err = hUnf1->GetBinError(i)/0.2;
			hUnf1->SetBinContent(i,nUnf1);
			hUnf1->SetBinError(i,nUnf1err);
			double nUnf2 = hUnf2->GetBinContent(i)/0.2;
			double nUnf2err = hUnf2->GetBinError(i)/0.2;
			hUnf2->SetBinContent(i,nUnf2);
			hUnf2->SetBinError(i,nUnf2err);
			double nUnf3 = hUnf3->GetBinContent(i)/0.2;
			double nUnf3err = hUnf3->GetBinError(i)/0.2;
			hUnf3->SetBinContent(i,nUnf3);
			hUnf3->SetBinError(i,nUnf3err);
		}
	}
	Mont->Divide(hDat,hGen,1.0,1.0,"b");
	Mont2->Divide(hUnf1,hGen,1.0,1.0,"b");
	Mont3->Divide(hUnf2,hGen,1.0,1.0,"b");
	Mont4->Divide(hUnf3,hGen,1.0,1.0,"b");

	TCanvas *c1 = new TCanvas("c1","",600,600);
	c1->Divide(1,2,0.01,0);
	c1->cd(1); 
	p11_1 = (TPad*)c1->GetPad(1);
	p11_1->SetTopMargin(0.01);
	p11_1->SetBottomMargin(0);
	p11_1->SetRightMargin(0.04);
	p11_1->SetLogy();
	colorIt(hDat,kRed);
	hDat->SetYTitle(Ytitle);
	//if(plot==1) hDat->GetXaxis()->SetRangeUser(1.,7.);
	hDat->Draw("eX0C");
	colorIt(hUnf1,kGreen);
	hUnf1->Draw("sameeX0C");
	colorIt(hUnf2,kBlue);
	hUnf2->Draw("sameeX0C");
	colorIt(hUnf3,kOrange+8);
	hUnf3->Draw("sameeX0C");
	hGen->SetLineColor(kBlack);
	hGen->Draw("hhistsames");
	leg1->Draw("same");

	c1->cd(2);
	p11_2 = (TPad*)c1->GetPad(2);
	p11_2->SetRightMargin(0.04);
	p11_2->SetTopMargin(0);
	p11_2->SetBottomMargin(0.2);
	p11_2->SetGridy();
	colorIt(Mont,kRed);
	Mont->SetMinimum(0.5);
	Mont->SetMaximum(1.5);
	Mont->SetXTitle(Xtitle); //Y jet
	Mont->SetYTitle("Data/Truth");
	//if(plot==1) Mont->GetXaxis()->SetRangeUser(1.,7.);
	Mont->Draw("peX0C");// rec/gen
	colorIt(Mont2,kGreen);
	Mont2->Draw("sameeX0C");
	colorIt(Mont3,kBlue);
	Mont3->Draw("sameeX0C");
	colorIt(Mont4,kOrange+8);
	Mont4->Draw("sameeX0C");

	double TomNjet[8]={1,60.6,13.1,2.66,0.516,0.102,0.0184,0.00335};
	if(plot==1){
		cout<<"\\begin{table}[hbtp]"<<endl;
		cout<<"\\begin{center}"<<endl;
		cout<<"\\caption{Number of jets electron channel}"<<endl;
		cout<<"\\begin{tabular}{|c|c|c|} \\hline"<<endl;
		cout<<"Njet & $d \\sigma / d(NJet)$ & stat  \\\\ \\hline "<<endl;
		for(int i=1;i<9;i++){
			cout<<i-1<<" &"<<hUnf1->GetBinContent(i)<<" &"<<hUnf1->GetBinError(i)<<"\\\\"<<endl;
		}
		cout<<"\\hline"<<endl;
		cout<<"\\end{tabular}" <<endl;
		cout<<"\\label{tab:njetunf}" <<endl;
		cout<<"\\end{center}" <<endl;
		cout<<"\\end{table}" <<endl;
	}
	if(plot==2){
		cout<<"\\begin{table}[hbtp]"<<endl;
		cout<<"\\begin{center}"<<endl;
		cout<<"\\caption{Leading jet pt electron channel}"<<endl;
		cout<<"\\begin{tabular}{|c|c|c|} \\hline"<<endl;
		cout<<"jet rapidity range & $d \\sigma / dPT$ & stat \\\\ \\hline "<<endl;
		for(int i=1;i<12;i++){
			cout<<jetPt_Zinc1jet[i-1]<<"-"<<jetPt_Zinc1jet[i]<<" &"<<hUnf1->GetBinContent(i)<<" &"<<hUnf1->GetBinError(i)<<"\\\\"<<endl;
		}
		cout<<"\\hline"<<endl;
		cout<<"\\end{tabular}" <<endl;
		cout<<"\\label{tab:jet1PTunf}" <<endl;
		cout<<"\\end{center}" <<endl;
		cout<<"\\end{table}" <<endl;
	}
	if(plot==3){
		cout<<"\\begin{table}[hbtp]"<<endl;
		cout<<"\\begin{center}"<<endl;
		cout<<"\\caption{Second leading jet pt electron channel}"<<endl;
		cout<<"\\begin{tabular}{|c|c|c|} \\hline"<<endl;
		cout<<"jet rapidity range & $d \\sigma / dPT$ & stat \\\\ \\hline "<<endl;
		for(int i=1;i<11;i++){
			cout<<jetPt_Zinc2jet[i-1]<<"-"<<jetPt_Zinc2jet[i]<<" &"<<hUnf1->GetBinContent(i)<<" &"<<hUnf1->GetBinError(i)<<"\\\\"<<endl;
		}
		cout<<"\\hline"<<endl;
		cout<<"\\end{tabular}" <<endl;
		cout<<"\\label{tab:jet2PTunf}" <<endl;
		cout<<"\\end{center}" <<endl;
		cout<<"\\end{table}" <<endl;

	}
	if(plot==4 || plot==5){
		cout<<"\\begin{table}[hbtp]"<<endl;
		cout<<"\\begin{center}"<<endl;
		if(plot==4)  cout<<"\\caption{Leading jet Y electron channel}"<<endl;
		if(plot==5)  cout<<"\\caption{Second leading jet Y electron channel}"<<endl;
		cout<<"\\begin{tabular}{|c|c|c|} \\hline"<<endl;
		cout<<"jet rapidity range & $d \\sigma / dY$ & stat \\\\ \\hline "<<endl;

		double low=0;
		double high=0;
		for(int i=1;i<25;i++){
			if(i<13){
				low=-1*(2.6-(i*0.2));
				high=low+0.2;
				if(i==12) high=0;
			}
			if(i>12){
				low=((i-13)*0.2);
				high=low+0.2;
			}
			cout<<low<<" - "<<high<<" & "<<hUnf1->GetBinContent(i)<<" &"<<hUnf1->GetBinError(i)<<"\\\\"<<endl;
		}
		cout<<"\\hline"<<endl;
		cout<<"\\end{tabular}" <<endl;
		if(plot==4)    cout<<"\\label{tab:jet1Yunf}" <<endl;
		if(plot==5)    cout<<"\\label{tab:jet2Yunf}" <<endl;
		cout<<"\\end{center}" <<endl;
		cout<<"\\end{table}" <<endl;
	}
	TCanvas *c2 = new TCanvas("c2","",600,600);
	char name[100];
	normalizeTH2D(hMat);
	sprintf(name,"Reco %s",Xtitle);
	hMat->SetXTitle(name);
	sprintf(name,"Gen %s",Xtitle);
	hMat->GetYaxis()->SetTitleOffset(1.5);
	hMat->SetYTitle(name);
	//if(plot==1) hMat->GetXaxis()->SetRangeUser(1.,7.);
	//if(plot==1) hMat->GetYaxis()->SetRangeUser(1.,7.); 
	//gPad->SetLogz();
	hMat->Draw("COLZTEXT");

	sprintf(name,"unfold_BGsub_070313/%s.root",varData.c_str());
	c1->Print(name);
	sprintf(name,"unfold_BGsub_070313/%s.pdf",varData.c_str());
	c1->Print(name);
	sprintf(name,"unfold_BGsub_070313/%s.png",varData.c_str());
	c1->Print(name);

	sprintf(name,"unfold_BGsub_070313/%s_res.root",varData.c_str());
	c2->Print(name);
	sprintf(name,"unfold_BGsub_070313/%s_res.pdf",varData.c_str());
	c2->Print(name);
	sprintf(name,"unfold_BGsub_070313/%s_res.png",varData.c_str());
	c2->Print(name);

}

#ifndef __CINT__
int main () { RooUnfoldExample(); return 0; }  // Main program when run stand-alone
#endif


