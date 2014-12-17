#include <TFile.h>
#include <TLatex.h>
#include <TProfile.h>
#include <THStack.h>
#include <TH2.h>
#include <TLine.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <iostream>
#include "../../../AnalyzerProcedures.h"
#include "../../../Plots/plot_procedures.h"
using namespace std;
void groomSpectrum(TH1* h, unsigned neta);
TH1D* makeNLOSpectrum(char* infile,char* title);
void groomSpecClosure(TH1* h, unsigned eta);
void drawClosure(TCanvas* c,TH1* hdt, TH1* hmc,TH1* hclo);
void spFillHisto(TH1D* hin, TH1D* hout);
void spFillHistoWithErr(TH1D* hin, TH1D* hout);
void applyNP(TH1* h,TF1* f);
void expuncDraw(TH1*  h, int uncColor=kRed+1);
void thyuncDraw(TH1*  h);
void pdfsDraw(TH1*  h,int col);
void datatheoryplots(TH1* h, unsigned eta);
TH1D* makePtRangeHisto(TH1D* h,unsigned eta);
//=============================================================================
//DRAW NLO AND RECO PT SPECTRUM------------------------------------------
void draw_datatheory_pdfs()
{
	const bool saveplots(true);
	if(saveplots)   gROOT->SetBatch(kTRUE); ///display mode off----------
	const unsigned nybinsInfo(6);
	const char* infile("UnfoldJetsData.root");
	const char* nloDir("/uscms_data/d3/tlibeiro/inclusivexsec/Generators/generateIncJets/inputXsecs");
	TFile* npfile= new TFile("npCorrections.root","read");
	TFile* thyfile= new TFile("TheoryUncertainties.root","read");
	TFile* expfile= new TFile("ExpUncertainties.root","read");
	///////////////////////////////////////////////////
	char* nloSpec[] = {"abm11","nnpdf2.1","ct10","hera1.5","mstw08"};
	char* nloSpecout[] = {"abm11","nnpdf21","ct10","hera15","mstw08"};
	const unsigned numNLO(sizeof(nloSpec)/sizeof(string));
	char* nloLabels[numNLO] = {"ABM 11","NNPDF 2.1","CT10","HERA1.5","MSTW2008"};
  unsigned markerSty[8]={20,24,21,25,22,26,23,27};
	unsigned markerCol[8]={kRed,kMagenta,kBlue+1,kBlue,
		kGreen,kOrange+7,kRed-4,kGreen-7};

	TFile* fileUnfold =new TFile(infile,"read");
  TH1D* hpt [nybinsInfo];//Unfolded data
  TH1D* hnlo[nybinsInfo][numNLO];//compare nlo spectrum
  TH1D* hthy[nybinsInfo][numNLO][2];//theory uncertainties
  TH1D* hexp[nybinsInfo][2];        //experimetnal uncertainties
  TH1D* hclo[nybinsInfo][numNLO];//dtmc ratio
  for(unsigned neta(0);neta<nybinsInfo;++neta) {
    //unfolded spectrum and uncertainties
		sprintf(name,"hUnf_eta%d",neta);
		hpt[neta] = (TH1D*)fileUnfold->Get(name);
    sprintf(name,"htotal_up_eta%d",neta);
    hexp[neta][1] = (TH1D*)expfile->Get(name);
    sprintf(name,"htotal_dn_eta%d",neta);
    hexp[neta][0] = (TH1D*)expfile->Get(name);
		//get nlo plots
		for(unsigned nspec(0);nspec<numNLO;++nspec) {
      //theory uncertainties
     sprintf(name,"htotalUp_eta%d_%s",neta,nloSpec[nspec]);
     hthy[neta][nspec][1]=(TH1D*)thyfile->Get(name);
     sprintf(name,"htotalDn_eta%d_%s",neta,nloSpec[nspec]);
     hthy[neta][nspec][0]=(TH1D*)thyfile->Get(name);
      //
			sprintf(name,"hclosure_%s_eta%d",nloSpec[nspec],neta);
			hclo[neta][nspec] = (TH1D*)hpt[neta]->Clone(name);
			hclo[neta][nspec]->Reset();
			char infile[500]; char title[500]; 
			sprintf(infile,"%s/%sy%i.cppInput",nloDir,nloSpec[nspec],neta);
			sprintf(title ,"nloy%i%s",neta,nloSpec[nspec]);
			hnlo[neta][nspec] = makeNLOSpectrum(infile,title);
			sprintf(name,"fratio_eta%d",neta);
			TF1* f = (TF1*)npfile->Get(name);
			applyNP(hnlo[neta][nspec],f);
		}
	}//eta
  cout<<"Finished getting histograms\n";
	////////////////////////////////////////////////////
	///draw data/theory histogram===============================
	TCanvas* cratio[numNLO][6];
	TCanvas* ccombpdf[numNLO];
	for(unsigned neta(0);neta<nybinsInfo;++neta) {
		for(unsigned nspec(0);nspec<numNLO ;++nspec) {
			//closure
			TH1D* htmp = (TH1D*)hpt[neta]->Clone("htmp");
			spFillHisto(hnlo[neta][nspec],htmp);
			hclo[neta][nspec]->Divide(hpt[neta],htmp);
			datatheoryplots(hclo[neta][nspec],neta);
			//draw 
			sprintf(name,"cratio_eta%d_%s",neta,nloSpecout[nspec]);
			cratio[nspec][neta]=new TCanvas(name,name,600,600);
      cratio[nspec][neta]->SetLogx();
      cratio[nspec][neta]->SetTickx(); cratio[nspec][neta]->SetTicky();
			hclo[neta][nspec]->Draw();//data theory ratio
			if(nspec==0) {
				expuncDraw(hexp[neta][0]);//experimentla uncertainty
				expuncDraw(hexp[neta][1]);
			}
			hexp[neta][0]->Draw("same");
			hexp[neta][1]->Draw("same");
      thyuncDraw(hthy[neta][nspec][0]);//theory uncertainty
      thyuncDraw(hthy[neta][nspec][1]);
      hthy[neta][nspec][0]->Draw("same");
      hthy[neta][nspec][1]->Draw("same");
      hclo[neta][nspec]->Draw("axissame");
    }//pdfs
	}//eta

	///draw compare pdfs histogram============================
  ///draw with hclo[][] made in data/theory part
	unsigned lineCol[8]={kRed,kMagenta,kCyan,kBlue,
		kGreen,kOrange+7,kRed-4,kGreen-7};
	TCanvas* cpdfs[numNLO][6];
  TH1D*    hpdf[6][numNLO][numNLO];
	for(unsigned neta(0);neta<nybinsInfo;++neta) {
		for(unsigned nspec(0);nspec<numNLO;++nspec) {
			//draw
			sprintf(name,"cpdfs_eta%d_%s",neta,nloSpecout[nspec]);
			cpdfs[nspec][neta]=new TCanvas(name,name,600,600);
			cpdfs[nspec][neta]->SetLogx();
			cpdfs[nspec][neta]->SetTickx(); cpdfs[nspec][neta]->SetTicky();
			hclo[neta][nspec]->Draw();//data theory ratio
       hexp[neta][0]->SetLineColor(kBlack);
       hexp[neta][1]->SetLineColor(kBlack);
			hexp[neta][0]->Draw("same");//experimental uncertainty
			hexp[neta][1]->Draw("same");
			//compare pdfs
			for(unsigned nsp(0);nsp<numNLO;++nsp) 
				if(nsp!=nspec)
				{
					sprintf(name,"cpdfs_eta%d_%s_ov_%s",neta,nloSpec[nsp],nloSpecout[nspec]);
					TH1D* htmp = (TH1D*)hnlo[neta][nsp]->Clone(name);
          htmp->Divide(hnlo[neta][nspec]);
          pdfsDraw(htmp,lineCol[nsp]);
          hpdf[neta][nspec][nsp] = htmp;
          hpdf[neta][nspec][nsp]->Draw("same");
				}
			hclo[neta][nspec]->Draw("sameaxis");
		}//pdfs
	}//neta

	///draw legends
	TLegend* leg = new TLegend(0.2,0.6,0.4,0.85); // xy coordinates
	leg->SetFillColor(0);  leg->SetTextFont(42);
	leg->SetBorderSize(0); leg->SetTextSize(0.04);
	//latex and legend
	for(unsigned neta(0);neta<nybinsInfo;++neta) {
		for(unsigned nspec(0);nspec<numNLO;++nspec) {
			cratio[nspec][neta]->cd();
			char lumilabel[100], etalabel[100];
			TLatex ltxclo1;ltxclo1.SetTextSize(0.04);ltxclo1.SetTextFont(42);
			TLatex ltxclo2;ltxclo2.SetTextSize(0.04);ltxclo2.SetTextFont(42);
//			sprintf(lumilabel,"#sqrt{s} = 2.76 TeV%*s L = %.2f pb^{-1}%*s CMS 2013",10,"",lumi,10,"");
		  sprintf(lumilabel,"%*s %.2f pb^{-1} (2.76 TeV)",50,"",lumi);
      ltxclo1.SetNDC();
			ltxclo1.DrawLatex(0.1,0.91,lumilabel);
			sprintf(etalabel,"%s",nloLabels[nspec]);
			ltxclo2.DrawTextNDC(0.2,0.30,etalabel);
			sprintf(etalabel,"%.1f < |y| < %.1f",ybins[neta],ybins[neta+1]);
			ltxclo2.DrawTextNDC(0.2,0.25,etalabel);
		    //CMS Logo
    extraCmsText="Preliminary";
		TLatex ltxclo3;ltxclo3.SetTextSize(0.04);ltxclo3.SetTextFont(42);
		ltxclo3.SetTextSize(0.04); ltxclo3.SetTextFont(cmsTextFont);
    ltxclo3.SetNDC();
		ltxclo3.DrawTextNDC(0.1,0.91,cmsText);
		ltxclo3.SetTextSize(0.04*0.76); ltxclo3.SetTextFont(extraTextFont);
		ltxclo3.DrawTextNDC(0.2,0.91,extraCmsText);
			////legend
			if(neta==0 && nspec==0)
			{
				leg->AddEntry(hclo[neta][nspec],"Data/Theory","pl");
				leg->AddEntry(hexp[neta][0]    ,"Exp. Uncertainty ","l");
				leg->AddEntry(hthy[neta][nspec][0],"Theo. Uncertainty","l");
			}
			leg->Draw("same");
			/////for pdfs
			cpdfs[nspec][neta]->cd();
      ltxclo1.SetNDC();
			ltxclo1.DrawLatex(0.1,0.91,lumilabel);
      sprintf(etalabel,"%s",nloLabels[nspec]);
			ltxclo2.DrawTextNDC(0.2,0.30,etalabel);
			sprintf(etalabel,"%.1f < |y| < %.1f",ybins[neta],ybins[neta+1]);
			ltxclo2.DrawTextNDC(0.2,0.25,etalabel);
    ltxclo3.SetNDC();
		ltxclo3.SetTextSize(0.04); ltxclo3.SetTextFont(cmsTextFont);
		ltxclo3.DrawTextNDC(0.1,0.91,cmsText);
		ltxclo3.SetTextSize(0.04*0.76); ltxclo3.SetTextFont(extraTextFont);
		ltxclo3.DrawTextNDC(0.2,0.91,extraCmsText);
			TLegend* lpdf = new TLegend(0.2,0.6,0.4,0.85); // xy coordinates
			lpdf->SetFillColor(0);  lpdf->SetTextFont(42);
			lpdf->SetBorderSize(0); lpdf->SetTextSize(0.03);
			lpdf->AddEntry(hclo[neta][nspec],"Data/Theory","pl");
			lpdf->AddEntry(hexp[neta][0]    ,"Exp. Uncertainty","l");
			for(unsigned nsp(0);nsp<numNLO;++nsp)
				if(nsp!=nspec)
				{
//					sprintf(name,"%s_over_%s",nloLabels[nsp],nloLabels[nspec]);
					sprintf(name,"%s",nloLabels[nsp]);
					lpdf->AddEntry( hpdf[neta][nspec][nsp],name,"l");
				}//loop over pdfs
			lpdf->Draw("same");
			//draw all eta bins together 
			if(neta==0)
			{
				sprintf(name,"ccombpdf_%s",nloSpecout[nspec]);
				ccombpdf[nspec] = new TCanvas(name,name,1000,700);
				ccombpdf[nspec]->Divide(3,2);
			}
			ccombpdf[nspec]->cd(neta+1);
			cratio[nspec][neta]->DrawClonePad();
		}//pdfs
	}//eta

///save plots
	if(saveplots)
		for(unsigned neta(0);neta<nybinsInfo;++neta) 
			for(unsigned nspec(0);nspec<numNLO;++nspec) {
				sprintf(name,"plots/%s.eps",cratio[nspec][neta]->GetName());
//				cratio[nspec][neta]->SaveAs(name);
				sprintf(name,"plots/%s.eps",cpdfs[nspec][neta]->GetName());
				cpdfs[nspec][neta]->SaveAs(name);
				if(neta==0)
				{
					sprintf(name,"plots/%s.eps",ccombpdf[nspec]->GetName());
					ccombpdf[nspec]->SaveAs(name);
				}
			}

	TFile* fout = new TFile("DataTheoryRatio.root","recreate");
	fout->cd();
	//write to file
	for(unsigned neta(0);neta<nybinsInfo;++neta) 
	{
		hexp[neta][0]->Write();
		hexp[neta][1]->Write();
		for(unsigned nspec(0);nspec<numNLO;++nspec) {
			hthy[neta][nspec][0]->Write();
			hthy[neta][nspec][1]->Write();
			hclo[neta][nspec]->Write();
			hnlo[neta][nspec]->Write();
		}
	}
	fout->Write();
	fout->Close();
}

void draw_nlopt_spectrum()
{
	const bool saveplots(true);
	if(saveplots)   gROOT->SetBatch(kTRUE); ///display mode off----------
	const unsigned nybinsInfo(7);
	const string infile("UnfoldJetsData.root");
	const char* nloDir("/uscms_data/d3/tlibeiro/inclusivexsec/Generators/generateIncJets/inputXsecs");
	TFile* npfile= new TFile("npCorrections.root","read");
	///////////////////////////////////////////////////
	char* nloSpec[] = {"abm11","nnpdf2.1","ct10","hera1.5","mstw08"};
	char* nloSpecout[] = {"abm11","nnpdf21","ct10","hera15","mstw08"};
	const unsigned numNLO(sizeof(nloSpec)/sizeof(string));
	char* nloLabels[numNLO] = {"ABM 11","NNPDF 2.1","CT10","HERA1.5","MSTW2008"};
	unsigned markerSty[8]={20,24,21,25,22,26,23,27};
	unsigned markerCol[8]={kRed,kMagenta,kBlue+1,kBlue,
		kGreen,kOrange+7,kRed-4,kGreen-7};

	TFile* fileUnfold;
	fileUnfold=new TFile((infile).c_str(),"read");
	TH1D* hpt [nybinsInfo-1];//Unfolded data 
	TH1D* hnlo[nybinsInfo-1][numNLO];//compare nlo spectrum
	TH1D* hclo[nybinsInfo-1][numNLO];//data mc closure
	for(unsigned neta(0);neta<nybinsInfo-1;++neta) {
		sprintf(name,"hUnf_eta%d",neta);
//		hpt[neta] = (TH1D*)fileUnfold->Get(name); //full range
    hpt[neta] = makePtRangeHisto((TH1D*)fileUnfold->Get(name),neta);
		//get nlo plots
		for(unsigned nspec(0);nspec<numNLO;++nspec) {
			sprintf(name,"hclosure_eta%.2f",ybins[neta+1]);
			hclo[neta][nspec] = (TH1D*)hpt[neta]->Clone(name);
			hclo[neta][nspec]->Reset();
			char infile[500]; char title[500]; 
			sprintf(infile,"%s/%sy%i.cppInput",nloDir,nloSpec[nspec],neta);
			sprintf(title ,"nloy%i%s",neta,nloSpec[nspec]);
			hnlo[neta][nspec] = makePtRangeHisto(makeNLOSpectrum(infile,title),neta);
			sprintf(name,"fratio_eta%d",neta);
			TF1* f = (TF1*)npfile->Get(name);
			applyNP(hnlo[neta][nspec],f);
		}
	}
	//draw histograms------------------------------------------
	TCanvas* ptsp[numNLO]; // combined spectrum 
	TCanvas* ptslice[nybinsInfo-1][numNLO]; // ratio data-theory
	TCanvas* cloComb[nybinsInfo-1];//closure combined
	gStyle->SetErrorX(0);
	//////draw legends
	TLegend* leg = new TLegend(0.55,0.65,0.85,0.88); // xy coordinates
	leg->SetFillColor(0);  leg->SetTextFont(42);
	leg->SetBorderSize(0); leg->SetTextSize(0.04);

	//	////draw data/mc histograms
	for(unsigned nspec(0);nspec<numNLO;++nspec) {
		sprintf(name,"c%s",nloSpecout[nspec]);
		ptsp[nspec] = new TCanvas(name,name,900,600);
		ptsp[nspec]->SetLogx(); ptsp[nspec]->SetLogy();
		ptsp[nspec]->SetTickx(); ptsp[nspec]->SetTicky();
		for(unsigned neta(0);neta<nybinsInfo-1-3;++neta) {
			if(nspec==0)
				groomSpectrum(hpt [neta],neta);
			groomSpectrum(hnlo[neta][nspec],neta);
			hnlo[neta][nspec]->SetLineColor(kRed);
//			hnlo[neta][nspec]->SetLineWidth(2);
			ptsp[nspec]->cd();
			if(neta==0) 			hpt [neta]->Draw("p");
			else              hpt [neta]->Draw("p same");
			hnlo[neta][nspec]->Draw("hist same");
			hpt [neta]->Draw("axissame");
			hpt [neta]->Draw("p same");
			if(nspec==0) {
				sprintf(name,"%.1f < |y|< %.1f (#times 10^{%.0f})",ybins[neta],ybins[neta+1],log10(yScales[neta]));
				leg->AddEntry(hpt[neta],name,"p");
			}
			//closure
			TH1D* htmp = (TH1D*)hpt[neta]->Clone("htmp");
			spFillHisto(hnlo[neta][nspec],htmp);
			hclo[neta][nspec]->Divide(hpt[neta],htmp);
			groomSpecClosure(hclo[neta][nspec],neta);
			//comment below to stop draw
			if(nloSpec[nspec]=="ct10")
			{
				sprintf(name,"data_theory_%s_eta%d",nloSpec[nspec],neta);
				ptslice[neta][nspec] = new TCanvas(name,name,600,600);
				ptslice[neta][nspec]->SetLogx();
				ptslice[neta][nspec]->SetTickx();
				ptslice[neta][nspec]->cd();
				hclo[neta][nspec]->GetYaxis()->SetRangeUser(0,2.);
				hclo[neta][nspec]->SetMarkerStyle(markerSty[nspec]);
				hclo[neta][nspec]->SetMarkerColor(markerCol[nspec]);
				hclo[neta][nspec]->SetLineColor(markerCol[nspec]);
				hclo[neta][nspec]->Draw(); //print latex
				TLatex ltxclo1;ltxclo1.SetTextSize(0.04);ltxclo1.SetTextFont(42);
				TLatex ltxclo2;ltxclo2.SetTextSize(0.04);ltxclo2.SetTextFont(42);
				sprintf(name,"#sqrt{s} = 2.76 TeV%*s   L = %.2f pb^{-1}%*s   CMS 2013",7,"",lumi,7,"");
				ltxclo1.DrawLatex(74,2.01,name);
				sprintf(name,"%s, %.1f<|y|<%.1f",nloLabels[nspec],ybins[neta],ybins[neta+1]);
				ltxclo2.DrawLatex(100,1.5,name);//draw line at ration 1
				TF1* fline = new TF1("fline","[0]",74,maxPt[neta]-1);
				fline->SetParameter(0,1);
				fline->SetLineColor(kBlack);
				fline->Draw("same");//////save the plot
				sprintf(name,"plots/%s.eps",ptslice[neta][nspec]->GetName());
				hclo[neta][nspec]->Draw("same"); //print latex
				if(saveplots)
					ptslice[neta][nspec]->SaveAs(name);
			}
		}
		ptsp[nspec]->cd();
    ptsp[nspec]->SetLeftMargin(0.12);
		leg->Draw("same");
		///print eta bin
//		sprintf(name,"#sqrt{s} = 2.76 TeV%*s   L = %.2f pb^{-1}%*s   CMS 2013",20,"",lumi,20,"");
		sprintf(name,"%*s %.2f pb^{-1} (2.76 TeV)",60,"",lumi);
		TLatex latex1;latex1.SetTextSize(0.05);latex1.SetTextFont(42);
		TLatex latex2;latex2.SetTextSize(0.04);latex2.SetTextFont(42);
		latex1.DrawLatex(74,2e11,name);
		sprintf(name,"%s NLO #otimes NP",nloLabels[nspec]);
		latex2.DrawLatex(80,2e-1,name);
		    //CMS Logo
    extraCmsText="Preliminary";
		latex1.SetTextSize(0.05); latex1.SetTextFont(cmsTextFont);
		latex1.DrawTextNDC(0.12,0.91,cmsText);
		latex1.SetTextSize(0.04); latex1.SetTextFont(extraTextFont);
		latex1.DrawTextNDC(0.22,0.91,extraCmsText);


    ///////save the plot
		sprintf(name,"plots/%s.eps",ptsp[nspec]->GetName());
		if(saveplots)
			ptsp[nspec]->SaveAs(name);

	}
	////combine closure for all spectrum	
	TLegend* legcl = new TLegend(0.6,0.6,0.8,0.85); // xy coordinates
	legcl->SetFillColor(0);  legcl->SetTextFont(62);
	legcl->SetBorderSize(0); legcl->SetTextSize(0.025);
	for(unsigned neta(0);neta<nybinsInfo-1;++neta) {
		sprintf(name,"data_theory_allTheory_eta%d",neta);
		cloComb[neta] = new TCanvas(name,name,600,600);
		cloComb[neta]->SetLogx();
		cloComb[neta]->SetTickx();
		cloComb[neta]->cd();
		for(unsigned nspec(0);nspec<numNLO;++nspec) {
			hclo[neta][nspec]->SetMarkerStyle(markerSty[nspec]);
			hclo[neta][nspec]->SetMarkerColor(markerCol[nspec]);
			hclo[neta][nspec]->SetLineColor(markerCol[nspec]);
			if(nspec==0)		hclo[neta][nspec]->Draw();
			else          hclo[neta][nspec]->Draw("same");
			//make legend
			if(neta==0) legcl->AddEntry(hclo[neta][nspec],nloLabels[nspec],"lp");
		}//spectrum loop
		TLatex ltxclo1;ltxclo1.SetTextSize(0.04);ltxclo1.SetTextFont(42);
		TLatex ltxclo2;ltxclo2.SetTextSize(0.04);ltxclo2.SetTextFont(42);
		sprintf(name,"#sqrt{s} = 2.76 TeV%*s   L = %.2f pb^{-1}%*s   CMS 2013",7,"",lumi,7,"");
		ltxclo1.DrawLatex(74,2.51,name);
		sprintf(name,"%.1f<|y|<%.1f",ybins[neta],ybins[neta+1]);
		ltxclo2.DrawTextNDC(0.4,0.8,name);
		legcl->Draw("same");//////////draw a line at ratio 1/////////////
		TF1* fline = new TF1("fline","[0]",74,maxPt[neta]-1);
		fline->SetParameter(0,1);
		fline->SetLineColor(kBlack);
		fline->Draw("same");///////save the plot
		sprintf(name,"plots/%s.eps",cloComb[neta]->GetName());
		if(saveplots)
			cloComb[neta]->SaveAs(name);
	}//eta loop
	////end combine closure for all spectrum


	//	/////drawing individual eta slices
	//	for(unsigned neta(0);neta<nybinsInfo-1;++neta) {
	//		sprintf(name,"ptmcdtslice_eta%d",neta);
	//		ptmcdtslice[0][neta] = new TCanvas(name,name,800,600);
	//		ptmcdtslice[0][neta]->SetLogx(1);
	//		ptmcdtslice[0][neta]->SetLogy(1);
	//		TH1D* hdt = (TH1D*)hptc[1][neta][0]->Clone("hdt");
	//		TH1D* hmc = (TH1D*)hptc[1][neta][1]->Clone("hmc");
	//		TH1D* hnl = (TH1D*)hnlo[neta]->Clone("hnlo");
	//		hmc->SetMarkerStyle(1); hmc->SetLineColor(kMagenta);
	//		hmc->SetMarkerSize(0.01);
	//		hnl->SetLineColor(kRed+1);
	//		drawClosure(ptmcdtslice[0][neta],hptc[1][neta][0],hnlo[neta],hclo[neta]);
	//		hdt->Draw("p");
	//		//	hmc->Draw("hist e same");
	//		hnl->Draw("hist same");
	//		hdt->Draw("p same");
	//		///print eta bin
	//		sprintf(name,"%.2f<|y|<%.2f",ybins[neta],ybins[neta+1]);
	//		TLatex latex;latex.SetTextSize(0.05);latex.SetTextFont(42);
	//		latex.DrawTextNDC(.5,.7,name);
	//
	//		///save to file
	//		sprintf(name,"plots/ptmcdtslice_eta%d.eps",neta);
	//		if(saveplots)
	//			ptmcdtslice[0][neta]->SaveAs(name);
	//	}
	//

	//	for(unsigned neta(0);neta<nybinsInfo-1;++neta) 
	//	{
	//		ostringstream legstr, scale;
	//		legstr<<fixed;
	//		legstr<<setprecision(1);
	//		legstr<<ybins[neta]<<"<|y|<";
	//		legstr<<ybins[neta+1];
	//		legstr<<scientific;
	//		legstr<<setprecision(1);
	//		legstr<<" (#times "<<1.0*yScales[neta]<<")";
	//	}
	//	ptmcdt->cd();
	//	leg->Draw("P same");
	//	if(saveplots) {
	//		sprintf(name,"plots/ptmcdtComdined.eps");
	//		ptmcdt->SaveAs(name);
	//	}
	//
}

void draw_datatheory_pdfs_testjec()
{
	const bool saveplots(true);
	if(saveplots)   gROOT->SetBatch(kTRUE); ///display mode off----------
	const unsigned nybinsInfo(6);
	const char* infile("UnfoldJetsData.root");
	const char* nloDir("/uscms_data/d3/tlibeiro/inclusivexsec/Generators/generateIncJets/inputXsecs");
	TFile* npfile= new TFile("npCorrections.root","read");
	TFile* thyfile= new TFile("TheoryUncertainties.root","read");
	TFile* expfile= new TFile("ExpUncertainties.root","read");
	///////////////////////////////////////////////////
	char* nloSpec[] = {"abm11","nnpdf2.1","ct10","hera1.5","mstw08"};
	const unsigned numNLO(sizeof(nloSpec)/sizeof(string));
	char* nloLabels[numNLO] = {"ABM 11","NNPDF 2.1","CT10","HERA1.5","MSTW2008"};
  unsigned markerSty[8]={20,24,21,25,22,26,23,27};
	unsigned markerCol[8]={kRed,kMagenta,kBlue+1,kBlue,
		kGreen,kOrange+7,kRed-4,kGreen-7};

	TFile* fileUnfold =new TFile(infile,"read");
  TH1D* hpt [nybinsInfo];//Unfolded data
  TH1D* hdat[nybinsInfo];//folded data
  TH1D* hnlo[nybinsInfo][numNLO];//compare nlo spectrum
  TH1D* hthy[nybinsInfo][numNLO][2];//theory uncertainties
  TH1D* hexp[nybinsInfo][2];        //experimetnal uncertainties
  TH1D* hclo[nybinsInfo][numNLO];//dtmc ratio
  TH1D* hclodat[nybinsInfo][numNLO];//dtmc ratio with folded data
  for(unsigned neta(0);neta<nybinsInfo;++neta) {
    //folded spectrum
    sprintf(name,"hDat_%d",neta);
    hdat[neta] = (TH1D*)fileUnfold->Get(name);  
    ScaleHistByBW(hdat[neta]);
    hdat[neta]->Scale(1./lumi);
    //unfolded spectrum and uncertainties
		sprintf(name,"hUnf_eta%d",neta);
		hpt[neta] = (TH1D*)fileUnfold->Get(name);
    sprintf(name,"htotal_up_eta%d",neta);
    hexp[neta][1] = (TH1D*)expfile->Get(name);
    sprintf(name,"htotal_dn_eta%d",neta);
    hexp[neta][0] = (TH1D*)expfile->Get(name);
		//get nlo plots
		for(unsigned nspec(0);nspec<numNLO;++nspec) {
      //theory uncertainties
     sprintf(name,"htotalUp_eta%d_%s",neta,nloSpec[nspec]);
     hthy[neta][nspec][1]=(TH1D*)thyfile->Get(name);
     sprintf(name,"htotalDn_eta%d_%s",neta,nloSpec[nspec]);
     hthy[neta][nspec][0]=(TH1D*)thyfile->Get(name);
      //
			sprintf(name,"hclosure_%s_eta%d",nloSpec[nspec],neta);
			hclo[neta][nspec] = (TH1D*)hpt[neta]->Clone(name);
			hclo[neta][nspec]->Reset();
      sprintf(name,"hclosuredat_%s_eta%d",nloSpec[nspec],neta);
			hclodat[neta][nspec] = (TH1D*)hdat[neta]->Clone(name);
			hclodat[neta][nspec]->Reset();
			char infile[500]; char title[500]; 
			sprintf(infile,"%s/%sy%i.cppInput",nloDir,nloSpec[nspec],neta);
			sprintf(title ,"nloy%i%s",neta,nloSpec[nspec]);
			hnlo[neta][nspec] = makeNLOSpectrum(infile,title);
			sprintf(name,"fratio_eta%d",neta);
			TF1* f = (TF1*)npfile->Get(name);
			applyNP(hnlo[neta][nspec],f);
		}
	}//eta
  cout<<"Finished getting histograms\n";
	////////////////////////////////////////////////////
	///draw data/theory histogram===============================
	TCanvas* cratio[numNLO][6];
	TCanvas* ccombpdf[numNLO];
	for(unsigned neta(0);neta<nybinsInfo;++neta) {
		for(unsigned nspec(0);nspec<numNLO ;++nspec) {
      //closure unfolded data
			TH1D* htmp = (TH1D*)hpt[neta]->Clone("htmp");
			spFillHisto(hnlo[neta][nspec],htmp);
			hclo[neta][nspec]->Divide(hpt[neta],htmp);
			datatheoryplots(hclo[neta][nspec],neta);
      //closure folded data
			TH1D* htmpdat = (TH1D*)hdat[neta]->Clone("htmpdat");
			spFillHisto(hnlo[neta][nspec],htmpdat);
			hclodat[neta][nspec]->Divide(hdat[neta],htmpdat);
			datatheoryplots(hclodat[neta][nspec],neta);
      hclodat[neta][nspec]->SetMarkerColor(kRed+1);
			//draw 
			sprintf(name,"cratio_eta%d_%s",neta,nloSpec[nspec]);
			cratio[nspec][neta]=new TCanvas(name,name,600,600);
      cratio[nspec][neta]->SetLogx();
      cratio[nspec][neta]->SetTickx(); cratio[nspec][neta]->SetTicky();
			hclo[neta][nspec]->Draw();//data theory ratio
      hclodat[neta][nspec]->Draw("same");
			if(nspec==0) {
				expuncDraw(hexp[neta][0]);//experimentla uncertainty
				expuncDraw(hexp[neta][1]);
			}
			hexp[neta][0]->Draw("same");
			hexp[neta][1]->Draw("same");
      thyuncDraw(hthy[neta][nspec][0]);//theory uncertainty
      thyuncDraw(hthy[neta][nspec][1]);
      hthy[neta][nspec][0]->Draw("same");
      hthy[neta][nspec][1]->Draw("same");
       hclo[neta][nspec]->Draw("same");
       hclodat[neta][nspec]->Draw("same");
       hclo[neta][nspec]->Draw("axissame");
    }//pdfs
	}//eta

	///draw compare pdfs histogram============================
  ///draw with hclo[][] made in data/theory part
	unsigned lineCol[8]={kRed,kMagenta,kCyan,kBlue,
		kGreen,kOrange+7,kRed-4,kGreen-7};
	TCanvas* cpdfs[numNLO][6];
  TH1D*    hpdf[6][numNLO][numNLO];
	for(unsigned neta(0);neta<nybinsInfo;++neta) {
		for(unsigned nspec(0);nspec<numNLO;++nspec) {
			//draw
			sprintf(name,"cpdfs_eta%d_%s",neta,nloSpec[nspec]);
			cpdfs[nspec][neta]=new TCanvas(name,name,600,600);
			cpdfs[nspec][neta]->SetLogx();
			cpdfs[nspec][neta]->SetTickx(); cpdfs[nspec][neta]->SetTicky();
			hclo[neta][nspec]->Draw();//data theory ratio
//			hclodat[neta][nspec]->Draw("same");//data theory ratio
			hexp[neta][0]->Draw("same");//experimental uncertainty
			hexp[neta][1]->Draw("same");
			//compare pdfs
			for(unsigned nsp(0);nsp<numNLO;++nsp) 
				if(nsp!=nspec)
				{
					sprintf(name,"cpdfs_eta%d_%s_ov_%s",neta,nloSpec[nsp],nloSpec[nspec]);
					TH1D* htmp = (TH1D*)hnlo[neta][nsp]->Clone(name);
          htmp->Divide(hnlo[neta][nspec]);
          pdfsDraw(htmp,lineCol[nsp]);
          hpdf[neta][nspec][nsp] = htmp;
          hpdf[neta][nspec][nsp]->Draw("same");
				}
			hclo[neta][nspec]->Draw("sameaxis");
		}//pdfs
	}//neta

	///draw legends
	TLegend* leg = new TLegend(0.55,0.6,0.75,0.85); // xy coordinates
	leg->SetFillColor(0);  leg->SetTextFont(42);
	leg->SetBorderSize(0); leg->SetTextSize(0.04);
	//latex and legend
	for(unsigned neta(0);neta<nybinsInfo;++neta) {
		for(unsigned nspec(0);nspec<numNLO;++nspec) {
			cratio[nspec][neta]->cd();
			char lumilabel[100], etalabel[100];
			TLatex ltxclo1;ltxclo1.SetTextSize(0.04);ltxclo1.SetTextFont(42);
			TLatex ltxclo2;ltxclo2.SetTextSize(0.04);ltxclo2.SetTextFont(42);
			sprintf(lumilabel,"#sqrt{s} = 2.76 TeV%*s L = %.2f pb^{-1}%*s CMS 2013",10,"",lumi,10,"");
			ltxclo1.DrawLatex(74,2.51,lumilabel);
			sprintf(etalabel,"%s  %.1f < |y| < %.1f",nloLabels[nspec],ybins[neta],ybins[neta+1]);
			ltxclo2.DrawTextNDC(0.3,0.2,etalabel);
			////legend
			if(neta==0 && nspec==0)
			{
				leg->AddEntry(hclo[neta][nspec],"Unfolded Data","pl");
				leg->AddEntry(hclodat[neta][nspec],"Raw Data","pl");
				leg->AddEntry(hexp[neta][0]    ,"Experiment ","l");
				leg->AddEntry(hthy[neta][nspec][0],"Theory","l");
			}
			leg->Draw("same");
			/////for pdfs
			cpdfs[nspec][neta]->cd();
			ltxclo1.DrawLatex(74,2.51,lumilabel);
			ltxclo2.DrawTextNDC(0.3,0.2,etalabel);
			TLegend* lpdf = new TLegend(0.25,0.6,0.55,0.85); // xy coordinates
			lpdf->SetFillColor(0);  lpdf->SetTextFont(42);
			lpdf->SetBorderSize(0); lpdf->SetTextSize(0.03);
			lpdf->AddEntry(hclo[neta][nspec],"Data/Theory","pl");
			lpdf->AddEntry(hexp[neta][0]    ,"Experimental Unc ","l");
			for(unsigned nsp(0);nsp<numNLO;++nsp)
				if(nsp!=nspec)
				{
					sprintf(name,"%s_over_%s",nloLabels[nsp],nloLabels[nspec]);
					lpdf->AddEntry( hpdf[neta][nspec][nsp],name,"l");
				}//loop over pdfs
			lpdf->Draw("same");
			//draw all eta bins together 
			if(neta==0)
			{
				sprintf(name,"ccombpdf_%s",nloSpec[nspec]);
				ccombpdf[nspec] = new TCanvas(name,name,1000,700);
				ccombpdf[nspec]->Divide(3,2);
			}
			ccombpdf[nspec]->cd(neta+1);
			cratio[nspec][neta]->DrawClonePad();
		}//pdfs
	}//eta

///save plots
	if(saveplots)
		for(unsigned neta(0);neta<nybinsInfo;++neta) 
			for(unsigned nspec(0);nspec<numNLO;++nspec) {
				sprintf(name,"plots/%s_testjec.eps",cratio[nspec][neta]->GetName());
				cratio[nspec][neta]->SaveAs(name);
				sprintf(name,"plots/%s_testjec.eps",cpdfs[nspec][neta]->GetName());
				cpdfs[nspec][neta]->SaveAs(name);
				if(neta==0)
				{
					sprintf(name,"plots/%s_testjec.eps",ccombpdf[nspec]->GetName());
					ccombpdf[nspec]->SaveAs(name);
				}
			}

	TFile* fout = new TFile("DataTheoryRatio.root","recreate");
	fout->cd();
	//write to file
	for(unsigned neta(0);neta<nybinsInfo;++neta) 
	{
		hexp[neta][0]->Write();
		hexp[neta][1]->Write();
		for(unsigned nspec(0);nspec<numNLO;++nspec) {
			hthy[neta][nspec][0]->Write();
			hthy[neta][nspec][1]->Write();
			hclo[neta][nspec]->Write();
			hnlo[neta][nspec]->Write();
		}
	}
	fout->Write();
	fout->Close();
}


void draw_datatheory_ResolutionFits()
{
	const bool saveplots(true);
	if(saveplots)   gROOT->SetBatch(kTRUE); ///display mode off----------
	const unsigned nybinsInfo(6);
	const char* infileQuantile("UnfoldJetsData_QuantileFit.root");
	const char* infileGauss   ("UnfoldJetsData_GaussFit.root");
	const char* nloDir("/uscms_data/d3/tlibeiro/inclusivexsec/Generators/generateIncJets/inputXsecs");
	TFile* npfile= new TFile("npCorrections.root","read");
	TFile* thyfile= new TFile("TheoryUncertainties.root","read");
	TFile* expfile= new TFile("ExpUncertainties.root","read");
	///////////////////////////////////////////////////
	char* nloSpec[] = {"abm11","nnpdf2.1","ct10","hera1.5","mstw08"};
	char* nloSpecout[] = {"abm11","nnpdf21","ct10","hera15","mstw08"};
	const unsigned numNLO(sizeof(nloSpec)/sizeof(string));
	char* nloLabels[numNLO] = {"ABM 11","NNPDF 2.1","CT10","HERA1.5","MSTW2008"};
  unsigned markerSty[8]={20,24,21,25,22,26,23,27};
	unsigned markerCol[8]={kRed,kMagenta,kBlue+1,kBlue,
		kGreen,kOrange+7,kRed-4,kGreen-7};

	TFile* fileQuantile =new TFile(infileQuantile,"read");
	TFile* fileGauss    =new TFile(infileGauss   ,"read");
  TH1D* hpt [nybinsInfo][2];//Unfolded data - 0- gauss and 1- quantile
  TH1D* hnlo[nybinsInfo][numNLO];//compare nlo spectrum
  TH1D* hthy[nybinsInfo][numNLO][2];//theory uncertainties
  TH1D* hexp[nybinsInfo][2];        //experimetnal uncertainties
  TH1D* hclo[nybinsInfo][numNLO][2];//dtmc ratio
  for(unsigned neta(0);neta<nybinsInfo;++neta) {
    //unfolded spectrum and uncertainties
		sprintf(name,"hUnf_eta%d",neta);
		hpt[neta][0] = (TH1D*)fileGauss   ->Get(name);
		hpt[neta][1] = (TH1D*)fileQuantile->Get(name);
    sprintf(name,"hUnfUnc_eta%d",neta);
    hexp[neta][1] = (TH1D*)expfile->Get(name);
    sprintf(name,"hUnfUnc_dn_eta%d",neta);
    hexp[neta][0] = (TH1D*)expfile->Get(name);
		//get nlo plots
		for(unsigned nspec(0);nspec<numNLO;++nspec) {
      //theory uncertainties
     sprintf(name,"htotalUp_eta%d_%s",neta,nloSpec[nspec]);
     hthy[neta][nspec][1]=(TH1D*)thyfile->Get(name);
     sprintf(name,"htotalDn_eta%d_%s",neta,nloSpec[nspec]);
     hthy[neta][nspec][0]=(TH1D*)thyfile->Get(name);
      //
			sprintf(name,"hclosure_gauss_%s_eta%d",nloSpec[nspec],neta);
			hclo[neta][nspec][0] = (TH1D*)hpt[neta][0]->Clone(name);
			hclo[neta][nspec][0]->Reset();
      sprintf(name,"hclosure_quant_%s_eta%d",nloSpec[nspec],neta);
			hclo[neta][nspec][1] = (TH1D*)hpt[neta][0]->Clone(name);
			hclo[neta][nspec][1]->Reset();

			char infile[500]; char title[500]; 
			sprintf(infile,"%s/%sy%i.cppInput",nloDir,nloSpec[nspec],neta);
			sprintf(title ,"nloy%i%s",neta,nloSpec[nspec]);
			hnlo[neta][nspec] = makeNLOSpectrum(infile,title);
			sprintf(name,"fratio_eta%d",neta);
			TF1* f = (TF1*)npfile->Get(name);
			applyNP(hnlo[neta][nspec],f);
		}
	}//eta
  cout<<"Finished getting histograms\n";
	////////////////////////////////////////////////////
	///draw data/theory histogram===============================
	TCanvas* cratio[numNLO][6];
	TCanvas* ccombpdf[numNLO];
	for(unsigned neta(0);neta<nybinsInfo;++neta) {
		for(unsigned nspec(0);nspec<numNLO ;++nspec) {
			//closure
			TH1D* htmp = (TH1D*)hpt[neta][0]->Clone("htmp");
			spFillHisto(hnlo[neta][nspec],htmp);
			hclo[neta][nspec][0]->Divide(hpt[neta][0],htmp);
			datatheoryplots(hclo[neta][nspec][0],neta);
			hclo[neta][nspec][1]->Divide(hpt[neta][1],htmp);
			datatheoryplots(hclo[neta][nspec][1],neta);
      hclo[neta][nspec][1]->SetMarkerStyle(25);
      hclo[neta][nspec][0]->SetMarkerColor(kBlue);
      hclo[neta][nspec][0]->SetLineColor(kBlue);
      hclo[neta][nspec][0]->GetYaxis()->SetRangeUser(0.9,1.1);
			//draw 
			sprintf(name,"cratio_eta%d_%s",neta,nloSpecout[nspec]);
			cratio[nspec][neta]=new TCanvas(name,name,600,600);
      cratio[nspec][neta]->SetLogx();
      cratio[nspec][neta]->SetTickx(); cratio[nspec][neta]->SetTicky();
      hclo[neta][nspec][0]->Divide( hclo[neta][nspec][1]);
		hclo[neta][nspec][0]->Draw();//data theory ratio
//			hclo[neta][nspec][1]->Draw("same");//data theory ratio
			if(nspec==0) {
				expuncDraw(hexp[neta][0]);//experimentla uncertainty
				expuncDraw(hexp[neta][1]);
			}
			hexp[neta][0]->Draw("same");
			hexp[neta][1]->Draw("same");
//      thyuncDraw(hthy[neta][nspec][0]);//theory uncertainty
//      thyuncDraw(hthy[neta][nspec][1]);
//      hthy[neta][nspec][0]->Draw("same");
//      hthy[neta][nspec][1]->Draw("same");
      hclo[neta][nspec][0]->Draw("axissame");
    }//pdfs
	}//eta

		///draw legends
	TLegend* leg = new TLegend(0.2,0.6,0.4,0.85); // xy coordinates
	leg->SetFillColor(0);  leg->SetTextFont(42);
	leg->SetBorderSize(0); leg->SetTextSize(0.04);
	//latex and legend
	for(unsigned neta(0);neta<nybinsInfo;++neta) {
		for(unsigned nspec(0);nspec<numNLO;++nspec) {
			cratio[nspec][neta]->cd();
			char lumilabel[100], etalabel[100];
			TLatex ltxclo1;ltxclo1.SetTextSize(0.04);ltxclo1.SetTextFont(42);
			TLatex ltxclo2;ltxclo2.SetTextSize(0.04);ltxclo2.SetTextFont(42);
			//			sprintf(lumilabel,"#sqrt{s} = 2.76 TeV%*s L = %.2f pb^{-1}%*s CMS 2013",10,"",lumi,10,"");
			sprintf(lumilabel,"%*s %.2f pb^{-1} (2.76 TeV)",50,"",lumi);
			ltxclo1.SetNDC();
			ltxclo1.DrawLatex(0.1,0.91,lumilabel);
			sprintf(etalabel,"%s",nloLabels[nspec]);
			ltxclo2.DrawTextNDC(0.2,0.30,etalabel);
			sprintf(etalabel,"%.1f < |y| < %.1f",ybins[neta],ybins[neta+1]);
			ltxclo2.DrawTextNDC(0.2,0.25,etalabel);
			//CMS Logo
			extraCmsText="Preliminary";
			TLatex ltxclo3;ltxclo3.SetTextSize(0.04);ltxclo3.SetTextFont(42);
			ltxclo3.SetTextSize(0.04); ltxclo3.SetTextFont(cmsTextFont);
			ltxclo3.SetNDC();
			ltxclo3.DrawTextNDC(0.1,0.91,cmsText);
			ltxclo3.SetTextSize(0.04*0.76); ltxclo3.SetTextFont(extraTextFont);
			ltxclo3.DrawTextNDC(0.2,0.91,extraCmsText);
			////legend
			if(neta==0 && nspec==0)
			{
				leg->AddEntry(hclo[neta][nspec][0],"Data/Theory Gaussian","pl");
				leg->AddEntry(hclo[neta][nspec][1],"Data/Theory Quantile","pl");
				leg->AddEntry(hexp[neta][0]    ,"Unfolding Uncertainty ","l");
//				leg->AddEntry(hthy[neta][nspec][0],"Theo. Uncertainty","l");
			}
			leg->Draw("same");
		}//pdfs
	}//eta

	///save plots
	if(saveplots)
		for(unsigned neta(0);neta<nybinsInfo;++neta) 
			for(unsigned nspec(0);nspec<numNLO;++nspec) {
				sprintf(name,"plots/%s_testResolutionFits.eps",cratio[nspec][neta]->GetName());
				cratio[nspec][neta]->SaveAs(name);
			}

	TFile* fout = new TFile("DataTheoryRatio_ResolutionFits.root","recreate");
	fout->cd();
	//write to file
	for(unsigned neta(0);neta<nybinsInfo;++neta) 
	{
		hexp[neta][0]->Write();
		hexp[neta][1]->Write();
		for(unsigned nspec(0);nspec<numNLO;++nspec) {
			hthy[neta][nspec][0]->Write();
			hthy[neta][nspec][1]->Write();
			hclo[neta][nspec][0]->Write();
			hclo[neta][nspec][1]->Write();
			hnlo[neta][nspec]->Write();
		}
	}
	fout->Write();
	fout->Close();
}


///HELPER FUNCTIONS----------------------------------------------------
void groomSpectrum(TH1* h, unsigned neta)
{
	float max (neta?maxPt[neta]-1:maxPt[neta]+10);
	h->SetTitle("");
	h->Scale(yScales[neta]);
	h->GetXaxis()->SetNoExponent(kTRUE);
	h->GetXaxis()->SetMoreLogLabels(kTRUE);
	h->GetYaxis()->SetRangeUser(1e-2,1e11);
	h->GetXaxis()->SetRangeUser(74,max);
	h->GetXaxis()->SetTitle("Jet P_{T} (GeV)");
	h->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{dp_{T}dy}(pb/GeV)");
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
	h->GetYaxis()->SetTitleSize(0.05);
	h->GetXaxis()->SetTitleSize(0.05);
	//	h->GetYaxis()->SetTitleOffset(1);
	h->GetXaxis()->SetLabelSize(0.04);
	h->GetYaxis()->SetTitleOffset(1.0);
	h->GetXaxis()->SetTitleOffset(0.9);
	h->GetYaxis()->SetTickLength(0.03*0.5);
	h->SetMarkerColor(kBlack);
	h->SetLineColor(kBlack);
	h->SetMarkerSize(0.9);
	h->SetMarkerStyle(yMarkers[neta]);
	h->SetStats(kFALSE);
}

void datatheoryplots(TH1* h, unsigned eta)
{
	h->SetTitle("");
	h->GetXaxis()->SetNoExponent(kTRUE);
	h->GetXaxis()->SetMoreLogLabels(kTRUE);
	h->GetYaxis()->SetRangeUser(0,2.5);
	h->GetXaxis()->SetRangeUser(74,maxPt[eta]-1);
	h->GetXaxis()->SetTitle("Jet P_{T} (GeV)");
	h->GetYaxis()->SetTitle("Data/Theory");
	h->GetYaxis()->SetTitleOffset(1.2);
	h->GetYaxis()->CenterTitle(kTRUE);
	h->GetXaxis()->SetTitleOffset(1.1);
	h->GetXaxis()->SetTitleSize(0.04);
	h->GetYaxis()->SetTitleSize(0.04);
	h->GetXaxis()->SetLabelSize(0.04);
	h->GetYaxis()->SetLabelSize(0.04);
	h->SetStats(kFALSE);
	h->SetLineColor(kBlack);
	h->SetMarkerStyle(20);
}



void groomSpecClosure(TH1* h, unsigned eta)
{
	h->SetTitle("");
	h->GetXaxis()->SetNoExponent(kTRUE);
	h->GetXaxis()->SetMoreLogLabels(kTRUE);
	h->GetYaxis()->SetRangeUser(0,2.5);
	h->GetXaxis()->SetRangeUser(74,maxPt[eta]-1);
	h->GetXaxis()->SetTitle("Jet Pt GeV/c");
	h->GetYaxis()->SetTitle("Ratio to theory");
	//	h->GetYaxis()->SetTitleOffset(0.5);
	h->GetYaxis()->SetTitleSize(0.04);
	h->GetXaxis()->SetLabelSize(0.04);
	h->GetYaxis()->SetLabelSize(0.04);
	h->SetStats(kFALSE);
	h->SetLineColor(kBlack);
	h->SetMarkerStyle(20);
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
}

void spFillHisto(TH1D* hin, TH1D* hout)
{
	hout->Reset();
	const  unsigned numbins(hout->GetNbinsX());
	const  unsigned numbinsIn(hin->GetNbinsX());
	cout<<"numbins "<<numbins<<' '<<hout->GetName()<<endl;
	for(unsigned i(1);i<=numbins;++i)
		for(unsigned j(1);j<=numbinsIn;++j)
			if(((hin->GetBinLowEdge(j)+hin->GetBinWidth(j))==
						(hout->GetBinLowEdge(i)+hout->GetBinWidth(i))	)	
					&& hin->GetBinLowEdge(j)==hout->GetBinLowEdge(i))
			{
				hout->SetBinContent(i,hin->GetBinContent(j));
				//				cout<<"Setting bin content "<<hin->GetBinContent(j)<<endl;
			}
}

void spFillHistoWithErr(TH1D* hin, TH1D* hout)
{
	hout->Reset();
	const  unsigned numbins(hout->GetNbinsX());
	const  unsigned numbinsIn(hin->GetNbinsX());
	cout<<"numbins "<<numbins<<' '<<hout->GetName()<<endl;
	for(unsigned i(1);i<=numbins;++i)
		for(unsigned j(1);j<=numbinsIn;++j)
			if(((hin->GetBinLowEdge(j)+hin->GetBinWidth(j))==
						(hout->GetBinLowEdge(i)+hout->GetBinWidth(i))	)	
					&& hin->GetBinLowEdge(j)==hout->GetBinLowEdge(i))
			{
				hout->SetBinContent(i,hin->GetBinContent(j));
				hout->SetBinError(i,hin->GetBinError(j));
				//				cout<<"Setting bin content "<<hin->GetBinContent(j)<<endl;
			}
}


void drawClosure(TCanvas* c,TH1* hdt, TH1* hmc,TH1* hclo)
{
	hclo->GetXaxis()->SetTitle(hmc->GetXaxis()->GetTitle());
	hclo->GetXaxis()->SetTitleSize(0.1);
	hdt->GetXaxis()->SetTitle("");
	c->cd();
	sprintf(name,"tpad_%s",hdt->GetName());
	double ytop(0.25);
	TPad* tpad = new TPad(name,name,0,ytop, 1.0, 1.0);
	tpad->SetLogx(); tpad->SetLogy();
	tpad->Draw();
	tpad->cd();
	tpad->SetBottomMargin(0.1);
	hdt->Draw();
	hmc->Draw("hist same");

	c->cd();
	double ybot(0.3);
	sprintf(name,"bpad_%s",hclo->GetName());
	TPad* bpad = new TPad("bpad", "bpad", 0, 0.03, 1.0, 0.3);
	bpad->SetLogx();
	bpad->Draw();
	bpad->cd();
	bpad->SetTopMargin(0.1);
	bpad->SetBottomMargin(0.2);
	bpad->SetGrid(1,1);
	hclo->Draw();
	tpad->cd();
}


void calculateErrorOnStdDev()
{
	const unsigned numIterations(10000);
	const double sdev(1);
	TF1* f1 =  new TF1("f1","gaus",10,20);
	f1->SetParameter(1,15);//mean
	f1->SetParameter(2,sdev);//std dev
	f1->SetParameter(0,1);//normalization
	f1->Draw();

	TH1D* h1 = new TH1D("h1","h1",1000,10,20);
	TH1D* h2 = new TH1D("h2","h2",500,-0.2,0.2);
	double xq[3]={0.1583,0.50,0.8417};
	double yq[3];
	for(unsigned i(0);i<numIterations;++i)
	{
		h1->FillRandom("f1",10000);
		h1->GetQuantiles(3,yq,xq);
		double stddev(yq[2]-yq[0]);
		h2->Fill(stddev/2-sdev);
		//  cout<<stddev<<endl;
		const unsigned bins(h1->GetNbinsX());
		if(i!=numIterations-1)
			for(unsigned bin(0);bin<bins;++bin)
			{	
				h1->SetBinContent(bin,0);
			}
	}
	TCanvas* c1 = new TCanvas();
	h2->Draw();
	h2->GetQuantiles(3,yq,xq);
	cout<<"Std Dev h2 "<<(yq[2]-yq[0])/2<<endl;
	TCanvas* c2 = new TCanvas();
	h1->Draw();
	f1->Draw("same");
}


void expuncDraw(TH1*  h, int uncColor)
{
	h->SetFillColor(0);
	h->SetLineColor(uncColor);
	h->SetLineStyle(1);
	h->SetLineWidth(1);
	h->GetYaxis()->SetRangeUser(0,2.5);
	h->GetXaxis()->SetTitle("Jet P_{T} (GeV)");

	unsigned bins(h->GetNbinsX());
	for(unsigned i(1);i<=bins;++i)
		h->SetBinContent(i,1.+h->GetBinContent(i));
}

void thyuncDraw(TH1*  h)
{
	h->SetFillColor(0);
	h->SetLineColor(kMagenta);
	h->GetYaxis()->SetRangeUser(0,2.5);
	h->GetXaxis()->SetTitle("Jet P_{T} (GeV)");
	h->SetLineStyle(5);
	h->SetLineWidth(1);
	h->GetYaxis()->SetRangeUser(0,2.5);
	unsigned bins(h->GetNbinsX());
	for(unsigned i(1);i<=bins;++i)
		h->SetBinContent(i,1.+h->GetBinContent(i));
}

void pdfsDraw(TH1*  h,int col)
{
	h->SetLineColor(col);
	h->SetLineStyle(9);
	h->SetLineWidth(1);
	h->GetYaxis()->SetRangeUser(0,2.5);
	h->GetYaxis()->SetRangeUser(0,2.5);
	h->GetXaxis()->SetTitle("Jet P_{T} (GeV)");
}


TH1D* makePtRangeHisto(TH1D* hin,unsigned eta)
{
	char* name = (char*)hin->GetName();
	sprintf(name,"%s_%s",name,"inrange");
	TH1D*  hout = new TH1D(name,name,n2x[eta],x2[eta]);
	spFillHistoWithErr(hin,hout);
	return hout;
}


void applyNP(TH1* h,TF1* f)
{
	const  unsigned numbins(h->GetNbinsX());
	for(unsigned i(1);i<=numbins;++i)
	{
		double np = f->Eval(h->GetBinCenter(i));
		h->SetBinContent(i,h->GetBinContent(i)*np);
		//		cout<<"NPCorr "<<h->GetBinCenter(i)<<" "<<np<<endl;
	}
}
