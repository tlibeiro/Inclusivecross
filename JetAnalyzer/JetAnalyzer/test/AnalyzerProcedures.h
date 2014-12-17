#ifndef ANALYZERPROCEDURES_HH
#define ANALYZERPROCEDURES_HH
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TProfile.h"
#include <TMath.h>
#include <iomanip>
#include <TColor.h>
#include <TH1F.h>
#include <TLorentzVector.h>
using namespace std;
double dR = 0.5;
double ptMin = 20;
//double ptbins[] = {74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1784,2116,2500};
double ybins[]  = {0.0,0.5,1.0,1.5,2.0,2.5,3.0};
//const int nptbins = sizeof(ptbins)/sizeof(double);
const int nybins  = sizeof(ybins) /sizeof(double);
const int nx[6] = {37,37,36,32,25,18};
const int n1x[6] = {47,47,46,42,35,28};
const double x[6][50]= {
	{74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1784,2116,2500},

	{74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1784,2116,2500},

	{74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1784,2116},

	{74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410},

	{74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905},

	{74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548}
};
//21 24 27 31 35 40 45 51 58 66
const int nxhlt[6] = {47,47,46,42,35,28};
const double xhlt[6][50]= {
	{21,24,27,31,35,40,45,51,58,66,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1784,2116,2500},

	{21,24,27,31,35,40,45,51,58,66,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1784,2116,2500},

	{21,24,27,31,35,40,45,51,58,66,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1784,2116},

	{21,24,27,31,35,40,45,51,58,66,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410},

	{21,24,27,31,35,40,45,51,58,66,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905},

	{21,24,27,31,35,40,45,51,58,66,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548}
};



const double x1[6][50]= {
	{18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1784,2116,2500},

	{18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1784,2116,2500},

	{18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1784,2116},

	{18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410},

	{18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905},

	{18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548}
};

//for unfolding test only
const int n2x[6] = {19,18,16,13,9,6};
const double x2[6][50]= {
	{74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592},

	{74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548},

	{74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468},

	{74,84,97,114,133,153,174,196,220,245,272,300,330,362},

	{74,84,97,114,133,153,174,196,220,245},

	{74,84,97,114,133,153,174}
};

const int n3x[6] = {27,26,24,21,17,14};
const double x3[6][50]= {
	{28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638},
                        
	{28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592},
                        
	{28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507},
                        
	{28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395},
                        
	{28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272},
                        
	{28,32,37,43,49,56,64,74,84,97,114,133,153,174,196}
};

const int n4x[6] = {23,22,20,17,13,10};
const double x4[6][50]= {
	{49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638},
            
	{49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592},
            
	{49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507},
            
	{49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395},
            
	{49,56,64,74,84,97,114,133,153,174,196,220,245,272},
            
	{49,56,64,74,84,97,114,133,153,174,196}
};


double maxPt[nybins] = {592,548,468,362,245,174};
double maxPtBins[nybins] = {19,18,16,13,9,6}; // bins from 74 gev to maxpt
//triggers 
string trignames[] = { "HLT_PAJet20","HLT_PAJet40","HLT_PAJet60","HLT_PAJet80"};
//string trignames[] = { "HLT_PAJet20","HLT_PAJet40", "HLT_PAJet60","HLT_PAJet80","HLT_PAJet100","HLT_PAJet120" } ;
string trignames8Tev[] = {"HLT_PFJet40","HLT_PFJet80","HLT_PFJet140",
	"HLT_PFJet200","HLT_PFJet260","HLT_PFJet320",
	"HLT_PFJet400"};
const unsigned numtrigs (sizeof(trignames)/sizeof(string));
const unsigned numtrigs8Tev (sizeof(trignames8Tev)/sizeof(string));

char name[200],title[200],trigtitle[200];
//bool hltPassj[numtrigs]; int hltPsc[numtrigs];
//bool l1Passj[numtrigs];  int l1Psc[numtrigs];
//int  ATrig[numtrigs] = {16,36,36,36,36,36};
//int  HLTJetPtN[numtrigs] = {20,40,60,80,100,120};
int  HLTJetPtN[numtrigs] = {20,40,60,80};
int  ATrig8Tev[numtrigs8Tev] = {16,36,68,92,128,128};
int  HLTJetPtN8Tev[numtrigs8Tev] = {40,80,140,200,260,320,400};

//int  HLTJetPtS[numtrigs] = {--,--,97,112,138,158}; //turn on thresholds
//int  HLTJetPtS[numtrigs+1+2] = {40,74,97,133,153,174,2500}; //turn on thresholds aligned with bins
//int HLTJetPtS[numtrigs+1] = {40,74,97,133,153,174,2500}; //turn on thresholds aligned with bins
int HLTJetPtS[numtrigs+1] = {40,74,97,133,2500}; //turn on thresholds aligned with bins
int HLTJetPtS8Tev[numtrigs8Tev+1] = {74,133,220,300,395,507,2500}; //turn on thresholds aligned with bins
int hltColors[numtrigs8Tev+2] = {kBlack,kBlack,kRed,kGreen,kBlue,kMagenta,kCyan,kBlack};
int yMarkers[nybins] = {21,22,23,24,25,26};
int yMarkersRes[nybins] = {20,21,22,23,33,34};
int yColors [nybins] = {kBlack,kRed,kBlue,kGreen,kMagenta,kCyan};
unsigned yScales[nybins] = {1e6,1e5,1e4,1e3,1e2,1e1};
int  mcSamplePt[] = {15,30,50,80,120,170,220,280,370,460,540};
const int nMCSamples(sizeof(mcSamplePt)/sizeof(int));
double mcSampleXsecmb[nMCSamples] = {
	2.034E-01, 1.075E-02, 1.025E-03,
  9.865E-05, 1.129E-05,	1.465E-06,
	2.837E-07, 5.323E-08, 5.934E-09,
  8.125E-10, 1.467E-10  
};

// double herwigXsecnb[11] = {
//	0.3003E+06, 15.08E+03,  1.401E+03,
//  0.1284E+03, 13.51E+0 ,	1.529E+0 ,
//	0.2815E+0 , 54.10E-03,  5.59E-03 ,
//  0.704E-3  , 0.1532E-3  
//};

 double herwigXsecnb[11] = {
	0.1522E+06, 7.71E+03,   0.737E+03,
  0.070E+03,  7.81E+0 ,	  0.951E+0 ,
  0.1870E+0,  38.23E-03,  4.187E-03 ,
  0.544E-3  , 0.1214E-3  
};



//temporary change from default 500000,
//some qcd jobs failed
double mcSampleEvents[nMCSamples] = {
550000,600000,500000,
600000,500000,550000,
550000,550000,650000,
750000,700000
};

//double mcSampleEvents[nMCSamples] = {
//500000, 500000, 500000,
//500000, 500000, 500000,
//500000, 500000, 500000,
//500000, 500000
//};


const double lumi(5.368);///
//const double lumi(5.07);///winter14 corr 
//const double lumi(4.689);///winter14 corr test
//const double lumi(0.370427);///summer13 corr test
//const double lumi(5.219);///ft53v10

//const double lumi(5.256);
//const double lumi(4.277);

///for CMS logo settings
char* cmsText        = "CMS";
char* extraCmsText   = "Preliminary";
float  cmsTextFont   = 61;  // default is helvetic-bold
float  extraTextFont = 52;  // default is helvetica-italics
// text sizes and text offsets with respect to the top frame
// in unit of the top margin size
float lumiTextSize     = 0.6;
float lumiTextOffset   = 0.2;
float cmsTextSize      = 0.75;
// ratio of "CMS" and extra text size
float extraCmsTextSize  = 0.76*cmsTextSize;

/////Jet Class 
class Jet {
	public : 
		double pt;              double neutralSum   ;
		double eta;             double neutralN     ;
		double phi;             double hcalSum      ;
		double mass;            double ecalSum      ;
		double chargedMax   ;   double nHCAL;
		double chargedSum   ;   double nECAL;
		double chargedN     ;   double eMax         ;
		double chargedHardSum;  double eSum         ;
		double chargedHardN ;   double eN           ;
		double photonMax    ;   double muMax;
		double photonSum    ;   double muSum;
		double photonN      ;   double muN;
		double photonHardSum;   double chf, nhf;
		double photonHardN  ;   double y;
		double neutralMax   ;
		//constructors
		Jet(double p, double e, double ph, double m) 
		{
			pt = p;              	neutralSum    = -1;
			eta = e;             	neutralN      = -1;
			phi = ph;            	hcalSum       = -1;
			mass = m;            	ecalSum       = -1;
			chargedMax    = -1;  	nHCAL = -1;
			chargedSum    = -1;  	nECAL = -1;
			chargedN      = -1;  	eMax          = -1;
			chargedHardSum = -1; 	eSum          = -1;
			chargedHardN  = -1;  	eN            = -1;
			photonMax     = -1;  	muMax = -1;
			photonSum     = -1;  	muSum = -1;
			photonN       = -1;  	muN = -1;
			photonHardSum = -1;   nhf=-1;
			photonHardN   = -1;   chf=-1;
			neutralMax    = -1;   y=-1;
		};

		Jet()
			: pt(1e-8),eta(1e-8),phi(1e-8),mass(1e-8)	
		{ 
			chargedMax    = -1;        			ecalSum       = -1;
			chargedSum    = -1;  			nHCAL = -1;
			chargedN      = -1;  			nECAL = -1;
			chargedHardSum = -1; 			eMax          = -1;
			chargedHardN  = -1;  			eSum          = -1;
			photonMax     = -1;  			eN            = -1;
			photonSum     = -1;  			muMax = -1;
			photonN       = -1;  			muSum = -1;
			photonHardSum = -1;  			muN = -1;
			photonHardN   = -1;       neutralN      = -1;
			neutralMax    = -1;       hcalSum       = -1;
			neutralSum    = -1;       nhf=-1; chf=-1;
			y=-1;
		};
		~Jet(){ };

		Jet(const TLorentzVector& vec ){
			if(vec.Pt()>0)
			{
				pt = vec.Pt();eta = vec.Eta();phi = vec.Phi();mass = vec.M(); 
			}
			else 
			{
				pt = 1e-8;eta = 1e-8;phi = 1e-8;mass = 1e-8; 
			}
		};
		inline TLorentzVector  getLorentzVector(){
			TLorentzVector vec;
			vec.SetPtEtaPhiM(pt,eta,phi,mass);
			return (vec);
		};
}; 

double deltaPhi(double phi1,double phi2)
{ 
	double deltaPhi = phi1 - phi2;
	if (deltaPhi < -M_PI)
		deltaPhi += 2.0*M_PI;
	if (deltaPhi > M_PI)
		deltaPhi -= 2.0*M_PI;
	return deltaPhi;
};

void expoBins(double from, double to, unsigned steps, double* bins)
{
  double delta=log(to)-log(from);
  delta /= steps-1;
  cout<<fixed;
  cout<<setprecision(2);
  for(unsigned i(0);i<steps;++i)
  {
//    cout<<i<<" "<<pow(2.71828,log(from)+delta*i)<<endl;
   bins[i]= pow(2.71828,log(from)+delta*i);
  }
}

unsigned  findBinsFromArray(double xlo,double xhi,const double* xarr, unsigned arrsize) //finds bins in the interval xlow  - xhigh 
{
	unsigned bins(0);
	for(unsigned i(0);i<arrsize;++i) 
 {  
		if(xarr[i]<xhi && xarr[i]>=xlo)
			++bins;
  }
   return unsigned(bins);
};

#endif
