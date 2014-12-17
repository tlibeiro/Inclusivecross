#ifndef  JETANALYZER_HH
#define  JETANALYZER_HH

namespace JA {
/////Jet Class 
class Jet {
	public : 
		double pt;       double mass;     
		double eta;      double y;
		double phi;      
		//constructors
		Jet(double p, double e, double ph, double m) 
		{
			pt = p;   	mass = m; 
			eta = e;  	y=-1;
			phi = ph; 
		};

		Jet()
			: pt(1e-8),eta(1e-8),phi(1e-8),mass(1e-8)	
		{ 
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

class DeltaRDistance
{
	public:
		template<class JetType>
			double operator()(const JetType& jet1, const JetType& jet2) const
			{
				const double deltaEta = jet1.eta - jet2.eta;
				double deltaPhi = jet1.phi - jet2.phi;
				if (deltaPhi < -M_PI)
					deltaPhi += 2.0*M_PI;
				if (deltaPhi > M_PI)
					deltaPhi -= 2.0*M_PI;
				const double distance = sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
				return (distance);
			}
};

}
#endif
