#define Analysis_cxx
#include "Analysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TF1.h>
#define PI 3.141592

const int NUM_ANAL_E = 12; // 4-15 GeV
const int NUM_FLAVOUR = 4; // NuEl, AntiNuEl, NuMu, AntiNuMu
const int NUM_CUT = 6;
const int NUM_ANAL = 7;
const double LEP_CUT = 30;
const double PROT_CUT = 50;
const double NEUT_CUT = 50;
const double CHPI_CUT = 100;
const double PHOT_CUT = 30;
const double CHKA_CUT = 50; // Charged Kaon Energy Threshold
const double NEKL_CUT = 50; // Neutral K-long Energy Threshold
const double THETA_L_CUT = PI/90; // Lower bounds of theta_l
const double THETA_J_CUT = PI/36; // Lower bounds of theta_j
const double THETA_CUT1 = PI/36; // Upper bounds of theta_plane
const double THETA_L_UPP = PI/15; // Upper bounds of theta_l
const double RES_E_LC = .05; // Resolution of the Energy of the Contained Lepton
const double RES_E_LE = .15; // Resolution of the Energy of the Exiting Lepton
const double RES_ANG_L = PI/180; // Angular Resolution of the Lepton
const double RES_E_J = .4; // Resolution of the Energy of the Jet
const double RES_ANG_J1 = PI/180; // First Angular Resolution of the Jet
const double RES_ANG_J3 = PI/60; // Second Angular Resolution of the Jet
const double RES_ANG_J5 = PI/36; // Third Angular Resolution of the Jet
const double TOTAL_ENERGY_THRESHOLD = 0; //ANAL_E*500;
const double TOTAL_ENERGY_UPP = 100000; //ANAL_E*1500;

const int NUM_HISTO = 39;
const int NOCUT = 21;

enum {SIGNAL=11, SIGNAL_C1=12, SIGNAL_C3=13, SIGNAL_C5=14, SIGNAL_E1=15, SIGNAL_E3=16, SIGNAL_E5=17};
enum {ELECTRON=11, ELECTRON_NEUTRINO=12, ANTI_ELECTRON_NEUTRINO=-12, MUON=13, MUON_NEUTRINO=14, ANTI_MUON_NEUTRINO=-14, PROTON=2212, NEUTRON=2112, CHARGED_PION=211, NEUTRAL_PION=111, PHOTON=22, CHARGED_KAON=321, NEUTRAL_K_LONG=130};

class Hadron {
private:
	TLorentzVector had;
	int pdgNum; // The PDG number of Hadron
public:
	Hadron() : had(0,0,0,0), pdgNum(0) {}
	Hadron(TLorentzVector par, int num) : had(par), pdgNum(num) {}

	inline int Num() const { return pdgNum; }
	inline TLorentzVector Had() const { return had; }
	inline void SetHadron(TLorentzVector par, int num) {
		had = par;
		pdgNum = num;
	}
	inline void ShowHadronInfo() const {
		cout<<"pdgNum: "<<pdgNum<<endl;
		cout<<"E: "<<had.E()<<endl;
		cout<<"theta: "<<had.Theta()/PI*180.<<endl<<endl;
	}
};

class Control {
private:
	TLorentzVector lep[NUM_ANAL], jet[NUM_ANAL]; // 4-vectors of the lepton and the jet
	Hadron* had[100]; // storage of outgoing hadrons' 4-vectors
	int hadNum; // The number of Hadrons
public:
	Control() : hadNum(0) {}
	~Control() { for(int i=0; i<hadNum; i++) delete had[i]; }

	inline void AddToLep(TLorentzVector par);
	inline void AddToLepNoCut(TLorentzVector par);
	inline void AddToJet(int pdgNum, TLorentzVector par, double KEThreshold);
	inline void MakeJet(int pdgNum, TLorentzVector par);
	inline void MakeJetNoCut(int pdgNum, TLorentzVector par);
	inline void CalculateCutPass(int cutNum);
	inline void RemoveNucleonMassFromTarget(Int_t pr, Int_t nr);
	inline void AnalyseBasics();
	inline bool DetermineIfPositiveEnergy();
	inline TLorentzVector Smear(TLorentzVector par, double energy_resolution, double theta_r, double phi_r);
	inline void CreateSmearedLepAndJet();
	inline void ProjectJet(TLorentzVector lep, TLorentzVector &jet);
	inline void Analyse();
	inline void AnalyseNoCut();
	inline void ShowEventInformation(const int INPUT_E, int anal) const;
};

inline void DetermineBinnedCutPass();
inline void PrintAllCutPass();
inline void WriteHistogramsToDataFiles();
inline void PrintMeansSDsSigmas(const int INPUT_E);
inline void PrintMeansSDsSigmasNoCut(const int INPUT_E, const int FLAVOUR);

TH1F * histo[NUM_HISTO];
TRandom3 random3;
int cutPass[NUM_ANAL][NUM_CUT+1]={0};

void Analysis::Loop(const int INPUT_E, const int FLAVOUR) {
//   In a ROOT session, you can do:
//      Root > gSystem->CompileMacro("Analysis.C");
//      Root > Analysis t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//
//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
	if (fChain == 0) return;

	const double NUM_BIN = 20*INPUT_E+1;
	const double END_RANGE = 100*NUM_BIN;
	histo[0] = new TH1F("E_l","Lepton Energy", NUM_BIN, 0, END_RANGE);
	histo[1] = new TH1F("theta_l","Lepton Angle with respect to the direction from the Sun",180,0,180);
	histo[2] = new TH1F("E_j","Jet Energy", NUM_BIN, 0, END_RANGE);
	histo[3] = new TH1F("theta_j","Jet Angle with respect to the direction from the Sun",180,0,180);
	histo[4] = new TH1F("delta_phi","Angular jet size",180,0,180);
	histo[5] = new TH1F("theta_plane_C1","Angle between the jet and the event plane with C1",90,0,90);
	histo[6] = new TH1F("theta_plane_C3","Angle between the jet and the event plane with C3",90,0,90);
	histo[7] = new TH1F("theta_plane_C5","Angle between the jet and the event plane with C5",90,0,90);
	histo[8] = new TH1F("theta_plane_E1","Angle between the jet and the event plane with E1",90,0,90);
	histo[9] = new TH1F("theta_plane_E3","Angle between the jet and the event plane with E3",90,0,90);
	histo[10] = new TH1F("theta_plane_E5","Angle between the jet and the event plane with E5",90,0,90);
	histo[SIGNAL] = new TH1F("E_nu","Reconstructed neutrino Energy without smearing", NUM_BIN, 0, END_RANGE);
	histo[SIGNAL_C1] = new TH1F("E_nu_C1","Reconstructed neutrino Energy with C1", NUM_BIN, 0, END_RANGE);
	histo[SIGNAL_C3] = new TH1F("E_nu_C3","Reconstructed neutrino Energy with C3", NUM_BIN, 0, END_RANGE);
	histo[SIGNAL_C5] = new TH1F("E_nu_C5","Reconstructed neutrino Energy with C5", NUM_BIN, 0, END_RANGE);
	histo[SIGNAL_E1] = new TH1F("E_nu_E1","Reconstructed neutrino Energy with E1", NUM_BIN, 0, END_RANGE);
	histo[SIGNAL_E3] = new TH1F("E_nu_E3","Reconstructed neutrino Energy with E3", NUM_BIN, 0, END_RANGE);
	histo[SIGNAL_E5] = new TH1F("E_nu_E5","Reconstructed neutrino Energy with E5", NUM_BIN, 0, END_RANGE);
	histo[18] = new TH1F("p_nu","Reconstructed neutrino Momentum without smearing", NUM_BIN, 0, END_RANGE);
	histo[19] = new TH1F("p_nu_C1","Reconstructed neutrino Momentum with C1", NUM_BIN, 0, END_RANGE);
	histo[20] = new TH1F("p_nu_C3","Reconstructed neutrino Momentum with C3", NUM_BIN, 0, END_RANGE);
	histo[21] = new TH1F("p_nu_C5","Reconstructed neutrino Momentum with C5", NUM_BIN, 0, END_RANGE);
	histo[22] = new TH1F("p_nu_E1","Reconstructed neutrino Momentum with E1", NUM_BIN, 0, END_RANGE);
	histo[23] = new TH1F("p_nu_E3","Reconstructed neutrino Momentum with E3", NUM_BIN, 0, END_RANGE);
	histo[24] = new TH1F("p_nu_E5","Reconstructed neutrino Momentum with E5", NUM_BIN, 0, END_RANGE);
	histo[25] = new TH1F("E_tot","Lepton energy plus jet energy", NUM_BIN, 0, END_RANGE);
	histo[26] = new TH1F("E_tot_C1","Lepton energy plus jet energy with C1", NUM_BIN, 0, END_RANGE);
	histo[27] = new TH1F("E_tot_C3","Lepton energy plus jet energy with C3", NUM_BIN, 0, END_RANGE);
	histo[28] = new TH1F("E_tot_C5","Lepton energy plus jet energy with C5", NUM_BIN, 0, END_RANGE);
	histo[29] = new TH1F("E_tot_E1","Lepton energy plus jet energy with E1", NUM_BIN, 0, END_RANGE);
	histo[30] = new TH1F("E_tot_E3","Lepton energy plus jet energy with E3", NUM_BIN, 0, END_RANGE);
	histo[31] = new TH1F("E_tot_E5","Lepton energy plus jet energy with E5", NUM_BIN, 0, END_RANGE);
	histo[SIGNAL+NOCUT] = new TH1F("E_nu_B_nocut", "Reconstructed background neutrino Energy without smearing nocut", NUM_BIN, 0, END_RANGE);
	histo[SIGNAL_C1+NOCUT] = new TH1F("E_nu_C1_nocut", "Reconstructed neutrino Energy with C1 nocut", NUM_BIN, 0, END_RANGE);
	histo[SIGNAL_C3+NOCUT] = new TH1F("E_nu_C3_nocut", "Reconstructed neutrino Energy with C3 nocut", NUM_BIN, 0, END_RANGE);
	histo[SIGNAL_C5+NOCUT] = new TH1F("E_nu_C5_nocut", "Reconstructed neutrino Energy with C5 nocut", NUM_BIN, 0, END_RANGE);
	histo[SIGNAL_E1+NOCUT] = new TH1F("E_nu_E1_nocut", "Reconstructed neutrino Energy with E1 nocut", NUM_BIN, 0, END_RANGE);
	histo[SIGNAL_E3+NOCUT] = new TH1F("E_nu_E3_nocut", "Reconstructed neutrino Energy with E3 nocut", NUM_BIN, 0, END_RANGE);
	histo[SIGNAL_E5+NOCUT] = new TH1F("E_nu_E5_nocut", "Reconstructed neutrino Energy with E5 nocut", NUM_BIN, 0, END_RANGE);

	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries; jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if(ientry < 0) break;
		nb = fChain->GetEntry(jentry); nbytes += nb;
		// if (Cut(ientry) < 0) continue;

	    Control control, nocut;
		if(post_>100) cout<<"Warning: The number of post particles exceeds 100."<<endl;

	//CYCLE THROUGH EACH "POST" PARTICLE IN THE EVENT.
		for(int outpar=0; outpar<post_; outpar++) {
			//Get the PDG particle-code, the energy, and the momentum of particles
			int pdgNum = post_pdg[outpar];
			TLorentzVector par(post_x[outpar], post_y[outpar], post_z[outpar], post_t[outpar]);

			if(abs(pdgNum) == MUON || abs(pdgNum) == ELECTRON) {
				control.AddToLep(par);
				nocut.AddToLepNoCut(par);
			}
			else {
				control.MakeJet(pdgNum, par);
				nocut.MakeJetNoCut(pdgNum, par);
			}
		} //THE END OF LOOP OVER THE PARTICLES IN AN EVENT

		control.CalculateCutPass(0); // Number of Raw Data Events
		control.RemoveNucleonMassFromTarget(pr, nr);
		if(control.DetermineIfPositiveEnergy()) {
			control.CalculateCutPass(1); // Number of Events above KE Thresholds
			control.CreateSmearedLepAndJet();
			control.AnalyseBasics();
			control.Analyse();
			control.ShowEventInformation(INPUT_E, 6); // E5 = 6
		}

		nocut.RemoveNucleonMassFromTarget(pr, nr);
		if(nocut.DetermineIfPositiveEnergy()) {
			nocut.CreateSmearedLepAndJet();
			nocut.AnalyseNoCut();
		}

	} // THE END OF LOOP OVER THE EVENTS

	DetermineBinnedCutPass();
	PrintAllCutPass();
	WriteHistogramsToDataFiles();
	PrintMeansSDsSigmas(INPUT_E);
	PrintMeansSDsSigmasNoCut(INPUT_E, FLAVOUR);
}





inline void Control::AddToLep(TLorentzVector par) { if(par.E()-par.M() > LEP_CUT) lep[0] += par; }

inline void Control::AddToLepNoCut(TLorentzVector par) { if(par.E()-par.M() > 0) lep[0] += par; }

inline void Control::AddToJet(int pdgNum, TLorentzVector par, double KEThreshold) {
    if(par.E()-par.M() > KEThreshold) {
        jet[0] += par;
        had[hadNum++] = new Hadron(par, pdgNum);
    }
}

inline void Control::MakeJet(int pdgNum, TLorentzVector par) {
    switch(abs(pdgNum)) {
        case PROTON:
            AddToJet(pdgNum, par, PROT_CUT);
            break;
        case NEUTRON:
            AddToJet(pdgNum, par, NEUT_CUT);
            break;
        case CHARGED_PION:
            AddToJet(pdgNum, par, CHPI_CUT);
            break;
        case PHOTON:
            AddToJet(pdgNum, par, PHOT_CUT);
            break;
        case CHARGED_KAON:
            AddToJet(pdgNum, par, CHKA_CUT);
            break;
        case NEUTRAL_K_LONG:
            AddToJet(pdgNum, par, NEKL_CUT);
            break;
        case NEUTRAL_PION:
            double E = par.M()/2;
            double X, Y, Z;
            random3.Sphere(X, Y, Z, E);
            TLorentzVector ph1(X, Y, Z, E), ph2(-X, -Y, -Z, E);
            TVector3 boostfactor = par.BoostVector();
            ph1.Boost(boostfactor);
            ph2.Boost(boostfactor);
            AddToJet(pdgNum, ph1, PHOT_CUT);
            AddToJet(pdgNum, ph2, PHOT_CUT);
    }
}

inline void Control::MakeJetNoCut(int pdgNum, TLorentzVector par) {
    switch(abs(pdgNum)) {
        case PROTON:
            AddToJet(pdgNum, par, 0);
            break;
        case NEUTRON:
            AddToJet(pdgNum, par, 0);
            break;
        case CHARGED_PION:
            AddToJet(pdgNum, par, 0);
            break;
        case PHOTON:
            AddToJet(pdgNum, par, 0);
            break;
        case CHARGED_KAON:
            AddToJet(pdgNum, par, 0);
            break;
        case NEUTRAL_K_LONG:
            AddToJet(pdgNum, par, 0);
            break;
        case NEUTRAL_PION:
            double E = par.M()/2;
            double X, Y, Z;
            random3.Sphere(X, Y, Z, E);
            TLorentzVector ph1(X, Y, Z, E), ph2(-X, -Y, -Z, E);
            TVector3 boostfactor = par.BoostVector();
            ph1.Boost(boostfactor);
            ph2.Boost(boostfactor);
            AddToJet(pdgNum, ph1, 0);
            AddToJet(pdgNum, ph2, 0);
    }
}

inline void Control::CalculateCutPass(int cutNum) { for(int anal=0; anal<NUM_ANAL; anal++) cutPass[anal][cutNum]++; }

inline void Control::RemoveNucleonMassFromTarget(Int_t pr, Int_t nr) { jet[0].SetE((pr-18.)*938.272+(nr-22.)*939.5654+jet[0].E()); }

inline void Control::AnalyseBasics() {
	histo[0]->Fill(lep[0].E());
	histo[1]->Fill(lep[0].Theta()*180./PI);
	histo[2]->Fill(jet[0].E());
	histo[3]->Fill(jet[0].Theta()*180./PI);
	double p_sep=0;
	for(int i=0; i<hadNum; i++) p_sep+=had[i]->Had().P()*sin(had[i]->Had().Angle(jet[0].Vect()));
	histo[4]->Fill(atan(p_sep/jet[0].P())*180./PI);
}

inline bool Control::DetermineIfPositiveEnergy() {
	if(lep[0].E()>0 && jet[0].E()>0) return true;
	return false;
}

inline TLorentzVector Control::Smear(TLorentzVector par, double energy_resolution, double theta_r, double phi_r) {
    double E_k = par.E()-par.M();
    double smearedE = par.M()+random3.Gaus(E_k, E_k*energy_resolution);
    double smearedP = sqrt(smearedE*smearedE-par.M2());
    double theta = par.Theta();
    double phi = par.Phi();
    double smearedX = smearedP*(sin(theta_r)*(cos(phi_r)*cos(theta)*cos(phi)-sin(phi_r)*sin(phi))+cos(theta_r)*sin(theta)*cos(phi));
	double smearedY = smearedP*(sin(theta_r)*(cos(phi_r)*cos(theta)*sin(phi)+sin(phi_r)*cos(phi))+cos(theta_r)*sin(theta)*sin(phi));
	double smearedZ = smearedP*(cos(theta_r)*cos(theta)-sin(theta_r)*cos(phi_r)*sin(theta));
	return TLorentzVector(smearedX, smearedY, smearedZ, smearedE);
}

inline void Control::CreateSmearedLepAndJet() {
	double theta_r = acos(1-abs(random3.Gaus(0, 1-cos(RES_ANG_L))));
	double phi_r = random3.Uniform(0, 2.*PI);
    TLorentzVector lepC = Smear(lep[0], RES_E_LC, theta_r, phi_r);
    TLorentzVector lepE = Smear(lep[0], RES_E_LE, theta_r, phi_r);
	theta_r = acos(1-abs(random3.Gaus(0, 1-cos(RES_ANG_J1))));
	phi_r = random3.Uniform(0, 2.*PI);
    TLorentzVector jet1 = Smear(jet[0], RES_E_J, theta_r, phi_r);
	theta_r = acos(1-abs(random3.Gaus(0, 1-cos(RES_ANG_J3))));
    TLorentzVector jet3 = Smear(jet[0], RES_E_J, theta_r, phi_r);
	theta_r = acos(1-abs(random3.Gaus(0, 1-cos(RES_ANG_J5))));
    TLorentzVector jet5 = Smear(jet[0], RES_E_J, theta_r, phi_r);
	for(int anal=1; anal<NUM_ANAL; anal++) {
		if(anal<4) lep[anal] = lepC;
		else lep[anal] = lepE;
		if(anal%3==1) jet[anal] = jet1;
		else if(anal%3==2) jet[anal] = jet3;
		else jet[anal] = jet5;
	}
}

inline void Control::ProjectJet(TLorentzVector lep, TLorentzVector &jet) {
	TVector3 z(0,0,1);
	TVector3 n = z.Cross(lep.Vect());
	TVector3 jetP = jet.Vect()-n.Dot(jet.Vect())/n.Mag2()*n;
	jet.SetVect(jetP);
}

inline void Control::Analyse() {
	for(int anal=0; anal<NUM_ANAL; anal++) {
		double theta_l = lep[anal].Theta();
		if(theta_l > THETA_L_CUT) {
			cutPass[anal][2]++;
			double theta_j = jet[anal].Theta();
			if(theta_j > THETA_J_CUT) {
				cutPass[anal][3]++;
				if(theta_l < THETA_L_UPP) {
					cutPass[anal][4]++;
					double theta_plane = asin(abs(lep[anal].Px()*jet[anal].Py()-lep[anal].Py()*jet[anal].Px())/lep[anal].Pt()/jet[anal].P());
					if(anal>0) histo[4+anal]->Fill(theta_plane*180/PI);
					if(theta_plane < THETA_CUT1) {
						cutPass[anal][5]++;
						ProjectJet(lep[anal], jet[anal]);
						double E_tot = lep[anal].E()+jet[anal].E();
						if(E_tot>TOTAL_ENERGY_THRESHOLD && E_tot<TOTAL_ENERGY_UPP) {
							double E_nu = lep[anal].E()/2*(1+cos(theta_l)+sin(theta_l)*(1.+cos(theta_j))/sin(theta_j));
							double p_nu = lep[anal].Pz()+lep[anal].Pt()/tan(theta_j);
							histo[SIGNAL+anal]->Fill(E_nu);
							histo[18+anal]->Fill(p_nu);
							histo[25+anal]->Fill(E_tot);
						}
					}
				}
			}
		}
	}
}

inline void Control::AnalyseNoCut() {
	for(int anal=0; anal<NUM_ANAL; anal++) {
		double theta_l = lep[anal].Theta();
		double theta_j = jet[anal].Theta();
		ProjectJet(lep[anal], jet[anal]);
		double E_nu = lep[anal].E()/2*(1+cos(theta_l)+sin(theta_l)*(1.+cos(theta_j))/sin(theta_j));
		histo[anal+SIGNAL+NOCUT]->Fill(E_nu);
	}
}

inline void Control::ShowEventInformation(const int INPUT_E, int anal) const {
	const double SEARCH_ENERGY_THRESHOLD = INPUT_E*5000; // Print Event information of energy over the threshold
	double theta_l = lep[anal].Theta();
	double theta_j = jet[anal].Theta();
	double theta_plane = asin(abs(lep[anal].Px()*jet[anal].Py()-lep[anal].Py()*jet[anal].Px())/lep[anal].Pt()/jet[anal].P());
	if(theta_l>THETA_L_CUT && theta_j>THETA_J_CUT && theta_plane<THETA_CUT1 && theta_l<THETA_L_UPP) {
		double E_nu = lep[anal].E()/2*(1+cos(theta_l)+sin(theta_l)*(1.+cos(theta_j))/sin(theta_j));
		if(E_nu > SEARCH_ENERGY_THRESHOLD) {
			cout<<"E_nu: "<<E_nu<<" ------------------------------------------------"<<endl;
			cout<<"E_l: "<<lep[anal].E()<<endl;
			cout<<"theta_l: "<<theta_l/PI*180.<<endl;
			cout<<"E_j: "<<jet[anal].E()<<endl;
			cout<<"theta_j: "<<theta_j/PI*180.<<endl;
			cout<<"theta_plane: "<<theta_plane/PI*180.<<endl;
			cout<<"hadNum: "<<hadNum<<endl<<endl;
			for(int num=0; num<hadNum; num++) had[num]->ShowHadronInfo();
			cout<<"----------------------------------------------------------------"<<endl;
		}	
	}
}





inline void DetermineBinnedCutPass() {
	for(int anal=0; anal<NUM_ANAL; anal++) {
		double N = 0;
		for(Int_t i=1; i<histo[anal+11]->GetSize()-1; i++) N += histo[anal+11]->GetBinContent(i);
		cutPass[anal][6] = N;
	}
}

inline void PrintAllCutPass() {
	cout<<"<Event Numbers passing Cuts>"<<endl;
	for(int anal=0; anal<NUM_ANAL; anal++) {
		for(int cut=0; cut<NUM_CUT+1; cut++) cout<<cutPass[anal][cut]<<'\t';
		cout<<endl;
	}
}

inline void WriteHistogramsToDataFiles() {
	ofstream histoFile;
	char filename[100];
		for(int histInd=0; histInd<NUM_HISTO; histInd++) {
		sprintf(filename, "hist_%s.dat", histo[histInd]->GetName());
		histoFile.open(filename, ofstream::trunc);
		for(Int_t idx=1; idx<(histo[histInd]->GetSize())-1; idx++) { histoFile << histo[histInd]->GetBinLowEdge(idx) << "\t" << histo[histInd]->GetBinContent(idx) << "\n"; }
		histoFile.close();
	}
}

inline void PrintMeansSDsSigmas(const int INPUT_E) {
	cout<<INPUT_E<<"GeV\tMean\tSD\t1sigma"<<endl;
	for(int histInd=SIGNAL; histInd<NUM_HISTO-7; histInd++) {
		double N=0, Sum_E=0, Sum_E2=0;
		for(Int_t idx=1; idx<histo[histInd]->GetSize()-2; idx++) {
			N += histo[histInd]->GetBinContent(idx);
			Sum_E += (histo[histInd]->GetBinCenter(idx))*(histo[histInd]->GetBinContent(idx));
			Sum_E2 += (histo[histInd]->GetBinCenter(idx))*(histo[histInd]->GetBinCenter(idx))*(histo[histInd]->GetBinContent(idx));
		}
		double Mean = Sum_E/N;
		double Centre = INPUT_E*1000;
		double SD = sqrt(Sum_E2/N-Sum_E*Sum_E/N/N);
		double N1s = histo[histInd]->GetBinContent((Centre-SD)/100+1)/100.*(histo[histInd]->GetBinLowEdge((Centre-SD)/100+2)-(Centre-SD))+histo[histInd]->GetBinContent((Centre+SD)/100+1)/100.*((Centre+SD)-histo[histInd]->GetBinLowEdge((Centre+SD)/100+1));
		for(Int_t idx=(Centre-SD)/100+2; idx<(Centre+SD)/100; idx++) N1s += histo[histInd]->GetBinContent(idx);
		cout<<histo[histInd]->GetName()<<'\t'<<Mean<<'\t'<<SD<<'\t'<<N1s<<endl;
	}
}

inline void PrintMeansSDsSigmasNoCut(const int INPUT_E, const int FLAVOUR) {
	cout<<INPUT_E<<"GeV_nocut\tMean\tSD\t1sigma"<<endl;
	double SDref[NUM_FLAVOUR][NUM_ANAL_E][NUM_ANAL] = {{{408.552, 469.242, 490.258, 520.742, 742.919, 757.693, 774.798},
		{473.803, 561.115, 598.562, 651.938, 909.923, 930.133, 964.058},
		{529.854, 644.041, 701.362, 783.789, 1071.36, 1104.59, 1160.38},
		{585.914, 728.679, 808.085, 919.624, 1238.61, 1284.14, 1358.16},
		{642.117, 812.721, 923.045, 1068.66, 1405.83, 1469.93, 1561.49},
		{695.688, 892.465, 1030.58, 1218.51, 1568.73, 1650.19, 1770.68},
		{741.774, 977.295, 1146.31, 1368.94, 1730.55, 1826.37, 1970.37},
		{791.832, 1058.24, 1271.98, 1539, 1892.61, 2014.85, 2185.83},
		{838.468, 1144.59, 1393.04, 1711.71, 2056.38, 2199.14, 2406.47},
		{887.892, 1231.1, 1521.28, 1891, 2219.89, 2384.28, 2625.84},
		{932.761, 1313.17, 1651.97, 2067.64, 2381.3, 2572.96, 2847.38},
		{978.308, 1397.46, 1792.67, 2257.84, 2553.37, 2777.32, 3087.54}},

		{{412.998, 474.169, 496.672, 525.052, 745.321, 757.945, 777.028},
		{463.911, 552.459, 590.972, 641.547, 908.392, 931.429, 964.564},
		{523.6, 634.697, 687.628, 765.438, 1073.32, 1104.15, 1151.74},
		{576.057, 717.98, 792.703, 894.455, 1233.24, 1278.08, 1340.61},
		{623.121, 796.84, 892.485, 1028.39, 1393.58, 1449.86, 1538.7},
		{667.771, 871.32, 995.175, 1163.29, 1555.83, 1624.98, 1731.8},
		{714.294, 958.129, 1114.52, 1317.49, 1715.61, 1804.22, 1935.99},
		{758.697, 1035.32, 1219.45, 1460.63, 1879.71, 1982.32, 2133.65},
		{798.576, 1114.21, 1333.31, 1613.54, 2038.9, 2161.04, 2340.18},
		{837.152, 1193.69, 1447.03, 1770.51, 2203.45, 2344.84, 2550.68},
		{876.266, 1271.19, 1562.05, 1932.17, 2366.07, 2528.78, 2763.39},
		{914.334, 1352.54, 1679.27, 2092.02, 2528.68, 2708.29, 2972.84}},
		
		{{410.101, 469.655, 491.865, 526.25, 727.452, 740.444, 762.711},
		{473.969, 557.438, 595.448, 651.22, 895.054, 916.792, 955.294},
		{534.671, 642.516, 703.729, 786.389, 1057.24, 1093.83, 1148.68},
		{586.859, 723.135, 808.388, 923.867, 1219.33, 1270.11, 1344.03},
		{643.702, 811.25, 919.322, 1081.18, 1379.2, 1444.21, 1546.97},
		{695.985, 895.009, 1035.57, 1225.4, 1545.52, 1627.74, 1750.9},
		{746.16, 979.68, 1158.18, 1392.09, 1710.3, 1815.21, 1973.5},
		{793.476, 1056.13, 1273.33, 1549.06, 1872.32, 1998.61, 2180.23},
		{837.395, 1139.5, 1397.68, 1728.58, 2033.35, 2184.29, 2402.34},
		{885.953, 1226.47, 1522.65, 1892.8, 2201.65, 2369.56, 2617.32},
		{925.549, 1308.52, 1656.85, 2079.39, 2363.2, 2561.38, 2842.61},
		{973.047, 1390.69, 1785.81, 2264.76, 2522.48, 2747.8, 3072.58}},

		{{407.499, 465.482, 488.95, 515.694, 722.679, 736.54, 756.034},
		{467.584, 555.552, 588.417, 640.815, 890.729, 911.959, 944.166},
		{516.477, 628.29, 682.263, 759.209, 1050.12, 1082.33, 1131.17},
		{574.641, 711.282, 785.383, 889.592, 1213.01, 1258.52, 1325.14},
		{621.903, 790.673, 888.842, 1030.22, 1370.92, 1427.83, 1519.84},
		{667.097, 871.733, 996.482, 1163.83, 1533.79, 1606.75, 1716.06},
		{712.251, 949.055, 1102.73, 1311.32, 1698.86, 1786.53, 1920.36},
		{755.501, 1027.65, 1209.81, 1456.64, 1854.6, 1956.81, 2114.57},
		{796.833, 1109.42, 1325.74, 1612.28, 2015.58, 2138.44, 2321.05},
		{836.487, 1188.81, 1443.22, 1766.08, 2177.74, 2321.28, 2528.16},
		{872.755, 1267.56, 1558.04, 1928.02, 2342.8, 2503.91, 2746.62},
		{908.821, 1343.15, 1670.53, 2088.96, 2503.52, 2685.08, 2953.45}}};
	for(int histInd=SIGNAL+NOCUT; histInd<NUM_HISTO; histInd++) {
		double N=0, Sum_E=0, Sum_E2=0;
		for(Int_t idx=1; idx<histo[histInd]->GetSize()-2; idx++) {
			N += histo[histInd]->GetBinContent(idx);
			Sum_E += (histo[histInd]->GetBinCenter(idx))*(histo[histInd]->GetBinContent(idx));
			Sum_E2 += (histo[histInd]->GetBinCenter(idx))*(histo[histInd]->GetBinCenter(idx))*(histo[histInd]->GetBinContent(idx));
		}
		double Mean = Sum_E/N;
		double Centre = INPUT_E*1000;
		double SD = 0;
		if(FLAVOUR == ELECTRON_NEUTRINO) SD = SDref[0][INPUT_E-4][histInd-SIGNAL-NOCUT];
		else if(FLAVOUR == ANTI_ELECTRON_NEUTRINO) SD = SDref[1][INPUT_E-4][histInd-SIGNAL-NOCUT];
		else if(FLAVOUR == MUON_NEUTRINO) SD = SDref[2][INPUT_E-4][histInd-SIGNAL-NOCUT];
		else if(FLAVOUR == ANTI_MUON_NEUTRINO) SD = SDref[3][INPUT_E-4][histInd-SIGNAL-NOCUT];
		double N1s = histo[histInd]->GetBinContent((Centre-SD)/100+1)/100.*(histo[histInd]->GetBinLowEdge((Centre-SD)/100+2)-(Centre-SD))+histo[histInd]->GetBinContent((Centre+SD)/100+1)/100.*((Centre+SD)-histo[histInd]->GetBinLowEdge((Centre+SD)/100+1));
		for(Int_t idx=(Centre-SD)/100+2; idx<(Centre+SD)/100; idx++) N1s += histo[histInd]->GetBinContent(idx);
		SD = sqrt(Sum_E2/N-Sum_E*Sum_E/N/N);
		cout<<histo[histInd]->GetName()<<'\t'<<Mean<<'\t'<<SD<<'\t'<<N1s<<endl;
	}
}