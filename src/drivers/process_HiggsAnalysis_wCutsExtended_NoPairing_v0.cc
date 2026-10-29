#include "Muon.h"
#include "Electron.h"
#include "TrigObj.h"
#include "Jet.h"
#include "FsrPhoton.h"
#include "Vertex.h"
#include "Event.h"
#include "MET.h"
#include "Flag.h"
#include "Constants.h"
#include "Streamer.h"
#include "MetaHiggs.h"
#include "defs/Macros.h"

//	ROOT headers
#include "TFile.h"
#include "TChain.h"
#include "TLeaf.h"
#include "TString.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "LumiReweightingStandAlone.h"
#include "HistogramSets.h"
#include "TSystem.h"
#include "correction.h"			// JSON reader for correction files - Not Golden JSON.

//#include "FWCore/FWLite/interface/FWLiteEnabler.h"

#include "boost/program_options.hpp"
#include <signal.h>
#include <Progbar.h>
#include "Corrections.h"
#define FILL_JETS(set) \
	set.hDiJetMass->Fill(dijetmass, weight); \
	set.hDiJetdeta->Fill(std::abs(deta), weight); \
	set.hDiMuonpt->Fill(p4dimuon.Pt(), weight); \
	set.hDiMuonMass->Fill(p4dimuon.M(), weight); \
	set.hDiMuoneta->Fill(p4dimuon.Eta(), weight); \
	set.hDiMuondphi->Fill(dphi, weight); \
	set.hMuonpt->Fill(p4m1.Pt(), weight); \
	set.hMuonpt->Fill(p4m2.Pt(), weight); \
	set.hMuoneta->Fill(p4m1.Eta(), weight); \
	set.hMuoneta->Fill(p4m2.Eta(), weight); \
	set.hMuonphi->Fill(p4m1.Phi(), weight); \
	set.hMuonphi->Fill(p4m2.Phi(), weight);

#define FILL_NOJETS(set) \
	set.hDiMuonpt->Fill(p4dimuon.Pt(), weight); \
	set.hDiMuonMass->Fill(p4dimuon.M(), weight); \
	set.hDiMuoneta->Fill(p4dimuon.Eta(), weight); \
	set.hDiMuondphi->Fill(dphi, weight); \
	set.hMuonpt->Fill(p4m1.Pt(), weight); \
	set.hMuonpt->Fill(p4m2.Pt(), weight); \
	set.hMuoneta->Fill(p4m1.Eta(), weight); \
	set.hMuoneta->Fill(p4m2.Eta(), weight); \
	set.hMuonphi->Fill(p4m1.Phi(), weight); \
	set.hMuonphi->Fill(p4m2.Phi(), weight);

/*	Declare/Define all the service globals
 */
std::string __inputfilename;
std::string __outputfilename;
std::string __efffilename;
bool __isMC;
bool __isCondor;
double __lumi;
double __sigma_bkg;
double __sumw;
bool __genPUMC;
std::string __puMCfilename;
std::string __puDATAfilename;
std::atomic<bool> __continueRunning = true;

/*
 * Define all the Constants
 * for Muon:
 * lead - id 0
 * sublead - is 1
 */
enum MuonType {kLead=0, kSubLead=1};
double _muonMatchedPt = 26.;
double _muonMatchedEta = 2.4;
double _leadmuonPt = 26.;
double _subleadmuonPt = 20.;
double _leadmuonEta = 2.4;
double _subleadmuonEta = 2.4;
double _leadmuonIso = 0.25;
double _subleadmuonIso = 0.25;
double _muondxy = 0.5;
double _muondz = 1.0;
double _muon_dR = 0.1;
double _leadJetPt = 35.;
double _subleadJetPt = 25.;
double _metPt = 40.;
double _dijetMass_VBFTight = 650;
double _dijetdEta_VBFTight = 3.5;
double _dijetMass_ggFTight = 400.;
double _dimuonPt_ggFTight = 50.;
double _dimuonPt_01JetsTight = 10.;

double _muon_dR2 = _muon_dR * _muon_dR;

//std::string const NTUPLEMAKER_NAME = "ntuplemaker_H2DiMuonMaker";
std::string const NTUPLEMAKER_NAME = "";

namespace po = boost::program_options;
using namespace analysis::core;
using namespace analysis::dimuon;
using namespace analysis::processing;

TH1D *hEventWeights = NULL;
DimuonSet setNoVertexCuts("NoVertexCuts");			//
DimuonSet setNoMuonCuts("NoMuonCuts");				//
DimuonSet setNoCats("NoCats");
DimuonSet setnFsrZero("nFsrZero");
DimuonSet setnFsrNonzero("nFsrNonzero");
DimuonSet setAllJets("AllJets");
DimuonSet set2Jets("2Jets");
DimuonSet setVBFTight("VBFTight");
DimuonSet setggFTight("ggFTight");
DimuonSet setggFLoose("ggFLoose");

DimuonSet set01Jets("01Jets");
DimuonSet set01JetsTight("01JetsTight");
DimuonSet set01JetsLoose("01JetsLoose");
DimuonSet set01JetsBB("01JetsBB");
DimuonSet set01JetsBO("01JetsBO");
DimuonSet set01JetsBE("01JetsBE");
DimuonSet set01JetsOO("01JetsOO");
DimuonSet set01JetsOE("01JetsOE");
DimuonSet set01JetsEE("01JetsEE");
DimuonSet set01JetsTightBB("01JetsTightBB");
DimuonSet set01JetsTightBO("01JetsTightBO");
DimuonSet set01JetsTightBE("01JetsTightBE");
DimuonSet set01JetsTightOO("01JetsTightOO");
DimuonSet set01JetsTightOE("01JetsTightOE");
DimuonSet set01JetsTightEE("01JetsTightEE");
DimuonSet set01JetsLooseBB("01JetsLooseBB");
DimuonSet set01JetsLooseBO("01JetsLooseBO");
DimuonSet set01JetsLooseBE("01JetsLooseBE");
DimuonSet set01JetsLooseOO("01JetsLooseOO");
DimuonSet set01JetsLooseOE("01JetsLooseOE");
DimuonSet set01JetsLooseEE("01JetsLooseEE");

const double ptBins[] = {15, 30, 60, 120, 250, 500, 1000};
const double etaBins[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5};
int nPtBins = sizeof(ptBins) / sizeof(ptBins[0]) - 1;
int nEtaBins = sizeof(etaBins) / sizeof(etaBins[0]) - 1;

TH2F *h_b_loose = NULL;
TH2F *h_c_loose = NULL;
TH2F *h_l_loose = NULL;
TH2F *h_b_medium = NULL;
TH2F *h_c_medium = NULL;
TH2F *h_l_medium = NULL;
TH2F *h_b_all = NULL;
TH2F *h_c_all = NULL;
TH2F *h_l_all = NULL;

pair<int,int> year_era(const std::string fName=""){
	std::string yr = fName;
	size_t posmc = yr.find("mc_");
	size_t posdata = yr.find("data_");
	size_t poseff = yr.find("eff_");
	if (posmc != std::string::npos){
		yr.erase(0,posmc+3);
	}
	else if(posdata != std::string::npos){
		yr.erase(0,posdata+5);
	}
	else if(poseff != std::string::npos){
		yr.erase(0,poseff+4);
	}
	yr.erase(6,-1);
	std::string era = yr;
	yr.erase(4,6);
	era.erase(0,5);
	try{
		return std::make_pair(stoi(yr),stoi(era));
	}
	catch(std::invalid_argument &e1){
		try{
			return std::make_pair(stoi(yr),1);
		}
		catch(std::invalid_argument &e2){
			throw std::runtime_error("Invalid year on --input filename.");
		}
	}
}

bool isBarrel(Muon const& m)
{
	return std::abs(m._eta)<0.8;
}

bool isOverlap(Muon const& m)
{
	return std::abs(m._eta)>=0.8 && std::abs(m._eta)<1.6;
}

bool isEndcap(Muon const& m)
{
	return std::abs(m._eta)>=1.6 && std::abs(m._eta)<2.1;
}

bool passVertex(Vertices* v)
{
/*	for (Vertices::const_iterator it=v->begin(); it!=v->end() && n<20; ++it){
		if (std::abs(it->_z)<24 && it->_ndf>4) return true;
	}
*/
	if (!v || v->empty()) return false;
	return v->front()._npvsGood >= 1;
}

bool passMuon(Muon const& m, MuonType id=kLead)
{
	float pt = m._ptcorr;
	float eta = std::abs(m._eta);
	float dxy = std::abs(m._dxy);
	float dz = std::abs(m._dz);

	if (id==kLead)
	{
		if ((m._isGlobal || m._isTracker) &&
				pt > _leadmuonPt && eta < _leadmuonEta &&
				m._isMedium && m._pfRelIso04_all <= _leadmuonIso &&
				dxy < _muondxy && dz < _muondz)
			return true;
	}
	else
	{
		if ((m._isGlobal || m._isTracker) &&
				pt > _subleadmuonPt && eta < _subleadmuonEta &&
				m._isMedium && m._pfRelIso04_all <= _subleadmuonIso &&
				dxy < _muondxy && dz < _muondz)
			return true;
	}
	return false;
}

bool useBeamConstraintPt(const Muon& m){
	if(!std::isfinite(m._bsConstrainedChi2)) return false;
	if(m._bsConstrainedChi2 >= 30.f) return false;

	if(!std::isfinite(m._bsConstrainedPt)) return false;
	if(m._bsConstrainedPt <= 0.f) return false;

	return true;
}

float baseMuonPtForCorrections(const Muon& m){
	return useBeamConstraintPt(m) ? m._bsConstrainedPt : m._pt;
}

inline float nLayers(const Muon& m){
	return (m._nL > 0 ? (float)m._nL : 12.0f);
}

float muCorrScale(const Muon& m, const vector<std::unique_ptr<correction::CorrectionSet>>& csets, bool isMC, uint32_t evt = 0, uint32_t lumi = 0){
	const double pt0 = baseMuonPtForCorrections(m);
	if(!std::isfinite(pt0) || pt0 <= 0.0) return 1.0f;
	const int q = (m._charge > 0 ? +1 : -1);
	const double eta = m._eta;
	const double phi = m._phi;
	const float nL = nLayers(m);

	float LowPtThreshold = 0.;

	double pt_scaled = pt0;
	for(auto &cset : csets){
		pt_scaled = pt_scale(!__isMC, pt_scaled, eta, phi, q, *cset, LowPtThreshold);
		if (!std::isfinite(pt_scaled) || pt_scaled <= 0.0) return 1.0f;
	}

	return static_cast<float>(pt_scaled / pt0);
}

double corrBSPt(Muon const& m, const vector<std::unique_ptr<correction::CorrectionSet>>& csets, bool isMC, uint32_t evt = 0, uint32_t lumi = 0){
	const float pt0 = baseMuonPtForCorrections(m);
	return pt0*muCorrScale(m, csets, isMC, evt, lumi);
}

bool passMuonHLT(TrigObjs* trigobjs, Muon const& m, /*const vector<std::unique_ptr<correction::CorrectionSet>>& csets,*/ bool isMC, uint32_t evt, uint32_t lumi){
//	if (/*(m._isHLTMatched[1] || m._isHLTMatched[0]) &&*/
//		m._pt>_muonMatchedPt && std::abs(m._eta)<_muonMatchedEta)
//		return true;
	if(!trigobjs) return false;
	if(trigobjs->size() == 0) return false;

//	double corrected_pt = corrBSPt(m, csets, isMC, evt, lumi);
	double corrected_pt = m._ptcorr;
	if(corrected_pt < _muonMatchedPt) return false;

	for(auto it = trigobjs->begin(); it != trigobjs->end(); ++it){
		if(std::abs(it->_id) != 13) continue;
		if((it->_filterBits & (1 << 3)) == 0 ) continue;
		float deta = m._eta - it->_eta;
		float dphi = std::fabs(m._phi - it->_phi);
		if (dphi > M_PI) dphi = 2.0f * M_PI - dphi;
		float dR2 = deta * deta + dphi * dphi;
		if (dR2 < _muon_dR2) return true;	// Within loops, sqrt is computationally intensive. Therefore, the comparison is done through dR^2.
	}
	return false;
}

bool passIsoMu24Leg(TrigObjs *trigobjs, const Muon& mu, bool isMC, uint32_t evt, uint32_t lumi){
	if(!mu._isMedium)
		return false;

	if(mu._pfRelIso04_all >= 0.20)
		return false;

	if(std::abs(mu._eta) >= _muonMatchedEta)
		return false;

	return passMuonHLT(trigobjs, mu, isMC, evt, lumi);
}

bool passMuons(Muon const& m1, Muon const& m2, TrigObjs* trigobjs, /*const vector<std::unique_ptr<correction::CorrectionSet>>& csets, */bool isMC, uint32_t evt, uint32_t lumi)
{
	float pt1_corr = m1._ptcorr; //corrBSPt(m1, csets, isMC, evt, lumi);
	float pt2_corr = m2._ptcorr; //corrBSPt(m2, csets, isMC, evt, lumi);
	Muon const& leadmuon = pt1_corr>=pt2_corr ? m1 : m2;
	Muon const& subleadmuon = pt1_corr<=pt2_corr ? m1 : m2;
	float lead_pt_corr      = std::max(pt1_corr, pt2_corr);
	float sublead_pt_corr   = std::min(pt1_corr, pt2_corr);

	bool OS = m1._charge*m2._charge == -1;	// Opposite signed muons
	bool passLead = passMuon(leadmuon, kLead);
	bool passSublead = passMuon(subleadmuon, kSubLead);
	bool passHLT = passIsoMu24Leg(trigobjs, m1, /*csets, */isMC, evt, lumi) || passIsoMu24Leg(trigobjs, m2, /*csets, */isMC, evt, lumi);

	return OS && passLead && passSublead && passHLT;
}

bool vetoExtraLeptons(Electrons* E, Muons* M, bool isMC, uint32_t evt = 0, uint32_t lumi = 0)		// exclude events with certain extra leptons
{
	if(!M || M->empty()) return true;	// This safeguard is a bit overkill, since the analysis strictly uses the SingleMuon trigger(s).
	int Nmu = 0;
	for(const auto& m : *M){
		if(m._ptcorr>20 && std::abs(m._eta)<2.4 && m._isMedium && m._pfRelIso04_all<0.25 && std::abs(m._dxy)<0.5 && std::abs(m._dz)<1) Nmu++;
	}
	if(Nmu != 2) return true;

	if(!E || E->empty()) return false;	// No electron, no reason to veto.
	for(const auto& e : *E) if(e._pt > 20 && std::abs(e._eta) < 2.5 && e._mvaIso_WP90) return true;
	return false;
}

float jetMuondR(float jeta,float jphi, float meta, float mphi)
{
	TLorentzVector p4j,p4m;
	p4j.SetPtEtaPhiM(10, jeta, jphi, 0);
	p4m.SetPtEtaPhiM(10, meta, mphi, 0);
	return p4j.DeltaR(p4m);
}

UChar_t getJetId2024(float eta, float chHEF, float neHEF, float chEmEF, float neEmEF, float muEF, int chMult, int neMult){
	const float aeta = std::abs(eta);
	const int nConstituents = chMult + neMult;
	bool tight = false;
	bool tightLepVeto = false;

	if (aeta <= 2.6) {
		tight =
		neHEF < 0.99 &&
		neEmEF < 0.90 &&
		nConstituents > 1 &&
		chHEF > 0.01 &&
		chMult > 0;

		tightLepVeto =
		tight &&
		muEF < 0.80 &&
		chEmEF < 0.80;
	}

	else if (aeta <= 2.7) {
		tight =
		neHEF < 0.90 &&
		neEmEF < 0.99;

		tightLepVeto =
		tight &&
		muEF < 0.80;
	}

	else if (aeta <= 3.0) {
		tight = neHEF < 0.99;
		tightLepVeto = tight;
	}

	else if (aeta < 5.0) {
		tight = neEmEF < 0.40;
		tightLepVeto = tight;
	}

	UChar_t jetId = 0;

	if (tight) jetId |= static_cast<UChar_t>(2);
	if (tightLepVeto) jetId |= static_cast<UChar_t>(4);

	return jetId;
}

void categorize(TrigObjs* trigobjs, Jets* jets, Muon const& mu1, Muon const& mu2, FsrPhotons* fsrphotons, Vertex const& vertex, MET const& met, Event const& event,
	TFile* effroot, bool isgenjet, const correction::CorrectionSet& cset_btag, pair<int,int> Y_ERA, bool isMC=false, float scale=1., int im=-1, int ievent=0, float genweight=1.)
{
	// Muons
	TLorentzVector p4m1, p4m2, p4fsr1, p4fsr2, p4mfsr1, p4mfsr2;
	float weight = scale * genweight;
	int year = Y_ERA.first;
	int era = Y_ERA.second;

	p4m1.SetPtEtaPhiM(mu1._ptcorr, mu1._eta, mu1._phi, PDG_MASS_Mu);
	p4m2.SetPtEtaPhiM(mu2._ptcorr, mu2._eta, mu2._phi, PDG_MASS_Mu);

	FsrPhoton fsr1;
	FsrPhoton fsr2;
	int idx1 = mu1._fsrPhotonIdx;
	int idx2 = mu2._fsrPhotonIdx;

	if(idx1 > -1) fsr1 = fsrphotons->at(idx1);
	if(idx2 > -1) fsr2 = fsrphotons->at(idx2);

	p4fsr1.SetPtEtaPhiM(fsr1._pt, fsr1._eta, fsr1._phi, 0);
	p4fsr2.SetPtEtaPhiM(fsr2._pt, fsr2._eta, fsr2._phi, 0);

	p4mfsr1 = p4m1 + p4fsr1;
	p4mfsr2 = p4m2 + p4fsr2;

	TLorentzVector p4dimuon = p4m1 + p4m2;
	TLorentzVector p4dimuonfsr = p4mfsr1 + p4mfsr2;

	// Jets
	std::vector<TLorentzVector> p4jets;
	std::vector<Jet> jets_selected;
	int nJet_selected = 0;
	for (Jets::const_iterator it=jets->begin(); it!=jets->end(); ++it)
	{
		if (it->_ptcorr>25 && std::abs(it->_eta)<4.7 && (it->_jetId & 2) != 0)
		{
			if (jetMuondR(it->_eta, it->_phi, mu1._eta, mu1._phi)>0.4 &&
				jetMuondR(it->_eta, it->_phi, mu2._eta, mu2._phi)>0.4)
			{
				TLorentzVector p4;
				p4.SetPtEtaPhiM(it->_ptcorr, it->_eta, it->_phi, it->_mass);
				p4jets.push_back(p4);
				jets_selected.push_back(*it);
				nJet_selected++;

				// Print out the jet's properties
				/* std::cout << "Jet passing selection: " << std::endl;
				 s *td::cout << "  pT:  " << it->_pt << " GeV" << std::endl;
				 std::cout << "  eta: " << it->_eta << std::endl;
				 std::cout << "  phi: " << it->_phi << std::endl;
				 std::cout << "  mass: " << it->_mass << " GeV" << std::endl;*/
			}
		}
	}

	double p4m1_pt = p4m1.Pt();
	double p4m2_pt = p4m2.Pt();
	double p4fsr1_pt = p4mfsr1.Pt();
	double p4fsr2_pt = p4mfsr2.Pt();
	double p4m1_eta = p4m1.Eta();
	double p4m2_eta = p4m2.Eta();
	double p4fsr1_eta = p4mfsr1.Eta();
	double p4fsr2_eta = p4mfsr2.Eta();
	double p4m1_phi = p4m1.Phi();
	double p4m2_phi = p4m2.Phi();
	double p4fsr1_phi = p4mfsr1.Phi();
	double p4fsr2_phi = p4mfsr2.Phi();

	double mass_dimuon = p4dimuon.M();
	double mass_dimuon_fsr = p4dimuonfsr.M();
	double y_dimuon = p4dimuon.Rapidity();
	double y_dimuon_fsr = p4dimuonfsr.Rapidity();
	double p4m1_div_m_mumu = p4m1_pt/mass_dimuon;
	double p4m2_div_m_mumu = p4m2_pt/mass_dimuon;
	double p4m1_div_m_mumu_fsr = p4fsr1_pt/mass_dimuon_fsr;
	double p4m2_div_m_mumu_fsr = p4fsr2_pt/mass_dimuon_fsr;

	double p4jet_pt_l = 0;
	double p4jet_eta_l = 0;
	double p4jet_phi_l = 0;
	double p4jet_pt_s = 0;
	double p4jet_eta_s = 0;
	double p4jet_phi_s = 0;
	double zepp = -1;

	if(nJet_selected>0){
		double yj1 = p4jets[0].Rapidity();
		p4jet_pt_l = p4jets[0].Pt();
		p4jet_eta_l = p4jets[0].Eta();
		p4jet_phi_l = p4jets[0].Phi();
	}
	else if(nJet_selected>1){
		double yj1 = p4jets[0].Rapidity();
		double yj2 = p4jets[1].Rapidity();
		p4jet_pt_s = p4jets[1].Pt();
		p4jet_eta_s = p4jets[1].Eta();
		p4jet_phi_s = p4jets[1].Phi();
		zepp = (y_dimuon - (yj1 + yj2)/2)/std::abs(yj1 - yj2);
	}

	// Collins - Sopper variables
	double dphi = p4m1.DeltaPhi(p4m2);
	bool mu1_isNegative = (mu1._charge < 0);
	TLorentzVector muMinus = mu1_isNegative ? p4m1 : p4m2;
	TLorentzVector muPlus = mu1_isNegative ? p4m2 : p4m1;

	TVector3 B = -p4dimuon.BoostVector();

	TLorentzVector muMinus_b = muMinus;
	TLorentzVector muPlus_b = muPlus;
	muMinus_b.Boost(B);
	muPlus_b.Boost(B);

	TLorentzVector beam1(0., 0. ,1., 1.);
	TLorentzVector beam2(0., 0. ,-1., 1.);
	beam1.Boost(B);
	beam2.Boost(B);

	TVector3 b1 = beam1.Vect().Unit();
	TVector3 b2 = beam2.Vect().Unit();
	TVector3 zcs = (b1 - b2).Unit();
	TVector3 xcs = (b1 + b2).Unit();
	TVector3 ycs = zcs.Cross(xcs).Unit();

	TVector3 muMinus_b_hat = muMinus_b.Vect().Unit();

	double cosThetaCS = muMinus_b_hat.Dot(zcs);
	double phiCS = std::atan2(muMinus_b_hat.Dot(ycs), muMinus_b_hat.Dot(xcs)); // [-π, π]

	//	Fill the No Categorization Set
	setNoCats.hDiMuonpt->Fill(p4dimuon.Pt(), weight);
	setNoCats.hDiMuonMass_full->Fill(mass_dimuon, weight);
	setNoCats.hDiMuonMass_Jpsi->Fill(mass_dimuon, weight);
	setNoCats.hDiMuonMass_Upsi->Fill(mass_dimuon, weight);
	setNoCats.hDiMuonMass->Fill(mass_dimuon, weight);
	setNoCats.hDiMuonMass_Higgs->Fill(mass_dimuon, weight);
	setNoCats.hDiMuonMass_Higgs2->Fill(mass_dimuon, weight);
	setNoCats.hDiMuonMass_Higgs25->Fill(mass_dimuon, weight);
	setNoCats.hDiMuonMass_Higgs3->Fill(mass_dimuon, weight);
	setNoCats.hDiMuonMass_Higgs35->Fill(mass_dimuon, weight);
	setNoCats.hDiMuonMassFSR->Fill(mass_dimuon_fsr, weight);
	setNoCats.hDiMuonMassFSR_Higgs->Fill(mass_dimuon_fsr, weight);
	setNoCats.hDiMuonMassFSR_Higgs2->Fill(mass_dimuon_fsr, weight);
	setNoCats.hDiMuonMassFSR_Higgs25->Fill(mass_dimuon_fsr, weight);
	setNoCats.hDiMuonMassFSR_Higgs3->Fill(mass_dimuon_fsr, weight);
	setNoCats.hDiMuonMassFSR_Higgs35->Fill(mass_dimuon_fsr, weight);
	setNoCats.hDiMuoneta->Fill(p4dimuon.Eta(), weight);
	setNoCats.hDiMuondphi->Fill(dphi, weight);
	setNoCats.hMuonpt_l->Fill(p4m1_pt, weight);
	setNoCats.hMuonpt_s->Fill(p4m2_pt, weight);
	setNoCats.hMuoneta->Fill(p4m1_eta, weight);
	setNoCats.hMuoneta->Fill(p4m2_eta, weight);
	setNoCats.hMuonphi->Fill(p4m1_phi, weight);
	setNoCats.hMuonphi->Fill(p4m1_phi, weight);
	setNoCats.hDiMuonpteta->Fill(p4dimuon.Eta(), p4dimuon.Pt(), weight);
	setNoCats.hMuonpteta->Fill(p4m1_eta, p4m1_pt, weight);
	setNoCats.hMuonpteta->Fill(p4m2_eta, p4m2_pt, weight);
	setNoCats.hDiMuonptM->Fill(mass_dimuon, p4dimuon.Pt(), weight);
	setNoCats.hYvsX->Fill(vertex._x,vertex._y,weight);
	setNoCats.hZ->Fill(vertex._z,weight);
	setNoCats.hMuonptpt->Fill(p4m1_pt, p4m2_pt, weight);

	if(fsrphotons->size()==0)
	{
		setnFsrZero.hDiMuonpt->Fill(p4dimuon.Pt(), weight);
		setnFsrZero.hDiMuonMass->Fill(mass_dimuon, weight);
		setnFsrZero.hDiMuonMass_full->Fill(mass_dimuon, weight);
		setnFsrZero.hDiMuonMass_Higgs->Fill(mass_dimuon, weight);
		setnFsrZero.hDiMuonMass_Higgs2->Fill(mass_dimuon, weight);
		setnFsrZero.hDiMuonMass_Higgs25->Fill(mass_dimuon, weight);
		setnFsrZero.hDiMuonMass_Higgs3->Fill(mass_dimuon, weight);
		setnFsrZero.hDiMuonMass_Higgs35->Fill(mass_dimuon, weight);
		setnFsrZero.hDiMuonMassFSR->Fill(mass_dimuon_fsr, weight);
		setnFsrZero.hDiMuonMassFSR_Higgs->Fill(mass_dimuon_fsr, weight);
		setnFsrZero.hDiMuonMassFSR_Higgs2->Fill(mass_dimuon_fsr, weight);
		setnFsrZero.hDiMuonMassFSR_Higgs25->Fill(mass_dimuon_fsr, weight);
		setnFsrZero.hDiMuonMassFSR_Higgs3->Fill(mass_dimuon_fsr, weight);
		setnFsrZero.hDiMuonMassFSR_Higgs35->Fill(mass_dimuon_fsr, weight);
		setnFsrZero.hDiMuoneta->Fill(p4dimuon.Eta(), weight);
		setnFsrZero.hDiMuondphi->Fill(dphi, weight);
		setnFsrZero.hMuonpt_l->Fill(p4m1_pt, weight);
		setnFsrZero.hMuonpt_s->Fill(p4m2_pt, weight);
		setnFsrZero.hMuoneta->Fill(p4m1_eta, weight);
		setnFsrZero.hMuoneta->Fill(p4m2_eta, weight);
		setnFsrZero.hMuonphi->Fill(p4m1_phi, weight);
		setnFsrZero.hMuonphi->Fill(p4m1_phi, weight);
		setnFsrZero.hDiMuonpteta->Fill(p4dimuon.Eta(), p4dimuon.Pt(), weight);
		setnFsrZero.hMuonpteta->Fill(p4m1_eta, p4m1_pt, weight);
		setnFsrZero.hMuonpteta->Fill(p4m2_eta, p4m2_pt, weight);
		setnFsrZero.hDiMuonptM->Fill(mass_dimuon, p4dimuon.Pt(), weight);
	}

	else
	{
		setnFsrNonzero.hDiMuonpt->Fill(p4dimuon.Pt(), weight);
		setnFsrNonzero.hDiMuonMass->Fill(mass_dimuon, weight);
		setnFsrNonzero.hDiMuonMass_full->Fill(mass_dimuon, weight);
		setnFsrNonzero.hDiMuonMass_Higgs->Fill(mass_dimuon, weight);
		setnFsrNonzero.hDiMuonMass_Higgs2->Fill(mass_dimuon, weight);
		setnFsrNonzero.hDiMuonMass_Higgs25->Fill(mass_dimuon, weight);
		setnFsrNonzero.hDiMuonMass_Higgs3->Fill(mass_dimuon, weight);
		setnFsrNonzero.hDiMuonMass_Higgs35->Fill(mass_dimuon, weight);
		setnFsrNonzero.hDiMuonMassFSR->Fill(mass_dimuon_fsr, weight);
		setnFsrNonzero.hDiMuonMassFSR_Higgs->Fill(mass_dimuon_fsr, weight);
		setnFsrNonzero.hDiMuonMassFSR_Higgs2->Fill(mass_dimuon_fsr, weight);
		setnFsrNonzero.hDiMuonMassFSR_Higgs25->Fill(mass_dimuon_fsr, weight);
		setnFsrNonzero.hDiMuonMassFSR_Higgs3->Fill(mass_dimuon_fsr, weight);
		setnFsrNonzero.hDiMuonMassFSR_Higgs35->Fill(mass_dimuon_fsr, weight);
		setnFsrNonzero.hDiMuoneta->Fill(p4dimuon.Eta(), weight);
		setnFsrNonzero.hDiMuondphi->Fill(dphi, weight);
		setnFsrNonzero.hMuonpt_l->Fill(p4m1_pt, weight);
		setnFsrNonzero.hMuonpt_s->Fill(p4m2_pt, weight);
		setnFsrNonzero.hMuoneta->Fill(p4m1_eta, weight);
		setnFsrNonzero.hMuoneta->Fill(p4m2_eta, weight);
		setnFsrNonzero.hMuonphi->Fill(p4m1_phi, weight);
		setnFsrNonzero.hMuonphi->Fill(p4m1_phi, weight);
		setnFsrNonzero.hDiMuonpteta->Fill(p4dimuon.Eta(), p4dimuon.Pt(), weight);
		setnFsrNonzero.hMuonpteta->Fill(p4m1_eta, p4m1_pt, weight);
		setnFsrNonzero.hMuonpteta->Fill(p4m2_eta, p4m2_pt, weight);
		setnFsrNonzero.hDiMuonptM->Fill(mass_dimuon, p4dimuon.Pt(), weight);
	}

//	ADDITIONAL THINGS TO PUT IN...
//	Isolation Cuts
//	Track/Calorimeter Isolation
//	Impact Paraimeter Cuts
//	Transverse Impact Parameter (d_0)
//	Event Topology Cuts
//	Number of Jets
//	MET
//	Muon ID and Trigger Cuts
//	Muon Quality Cuts --> Ex. # of hits in tracker/muon system
//	Triggers --> Ex. Single-Muon Triggers
//	Background Rejection Cuts
//	Drell-Yan Processes
//	Pileup Mitigation

//	Filling efficiency maps
	double btag_m;
	double btag_l;
	if(year==2022){
		if(era==1){
			btag_l = 0.0470f;			// 2022
			btag_m = 0.2450f;			// 2022
		}
		else if(era==2){
			btag_l = 0.0499f;			// 2022EE
			btag_m = 0.2605f;			// 2022EE
		}
		else{
			throw std::runtime_error("Invalid era (B-tag WPs).");
		}
	}
	else if(year==2023){
		if(era==1){
			btag_l = 0.0358f;			// 2023
			btag_m = 0.1917f;			// 2023
		}
		else if(era==2){
			btag_l = 0.0359f;			// 2023 BPix
			btag_m = 0.1919f;			// 2023 BPix
		}
		else{
			throw std::runtime_error("Invalid era (B-tag WPs).");
		}
	}
	else if(year>2023){
		if(era==1){
			btag_l = 0.0246f;			// 2024 || 2025
			btag_m = 0.1272f;			// 2024 || 2025
		}
		else{
			throw std::runtime_error("Invalid era (B-tag WPs).");
		}
	}

	if(!effroot && __isMC){
		for(auto jet : jets_selected){
			if(std::abs(jet._eta)<2.5){
				if(jet._hadronFlavour == 5){
					h_b_all->Fill(jet._ptcorr, std::abs(jet._eta), weight);
					if(jet._btag>=btag_m) h_b_medium->Fill(jet._ptcorr, std::abs(jet._eta), weight);
					if(jet._btag>=btag_l) h_b_loose->Fill(jet._ptcorr, std::abs(jet._eta), weight);
				}
				if(jet._hadronFlavour == 4){
					h_c_all->Fill(jet._ptcorr, std::abs(jet._eta), weight);
					if(jet._btag>=btag_m) h_c_medium->Fill(jet._ptcorr, std::abs(jet._eta), weight);
					if(jet._btag>=btag_l) h_c_loose->Fill(jet._ptcorr, std::abs(jet._eta), weight);
				}
				if(jet._hadronFlavour == 0){
					h_l_all->Fill(jet._ptcorr, std::abs(jet._eta), weight);
					if(jet._btag>=btag_m) h_l_medium->Fill(jet._ptcorr, std::abs(jet._eta), weight);
					if(jet._btag>=btag_l) h_l_loose->Fill(jet._ptcorr, std::abs(jet._eta), weight);
				}
			}
		}
		return;
	}

//	JET B-TAG CUTS
	int nb_loose=0;
	for(auto jet : jets_selected){
		if(std::abs(jet._eta)<2.5){
			if(jet._btag>=btag_m) return;			// Medium b-tag cut
			if(jet._btag>=btag_l) nb_loose++;		// Loose b-tag cut
			if(nb_loose>1) return;
		}
	}

	if(effroot && __isMC){
		// Initialize corrections here, if necessary
//		for(auto flav : {0, 4, 5, 3}) for(int i=-1;i<2;i++) for(auto wp : {"L", "M"}) for(auto pt : {30., 98.}) cout<<flav<<" "<<i<<" "<<wp<<" "<<pt<<" "<<BTAG(wp, flav, pt, eta, cset_btag)[i+1]<<endl;
		auto efficiencies = EffFinder(effroot);
		double P0_MC=1., P0_Data=1.;
		double P1_MC=0., P1_Data=0.;
		for(auto jet : jets_selected){
			if(std::abs(jet._eta)<2.5){
				double eff_l, eff_m, sf_l, sf_m, p0_mc, p0_data, p1_mc, p1_data;
				if(jet._hadronFlavour==5){
					int bin_l = efficiencies[0]->FindBin(jet._ptcorr, std::abs(jet._eta));
					int bin_m = efficiencies[1]->FindBin(jet._ptcorr, std::abs(jet._eta));
					eff_l = efficiencies[0]->GetBinContent(bin_l);
					eff_m = efficiencies[1]->GetBinContent(bin_m);
					sf_l = BTAG("L", 5, jet._ptcorr, jet._eta, cset_btag, year)[0];	// [0]: Nominal, [1] and [2] are up and down, respectively.
					sf_m = BTAG("M", 5, jet._ptcorr, jet._eta, cset_btag, year)[0];	// [0]: Nominal, [1] and [2] are up and down, respectively.
				}
				else if(jet._hadronFlavour==4){
					int bin_l = efficiencies[2]->FindBin(jet._ptcorr, std::abs(jet._eta));
					int bin_m = efficiencies[3]->FindBin(jet._ptcorr, std::abs(jet._eta));
					eff_l = efficiencies[2]->GetBinContent(bin_l);
					eff_m = efficiencies[3]->GetBinContent(bin_m);
					sf_l = BTAG("L", 4, jet._ptcorr, jet._eta, cset_btag, year)[0];	// [0]: Nominal, [1] and [2] are up and down, respectively.
					sf_m = BTAG("M", 4, jet._ptcorr, jet._eta, cset_btag, year)[0];	// [0]: Nominal, [1] and [2] are up and down, respectively.
				}
				else if(jet._hadronFlavour==0){
					int bin_l = efficiencies[4]->FindBin(jet._ptcorr, std::abs(jet._eta));
					int bin_m = efficiencies[5]->FindBin(jet._ptcorr, std::abs(jet._eta));
					eff_l = efficiencies[4]->GetBinContent(bin_l);
					eff_m = efficiencies[5]->GetBinContent(bin_m);
					sf_l = BTAG("L", 0, jet._ptcorr, jet._eta, cset_btag, year)[0];	// [0]: Nominal, [1] and [2] are up and down, respectively.
					sf_m = BTAG("M", 0, jet._ptcorr, jet._eta, cset_btag, year)[0];	// [0]: Nominal, [1] and [2] are up and down, respectively.
				}
				else{
					return;
				}
				eff_l = std::clamp(eff_l, 0., 1.);
				eff_m = std::clamp(eff_m, 0., eff_l);
				p0_mc = 1 - eff_l;					// 0 tagged, MC
				p0_data = 1 - sf_l*eff_l;			// 0 tagged, Data
				p1_mc = eff_l - eff_m;				// 1 tagged, MC
				p1_data = eff_l*sf_l - eff_m*sf_m;	// 1 tagged, Data

				P1_MC = P0_MC*p1_mc + P1_MC*p0_mc;
				P1_Data = P0_Data*p1_data + P1_Data*p0_data;
				P0_MC *= p0_mc;
				P0_Data *= p0_data;
			}
		}
		double w_btag = (P0_Data + P1_Data)/(P0_MC + P1_MC);
		if(finite(w_btag)) weight *= w_btag;
	}

	if(p4jets.size()<2)
	{
		for(auto hSet : {set01Jets, setAllJets}){
			hSet.hDiMuonpt->Fill(p4dimuonfsr.Pt(), weight);
			hSet.hDiMuonMass->Fill(mass_dimuon_fsr, weight);
			hSet.hDiMuonMass_full->Fill(mass_dimuon_fsr, weight);
			hSet.hDiMuonMass_Higgs->Fill(mass_dimuon_fsr, weight);
			hSet.hDiMuoneta->Fill(p4dimuonfsr.Eta(), weight);
			hSet.hDiMuondphi->Fill(dphi, weight);
			hSet.hMuonpt_l->Fill(p4m1_pt, weight);
			hSet.hMuonpt_s->Fill(p4m2_pt, weight);
			hSet.hMuoneta->Fill(p4m1_eta, weight);
			hSet.hMuoneta->Fill(p4m2_eta, weight);
			hSet.hMuonphi->Fill(p4m1_phi, weight);
			hSet.hMuonphi->Fill(p4m1_phi, weight);
			hSet.hMuonptpt->Fill(p4m1_pt, p4m2_pt, weight);
		}
	}

	else if (p4jets.size()>=2)
	{
		TLorentzVector p4lead; //= p4jets[0];
		TLorentzVector p4sub; //= p4jets[1];
		for(auto p4jet : p4jets){
			if(p4jet.Pt() > p4lead.Pt()){
				p4sub = p4lead;
				p4lead = p4jet;
			}
			else if(p4jet.Pt() > p4sub.Pt()){
				p4sub = p4jet;
			}
		}

		for(auto p4jet : p4jets) if(p4sub.Pt() <= p4jet.Pt() && p4jet.Pt() < p4lead.Pt()) p4sub = p4jet;

		TLorentzVector dijet = p4lead + p4sub;
		float detajj = std::abs(p4lead.Eta() - p4sub.Eta());
		float dphijj = std::abs(p4lead.DeltaPhi(p4sub));
		float dijetmass = dijet.M();
		float dijeteta = dijet.Eta();
		double deta_mm_jj = std::abs(p4dimuon.Eta() - dijet.Eta());
		double dphi_mm_jj = std::abs(p4dimuon.DeltaPhi(dijet));
		if (p4lead.Pt()<35 || p4sub.Pt()<25 || dijetmass<400 || detajj<2.5)
		{
			for(auto hSet : {set2Jets, setAllJets}){
				hSet.hDiJetMass->Fill(dijetmass, weight);
				hSet.hDiJeteta->Fill(dijeteta, weight);
				hSet.hDiJetdeta->Fill(detajj, weight);
				hSet.hDiJetPt->Fill(dijet.Pt(), weight);
				hSet.hDiMuonpt->Fill(p4dimuonfsr.Pt(), weight);
				hSet.hDiMuonMass->Fill(mass_dimuon_fsr, weight);
				hSet.hDiMuonMass_full->Fill(mass_dimuon_fsr, weight);
				hSet.hDiMuonMass_Higgs->Fill(mass_dimuon_fsr, weight);
				hSet.hDiMuoneta->Fill(p4dimuonfsr.Eta(), weight);
				hSet.hDiMuondphi->Fill(dphi, weight);
				hSet.hMuonpt_l->Fill(p4m1_pt, weight);
				hSet.hMuonpt_s->Fill(p4m2_pt, weight);
				hSet.hMuoneta->Fill(p4m1_eta, weight);
				hSet.hMuoneta->Fill(p4m2_eta, weight);
				hSet.hMuonphi->Fill(p4m1_phi, weight);
				hSet.hMuonphi->Fill(p4m1_phi, weight);
				hSet.hMuonptpt->Fill(p4m1_pt, p4m2_pt, weight);
			}
		}
		//else cout<<"nJ>=2 Fail pTlead, Mjj, |deta|:"<<p4lead.Pt()<<", "<<dijetmass<<", "<<detajj<<endl;
	}
	return;
}

float sampleinfo(std::string const& inputname)
{
	Streamer s(inputname, NTUPLEMAKER_NAME+"/Meta");
	s.chainup();

	using namespace analysis::dimuon;
	MetaHiggs *meta=NULL;
	s._chain->SetBranchAddress("Meta", &meta);

	long long int numEvents = 0;
	long long int numEventsWeighted = 0;
	for (int i=0; i<s._chain->GetEntries(); i++)
	{
		s._chain->GetEntry(i);
		numEvents+=meta->_nEventsProcessed;
		numEventsWeighted+=meta->_sumEventWeights;
	}
	std::cout 
		<< "#events processed total = " << numEvents << std::endl
		<< "#events weighted total = " << numEventsWeighted << std::endl;

	return numEventsWeighted;
}

void generatePUMC()
{
	std::cout << "### Generate PU MC file...." << std::endl;
	TFile *pufile = new TFile(__puMCfilename.c_str(), "recreate");
	//TH1D *h = new TH1D("pileup", "pileup", 50, 0, 50);
	TH1D *h = new TH1D("pileup", "pileup", 99, 0, 99);

	Streamer s(__inputfilename, NTUPLEMAKER_NAME+"/Events");
	s.chainup();

	EventAuxiliary *aux=NULL;
	s._chain->SetBranchAddress("EventAuxiliary", &aux);
	int numEvents = s._chain->GetEntries();
	for (uint32_t i=0; i<numEvents; i++)
	{
		s._chain->GetEntry(i);
		h->Fill(aux->_nPU, aux->_genWeight);
	}

	pufile->Write();
	pufile->Close();
}

void process()
{
	double Norm = __lumi*__sigma_bkg/__sumw;
	TFile *outroot = new TFile(__outputfilename.c_str(), "recreate");
	TFile *effroot = TFile::Open(__efffilename.c_str(), "r");
	outroot->cd();
	setNoVertexCuts.init();
	setNoMuonCuts.init();
	setNoCats.init();
	setnFsrZero.init();
	setnFsrNonzero.init();
	setAllJets.init();
	set2Jets.init();
	set01Jets.init();
	setVBFTight.init();
	setggFTight.init();
	setggFLoose.init();
	set01JetsTight.init();
	set01JetsLoose.init();
	set01JetsBB.init();
	set01JetsBO.init();
	set01JetsBE.init();
	set01JetsOO.init();
	set01JetsOE.init();
	set01JetsEE.init();
	set01JetsTightBB.init();
	set01JetsTightBO.init();
	set01JetsTightBE.init();
	set01JetsTightOO.init();
	set01JetsTightOE.init();
	set01JetsTightEE.init();
	set01JetsLooseBB.init();
	set01JetsLooseBO.init();
	set01JetsLooseBE.init();
	set01JetsLooseOO.init();
	set01JetsLooseOE.init();
	set01JetsLooseEE.init();

	gDirectory->mkdir("Efficiencies");
	gDirectory->cd("Efficiencies");
	h_b_loose = new TH2F("h_b_loose", "h_b_loose", nPtBins, ptBins, nEtaBins, etaBins);		// Loose
	h_c_loose = new TH2F("h_c_loose", "h_c_loose", nPtBins, ptBins, nEtaBins, etaBins);		// Loose
	h_l_loose = new TH2F("h_l_loose", "h_l_loose", nPtBins, ptBins, nEtaBins, etaBins);		// Loose
	h_b_medium = new TH2F("h_b_medium", "h_b_medium", nPtBins, ptBins, nEtaBins, etaBins);	// Medium
	h_c_medium = new TH2F("h_c_medium", "h_c_medium", nPtBins, ptBins, nEtaBins, etaBins);	// Medium
	h_l_medium = new TH2F("h_l_medium", "h_l_medium", nPtBins, ptBins, nEtaBins, etaBins);	// Medium
	h_b_all = new TH2F("h_b_all", "h_b_all", nPtBins, ptBins, nEtaBins, etaBins);				// All
	h_c_all = new TH2F("h_c_all", "h_c_all", nPtBins, ptBins, nEtaBins, etaBins);				// All
	h_l_all = new TH2F("h_l_all", "h_l_all", nPtBins, ptBins, nEtaBins, etaBins);				// All
	h_b_loose->Sumw2();
	h_c_loose->Sumw2();
	h_l_loose->Sumw2();
	h_b_medium->Sumw2();
	h_c_medium->Sumw2();
	h_l_medium->Sumw2();
	h_b_all->Sumw2();
	h_c_all->Sumw2();
	h_l_all->Sumw2();
	gDirectory->cd("../");

	hEventWeights = new TH1D("eventWeights", "eventWeights", 100, -1., 1.);

	//	get the total events, etc...

	#if 0
	long long int numEventsWeighted = sampleinfo(__inputfilename);
	#endif

	/// TODO understand what this is
	/// how was filled
	/// how used downstream
	auto numEventsWeighted = 100;
	hEventWeights->Fill(0.5, numEventsWeighted);

	//	generate the MC Pileup histogram
	if (__genPUMC && __isMC)
		generatePUMC();

	Streamer streamer(__inputfilename, "Events");
	streamer.chainup();
	pair<int,int> Y = year_era(__inputfilename);
	int year = Y.first;
	int era = Y.second;

	Electrons *electrons=NULL;
	Muons *muons=NULL;
	Muons muons1;
	Muons muons2;
	TrigObjs *trigobjs=NULL;
	Jets *jets=NULL;
	GenJets *genjets=NULL;
	FsrPhotons *fsrphotons=NULL;
	Vertices *vertices=NULL;
	Event *event_=NULL;
	EventAuxiliary *aux=NULL;
	MET *met=NULL;

	#define ROOT_CHAIN streamer._chain

	bool intt_e=false, intt_m=false, intt_j=false, intt_gj=false, intt_f=false, intt_p=false, intt_t=false, intt_npvs=false, intt_npvsGood=false;
	Int_t nElectron_s, nMuon_s, nJet_s, nGenJet_s, nFsrPhoton_s, nTrigObj_s, PV_npvs_s, PV_npvsGood_s;
	UInt_t nElectron_u, nMuon_u, nJet_u, nGenJet_u, nFsrPhoton_u, nTrigObj_u;
	UChar_t PV_npvs_u, PV_npvsGood_u;

	if(strcmp(streamer._chain->GetLeaf("nElectron")->GetTypeName(),"Int_t")==0){
		streamer._chain->SetBranchAddress("nElectron", &nElectron_s);
		intt_e=true;
	}
	else{
		streamer._chain->SetBranchAddress("nElectron", &nElectron_u);
	}

	if(strcmp(streamer._chain->GetLeaf("nMuon")->GetTypeName(),"Int_t")==0){
		streamer._chain->SetBranchAddress("nMuon", &nMuon_s);
		intt_m=true;
	}
	else{
		streamer._chain->SetBranchAddress("nMuon", &nMuon_u);
	}

	if(strcmp(streamer._chain->GetLeaf("nJet")->GetTypeName(),"Int_t")==0){
		streamer._chain->SetBranchAddress("nJet", &nJet_s);
		intt_j=true;
	}
	else{
		streamer._chain->SetBranchAddress("nJet", &nJet_u);
	}

	if(strcmp(streamer._chain->GetLeaf("nFsrPhoton")->GetTypeName(),"Int_t")==0){
		streamer._chain->SetBranchAddress("nFsrPhoton", &nFsrPhoton_s);
		intt_f=true;
	}
	else{
		streamer._chain->SetBranchAddress("nFsrPhoton", &nFsrPhoton_u);
	}

	if(strcmp(streamer._chain->GetLeaf("nTrigObj")->GetTypeName(),"Int_t")==0){
		streamer._chain->SetBranchAddress("nTrigObj", &nTrigObj_s);
		intt_t=true;
	}
	else{
		streamer._chain->SetBranchAddress("nTrigObj", &nTrigObj_u);
	}

	if(strcmp(streamer._chain->GetLeaf("PV_npvs")->GetTypeName(),"Int_t")==0){
		streamer._chain->SetBranchAddress("PV_npvs", &PV_npvs_s);
		intt_npvs=true;
	}
	else{
		streamer._chain->SetBranchAddress("PV_npvs", &PV_npvs_u);
	}

	if(strcmp(streamer._chain->GetLeaf("PV_npvsGood")->GetTypeName(),"Int_t")==0){
		streamer._chain->SetBranchAddress("PV_npvsGood", &PV_npvsGood_s);
		intt_npvsGood=true;
	}
	else{
		streamer._chain->SetBranchAddress("PV_npvsGood", &PV_npvsGood_u);
	}

	if(streamer._chain->GetBranch("nGenJet")){
		if(strcmp(streamer._chain->GetLeaf("nGenJet")->GetTypeName(),"Int_t")==0){
			streamer._chain->SetBranchAddress("nGenJet", &nGenJet_s);
			intt_gj=true;
		}
		else{
			streamer._chain->SetBranchAddress("nGenJet", &nGenJet_u);
		}
	}

	bool isGenJet = true;
	UChar_t Jet_jetId[100];
	Short_t Jet_genJetIdx[100];
	Int_t nGenJet = 0;
	Float_t genWeight=1., Pileup_nTrueInt=1.;
	Float_t GenJet_pt[100], GenJet_eta[100], GenJet_phi[100], GenJet_mass[100], Jet_btagUParTAK4B[100], Jet_btagPNetB[100];
	UChar_t Jet_hadronFlavour[100];

	if(streamer._chain->GetBranch("Jet_btagPNetB")){
		streamer._chain->SetBranchAddress("Jet_btagPNetB", &Jet_btagPNetB);
	}

	if(streamer._chain->GetBranch("Jet_btagUParTAK4B")){
		streamer._chain->SetBranchAddress("Jet_btagUParTAK4B", &Jet_btagUParTAK4B);
	}

	if(streamer._chain->GetBranch("Jet_jetId")){
		streamer._chain->SetBranchAddress("Jet_jetId", &Jet_jetId);
	}

	if(__isMC){
		streamer._chain->SetBranchAddress("genWeight", &genWeight);
		streamer._chain->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt);

		streamer._chain->SetBranchAddress("Jet_hadronFlavour", &Jet_hadronFlavour);
		if(streamer._chain->GetBranch("Jet_genJetIdx")){
			streamer._chain->SetBranchAddress("Jet_genJetIdx", &Jet_genJetIdx);
		}
		else{
			for(Short_t& j : Jet_genJetIdx) j=-2;
			isGenJet = false;
		}

		streamer._chain->SetBranchAddress("GenJet_pt", &GenJet_pt);
		streamer._chain->SetBranchAddress("GenJet_eta", &GenJet_eta);
		streamer._chain->SetBranchAddress("GenJet_phi", &GenJet_phi);
		streamer._chain->SetBranchAddress("GenJet_mass", &GenJet_mass);
	}
//	Jet_btagDeepFlavB
//	Jet_btagPNetB

	SET_BRANCH_ULONG64(event);
	SET_BRANCH_UINT(luminosityBlock);
	SET_BRANCH_UINT(run);

	SET_BRANCH_FLOAT(PV_ndof);
	SET_BRANCH_FLOAT(PV_x);
	SET_BRANCH_FLOAT(PV_y);
	SET_BRANCH_FLOAT(PV_z);

	SET_BRANCH_BOOL(Flag_goodVertices);
	SET_BRANCH_BOOL(Flag_globalSuperTightHalo2016Filter);
	SET_BRANCH_BOOL(Flag_EcalDeadCellTriggerPrimitiveFilter);
	SET_BRANCH_BOOL(Flag_BadPFMuonFilter);
	SET_BRANCH_BOOL(Flag_BadPFMuonDzFilter);
	SET_BRANCH_BOOL(Flag_hfNoisyHitsFilter);
	SET_BRANCH_BOOL(Flag_eeBadScFilter);
	SET_BRANCH_BOOL(Flag_ecalBadCalibFilter);

//	SET_BRANCH_FLOAT_ARRAY(TrigObj_pt);
	SET_BRANCH_FLOAT_ARRAY(TrigObj_eta);
	SET_BRANCH_FLOAT_ARRAY(TrigObj_phi);
	SET_BRANCH_USHORT_ARRAY(TrigObj_id);
	SET_BRANCH_INT_ARRAY(TrigObj_filterBits);

	SET_BRANCH_FLOAT_ARRAY(Electron_pt);
	SET_BRANCH_FLOAT_ARRAY(Electron_eta);
	//SET_BRANCH_FLOAT_ARRAY(Electron_phi);
	SET_BRANCH_BOOL_ARRAY(Electron_mvaIso_WP90);
	//SET_BRANCH_UCHAR_ARRAY(Electron_cutBased);

	SET_BRANCH_INT_ARRAY(Muon_nStations);
	SET_BRANCH_FLOAT_ARRAY(Muon_pt);
	SET_BRANCH_FLOAT_ARRAY(Muon_bsConstrainedPt);
	SET_BRANCH_FLOAT_ARRAY(Muon_bsConstrainedChi2);
	SET_BRANCH_FLOAT_ARRAY(Muon_eta);
	SET_BRANCH_FLOAT_ARRAY(Muon_phi);
	SET_BRANCH_FLOAT_ARRAY(Muon_dz);
	SET_BRANCH_FLOAT_ARRAY(Muon_dxy);
	SET_BRANCH_INT_ARRAY(Muon_charge);
	SET_BRANCH_UCHAR_ARRAY(Muon_nTrackerLayers);
	SET_BRANCH_BOOL_ARRAY(Muon_isGlobal);
	SET_BRANCH_BOOL_ARRAY(Muon_isPFcand);
	SET_BRANCH_BOOL_ARRAY(Muon_isTracker);
	SET_BRANCH_BOOL_ARRAY(Muon_isStandalone);
	SET_BRANCH_BOOL_ARRAY(Muon_mediumId);
	SET_BRANCH_SHORT_ARRAY(Muon_fsrPhotonIdx);
	SET_BRANCH_FLOAT_ARRAY(Muon_pfRelIso04_all);
	SET_BRANCH_BOOL(HLT_IsoMu24);
//	Float_t Muon_bsConstrainedPt[100], Muon_bsConstrainedChi2[100];
//	if(year > 2023){
//		streamer._chain->SetBranchAddress("Muon_bsConstrainedPt", &Muon_bsConstrainedPt);
//		streamer._chain->SetBranchAddress("Muon_bsConstrainedChi2", &Muon_bsConstrainedChi2);
//	}
//	else{
//		for(int i=0;i<100;i++){
//			Muon_bsConstrainedPt[i] = -1.;
//			Muon_bsConstrainedChi2[i] = 256.;
//		}
//	}

	//	SET_BRANCH_UCHAR_ARRAY(Muon_miniIsoId);
	//	SET_BRANCH_BOOL_ARRAY(Muon_tightId);

	SET_BRANCH_FLOAT(MET_pt);
	SET_BRANCH_FLOAT(MET_phi);
	SET_BRANCH_FLOAT(PuppiMET_pt);
	SET_BRANCH_FLOAT(PuppiMET_phi);

	SET_BRANCH_FLOAT_ARRAY(Jet_eta);
	SET_BRANCH_FLOAT_ARRAY(Jet_mass);
	SET_BRANCH_FLOAT_ARRAY(Jet_pt);
	SET_BRANCH_FLOAT_ARRAY(Jet_phi);
	SET_BRANCH_FLOAT_ARRAY(Jet_area);
	SET_BRANCH_FLOAT_ARRAY(Jet_rawFactor);
	SET_BRANCH_FLOAT_ARRAY(Jet_chHEF);
	SET_BRANCH_FLOAT_ARRAY(Jet_neHEF);
	SET_BRANCH_FLOAT_ARRAY(Jet_chEmEF);
	SET_BRANCH_FLOAT_ARRAY(Jet_neEmEF);
	SET_BRANCH_FLOAT_ARRAY(Jet_muEF);
	SET_BRANCH_INT_ARRAY(Jet_chMultiplicity);
	SET_BRANCH_INT_ARRAY(Jet_neMultiplicity);
	SET_BRANCH_FLOAT(Rho_fixedGridRhoFastjetAll);
//	SET_BRANCH_FLOAT_ARRAY(Jet_btagDeepFlavB);
//	SET_BRANCH_FLOAT_ARRAY(Jet_btagPNetB);
//	SET_BRANCH_FLOAT_ARRAY(Jet_btagUParTAK4B);

	SET_BRANCH_FLOAT_ARRAY(FsrPhoton_pt);
	SET_BRANCH_FLOAT_ARRAY(FsrPhoton_phi);
	SET_BRANCH_FLOAT_ARRAY(FsrPhoton_eta);
	SET_BRANCH_INT_ARRAY(FsrPhoton_muonIdx);

	trigobjs = new std::vector<TrigObj>();
	electrons = new std::vector<Electron>();
	muons = new std::vector<Muon>();
	jets = new std::vector<Jet>();
	genjets = new std::vector<GenJet>();
	fsrphotons = new std::vector<FsrPhoton>();
	vertices = new std::vector<Vertex>();
	met = new MET();

#if 0
	streamer._chain->SetBranchAddress("Muons", &muons);
	streamer._chain->SetBranchAddress("Jets", &jets);
	streamer._chain->SetBranchAddress("FsrPhotons", &fsrphotons);
	streamer._chain->SetBranchAddress("Photons", &photons);
	streamer._chain->SetBranchAddress("Vertices", &vertices);
	streamer._chain->SetBranchAddress("Event", &event_);
	streamer._chain->SetBranchAddress("EventAuxiliary", &aux);
	streamer._chain->SetBranchAddress("MET", &met);
#endif

	//	init the PU reweighter
	reweight::LumiReWeighting *weighter = NULL;
/*	if (__isMC)	// CAUSES SEGFAULT
	{
		std::cout << "mcPU=" << __puMCfilename << "  dataPU="
			<< __puDATAfilename << std::endl;
		TString mc_pileupfile = __puMCfilename.c_str();
		TString data_pileupfile = __puDATAfilename.c_str();
		weighter = new reweight::LumiReWeighting(
		mc_pileupfile.Data(), data_pileupfile.Data(), "pileup", "pileup");
	}
*/
	std::string epoch = to_string(year) + "_" + to_string(era) + "_";
	auto cset_BSC = correction::CorrectionSet::from_file(epoch + "schemaV2_BSC.json");			// Mixed
	auto cset_Regular = correction::CorrectionSet::from_file(epoch + "schemaV2_Regular.json");		// Mixed
	auto cset_muonZ = correction::CorrectionSet::from_file(epoch + "muon_Z.json.gz");				// Mixed
	auto cset_jerc = correction::CorrectionSet::from_file(epoch + "jet_jerc.json.gz");			// MC and Data
	//auto cset_jer = correction::CorrectionSet::from_file(epoch + "jer_smear.json.gz");			// MC only
	auto cset_puweight = correction::CorrectionSet::from_file(epoch + "puWeights.json.gz");			// MC only
	auto cset_btag = correction::CorrectionSet::from_file(epoch + "btagging.json.gz");			// MC only

	Json json_golden;
	std::vector<std::string> json_files;
	if(!__isMC){
		if(year==2022){
			if(era==1){
				json_files.push_back("2022_1_B_355100_355769_Muon.json");
				json_files.push_back("2022_1_C_355862_357482_Muon.json");
				json_files.push_back("2022_1_D_357538_357900_Muon.json");
			}
			else if(era==2){
				json_files.push_back("2022_2_E_359022_360331_Muon.json");
				json_files.push_back("2022_2_F_360390_362167_Muon.json");
				json_files.push_back("2022_2_G_362433_362760_Muon.json");
			}
		}
		else if(year==2023){
			if(era==1){
				json_files.push_back("2023_1_C_367095_368823_Muon.json");
			}
			else if(era==2){
				json_files.push_back("2023_2_D_369803_370790_Muon.json");
			}
		}
		else if(year==2024){
			json_files.push_back("2024_1_C_Muon.json");
			json_files.push_back("2024_1_D_Muon.json");
			json_files.push_back("2024_1_E_Muon.json");
			json_files.push_back("2024_1_F_Muon.json");
			json_files.push_back("2024_1_G_Muon.json");
			json_files.push_back("2024_1_H_Muon.json");
			json_files.push_back("2024_1_I_Muon.json");
		}

		for (const auto& fName : json_files){
			std::ifstream file(fName);

			if (!file.is_open()){
				std::cerr<<"Could not open "<<fName<<endl;
				continue;
			}

			Json temp;
			file >> temp;
			json_golden.update(temp);
		}
	}
	//	Event loop
//	std::vector<std::unique_ptr<correction::CorrectionSet>> csets;
//	csets.push_back(std::move(cset_BSC));
//	csets.push_back(std::move(cset_Regular));
	Int_t nElectron, nMuon, nJet, nFsrPhoton, nPhoton, nTrigObj, PV_npvs, PV_npvsGood;

	TFile *fNtuple = new TFile(("ntuple_"+__outputfilename).c_str(),"recreate");
	TTree *ntuple = new TTree;
	ntuple->SetName("ntuple");
	ntuple->Branch("event", &event, "event/l");
	ntuple->Branch("luminosityBlock", &luminosityBlock, "luminosityBlock/i");
	ntuple->Branch("run", &run, "run/i");

	ntuple->Branch("nMuon", &nMuon, "nMuon/I");
	ntuple->Branch("Muon_pt", Muon_pt, "Muon_pt[nMuon]/F");
	ntuple->Branch("Muon_eta", Muon_eta, "Muon_eta[nMuon]/F");
	ntuple->Branch("Muon_phi", Muon_phi, "Muon_phi[nMuon]/F");
	ntuple->Branch("Muon_charge", Muon_charge, "Muon_charge[nMuon]/I");
	ntuple->Branch("Muon_pfRelIso04_all", Muon_pfRelIso04_all, "Muon_pfRelIso04_all[nMuon]/I");
	ntuple->Branch("Muon_mediumId", Muon_mediumId, "Muon_mediumId[nMuon]/O");
	ntuple->Branch("HLT_IsoMu24", &HLT_IsoMu24, "HLT_IsoMu24/O");

	ntuple->Branch("nJet", &nJet, "nJet/i");
	ntuple->Branch("Jet_pt", Jet_pt, "Jet_pt[nJet]/F");
	ntuple->Branch("Jet_eta", Jet_eta, "Jet_eta[nJet]/F");
	ntuple->Branch("Jet_phi", Jet_phi, "Jet_phi[nJet]/F");
	ntuple->Branch("Jet_btagUParTAK4B", Jet_btagUParTAK4B, "Jet_btagUParTAK4B[nJet]/F");

	ntuple->Branch("MET_pt", &MET_pt, "MET_pt/F");
	ntuple->Branch("MET_phi", &MET_phi, "MET_phi/F");
	ntuple->Branch("PuppiMET_pt", &PuppiMET_pt, "PuppiMET_pt/F");
	ntuple->Branch("PuppiMET_phi", &PuppiMET_phi, "PuppiMET_phi/F");

	float rho = 0.;
	uint32_t numEntries = streamer._chain->GetEntries();
	for (uint32_t i=0; i<numEntries && __continueRunning; i++)
	{
//		if(i>100000) break;		// Limit number of events for testing purposes
		trigobjs->clear();
		electrons->clear();
		muons->clear();
		jets->clear();
		genjets->clear();
		fsrphotons->clear();
		vertices->clear();
		muons1.clear(); muons2.clear();
		streamer._chain->GetEntry(i);
		cout<<HLT_IsoMu24<<endl;
		ntuple->Fill();

//		Golden JSON
		if(!__isMC){
			const bool goldenPass = passGoldenJSON(json_golden, run, luminosityBlock);
			if (!goldenPass) continue;
		}
//		HOW TO RUN THROUGH THE ENTIRE GOLDEN JSON FILE?
//		for(auto json : {json_b,json_c}) for(auto run : json.items()) for(auto lumi : run.value()) cout<<run.key()<<" "<<lumi[0]<<" - "<<lumi[1]<<endl;

		rho = Rho_fixedGridRhoFastjetAll;

		if(intt_e){
			nElectron=nElectron_s;
		}
		else{
			nElectron=nElectron_u;
		}
		if(intt_m){
			nMuon=nMuon_s;
		}
		else{
			nMuon=nMuon_u;
		}
		if(intt_j){
			nJet=nJet_s;
		}
		else{
			nJet=nJet_u;
		}
		if(intt_f){
			nFsrPhoton=nFsrPhoton_s;
		}
		else{
			nFsrPhoton=nFsrPhoton_u;
		}
		if(intt_t){
			nTrigObj=nTrigObj_s;
		}
		else{
			nTrigObj=nTrigObj_u;
		}
		if(intt_npvs){
			PV_npvs=PV_npvs_s;
		}
		else{
			PV_npvs=PV_npvs_u;
		}
		if(intt_npvsGood){
			PV_npvsGood=PV_npvsGood_s;
		}
		else{
			PV_npvsGood=PV_npvsGood_u;
		}
		if(intt_gj){
			nGenJet=nGenJet_s;
		}
		else{
			nGenJet=nGenJet_u;
		}
		double sf=1., pu = 1.;

		for (size_t i=0; i<nElectron; i++) {
			Electron electron;
			electron._pt = Electron_pt[i];
			electron._eta = Electron_eta[i];
//			electron._phi = Electron_phi[i];
			electron._mvaIso_WP90 = Electron_mvaIso_WP90[i];
//			electron._cutBased = Electron_cutBased[i];
			electrons->push_back(std::move(electron));
		}

		for (size_t i=0; i<nMuon; i++) {
			Muon mu;
			mu._pt = Muon_pt[i];
			mu._ptcorr = mu._pt;
			mu._eta = Muon_eta[i];
			mu._phi = Muon_phi[i];
			mu._dz = Muon_dz[i];
			mu._dxy = Muon_dxy[i];
			mu._charge = Muon_charge[i];
			mu._nL = Muon_nTrackerLayers[i];
			mu._isGlobal = Muon_isGlobal[i];
			mu._isPFcand = Muon_isPFcand[i];
			mu._isTracker = Muon_isTracker[i];
			mu._isStandalone = Muon_isStandalone[i];
			mu._isMedium = Muon_mediumId[i];
			mu._nStations = Muon_nStations[i];
			mu._fsrPhotonIdx = Muon_fsrPhotonIdx[i];
			mu._pfRelIso04_all = Muon_pfRelIso04_all[i];
			mu._bsConstrainedPt = Muon_bsConstrainedPt[i];
			mu._bsConstrainedChi2 = Muon_bsConstrainedChi2[i];

			// KIT Corrections
			bool isBSconstr =
			std::isfinite(mu._bsConstrainedChi2) &&
			mu._bsConstrainedChi2 < 30.f &&
			std::isfinite(mu._bsConstrainedPt) &&
			mu._bsConstrainedPt > 0.f;

			if(isBSconstr){
				mu._ptcorr=pt_scale(!__isMC, mu._bsConstrainedPt, mu._eta, mu._phi, mu._charge, *cset_BSC, 0);
				if(__isMC) mu._ptcorr=pt_resol(mu._ptcorr, mu._eta, mu._phi, mu._nL, event, luminosityBlock, *cset_BSC, 0);
			}

			else{
				mu._ptcorr=pt_scale(!__isMC, mu._pt, mu._eta, mu._phi, mu._charge, *cset_Regular, 0);
				if(__isMC) mu._ptcorr=pt_resol(mu._ptcorr, mu._eta, mu._phi, mu._nL, event, luminosityBlock, *cset_Regular, 0);
			}
			muons->push_back(std::move(mu));
		}

		for (size_t i=0; i<nTrigObj; i++) {
			TrigObj trigobj;
//			trigobj._pt = TrigObj_pt[i];
			trigobj._eta = TrigObj_eta[i];
			trigobj._phi = TrigObj_phi[i];
			trigobj._id = TrigObj_id[i];
			trigobj._filterBits = TrigObj_filterBits[i];
			trigobjs->push_back(std::move(trigobj));
		}

		if(__isMC){
			for (size_t i=0; i<nGenJet; i++) {
				GenJet genjet;
				genjet._pt = GenJet_pt[i];
				genjet._eta = GenJet_eta[i];
				genjet._phi = GenJet_phi[i];
				genjet._mass = GenJet_mass[i];
				genjets->push_back(std::move(genjet));
			}
		}

		for (size_t i=0; i<nJet; i++) {
			Jet jet;
			pair<double, double> jer;

			jet._pt = Jet_pt[i];
			jet._ptcorr = jet._pt;
			jet._eta = Jet_eta[i];
			jet._phi = Jet_phi[i];
			jet._mass = Jet_mass[i];
			jet._area = Jet_area[i];
//			jet._btagDeepFlavB = Jet_btagDeepFlavB[i];
			jet._genJetIdx = Jet_genJetIdx[i];
			jet._hadronFlavour = Jet_hadronFlavour[i];
			jet._rawFactor = Jet_rawFactor[i];
			jet._chHEF = Jet_chHEF[i];
			jet._neHEF = Jet_neHEF[i];
			jet._chEmEF = Jet_chEmEF[i];
			jet._neEmEF = Jet_neEmEF[i];
			jet._muEF = Jet_muEF[i];
			jet._chMult = Jet_chMultiplicity[i];
			jet._neMult = Jet_neMultiplicity[i];
			if(year<2024){
				jet._btag = Jet_btagPNetB[i];
				jet._jetId = Jet_jetId[i];
			}
			else{
				jet._btag = Jet_btagUParTAK4B[i];
				jet._jetId = getJetId2024(jet._eta, jet._chHEF, jet._neHEF, jet._chEmEF, jet._neEmEF, jet._muEF, jet._chMult, jet._neMult);
			}

			std::string jecName;
			std::string jerScaleName;
			std::string jerResName;
			if(year==2022){
				if(era==1){
					jecName = __isMC ? "Summer22_22Sep2023_V4_MC_L1L2L3Res_AK4PFPuppi" : "Summer22_22Sep2023_V4_DATA_L1L2L3Res_AK4PFPuppi";
					jerScaleName = "Summer22_22Sep2023_JRV2_MC_ScaleFactor_AK4PFPuppi";
					jerResName = "Summer22_22Sep2023_JRV2_MC_PtResolution_AK4PFPuppi";
				}
				else if(era==2){
					jecName = __isMC ? "Summer22EE_22Sep2023_V4_MC_L1L2L3Res_AK4PFPuppi" : "Summer22EE_22Sep2023_V4_DATA_L1L2L3Res_AK4PFPuppi";
					jerScaleName = "Summer22EE_22Sep2023_JRV2_MC_ScaleFactor_AK4PFPuppi";
					jerResName = "Summer22EE_22Sep2023_JRV2_MC_PtResolution_AK4PFPuppi";
				}
			}
			else if(year==2023){
				if(era==1){
					jecName = __isMC ? "Summer23Prompt23_V4_MC_L1L2L3Res_AK4PFPuppi" : "Summer23Prompt23_V4_DATA_L1L2L3Res_AK4PFPuppi";
					jerScaleName = "Summer23Prompt23_RunCv1234_JRV2_MC_ScaleFactor_AK4PFPuppi";
					jerResName = "Summer23Prompt23_RunCv1234_JRV2_MC_PtResolution_AK4PFPuppi";
				}
				else if(era==2){
					jecName = __isMC ? "Summer23BPixPrompt23_V4_MC_L1L2L3Res_AK4PFPuppi" : "Summer23BPixPrompt23_V4_DATA_L1L2L3Res_AK4PFPuppi";
					jerScaleName = "Summer23BPixPrompt23_RunD_JRV2_MC_ScaleFactor_AK4PFPuppi";
					jerResName =   "Summer23BPixPrompt23_RunD_JRV2_MC_PtResolution_AK4PFPuppi";
				}
			}
			else if(year==2024){
				jecName = __isMC ? "Summer24Prompt24_V5_MC_L1L2L3Res_AK4PFPuppi" : "Summer24Prompt24_V5_DATA_L1L2L3Res_AK4PFPuppi";
				jerScaleName = "Summer24Prompt24_JRV2_MC_ScaleFactor_AK4PFPuppi";
				jerResName = "Summer24Prompt24_JRV2_MC_PtResolution_AK4PFPuppi";
			}

			double rawscale = (1 - jet._rawFactor);
			jet._mass *= rawscale;
			jet._ptcorr *=  rawscale;
			double jecfactor;
			try{
				jecfactor = __isMC ? JEC(jet._area, jet._eta, jet._ptcorr, rho, jet._phi, *cset_jerc, jecName, year, era) : JEC(jet._area, jet._eta, jet._ptcorr, rho, jet._phi, *cset_jerc, jecName, year, era, run);
			}
			catch(...){
				throw std::runtime_error("JEC error.");
			}
			jet._mass *= jecfactor;
			jet._ptcorr *= jecfactor;
			if(__isMC){
				try{
					jer = JER(jet._eta, jet._ptcorr, rho, *cset_jerc, jerScaleName, jerResName);
				}
				catch(...){
					throw std::runtime_error("JER error.");
				}
				double scale = jer.first;
				double resol = jer.second;
				bool hybridMatch = false;

				double pt_beforejer = jet._ptcorr;
				if(jet._genJetIdx>=0 && static_cast<std::size_t>(jet._genJetIdx) < genjets->size()){
					const auto& gen = genjets->at(jet._genJetIdx);

					const double dphi = std::remainder(jet._phi - gen._phi, 2.0*M_PI);
					const double deta = jet._eta - gen._eta;
					const double dR = std::hypot(deta, dphi);
					const double dPt = std::abs(jet._ptcorr - gen._pt);

					hybridMatch = dR < 0.2 && dPt < 3.0 * resol * jet._ptcorr;

					if(hybridMatch){
						const double smearFactor = std::max(0.0, 1.0 + (scale - 1) * (pt_beforejer - gen._pt) / pt_beforejer);
						jet._ptcorr = pt_beforejer * smearFactor;
					}
				}
				if(!hybridMatch){
					double r = deterministicGaussian(event, luminosityBlock, run, i);
					jet._ptcorr *= std::max(0.0, (1. + r*resol*TMath::Sqrt(std::max(scale*scale-1., 0.))));
				}
				jet._mass *= jet._ptcorr/pt_beforejer;
			}
			jets->push_back(std::move(jet));
		}

		for (size_t i=0; i<nFsrPhoton; i++){
			FsrPhoton fsrphoton;
			fsrphoton._pt = FsrPhoton_pt[i];
			fsrphoton._eta = FsrPhoton_eta[i];
			fsrphoton._phi = FsrPhoton_phi[i];
			fsrphoton._muonIdx = FsrPhoton_muonIdx[i];
			fsrphotons->push_back(std::move(fsrphoton));
		}

		Vertex vtx;
		vtx._ndf = PV_ndof;
		vtx._x = PV_x;
		vtx._y = PV_y;
		vtx._z = PV_z;
		vtx._npvs = PV_npvs;
		vtx._npvsGood = PV_npvsGood;
		vertices->push_back(std::move(vtx));

		Flag flg;
		flg._goodVertices = Flag_goodVertices;
		flg._globalSuperTightHalo2016Filter = Flag_globalSuperTightHalo2016Filter;
		flg._EcalDeadCellTriggerPrimitiveFilter = Flag_EcalDeadCellTriggerPrimitiveFilter;
		flg._BadPFMuonFilter = Flag_BadPFMuonFilter;
		flg._BadPFMuonDzFilter = Flag_BadPFMuonDzFilter;
		flg._hfNoisyHitsFilter = Flag_hfNoisyHitsFilter;
		flg._eeBadScFilter = Flag_eeBadScFilter;
		flg._ecalBadCalibFilter = Flag_ecalBadCalibFilter;

		Event event_;
		event_._run = run;
		event_._lumi = luminosityBlock;
		event_._event = event;

		met->_pt = MET_pt;

		std::string puName;
		if(year==2022){
			if(era==1){
				puName = "Collisions2022_355100_357900_eraBCD_GoldenJson";
			}
			else if(era==2){
				puName = "Collisions2022_359022_362760_eraEFG_GoldenJson";
			}
		}
		else if(year==2023){
			if(era==1){
				puName = "Collisions2023_366403_369802_eraBC_GoldenJson";
			}
			else if(era==2){
				puName = "Collisions2023_369803_370790_eraD_GoldenJson";
			}
		}
		else if(year==2024){
			puName = "Collisions24_CDEFGHI_goldenJSON";	// What?
		}

		try{
			if(__isMC) pu *= PU(Pileup_nTrueInt, *cset_puweight, puName);
		}
		catch(...){
			throw std::runtime_error("PU cset error");
		}

//		Old progress display
//		if (i%1000==0) std::cout << "### Event " << i << " / " << numEntries << std::endl;

		if(i==0) Progbar(0, numEntries, __inputfilename);
		if(__isCondor){
			if(((i+1)*8)%numEntries<8) Progbar(i, numEntries, __inputfilename);
		}
		else{
		if((i+1)%1000==0) Progbar(i, numEntries, __inputfilename);
		}
		if(i+1==numEntries) Progbar(numEntries-1, numEntries, __inputfilename);


//		Selections
//

		if(!HLT_IsoMu24) continue;

		if(!flg._goodVertices) continue;
		if(!flg._globalSuperTightHalo2016Filter) continue;
		if(!flg._EcalDeadCellTriggerPrimitiveFilter) continue;
		if(!flg._BadPFMuonFilter) continue;
		if(!flg._BadPFMuonDzFilter) continue;
		if(!flg._hfNoisyHitsFilter) continue;
		if(!flg._eeBadScFilter) continue;
		if(year>2023 && !flg._ecalBadCalibFilter) continue;

		int nFlag = 0;
		if(year<2024){
			if(run>362433 && run<367144){
				for(auto& jet : *jets){
					if(PuppiMET_pt<100) continue;
					if(jet._pt<50) continue;
					if(-0.5 > jet._eta || jet._eta > -0.1) continue;
					if(-2.1 > jet._phi || jet._phi > -1.8) continue;
					if(jet._neEmEF < 0.9 && jet._chEmEF < 0.9) continue;
					if(std::abs(PuppiMET_phi-jet._phi)) continue;
					nFlag++;
				}
			}
		}
		if(nFlag>0) continue;

//		setNoVertexCuts.hFsrPhotonpt->Fill(pt);
//		setNoVertexCuts.hFsrPhotoneta->Fill(eta);
//		setNoVertexCuts.hFsrPhotonphi->Fill(phi);
//		setNoVertexCuts.hFsrPhotonpteta->Fill(eta, pt);

		if(!passVertex(vertices)) continue;

//		setNoMuonCuts.hFsrPhotonpt->Fill(pt);
//		setNoMuonCuts.hFsrPhotoneta->Fill(eta);
//		setNoMuonCuts.hFsrPhotonphi->Fill(phi);
//		setNoMuonCuts.hFsrPhotonpteta->Fill(eta, pt);

		if(vetoExtraLeptons(electrons, muons, __isMC, event, luminosityBlock)) continue;


//		prepare the pairs of muons
//

		bool dve_dlinnyx_spichki = false;	// Писатель: Что я сам выбрал? ОДНУ ДЛИННУЮ СРИЧКУ ИЗ ДВУХ ДЛИННЫХ!?
		for (analysis::core::Muons::const_iterator it=muons->begin(); it!=muons->end(); ++it){
			for (analysis::core::Muons::const_iterator jt=(it+1); jt!=muons->end(); ++jt)
			{
				if(!passMuons(*it, *jt, trigobjs, /*csets,*/ __isMC, event, luminosityBlock)) continue;
				dve_dlinnyx_spichki = true;
				muons1.push_back(*it); muons2.push_back(*jt);
				break;
			}
			if(dve_dlinnyx_spichki) break;
		}
		if(!dve_dlinnyx_spichki) continue;

//		USE THIS FOR fsr photon identification
//		int idx = 0;
//		for(auto it = trigobjs->begin(); it != trigobjs->end(); ++it, ++idx){

		for (uint32_t im=0; im<muons1.size(); im++)
		{
			if(__isMC){
				bool match1 = passIsoMu24Leg(trigobjs, muons1.at(im), __isMC, event, luminosityBlock);
				bool match2 = passIsoMu24Leg(trigobjs, muons2.at(im), __isMC, event, luminosityBlock);
				const Muon *triggerMuon = nullptr;
				if(match1 && match2){
					triggerMuon = (muons1.at(im)._ptcorr >= muons2.at(im)._ptcorr) ? &muons1.at(im) : &muons2.at(im);
				}
				else if(match1){
					triggerMuon = &muons1.at(im);
				}
				else if(match2){
					triggerMuon = &muons2.at(im);
				}
				if(!triggerMuon) continue;
				try{
					sf *= muon_ID_ISO_SF(muons1.at(im)._ptcorr, muons1.at(im)._eta, "nominal", *cset_muonZ);
					sf *= muon_ID_ISO_SF(muons2.at(im)._ptcorr, muons2.at(im)._eta, "nominal", *cset_muonZ);
					sf *= muon_HLT_SF(triggerMuon->_ptcorr, triggerMuon->_eta, "nominal", *cset_muonZ);
				}
				catch(...){
					throw std::runtime_error("Muon Z error.");
				}
			}
			categorize(trigobjs, jets, muons1.at(im), muons2.at(im), fsrphotons, vtx, *met, event_,
				effroot, isGenJet, *cset_btag, Y, __isMC, pu*sf, im, i, genWeight);
		}
	} // Event loop ends here.
	ntuple->Write();
	fNtuple->Close();

//	MC Scaling
	for(auto set : {setNoCats, setnFsrZero, setnFsrNonzero, set01Jets, set2Jets, setAllJets}) set.scale(Norm);

// Closing the output file
	outroot->Write();
	outroot->Close();

	return;
}	// process() ends here.

void sigHandler(int sig)
{
	cout << "\n### Signal: " << sig << " caughter. Exiting...\n\n";
	__continueRunning = false;
}

void printCuts()
{
	std::cout << std::endl << "\033[4mCuts:				" << std::endl
		<< "_muonMatchedPt		= " << _muonMatchedPt << "	" << std::endl
		<< "_muonMatchedEta		= " << _muonMatchedEta << "	" << std::endl
		<< "_leadmuonPt		= " << _leadmuonPt << "	" << std::endl
		<< "_subleadmuonPt		= " << _subleadmuonPt << "	" << std::endl
		<< "_leadmuonEta		= " << _leadmuonEta << "	" << std::endl
		<< "_subleadmuonEta		= " << _subleadmuonEta << "	" << std::endl
		<< "_leadmuonIso		= " << _leadmuonIso << "	" << std::endl
		<< "_subleadmuonIso		= " << _subleadmuonIso << "	" << std::endl
		<< "_muondxy		= " << _muondxy << "	" << std::endl
		<< "_muondz			= " << _muondz << "	" << std::endl
		<< "_muon_dR		= " << _muon_dR << "	" << std::endl
		<< "_leadJetPt		= " << _leadJetPt << "	" << std::endl
		<< "_subleadJetPt		= " << _subleadJetPt << "	" << std::endl
		<< "_metPt			= " << _metPt << "	" << std::endl
		<< "_dijetMass_VBFTight	= " << _dijetMass_VBFTight << "	" << std::endl
		<< "_dijetdEta_VBFTight	= " << _dijetdEta_VBFTight << "	" << std::endl
		<< "_dijetMass_ggFTight	= " << _dijetMass_ggFTight << "	" << std::endl
		<< "_dimuonPt_ggFTight	= " << _dimuonPt_ggFTight << "	" << std::endl
		<< "_dimuonPt_01JetsTight	= " << _dimuonPt_01JetsTight << "	\033[0m" << std::endl << std::endl;
}

int main(int argc, char** argv)
{
	//gSystem->Load("libFWCoreFWLite");
	//FWLiteEnabler::enable();
	/*
	 *	Register signals
	 */
	signal(SIGABRT, &sigHandler);
	signal(SIGTERM, &sigHandler);
	signal(SIGINT, &sigHandler);

	std::string none;
	bool genPUMC = false;

	/*
	 *	Pare Options
	 */
	po::options_description desc("Allowed Program Options");
	desc.add_options()
		("help", "produce help messages")
		("input", po::value<std::string>(), "a file specifying all the ROOT files to process")
		("isMC", po::value<bool>(), "type of data: DATA vs MC")
		("isCondor", po::value<bool>()->default_value(true), "setting this false, limits the update amount of the progress bar to 10, preventing condor log files from bloating")
		("lumi", po::value<double>()->default_value(1.), "Luminosity value")
		("sigma_bkg", po::value<double>()->default_value(1.), "Cross section of the background process")
		("sumw", po::value<double>()->default_value(1.), "Sum of the weights")
		("output", po::value<std::string>(), "output file name")
		("effput", po::value<std::string>()->default_value(""), "b-tag efficiency input file name")
		("genPUMC", po::value<bool>(&genPUMC)->default_value(false), "true if should generate the MC PU file")
		("puMC", po::value<std::string>(&none)->default_value("None"), "MC PU Reweight file")
		("puDATA", po::value<std::string>(&none)->default_value("None"), "DATA PU Reweight file")
		("muonMatchedPt", po::value<double>(&_muonMatchedPt)->default_value(_muonMatchedPt), "Muon Matched Pt Cut")
		("muonMatchedEta", po::value<double>(&_muonMatchedEta)->default_value(_muonMatchedEta), "Muon Matched Eta Cut")
		("leadmuonPt", po::value<double>(&_leadmuonPt)->default_value(_leadmuonPt), "Muon Pt Cut")
		("subleadmuonPt", po::value<double>(&_subleadmuonPt)->default_value(_subleadmuonPt), "Muon Pt Cut")
		("leadmuonEta", po::value<double>(&_leadmuonEta)->default_value(_leadmuonEta), "Muon Eta Cut")
		("subleadmuonEta", po::value<double>(&_subleadmuonEta)->default_value(_subleadmuonEta), "Muon Eta Cut")
		("leadmuonIso", po::value<double>(&_leadmuonIso)->default_value(_leadmuonIso), "Muon Isolation Cut")
		("subleadmuonIso", po::value<double>(&_subleadmuonIso)->default_value(_subleadmuonIso), "Muon Isolation Cut")
		("leadJetPt", po::value<double>(&_leadJetPt)->default_value(_leadJetPt), "Lead Jet Pt Cut")
		("subleadJetPt", po::value<double>(&_subleadJetPt)->default_value(_subleadJetPt), "SubLeading Jet Pt Cut")
		("metPt", po::value<double>(&_metPt)->default_value(_metPt), "MET Pt Cut")
		("dijetMass_VBFTight", po::value<double>(&_dijetMass_VBFTight)->default_value(_dijetMass_VBFTight), "DiJet Mass VBFTight-Category Cut")
		("dijetdEta_VBFTight", po::value<double>(&_dijetdEta_VBFTight)->default_value(_dijetdEta_VBFTight), "DiJet deta VBFTight-Category Cut")
		("dijetMass_ggFTight", po::value<double>(&_dijetMass_ggFTight)->default_value(_dijetMass_ggFTight), "DiJet Mass ggFTight-Category Cut")
		("dimuonPt_ggFTight", po::value<double>(&_dimuonPt_ggFTight)->default_value(_dimuonPt_ggFTight), "DiMuon Pt ggFTight-Category Cut")
		("dimuonPt_01JetsTight", po::value<double>(&_dimuonPt_01JetsTight)->default_value(_dimuonPt_01JetsTight), "DiMuon Pt 01JetsTight-Category Cut")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help") || argc<2)
	{
		std::cout << desc << std::endl;
		return 1;
	}

	//	Assign globals
	__inputfilename = vm["input"].as<std::string>();
	__isMC = vm["isMC"].as<bool>();
	__isCondor = vm["isCondor"].as<bool>();
	__lumi = vm["lumi"].as<double>();
	__sigma_bkg = vm["sigma_bkg"].as<double>();
	__sumw = vm["sumw"].as<double>();
	__outputfilename = vm["output"].as<std::string>();
	__efffilename = vm["effput"].as<std::string>();
	__genPUMC = vm["genPUMC"].as<bool>();
	__puMCfilename = vm["puMC"].as<std::string>();
	__puDATAfilename = vm["puDATA"].as<std::string>();
	_muonMatchedPt = vm["muonMatchedPt"].as<double>();
	_muonMatchedEta = vm["muonMatchedEta"].as<double>();
	_leadmuonPt = vm["leadmuonPt"].as<double>();
	_subleadmuonPt = vm["subleadmuonPt"].as<double>();
	_leadmuonEta = vm["leadmuonEta"].as<double>();
	_subleadmuonEta = vm["subleadmuonEta"].as<double>();
	_leadmuonIso = vm["leadmuonIso"].as<double>();
	_subleadmuonIso = vm["subleadmuonIso"].as<double>();
	_leadJetPt = vm["leadJetPt"].as<double>();
	_subleadJetPt = vm["subleadJetPt"].as<double>();
	_metPt = vm["metPt"].as<double>();
	_dijetMass_VBFTight = vm["dijetMass_VBFTight"].as<double>();
	_dijetdEta_VBFTight = vm["dijetdEta_VBFTight"].as<double>();
	_dijetMass_ggFTight = vm["dijetMass_ggFTight"].as<double>();
	_dimuonPt_ggFTight = vm["dimuonPt_ggFTight"].as<double>();
	_dimuonPt_01JetsTight = vm["dimuonPt_01JetsTight"].as<double>();
	printCuts();

	//	start processing
	process();
	return 0;
}

// VIM DID NOT DELETE ANYTHING FROM THIS FILE SO FAR!!!!!!! 
