#include <array>
#include <memory>
#include <stdexcept>
#include <string>
#include <boost/math/special_functions/erf.hpp>
#include <random>
#include <cstdint>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"
#include "correction.h"
#include <nlohmann/json.hpp>	// JSON reader for correction files - Golden JSON files aren't formatted as correctionlib for CMS Run 3.

using Json = nlohmann::json;

//  Used by JER correction
std::uint64_t hashCombine(std::uint64_t seed, std::uint64_t value){
    // Similar to the common boost hash-combine operation.
    seed ^= value + 0x9e3779b97f4a7c15ULL
    + (seed << 6)
    + (seed >> 2);
    return seed;
}

// Used by JER correction
double deterministicGaussian(std::uint64_t event, std::uint32_t luminosityBlock, std::uint32_t run, std::uint32_t jetIndex){
    std::uint64_t seed = 0;

    seed = hashCombine(seed, run);
    seed = hashCombine(seed, luminosityBlock);
    seed = hashCombine(seed, event);
    seed = hashCombine(seed, jetIndex);

    std::mt19937_64 generator(seed);
    std::normal_distribution<double> gaussian(0.0, 1.0);

    return gaussian(generator);
}

// Used by KIT corrections
struct CrystalBall{     // Crystal Ball Function
    double pi=3.14159;
    double sqrtPiOver2=sqrt(pi/2.0);
    double sqrt2=sqrt(2.0);
    double m;
    double s;
    double a;
    double n;
    double B;
    double C;
    double D;
    double N;
    double NA;
    double Ns;
    double NC;
    double F;
    double G;
    double k;
    double cdfMa;
    double cdfPa;
CrystalBall():m(0),s(1),a(10),n(10){
    init();
}
CrystalBall(double mean, double sigma, double alpha, double n)
    :m(mean),s(sigma),a(alpha),n(n){
    init();
}
void init(){
    double fa = fabs(a);
    double ex = exp(-fa*fa/2);
    double A  = pow(n/fa, n) * ex;
    double C1 = n/fa/(n-1) * ex;
    double D1 = 2 * sqrtPiOver2 * erf(fa/sqrt2);
    B = n/fa-fa;
    C = (D1+2*C1)/C1;
    D = (D1+2*C1)/2;
    N = 1.0/s/(D1+2*C1);
    k = 1.0/(n-1);
    NA = N*A;
    Ns = N*s;
    NC = Ns*C1;
    F = 1-fa*fa/n;
    G = s*n/fa;
    cdfMa = cdf(m-a*s);
    cdfPa = cdf(m+a*s);
}
double pdf(double x) const{
    double d=(x-m)/s;
    if(d<-a) return NA*pow(B-d, -n);
    if(d>a) return NA*pow(B+d, -n);
    return N*exp(-d*d/2);
}
double pdf(double x, double ks, double dm) const{
    double d=(x-m-dm)/(s*ks);
    if(d<-a) return NA/ks*pow(B-d, -n);
    if(d>a) return NA/ks*pow(B+d, -n);
    return N/ks*exp(-d*d/2);

}
double cdf(double x) const{
    double d = (x-m)/s;
    if(d<-a) return NC / pow(F-s*d/G, n-1);
    if(d>a) return NC * (C - pow(F+s*d/G, 1-n) );
    return Ns * (D - sqrtPiOver2 * erf(-d/sqrt2));
}
double invcdf(double u) const{
    if(u<cdfMa) return m + G*(F - pow(NC/u, k));
    if(u>cdfPa) return m - G*(F - pow(C-u/NC, -k) );
    return m - sqrt2 * s * boost::math::erf_inv((D - u/Ns )/sqrtPiOver2);
}
};

double get_rndm(double eta, double phi, float nL, std::uint64_t evtNumber, std::uint32_t lumiNumber, const correction::CorrectionSet& cset) {
    // obtain parameters from correctionlib
    double mean = cset.at("cb_params")->evaluate({abs(eta), nL, 0});
    double sigma = cset.at("cb_params")->evaluate({abs(eta), nL, 1});
    double n = cset.at("cb_params")->evaluate({abs(eta), nL, 2});
    double alpha = cset.at("cb_params")->evaluate({abs(eta), nL, 3});

    // instantiate CB and get random number following the CB
    CrystalBall cb(mean, sigma, alpha, n);
    const int evtSeed = static_cast<int>(evtNumber ^ (evtNumber >> 32));
    double rndm = cset.at("RandomSmearing")->evaluate({evtSeed, (int)lumiNumber, phi});

    return cb.invcdf(rndm);
}


double get_std(double pt, double eta, float nL, const correction::CorrectionSet& cset) {

    // obtain paramters from correctionlib
    double param_0 = cset.at("poly_params")->evaluate({abs(eta), nL, 0});
    double param_1 = cset.at("poly_params")->evaluate({abs(eta), nL, 1});
    double param_2 = cset.at("poly_params")->evaluate({abs(eta), nL, 2});

    // calculate value and return max(0, val)
    double sigma = param_0 + param_1 * pt + param_2 * pt*pt;
    if (sigma < 0) sigma = 0;
    return sigma;
}


double get_k(double eta, string var, const correction::CorrectionSet& cset) {

    // obtain parameters from correctionlib
    double k_data = cset.at("k_data")->evaluate({abs(eta), var});
    double k_mc = cset.at("k_mc")->evaluate({abs(eta), var});

    // calculate residual smearing factor
    // return 0 if smearing in MC already larger than in data
    double k = 0;
    if (k_mc < k_data) k = sqrt(k_data*k_data - k_mc*k_mc);
    return k;
}


double pt_resol(double pt, double eta, double phi, float nL, std::int64_t evtNumber, std::int32_t lumiNumber, const correction::CorrectionSet& cset, double low_pt_threshold = 26) {
    // load correction values
    double rndm = (double) get_rndm(eta, phi, nL, evtNumber, lumiNumber, cset);
    double std = (double) get_std(pt, eta, nL, cset);
    double k = (double) get_k(eta, "nom", cset);

    // calculate corrected value and return original value if a parameter is nan
    double ptc = pt * ( 1 + k * std * rndm);
    if (isnan(ptc)) ptc = pt;
    if(ptc / pt > 2 || ptc / pt < 0.1 || ptc < 0 || pt < low_pt_threshold || pt > 200){
        ptc = pt;
    }
    // TODO: Understand why for evts with pT < threshold the pt_corr is set to one
    return ptc;
}

double pt_resol_var(double pt_woresol, double pt_wresol, double eta, string updn, const correction::CorrectionSet& cset){

    double k = (double) get_k(eta, "nom", cset);

    if (k==0) return pt_wresol;

    double k_unc = cset.at("k_mc")->evaluate({abs(eta), "stat"});

    double std_x_rndm = (pt_wresol / pt_woresol - 1) / k;

    double pt_var = pt_wresol;

    if (updn=="up"){
        pt_var = pt_woresol * (1 + (k+k_unc) * std_x_rndm);
    }
    else if (updn=="dn"){
        pt_var = pt_woresol * (1 + (k-k_unc) * std_x_rndm);
    }
    else {
        cout << "ERROR: updn must be 'up' or 'dn'" << endl;
    }
    if(pt_var / pt_woresol > 2 || pt_var / pt_woresol < 0.1 || pt_var < 0){
        pt_var = pt_woresol;
    }

    return pt_var;
}

double pt_scale(bool is_data, double pt, double eta, double phi, int charge, const correction::CorrectionSet& cset, double low_pt_threshold = 0) {

    // use right correction
    string dtmc = is_data ? "data" : "mc";

    double a = cset.at("a_"+dtmc)->evaluate({eta, phi, "nom"});
    double m = cset.at("m_"+dtmc)->evaluate({eta, phi, "nom"});
    if(pt < low_pt_threshold)
        return pt;

    return 1. / (m/pt + charge * a);
}


double pt_scale_var(double pt, double eta, double phi, int charge, string updn, const correction::CorrectionSet& cset) {

    double stat_a = cset.at("a_mc")->evaluate({eta, phi, "stat"});
    double stat_m = cset.at("m_mc")->evaluate({eta, phi, "stat"});
    double stat_rho = cset.at("m_mc")->evaluate({eta, phi, "rho_stat"});

    double unc = pt*pt*sqrt(stat_m*stat_m / (pt*pt) + stat_a*stat_a + 2*charge*stat_rho*stat_m/pt*stat_a);

    double pt_var = pt;

    if (updn=="up"){
        pt_var = pt + unc;
    }
    else if (updn=="dn"){
        pt_var = pt - unc;
    }

    return pt_var;
}


double muon_ID_ISO_SF(double pt, double eta, const char *c, const correction::CorrectionSet& cset) {
    double sf = 1.;
    sf *= cset.at("NUM_MediumID_DEN_TrackerMuons")->evaluate({std::abs(eta), pt, c});
    sf *= cset.at("NUM_LoosePFIso_DEN_MediumID")->evaluate({std::abs(eta), pt, c});
    return sf;
}

double muon_HLT_SF(double pt, double eta, const char *variation, const correction::CorrectionSet& cset){
    return cset.at("NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium")->evaluate({std::abs(eta), pt, variation});
}

double JEC(double JetA, double eta, double pt, double rho, double phi, const correction::CorrectionSet& cset, const std::string& correctionName, int year, int era, double run=-1.){
    if(year>2023 || (year==2023 && era==2)){
        if(run<0) return cset.compound().at(correctionName)->evaluate({JetA, eta, pt, rho, phi});
        return cset.compound().at(correctionName)->evaluate({JetA, eta, pt, rho, phi, run});
    }
    else{
        if(run<0) return cset.compound().at(correctionName)->evaluate({JetA, eta, pt, rho});
        return cset.compound().at(correctionName)->evaluate({JetA, eta, pt, rho, run});
    }

//  double jec = 1.;
//  jec *= cset.at("Summer22_22Sep2023_RunCD_V4_DATA_L1L2L3Res_AK4PFPuppi")->evaluate({JetA, eta, pt, rho});
//  jec *= cset.at("Summer22_22Sep2023_V4_MC_L1FastJet_AK4PFPuppi")->evaluate({JetA, eta, pt, rho});
//  jec *= cset.at("Summer22_22Sep2023_V4_MC_L2Relative_AK4PFPuppi")->evaluate({eta, pt});     // 2 arg
//  jec *= cset.at("Summer22_22Sep2023_V4_MC_L3Absolute_AK4PFPuppi")->evaluate({eta, pt});     // 2 arg
//  return jec;
}

std::pair<double, double> JER(double eta, double pt, double rho, const correction::CorrectionSet& cset, const std::string& scaleName="Summer22_22Sep2023_JRV2_MC_ScaleFactor_AK4PFPuppi", const std::string& resName="Summer22_22Sep2023_JRV2_MC_PtResolution_AK4PFPuppi"){
    pair<double,double> jer;
    jer.first = cset.at(scaleName)->evaluate({eta, pt});          // Scale - 2 arg
    jer.second = cset.at(resName)->evaluate({eta, pt, rho});   // Resolution - 3 arg
    return jer;
}

double PU(double nTrueInt, const correction::CorrectionSet& cset, const std::string& correctionName="Collisions2022_355100_357900_eraBCD_GoldenJson"){
    double pu = 1.;
    pu *= cset.at(correctionName)->evaluate({nTrueInt, "nominal"});
    return pu;
}

vector<double> BTAG(const char* wp, const int flav, const double pt, const double eta, const correction::CorrectionSet& cset, const int year=2022){
    vector<double> sf;
    const char* key = "";
    double absEta = std::abs(eta);
    if(flav==4 || flav==5){
        key = (year>2023) ? "UParTAK4_comb" : "particleNet_comb";
    }
    else if(flav==0){
        key = (year>2023) ? "UParTAK4_light" : "particleNet_light";
    }
    else{
        return {1., 1., 1.};
    }
    sf.push_back(cset.at(key)->evaluate({"central", wp, flav, absEta, pt}));    // "central" or "up" or "down", flav: hadronFlavour, wp: working point
    sf.push_back(cset.at(key)->evaluate({"up", wp, flav, absEta, pt}));         // "central" or "up" or "down", flav: hadronFlavour, wp: working point
    sf.push_back(cset.at(key)->evaluate({"down", wp, flav, absEta, pt}));       // "central" or "up" or "down", flav: hadronFlavour, wp: working point
//  sf.push_back(cset.at("particleNet_comb")->evaluate({"central", wp, flav, absEta, pt}));     // "central" or "up" or "down", flav: hadronFlavour, wp: working point
//  sf.push_back(cset.at("particleNet_light")->evaluate({"central", wp, flav, absEta, pt}));     // "central" or "up" or "down", flav: hadronFlavour, wp: working point
//  sf.push_back(cset.at("particleNet_wp_values")->evaluate({"central", wp, flav, absEta, pt}));  // "central" or "up" or "down", flav: hadronFlavour, wp: working point
    return sf;
}

using HistPtr = std::unique_ptr<TH1>;
using EfficiencyArray = std::array<HistPtr, 6>;

EfficiencyArray EffFinder(TFile* effroot)
{
    if (!effroot || effroot->IsZombie()) {
        throw std::runtime_error("EffFinder: invalid ROOT file");
    }

    TDirectory* dir = effroot->GetDirectory("Efficiencies");

    if (!dir) {
        throw std::runtime_error(
            "EffFinder: directory 'Efficiencies' was not found"
        );
    }

    // Replace these with the actual names in your ROOT file.
    const std::array<std::string, 6> numeratorNames = {
        "h_b_loose",
        "h_b_medium",
        "h_c_loose",
        "h_c_medium",
        "h_l_loose",
        "h_l_medium"
    };

    const std::array<std::string, 6> denominatorNames = {
        "h_b_all",
        "h_b_all",
        "h_c_all",
        "h_c_all",
        "h_l_all",
        "h_l_all"
    };

    EfficiencyArray efficiencies;

    for (std::size_t i = 0; i < efficiencies.size(); ++i) {
        TH1* numerator =
        dynamic_cast<TH1*>(dir->Get(numeratorNames[i].c_str()));

        TH1* denominator =
        dynamic_cast<TH1*>(dir->Get(denominatorNames[i].c_str()));

        if (!numerator) {
            throw std::runtime_error(
                "EffFinder: histogram not found: " + numeratorNames[i]
            );
        }

        if (!denominator) {
            throw std::runtime_error(
                "EffFinder: histogram not found: " + denominatorNames[i]
            );
        }

        const std::string efficiencyName =
        "eff_" + numeratorNames[i];

        efficiencies[i].reset(
            dynamic_cast<TH1*>(numerator->Clone(efficiencyName.c_str()))
        );

        if (!efficiencies[i]) {
            throw std::runtime_error(
                "EffFinder: could not clone " + numeratorNames[i]
            );
        }

        // Prevent the histogram from being owned by the ROOT file.
        efficiencies[i]->SetDirectory(nullptr);
        efficiencies[i]->Divide(denominator);
    }

    return efficiencies;
}

bool passGoldenJSON(const Json& golden,
                    UInt_t run,
                    UInt_t luminosityBlock)
{
    const auto runIt = golden.find(std::to_string(run));

    if (runIt == golden.end())
        return false;

    for (const auto& range : *runIt) {
        if (range.size() != 2) throw std::runtime_error("Malformed Golden JSON lumi range");

        const UInt_t firstLumi = range[0].get<UInt_t>();
        const UInt_t lastLumi  = range[1].get<UInt_t>();

        if (firstLumi <= luminosityBlock && luminosityBlock <= lastLumi) return true;
    }

    return false;
}
