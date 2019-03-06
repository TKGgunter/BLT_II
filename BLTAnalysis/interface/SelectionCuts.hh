#ifndef SelectionCuts_hh
#define SelectionCuts_hh

#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TGenJet.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"

#include "BLT_II/BLTAnalysis/interface/TGPhysicsObject.hh"

//CMSSW libraries
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"


//https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_Isolation
//==============Good Muon Selection Cuts ==========
bool muIsGood									= true; 
float muEta										= 2.4;
float muPt										= 10;
float muNormChi2							= 10;
int muNumbValidHits						= 0;
float muSumPFIso							= 0.15;
bool muIsGLB									=	true;
bool muIsPF										=	true;
int muNumbMatchedStat				  =	1;
int muTrackLayersMeasurements = 5;
float muDxy										=	0.2;
float muDz										=	0.5;


//https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2 
//x80
//==============Good Tight Electron Selection Cuts=========
//bool elIsGood								= true;
//EB
float elEtaEB									= 1.479;
float elSigmaIEtaIEtaEB				= 0.00998;
float elEtaInSeedEB						= 0.00308;
float elDeltaPhiInEB					=	0.0816;
float elHadOverEmEB						= 0.0414;
float elSumPFIsoEB						= 0.0588;
float elInvE_InvPEB						= 0.0129;
float elMissingHitsEB					= 1;
float elDxyEB									= 0.05;
float elDzEB									= 0.10;


//EE
float elEtaEE									=	1.479;
float elSigmaIEtaIEtaEE       = 0.0292; 
float elSCDeltaEtaInEE        = 0.00605; 
float elSCDeltaPhiInEE        = 0.0394;
float elHadOverEmEE           = 0.0641;  
float elSumPFIsoEE            = 0.0571;
float elInvE_InvPEE						= 0.0129; 
float elMissingHitsEE					= 1;
float elDxyEE									= 0.1;
float elDzEE									= 0.2;

//==============Loose Electron Selection Cuts=========
//bool elIsGood								= true;
//EB
float l_elSigmaIEtaIEtaEB				= 0.011;
float l_elEtaInSeedEB						= 0.00477;
float l_elDeltaPhiInEB					=	0.222;
float l_elHadOverEmEB						= 0.298;
float l_elSumPFIsoEB						= 0.0994;
float l_elInvE_InvPEB						= 0.241;
float l_elMissingHitsEB					= 1;


//EE
float l_elSigmaIEtaIEtaEE       = 0.0314; 
float l_elSCDeltaEtaInEE        = 0.00868; 
float l_elSCDeltaPhiInEE        = 0.213;
float l_elHadOverEmEE           = 0.101;  
float l_elSumPFIsoEE            = 0.107;
float l_elInvE_InvPEE						= 0.14; 
float l_elMissingHitsEE					= 1;

//General
float elPt										=	10.;
float elEta										= 2.5;

//https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
//==============Good Jet Selection Cuts==============


//ETA < 2.4
float jetChEmFrac							=	.99;
float jetChHadFrac						=	0;
unsigned int jetChMult				= 0;
//charged multiplicity > 0

//ETA < 2.7
float jetNeuHadFrac						=	0.99;
float jetNeuEmFrac						= 0.99;
float jetNumConstit						=	1.;
float jetPt										=	20.;
float jetDR_mu_el							=	0.4;
float jetDR_me								= 0.3;

//ETA < 3.0
float jetNeuEmFrac_eta_3			=	0.01;
float jetNeuHadFrac_eta_3			=	0.98;
unsigned int jetnNeutrals_eta_3    = 2;

//ETA > 3.0
float jetNeuEmFrac_eta_4						=	.90;
unsigned int	jetnNeutrals_eta_4		= 10;

float jetEta									=	4.7;
float vtxCut									=	.3;

bool test_bits_(unsigned int bits, unsigned int test);

void jetSelection( std::vector<TGPhysObject> &jetList , TClonesArray* jets, std::vector<TGPhysObject> &leptonList);
void muonSelection( std::vector<TGPhysObject> &muonList, TClonesArray* muons);
void electronSelection( std::vector<TGPhysObject> &elecList, TClonesArray* electrons, float rhoFactor);
void loosemuonSelection( std::vector<TGPhysObject> &muonList, TClonesArray* muons);
void looseelectronSelection( std::vector<TGPhysObject> &elecList, TClonesArray* electrons, float rhoFactor);

void wwSelection(std::vector<TGPhysObject> &wwList, TClonesArray* leptons);
int genJetSelection(TClonesArray* gen_particles, TClonesArray* reco_jets, float pt_cut, float eta_cut, std::vector<TGPhysObject> &genJetList);
void wwGenLepSelection(std::vector<TGPhysObject> &genLepList, TClonesArray* gen_particles);
void wwGenNeutrinoSelection(std::vector<TGPhysObject> &genNeutrinoList, TClonesArray* gen_particles);

#endif
