#ifndef TGPHYSOBJECTS_HH
#define TGPHYSOBJECTS_HH

#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TGenJet.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include <TLorentzVector.h>



struct _muon{
      float          pt, eta, phi, ptErr;                   // kinematics
      float          staPt, staEta, staPhi;                 // STA track kinematics
      float          pfPt, pfEta, pfPhi;                    // matched PFCandidate
      float          trkIso, ecalIso, hcalIso;              // detector isolation (R=0.3)
      float          chHadIso, gammaIso, neuHadIso, puIso;  // PF isolation variables (R=0.4)
      float          puppiChHadIso,      puppiGammaIso,      puppiNeuHadIso;  // Puppi Isolation R=0.4
      float          puppiChHadIsoNoLep, puppiGammaIsoNoLep, puppiNeuHadIsoNoLep; // Puppi Isolation R=0.4 no lep
      float          d0, dz, sip3d;                         // impact parameter
      float          tkNchi2, muNchi2;                      // track fit normalized chi-square
      float          trkKink, glbKink;                      // track kink
      float          trkHitFrac;                            // fraction of valid tracker hits
      float          chi2LocPos;                            // TRK-STA position match
      float          segComp;                               // compatibility of tracker track with muon segment
      float          caloComp;                              // muon hypothesis compatibility with calo energy
      int            q;                                     // charge
      int            nValidHits;                            // number of valid muon hits in global fit
      unsigned int   typeBits;                              // muon type bits
      unsigned int   selectorBits;                          // MuonSelector bits
      unsigned int   pogIDBits;                             // POG muon IDs from CMSSW
      unsigned int   nTkHits, nPixHits;                     // number of hits in tracker
      unsigned int   nTkLayers, nPixLayers;                 // number of hit layers in tracker
      unsigned int   nMatchStn;                             // number of stations with muon segments
      int            trkID;                                 // tracker track ID (unique per event)
};


struct _electron{

      float          pt, eta, phi;                             // kinematics
      float          scEt, scEta, scPhi;                       // supercluster kinematics
      float          ecalEnergy;                               // ECAL energy
      float          pfPt, pfEta, pfPhi;                       // matching PF-candidate kinematics
      float          trkIso, ecalIso, hcalIso, hcalDepth1Iso;  // detector isolation
      float          chHadIso, gammaIso, neuHadIso, puIso;     // PF isolation variables
      float          ecalPFClusIso, hcalPFClusIso;             // PF cluster isolation variables
      float          puppiChHadIso,      puppiGammaIso,      puppiNeuHadIso;  // Puppi Isolation R=0.4
      float          puppiChHadIsoNoLep, puppiGammaIsoNoLep, puppiNeuHadIsoNoLep; // Puppi Isolation R=0.4 no lep
      float          d0, dz, sip3d;                            // impact parameter
      float          sieie, e1x5, e2x5, e5x5, r9;              // shower shape
      float          eoverp;                                   // E/p
      float          hovere;                                   // H/E
      float          fbrem;                                    // brem fraction
      float          dEtaInSeed, dEtaIn, dPhiIn;               // track-supercluster matching
      float          mva;                                      // electron ID MVA value
      int            q;                                        // charge
      bool           isConv;                                   // identified by track fit based conversion finder?
      unsigned int   nMissingHits;                             // number of missing expected inner hits 
      unsigned int   typeBits;                                 // electron type
      unsigned int   fiducialBits;                             // ECAL fiducial region bits
      int            classification;                           // electron classification
      int            scID;                                     // supercluster ID number (unique per event)
      int            trkID;                                    // track ID number (unique per event)
};



struct _jet{
      float          pt, eta, phi, mass, ptRaw, unc;                          // kinematics
      float          area;                                                    // jet area (from FastJet)
      float          d0, dz;                                                  // impact parameter of leading charged constituent
      float          csv,bmva,cvb,cvl;                                        // CSV b-taggers/c-tagger for the jet
      float          qgid, axis2, ptD;                                        // q/g discriminator and input variables
      int            mult;
      float          q;                                                       // Charge for jet and subjets
      float          mva;                                                     // PU discriminator MVA
      float          beta, betaStar, dR2Mean;                                 // input variables for PU and q/g discriminators
      float          pullY, pullPhi;                                          // Jet pull
      float          chPullY, chPullPhi, neuPullY, neuPullPhi;
      float          chEmFrac, neuEmFrac, chHadFrac, neuHadFrac, muonFrac;    // fractional energy contribution by type
      float          genpt, geneta, genphi, genm;                             // Matched GenJet
      int            partonFlavor, hadronFlavor;                              // Flavor
      unsigned int   nCharged, nNeutrals, nParticles;                         // constituent multiplicity

};

struct _genjet{
    float pdgId; 
    float pt   ;
    float eta  ;
    float phi  ;
    float mass ;

    float elept       ;
    float eleeta      ;
    float elephi      ;
    float elem        ;

    float mupt        ;
    float mueta       ;
    float muphi       ;
    float mum         ;

    float gapt        ;
    float gaeta       ;
    float gaphi       ;
    float gam         ;

    float totpt       ;
    float toteta      ;
    float totphi      ;
    float totm        ;

    float iso03       ;
    float iso04       ;
    float iso05       ;

    float mtrim       ;
    float tau1        ;
    float tau2        ;
};

struct _genparticle{
      int   parent;
      int   pdgId;
      int   status;
      float pt, eta, phi, mass, y;
};


union PhysicsObject{
		_muon					muon;
		_electron 		electron;
		_jet					jet;
		_genjet				genjet;	
		_genparticle	genparticle;	
};


enum TGPhysType{
		MUON,
		ELECTRON,
		TAU,   //NOT IMPLEMENTED
		GAMMA, //NOT IMPLEMENTED 
		JET,
		GENJET,
		GENPARTICLE,
		NONE
};


struct TGPhysObject{
		TGPhysObject();	
		TGPhysObject(baconhep::TMuon* muon);
		TGPhysObject(baconhep::TElectron* electron);
		TGPhysObject(baconhep::TJet* jet);
		TGPhysObject(baconhep::TGenParticle* genparticle);
		TGPhysObject(baconhep::TGenJet* genjet);


		TLorentzVector p4();
		float pt();
		float eta();
		float phi();
		int		q();

    void set_p4(TLorentzVector p4);

		TGPhysType  type_flag;
		PhysicsObject particle;
		TriggerObjects hltMatchBits;


};




#endif  
