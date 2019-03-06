#include "BLT_II/BLTAnalysis/interface/TGPhysicsObject.hh"
#include "BLT_II/BLTAnalysis/interface/BLTHelper.hh"
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include<cmath>

TGPhysObject::TGPhysObject(){
		type_flag = NONE;
}

TGPhysObject::TGPhysObject(baconhep::TMuon* muon){
		type_flag = MUON;

		particle.muon.pt   = muon->pt;
		particle.muon.eta  = muon->eta;
		particle.muon.phi  = muon->phi;
		particle.muon.ptErr= muon->ptErr;

		particle.muon.staPt  = muon->staPt;
		particle.muon.staEta = muon->staEta;
		particle.muon.staPhi = muon->staPhi;

		particle.muon.pfPt  = muon->pfPt;
		particle.muon.pfEta = muon->pfEta;
		particle.muon.pfPhi = muon->pfPhi;

		particle.muon.trkIso	= muon->trkIso; 
		particle.muon.ecalIso = muon->ecalIso;
		particle.muon.hcalIso = muon->hcalIso;

		particle.muon.chHadIso = muon->chHadIso;
		particle.muon.gammaIso = muon->gammaIso;
		particle.muon.neuHadIso= muon->neuHadIso;
		particle.muon.puIso	   = muon->puIso;

		particle.muon.puppiChHadIso = muon->puppiChHadIso;
		particle.muon.puppiGammaIso = muon->puppiGammaIso;
		particle.muon.puppiNeuHadIso= muon->puppiNeuHadIso;

		particle.muon.puppiChHadIsoNoLep = muon->puppiChHadIsoNoLep;
		particle.muon.puppiGammaIsoNoLep = muon->puppiGammaIsoNoLep;
		particle.muon.puppiNeuHadIsoNoLep= muon->puppiNeuHadIsoNoLep;

		particle.muon.d0   = muon->d0;
		particle.muon.dz   = muon->dz;
		particle.muon.sip3d= muon->sip3d;

		particle.muon.tkNchi2   = muon->tkNchi2;
		particle.muon.muNchi2   = muon->muNchi2;

		particle.muon.trkKink      = muon->trkKink;
		particle.muon.glbKink      = muon->glbKink;
		particle.muon.trkHitFrac   = muon->trkHitFrac;

		particle.muon.chi2LocPos   = muon->chi2LocPos;
		particle.muon.segComp      = muon->segComp;

		particle.muon.caloComp      = muon->caloComp;
		particle.muon.q					    = muon->q;
		particle.muon.nValidHits    = muon->nValidHits;
		particle.muon.typeBits	    = muon->typeBits;
		particle.muon.selectorBits  = muon->selectorBits;
		particle.muon.pogIDBits     = muon->pogIDBits;

		particle.muon.nTkHits     = muon->nTkHits;
		particle.muon.nPixHits    = muon->nPixHits;

		particle.muon.nTkLayers     = muon->nTkLayers;
		particle.muon.nPixLayers    = muon->nPixLayers;


		particle.muon.nMatchStn = muon->nMatchStn;
		particle.muon.trkID     = muon->trkID;

		hltMatchBits = muon->hltMatchBits;
}

TGPhysObject::TGPhysObject(baconhep::TElectron* electron){
		type_flag = ELECTRON;

		particle.electron.pt   = electron->pt;
		particle.electron.eta  = electron->eta;
		particle.electron.phi  = electron->phi;

		particle.electron.scEt  = electron->scEt;
		particle.electron.scEta = electron->scEta;
		particle.electron.scPhi = electron->scPhi;

		particle.electron.ecalEnergy= electron->ecalEnergy;

		particle.electron.pfPt  = electron->pfPt;
		particle.electron.pfEta = electron->pfEta;
		particle.electron.pfPhi = electron->pfPhi;

		particle.electron.trkIso	= electron->trkIso; 
		particle.electron.ecalIso = electron->ecalIso;
		particle.electron.hcalIso = electron->hcalIso;
		particle.electron.hcalDepth1Iso = electron->hcalDepth1Iso;

		particle.electron.chHadIso = electron->chHadIso;
		particle.electron.gammaIso = electron->gammaIso;
		particle.electron.neuHadIso= electron->neuHadIso;
		particle.electron.puIso	   = electron->puIso;

		particle.electron.puppiChHadIso = electron->puppiChHadIso;
		particle.electron.puppiGammaIso = electron->puppiGammaIso;
		particle.electron.puppiNeuHadIso= electron->puppiNeuHadIso;

		particle.electron.puppiChHadIsoNoLep = electron->puppiChHadIsoNoLep;
		particle.electron.puppiGammaIsoNoLep = electron->puppiGammaIsoNoLep;
		particle.electron.puppiNeuHadIsoNoLep= electron->puppiNeuHadIsoNoLep;

		particle.electron.d0   = electron->d0;
		particle.electron.dz   = electron->dz;
		particle.electron.sip3d= electron->sip3d;

		//TODO
		/**
      float          sieie, e1x5, e2x5, e5x5, r9;              // shower shape
      float          eoverp;                                   // E/p
      float          hovere;                                   // H/E
      float          fbrem;                                    // brem fraction
      float          dEtaInSeed, dEtaIn, dPhiIn;               // track-supercluster matching
      float          mva;                                      // electron ID MVA value
			*/

			particle.electron.q = electron->q;                                        // charge
			/*
      bool           isConv;                                   // identified by track fit based conversion finder?
      unsigned int   nMissingHits;                             // number of missing expected inner hits 
      unsigned int   typeBits;                                 // electron type
      unsigned int   fiducialBits;                             // ECAL fiducial region bits
      int            classification;                           // electron classification
      int            scID;                                     // supercluster ID number (unique per event)
      int            trkID;                                    // track ID number (unique per event)

		**/


		hltMatchBits = electron->hltMatchBits;
		
}

TGPhysObject::TGPhysObject(baconhep::TJet* jet){
		type_flag = JET;

		particle.jet.pt   = jet->pt;
		particle.jet.eta  = jet->eta;
		particle.jet.phi  = jet->phi;
		particle.jet.mass = jet->mass;
		particle.jet.ptRaw= jet->ptRaw;
		particle.jet.unc  = jet->unc;

		particle.jet.area   = jet->area;

		particle.jet.d0   = jet->d0;
		particle.jet.dz   = jet->dz;
    particle.jet.csv  = jet->csv;
		particle.jet.bmva = jet->bmva;
		particle.jet.cvb  = jet->cvb;
		particle.jet.cvl  = jet->cvl;                                        // CSV b-taggers/c-tagger for the jet

		//TODO:
		/**
      float          qgid, axis2, ptD;                                        // q/g discriminator and input variables
      int            mult;
		**/
    particle.jet.q = jet->q;                                                       // Charge for jet and subjets
/**
      float          mva;                                                     // PU discriminator MVA
      float          beta, betaStar, dR2Mean;                                 // input variables for PU and q/g discriminators
      float          pullY, pullPhi;                                          // Jet pull
      float          chPullY, chPullPhi, neuPullY, neuPullPhi;
      float          chEmFrac, neuEmFrac, chHadFrac, neuHadFrac, muonFrac;    // fractional energy contribution by type
      float          genpt, geneta, genphi, genm;                             // Matched GenJet
		**/
      particle.jet.partonFlavor = jet->partonFlavor;
			particle.jet.hadronFlavor = jet->hadronFlavor;                          // Flavor
      
			particle.jet.nCharged   = jet->nCharged;
			particle.jet.nNeutrals	= jet->nNeutrals; 
			particle.jet.nParticles = jet->nParticles;           // constituent multiplicity

		hltMatchBits = jet->hltMatchBits;

}

TGPhysObject::TGPhysObject(baconhep::TGenParticle* genparticle ){
		type_flag = GENPARTICLE;

		particle.genparticle.pdgId   = genparticle->pdgId;
		particle.genparticle.pt   = genparticle->pt;
		particle.genparticle.eta  = genparticle->eta;
		particle.genparticle.phi  = genparticle->phi;
    particle.genparticle.mass = genparticle->mass;
}

TGPhysObject::TGPhysObject(baconhep::TGenJet* genjet){
		type_flag = GENJET;

		particle.genjet.pt   = genjet->pt;
		particle.genjet.eta  = genjet->eta;
		particle.genjet.phi  = genjet->phi;
}

TLorentzVector TGPhysObject::p4(){
		TLorentzVector p4;

		switch (type_flag){
	
				case MUON			: p4.SetPtEtaPhiM(particle.muon.pt, particle.muon.eta, particle.muon.phi, MUON_MASS); break;
				case ELECTRON : p4.SetPtEtaPhiM(particle.electron.pt, particle.electron.eta, particle.electron.phi, ELE_MASS); break;
				case TAU			: printf("Tau Type Not implemented"); break;
				case GAMMA		: printf("GAMMA TypeNot implemented"); break;
				case JET			: p4.SetPtEtaPhiM(particle.jet.pt, particle.jet.eta, particle.jet.phi, particle.jet.mass); break;
				case GENJET		: p4.SetPtEtaPhiM(particle.jet.genpt, particle.genjet.eta, particle.genjet.phi, particle.genjet.mass); break;
				case GENPARTICLE : p4.SetPtEtaPhiM(particle.genparticle.pt, particle.genparticle.eta, particle.genparticle.phi, particle.genparticle.mass); break;
				case NONE : printf("Type NONE");
		}
		return p4;
}
float TGPhysObject::pt(){
		switch (type_flag){
				case MUON			: return particle.muon.pt;
				case ELECTRON : return particle.electron.pt;
				case TAU			: printf("Tau Type Not implemented");  return -999.0;
				case GAMMA		: printf("GAMMA TypeNot implemented"); return -999.0;
				case JET			: return particle.jet.pt;
				case GENJET		: return particle.genjet.pt;
				case GENPARTICLE : return particle.genparticle.pt;
				case NONE : printf("Type NONE"); return -999.0;
		}
		return -888.0;
}
float TGPhysObject::eta(){
		switch (type_flag){
				case MUON			: return particle.muon.eta;
				case ELECTRON : return particle.electron.eta;
				case TAU			: printf("Tau Type Not implemented");  return -999.0;
				case GAMMA		: printf("GAMMA TypeNot implemented"); return -999.0;
				case JET			: return particle.jet.eta;
				case GENJET		: return particle.genjet.eta;
				case GENPARTICLE : return particle.genparticle.eta;
				case NONE : printf("Type NONE"); return -999.0;
		}
		return -888.0;
}
float TGPhysObject::phi(){
		switch (type_flag){
				case MUON			: return particle.muon.phi;
				case ELECTRON : return particle.electron.phi;
				case TAU			: printf("Tau Type Not implemented");  return -999.0;
				case GAMMA		: printf("GAMMA TypeNot implemented"); return -999.0;
				case JET			: return particle.jet.phi;
				case GENJET		: return particle.genjet.phi;
				case GENPARTICLE : return particle.genparticle.phi;
				case NONE : printf("Type NONE"); return -999.0;
		}
		return -888.0;
}
int TGPhysObject::q(){
		switch (type_flag){
				case MUON			: return particle.muon.q;
				case ELECTRON : return particle.electron.q;
				case TAU			: printf("Tau Type Not implemented");  return -999;
				case GAMMA		: printf("GAMMA TypeNot implemented"); return -999;
				case JET			: return (int)particle.jet.q;
				case GENJET		: return particle.genjet.pdgId / fabs(particle.genjet.pdgId);
				case GENPARTICLE : return particle.genparticle.pdgId / abs(particle.genparticle.pdgId);
				case NONE : printf("Type NONE"); return -999;
		}
		return -888;
}

void TGPhysObject::set_p4(TLorentzVector p4){
		switch (type_flag){
				case MUON			: 
        {
            particle.muon.pt  = p4.Pt(); 
            particle.muon.eta = p4.Eta();
            particle.muon.phi = p4.Phi();
        }break;
				case ELECTRON :
        {
            particle.electron.pt  = p4.Pt(); 
            particle.electron.eta = p4.Eta();
            particle.electron.phi = p4.Phi();
        }break;
				case TAU			: printf("Tau Type Not implemented"); break;
				case GAMMA		: printf("GAMMA Type Not implemented"); break;
				case JET			: printf("JET Type Not implemented"); break;
				case GENJET		: printf("GENJET Type Not implemented"); break;
				case GENPARTICLE :
        {
            particle.genparticle.pt  = p4.Pt(); 
            particle.genparticle.eta = p4.Eta();
            particle.genparticle.phi = p4.Phi();
        } break;
				case NONE : printf("Type NONE");
		}
}

void print_tgphysobj(TGPhysObject* obj){
    printf("TGPhysObject: \nid: %d pt: %f eta: %f phi: %f charge: %d\n", 
          (int) obj->type_flag, (double)obj->pt(), (double)obj->eta(), (double)obj->phi(), obj->q());
}







