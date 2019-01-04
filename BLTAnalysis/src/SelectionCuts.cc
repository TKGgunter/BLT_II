#include "BLT_II/BLTAnalysis/interface/SelectionCuts.hh"
#include "BLT_II/BLTAnalysis/interface/tgUtility.hh"
#include <cmath>
#include <TClonesArray.h>


bool test_bits_(unsigned int bits, unsigned int test) {
    return (bits & test) == test;
}

void jetSelection( std::vector<TGPhysObject> &jetList , TClonesArray* jets, std::vector<TGPhysObject> &leptonList){

	
	for( Int_t i = 0; i < jets->GetEntries(); i++){
		baconhep::TJet* jet = (baconhep::TJet*) jets->At(i);
		bool isJetGood = false;
		if( (fabs(jet->eta) < jetEta) &&
				jet->pt > jetPt )
		{
			if( (fabs(jet->eta) < 2.7)  &&
					(jet->neuHadFrac < jetNeuHadFrac) &&
					(jet->neuEmFrac < jetNeuEmFrac) &&
					(jet->nParticles > jetNumConstit) ){
				isJetGood = true;
			}
			if( (fabs(jet->eta) < 2.4) &&
					(jet->chEmFrac < jetChEmFrac) &&
					(jet->neuHadFrac < jetNeuHadFrac) &&
					(jet->neuEmFrac < jetNeuEmFrac) &&
					(jet->nParticles > jetNumConstit) &&
					(jet->chHadFrac > jetChHadFrac) &&
					(jet->nCharged > jetChMult)) {
				isJetGood = true;
			}
			if( (fabs(jet->eta) < 3.) && (fabs(jet->eta) > 2.7) &&
					(jet->neuHadFrac < jetNeuHadFrac_eta_3) &&
					(jet->neuEmFrac > jetNeuEmFrac_eta_3) &&
					(jet->nNeutrals > jetnNeutrals_eta_3) ){
				isJetGood = true;
			}
			if( (fabs(jet->eta) > 3.) && 
					(jet->neuEmFrac < jetNeuEmFrac_eta_4) &&
					(jet->nNeutrals > jetnNeutrals_eta_4) ) isJetGood = true;


			/////////////////////////////
			if ( isJetGood == true){
				for(  auto lepton = leptonList.begin(); lepton != leptonList.end(); ++lepton){
					float dR_muon = dR(&*lepton, jet);
					if(dR_muon < jetDR_mu_el){
						isJetGood = false;
					}
				}
			}
			if( isJetGood == true){
				jetList.push_back(TGPhysObject(jet));
			}
		}
	}
}


//Muon Selection
void muonSelection(std::vector<TGPhysObject> &muonList, TClonesArray* muons){
	std::vector<baconhep::TMuon*> tmp_muons;
	float pfIso;
	for(Int_t i = 0; i < muons->GetEntries(); i++){
		baconhep::TMuon* muon = (baconhep::TMuon*) muons->At(i);
		pfIso = (muon->chHadIso + std::max( 0., muon->neuHadIso + muon->gammaIso - 0.5*muon->puIso)) / muon->pfPt;

		if (
				fabs(muon->pfEta) < muEta   &&
				test_bits_(muon->typeBits, baconhep::kPFMuon) == muIsGood &&
				test_bits_(muon->typeBits, baconhep::kGlobal) == muIsGLB && 
				muon->muNchi2 < muNormChi2 &&
				muon->nValidHits > muNumbValidHits &&
				int(muon->nMatchStn) > muNumbMatchedStat &&
				fabs(muon->d0) < muDxy &&
				fabs(muon->dz) < muDz &&
				int(muon->nTkLayers)  > muTrackLayersMeasurements &&
				muon->pfPt > 10 &&
				pfIso < muSumPFIso
			)
		{
			
			muonList.push_back(muon);

		}
	}
} 

//Note
//Oct 3, 2018 
//taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SMP18004QA
//the IFCA/MIT description
void loosemuonSelection( std::vector<TGPhysObject> &muonList, TClonesArray* muons){
	std::vector<baconhep::TMuon*> tmp_muons;
	float pfIso;
	for(Int_t i = 0; i < muons->GetEntries(); i++){
		baconhep::TMuon* muon = (baconhep::TMuon*) muons->At(i);
		pfIso = (muon->chHadIso + std::max( 0., muon->neuHadIso + muon->gammaIso - 0.5*muon->puIso)) / muon->pfPt;

		if (
				fabs(muon->pfEta) < muEta   &&
				test_bits_(muon->typeBits, baconhep::kPFMuon) == muIsGood &&
				test_bits_(muon->typeBits, baconhep::kGlobal) == muIsGLB && 
				muon->muNchi2 < muNormChi2 &&
				muon->nValidHits > muNumbValidHits &&
				int(muon->nMatchStn) > muNumbMatchedStat &&
				fabs(muon->d0) < muDxy &&
				fabs(muon->dz) < muDz &&
				int(muon->nTkLayers)  > muTrackLayersMeasurements &&
				muon->pfPt > 10 &&
				pfIso < 0.25 //muSumPFIso
			)
		{
			
			muonList.push_back(muon);

		}
	}
}




//Electron Selection
void electronSelection(std::vector<TGPhysObject> &elecList, TClonesArray* electrons, float rhoFactor){

	int iEta = 0;
	float etaBins[8] = {0., 1., 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};
//TODO
	float EAEL[7] = {0.13, 0.14, 0.07, 0.09, 0.11, 0.11, 0.14}; //NEEDS TO BE UPDATED
//https://github.com/ikrav/cmssw/blob/egm_id_80X_v1/RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt

	for(Int_t i = 0; i < electrons->GetEntries(); i++){
		baconhep::TElectron* electron = (baconhep::TElectron*) electrons->At(i);

		for( unsigned int j = 0; j < 8; ++j ){
			if( fabs(electron->scEta) > etaBins[j] && fabs(electron->scEta) < etaBins[j+1] ){
				iEta = j;
			}
		}


		float combIso = electron->chHadIso + std::max(0., (double)electron->neuHadIso + electron->gammaIso - rhoFactor * EAEL[iEta]);
		float sumPfIso = combIso / electron->pt;
		float invE_invP = fabs(1./ (electron->eoverp * electron->pt) - electron->eoverp * 1./ (electron->eoverp * electron->pt));
	
		//Barrel Region////////////////
		if(fabs(electron->scEta) < elEtaEB ){

			if(	fabs(electron->eta) < elEta &&
					electron->pt > elPt &&
					electron->sieie								< elSigmaIEtaIEtaEB &&
					electron->dEtaInSeed					< elEtaInSeedEB  &&
					fabs(electron->dPhiIn)				< elDeltaPhiInEB &&
          electron->hovere							< elHadOverEmEB &&
					sumPfIso											< elSumPFIsoEB &&
					invE_invP											< elInvE_InvPEB &&
					electron->nMissingHits				<= elMissingHitsEB && 
					!electron->isConv && 
					fabs(electron->d0)						< elDxyEB &&
					fabs(electron->dz)						< elDzEB //&&

					//energyInversion < 0.05 && //new
				)      
      {                                                                                 
        elecList.push_back(TGPhysObject(electron));
      } 
    }
		//////////////////////////
		//End Cap region//////////
		else if(fabs(electron->scEta) > elEtaEE){

      if(	fabs(electron->eta) < elEta && 
					electron->pt > elPt &&
					electron->sieie < elSigmaIEtaIEtaEE &&
					fabs(electron->dEtaIn) < elSCDeltaEtaInEE && 
          fabs(electron->dPhiIn) < elSCDeltaPhiInEE &&
          electron->hovere  < elHadOverEmEE &&
          sumPfIso < elSumPFIsoEE &&
					invE_invP < elInvE_InvPEE &&
					electron->nMissingHits <= elMissingHitsEE &&			
					!electron->isConv && 
        	fabs(electron->d0) < elDxyEE &&
        	fabs(electron->dz) < elDzEE //&&
				)         
      {
        elecList.push_back(TGPhysObject(electron));
      }
    }
    else;
	}	
}

void looseelectronSelection( std::vector<TGPhysObject> &elecList, TClonesArray* electrons, float rhoFactor){
	int iEta = 0;
	float etaBins[8] = {0., 1., 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};
	float EAEL[7] = {0.13, 0.14, 0.07, 0.09, 0.11, 0.11, 0.14};

	for(Int_t i = 0; i < electrons->GetEntries(); i++){
		baconhep::TElectron* electron = (baconhep::TElectron*) electrons->At(i);

		for( unsigned int j = 0; j < 8; ++j ){
			if( fabs(electron->scEta) > etaBins[j] && fabs(electron->scEta) < etaBins[j+1] ){
				iEta = j;
			}
		}


		float combIso = electron->chHadIso + std::max(0., (double)electron->neuHadIso + electron->gammaIso - rhoFactor * EAEL[iEta]);
		float sumPfIso = combIso / electron->pt;
		float invE_invP = fabs(1./ (electron->eoverp * electron->pt) - electron->eoverp * 1./ (electron->eoverp * electron->pt));
	
		//Barrel Region////////////////
		if(fabs(electron->scEta) < elEtaEB ){

			if(	fabs(electron->eta) < elEta &&
					electron->pt > 10 &&
					electron->sieie								<  l_elSigmaIEtaIEtaEB &&
					electron->dEtaInSeed					<  l_elEtaInSeedEB  &&
					fabs(electron->dPhiIn)				<  l_elDeltaPhiInEB &&
          electron->hovere							<  l_elHadOverEmEB &&
					sumPfIso											<  l_elSumPFIsoEB &&
					invE_invP											<  l_elInvE_InvPEB &&
					electron->nMissingHits				<= l_elMissingHitsEB && 
					!electron->isConv && 
					fabs(electron->d0)						<  elDxyEB &&
					fabs(electron->dz)						<  elDzEB //&&

					//energyInversion < 0.05 && //new
				)      
      {                                                                                 
        elecList.push_back(TGPhysObject(electron));
      } 
    }
		//////////////////////////
		//End Cap region//////////
		else if(fabs(electron->scEta) > elEtaEE){

      if(	fabs(electron->eta) < elEta && 
					electron->pt > 10 &&
					electron->sieie						<  l_elSigmaIEtaIEtaEE &&
					fabs(electron->dEtaIn)		<  l_elSCDeltaEtaInEE && 
          fabs(electron->dPhiIn) 		<  l_elSCDeltaPhiInEE &&
          electron->hovere					<  l_elHadOverEmEE &&
          sumPfIso									<  l_elSumPFIsoEE &&
					invE_invP 								<  l_elInvE_InvPEE &&
					electron->nMissingHits    <= l_elMissingHitsEE &&			
					!electron->isConv && 
        	fabs(electron->d0)				< elDxyEE &&
        	fabs(electron->dz) 				< elDzEE //&&
				)         
      {
        elecList.push_back(TGPhysObject(electron));
      }
    }
    else;
	}	
}




void wwSelection(std::vector<TGPhysObject> &wwList, TClonesArray* gen_particles){
		int id = 24;


		FOR_IN_fARR(particle, gen_particles, baconhep::TGenParticle){
				if( particle == NULL){
					continue;
				}

				if( abs(particle->pdgId ) == id ){
					int pdgId = abs( ((baconhep::TGenParticle*)gen_particles->At(particle->parent))->pdgId);
					//NOTE
					//quarks have ids <= 5
					//gluons are 21
					if( pdgId <= 5 || pdgId == 21){
						wwList.push_back(TGPhysObject(particle));
					}
				}
		}
}




int genJetSelection(TClonesArray* gen_particles, TClonesArray* gen_jets, float pt_cut, float eta_cut, std::vector<TGPhysObject> &genJetList){
	int n_gen_jets = 0;

	for(Int_t i = 0; i < gen_jets->GetSize(); i++){

		baconhep::TGenJet* jet = (baconhep::TGenJet*) gen_jets->At(i);
		if( jet == NULL){
			continue;
		}

		if(	jet->pt > pt_cut &&
				fabs(jet->eta) < eta_cut ){

			bool pass_dR = true;
			for(Int_t r_i = 0; r_i < gen_particles->GetSize(); r_i++){

				baconhep::TGenParticle* particle = (baconhep::TGenParticle*) gen_particles->At(r_i);
				if( particle == NULL ){ continue;}

				if(abs(particle->pdgId)  == 11 || abs(particle->pdgId) == 13 || abs(particle->pdgId) == 15){
					int mother_pdgId = -1;
					if (particle->parent >= 0){
						mother_pdgId = abs(((baconhep::TGenParticle*)gen_particles->At(particle->parent))->pdgId);
					}
					if(mother_pdgId != 24 ) continue;



					float dphi = particle->phi - jet->phi;
					if(dphi > M_PI){
						dphi -= 2.0*M_PI;
					}
					else if(dphi <= -M_PI){
						dphi += 2.0*M_PI;
					}
					float dR = pow(pow(particle->eta - jet->eta, 2.0) + pow(dphi, 2.0), 0.5);

					if(  dR < 0.4 && particle->pt > 10.0 ){
						pass_dR = false; 
					}
				}
			}	
			//std::cout << "##########################" << std::endl;

			if(pass_dR == true){ 
				n_gen_jets++;
				genJetList.push_back(TGPhysObject(jet));
				//std::cout << "eta " << eta_cut << " " << (eta_cut < 2.9) << " " << 2.4<< std::endl;
				if(eta_cut < 5.0 && jet->eta > 2.4) std::cout << "Eta_cut " << abs(2.7) << " Gen Jet summary: " << jet->pdgId << " " << jet->pt << " eta: " << jet->eta << " " << jet->phi << std::endl;
			}
		}
	}
	//std::cout << "////////////////////////////////////////////////" << std::endl;
	return n_gen_jets;
}




