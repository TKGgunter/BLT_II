#include "BLT_II/BLTAnalysis/interface/tgUtility.hh"
#include <cmath>



void progress_print(ProgressBar * progressbar){
		std::string bar = "";
		int current_bar_length = (int) ((float)progressbar->epoch_completed/ (float)progressbar->max_epochs * 50.0);
		for (int i = 0; i < 50; i++){
				if (i < current_bar_length) bar += "=";
				else if (i == current_bar_length) bar += ">";
				else bar += " " ;
		}
		printf("\rSeconds to completion %f   Avg Time: %f    Percent complete %f  %ld / %ld  [%s]", 
																																		progressbar->average_time_per_epoch * (progressbar->max_epochs - progressbar->epoch_completed), 
																																		progressbar->average_time_per_epoch , 
																																	 (float)progressbar->epoch_completed/ (float)progressbar->max_epochs * 100.0,
																																	 progressbar->epoch_completed, progressbar->max_epochs,
																																	 bar.c_str() );
		fflush(stdout);
}

void progress_update( ProgressBar * progressbar){
		progressbar->completed_events++;
		if ( progressbar->completed_events >= progressbar->prev_completed_events + progressbar->epoch){
				progressbar->prev_completed_events = progressbar->completed_events;

				double current_time = progressbar->timer.RealTime();
				progressbar->average_time_per_epoch += current_time;
				progressbar->average_time_per_epoch =  progressbar->average_time_per_epoch  / 2.0;


				progressbar->epoch_completed++;
				progressbar->timer.Reset();
				progressbar->timer.Start();
				progress_print(progressbar);
		}
}

void testFunc(std::string testStr, float testNumb){
	if( testNumb == -99) std::cout << "Is " << testStr << " working." << std::endl;

	else std::cout << "Is " << testStr << " working. Are you expecting the folling value: " << testNumb << std::endl; 
}

//======================================

float qtCalc( TGPhysObject* obj1, TGPhysObject* obj2 ){
	float qt = 0;
	qt = sqrt(pow(obj1->pt()*sin(obj1->phi()) + obj2->pt()*sin(obj2->phi()), 2) + pow(obj1->pt()*cos(obj1->phi()) + obj2->pt()*cos(obj2->phi()), 2));
	return qt; 
}


//=====================================

float deltaPhi(float dPhi){
	if(dPhi != dPhi){
		printf("Phi%f\n", dPhi);
		dPhi = 0;
	}
	dPhi = TVector2::Phi_mpi_pi(dPhi);
	//if( dPhi > M_PI) dPhi = abs(dPhi - M_PI);
	
	return dPhi;
}

float dR(baconhep::TMuon* lep1, baconhep::TJet* jet){
	float dR_ = 0.;
	float dPhi = deltaPhi(lep1->pfPhi-jet->phi);
	float dEta = lep1->pfEta - jet->eta;
	dR_ = sqrt(pow(dPhi, 2) + pow(dEta, 2));
	return dR_;
}

float dR(baconhep::TElectron* lep1, baconhep::TJet* jet){  
	float dR_ = 0.; 
	float dPhi = deltaPhi(lep1->pfPhi-jet->phi);   
	float dEta = lep1->pfEta - jet->eta; 
	dR_ = sqrt(pow(dPhi, 2) + pow(dEta, 2));
	return dR_;
}

float dR(TGPhysObject* obj1, TGPhysObject* obj2){  
	float dR_ = 0.; 
	float dPhi = deltaPhi(obj1->phi() - obj2->phi());   
	float dEta = obj1->eta() - obj2->eta(); 
	dR_ = sqrt(pow(dPhi, 2) + pow(dEta, 2));
	return dR_;
}
float dR(TGPhysObject* obj1, baconhep::TJet* obj2){  
	float dR_ = 0.; 
	float dPhi = deltaPhi(obj1->phi() - obj2->phi);   
	float dEta = obj1->eta() - obj2->eta; 
	dR_ = sqrt(pow(dPhi, 2) + pow(dEta, 2));
	return dR_;
}

float dR( float phi1, float phi2, float eta1, float eta2){
	float dR_  = 0.;
	float dPhi = deltaPhi( phi1-phi2);
	float dEta = eta1-eta2;
	dR_ = sqrt( pow(dPhi, 2.) + pow(dEta, 2) ); 
	return dR_;
}

/*
float dR(baconhep::TMuon* lep1, TGenJet* jet){                                                     
  float dR_ = 0.;
  float dPhi = deltaPhi(lep1->pfPhi-jet->pfPhi);                                        
  float dEta = lep1->pfEta - jet->pfEta;
  dR_ = sqrt(pow(dPhi, 2) + pow(dEta, 2));                                              
  return dR_;                                                                           
}                                                                                       

float dR(baconhep::TElectron* lep1, TGenJet* jet){                                                 
  float dR_ = 0.; 
  float dPhi = deltaPhi(lep1->pfPhi-jet->pfPhi);                                        
  float dEta = lep1->pfEta - jet->pfEta; 
  dR_ = sqrt(pow(dPhi, 2) + pow(dEta, 2));                                              
  return dR_;                                                                           
}

float dR(TGenParticle* gParticle, TGenJet* jet){
  float dR_ = 0.;
  float dPhi = deltaPhi(gParticle->pfPhi-jet->pfPhi);
  float dEta = gParticle->pfEta - jet->pfEta;
  dR_ = sqrt(pow(dPhi, 2) + pow(dEta, 2));
  return dR_;
}

float dR(TJet* jet, TGenJet* genJet){
  float dR_ = 0.;
  float dPhi = deltaPhi(genJet->pfPhi-jet->pfPhi);                                    
  float dEta = genJet->pfEta - jet->pfEta;                                            
  dR_ = sqrt(pow(dPhi, 2) + pow(dEta, 2));                                               
  return dR_;  
}
*/
/*
float dR(baconhep::TMuon* lep1, TMET* met){                                                     
  float dR_ = 0.;
	//dR_ = met->DeltaR(*lep1);
  return dR_;                                                                           
}                                                                                       

float dR(baconhep::TElectron* lep1, TMET* met){                                                 
  float dR_ = 0.; 
  //dR_ = met->DeltaR(*lep1);                                        
  return dR_;                                                                           
}                                                                                       
*/
//===================================
/*
float dotLepMet(baconhep::TMuon* lep, TMET* met ){
	float dotProduct = 0;
	dotProduct = lep->pfPt*sin(lep->pfPhi) * met->pfPt*sin(met->pfPhi) + lep->pfPt*sin(lep->pfPhi) * met->pfPt*sin(met->pfPhi);
	return dotProduct;
}

float dotLepMet(baconhep::TElectron* lep, TMET* met){
	float dotProduct = 0;
	dotProduct = lep->pfPt*sin(lep->pfPhi) * met->pfPt*sin(met->pfPhi) + lep->pfPt*sin(lep->pfPhi) * met->pfPt*sin(met->pfPhi);  
	return dotProduct;
}
*/  
//====================================




void tgSort(std::vector<TGPhysObject>& objList){
  float pfPt_1 = 0.;
  float pfPt_2 = 0.;
  for(unsigned int i=0; i < objList.size(); i++){

    auto obj_1 = objList.at(i);
    if (i == 0) continue;
    for(unsigned int j = 0; j < i ; j++ ){
      auto obj_2 = objList.at(j);
      pfPt_1 = obj_1.pt();
      pfPt_2 = obj_2.pt();
      if(pfPt_2 < pfPt_1 ){
        objList.insert(objList.begin()+j, obj_1);
        objList.erase(objList.begin()+i+1);
        continue;
      }
    }
  }
}


//================================================

void tgCleanVector(std::vector<TGPhysObject>& objList){
		for(unsigned int i = 0 ; i < objList.size()-1; i++ ){
				if(objList.size() <= i) break;
				TGPhysObject obj1 = (TGPhysObject) objList.at(i);

				for(unsigned int j = 1+i ; j < objList.size(); j++ ){
						if( objList.size() <= i || objList.size() <= j)  break;
						TGPhysObject obj2 = (TGPhysObject) objList.at(j);


						if (obj1.pt() == obj2.pt()){
								objList.erase(objList.begin() + i);
								i=i-1;
						}
				}
		}
}


