#include "recHitFilter.h"




 

//! ctor
recHitFilter::recHitFilter(const edm::ParameterSet& iConfig)
{}

// ----------------------------------------------------------------






//! dtor
recHitFilter::~recHitFilter()
{}

// ----------------------------------------------------------------






void recHitFilter::beginJob() 
{}

// ----------------------------------------------------------------






void recHitFilter::endJob() 
{}

// ----------------------------------------------------------------






//! check the conversionsCollection size
bool recHitFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
//   int run = iEvent.id().run();
//   if (run != 124009 &&  run != 124020 && run != 124022 && run != 124023 && run != 124024) 
//     return false;

  edm::Handle<EcalRecHitCollection> pBarrelEcalRecHits ;
  iEvent.getByLabel ("ecalRecHit","EcalRecHitsEB", pBarrelEcalRecHits) ;
  const EcalRecHitCollection* theBarrelEcalRecHits = pBarrelEcalRecHits.product () ;

  edm::Handle<EcalRecHitCollection> pEndcapEcalRecHits ;
  iEvent.getByLabel ("ecalRecHit","EcalRecHitsEE", pEndcapEcalRecHits) ;
  const EcalRecHitCollection* theEndcapEcalRecHits = pEndcapEcalRecHits.product () ;

  
  if((theBarrelEcalRecHits->size() == 0 && theEndcapEcalRecHits->size() != 0) ||
     (theBarrelEcalRecHits->size() != 0 && theEndcapEcalRecHits->size() == 0))
    {
      //std::cout << "false" << std::endl;
      return false;
    }
  
  else
    {
      //std::cout << "true" << std::endl;
      return true;  
    }
}




DEFINE_FWK_MODULE(recHitFilter);
