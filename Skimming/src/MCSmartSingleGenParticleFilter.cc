#include "PxlSkimmer/Skimming/interface/MCSmartSingleGenParticleFilter.h"

#include <iostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

using namespace edm;
using namespace std;


MCSmartSingleGenParticleFilter::MCSmartSingleGenParticleFilter( const ParameterSet& iConfig ) :
   label_( iConfig.getParameter< InputTag >( "genParSource" ) )
{
   // here do whatever other initialization is needed
   vector< int > defpid;
   defpid.push_back( 0 );
   particleID = iConfig.getUntrackedParameter< vector< int > >( "ParticleID", defpid );

   vector< double > defpmin;
   defpmin.push_back( 0. );
   pMin = iConfig.getUntrackedParameter< vector< double > >( "MinP", defpmin );

   vector< double > defptmin;
   defptmin.push_back( 0. );
   ptMin = iConfig.getUntrackedParameter< vector< double > >( "MinPt", defptmin );

   vector< double > defetamin;
   defetamin.push_back( -10. );
   etaMin = iConfig.getUntrackedParameter< vector< double > >( "MinEta", defetamin );

   vector< double > defetamax;
   defetamax.push_back( 10. );
   etaMax = iConfig.getUntrackedParameter< vector< double > >( "MaxEta", defetamax );

   vector< int > defstat;
   defstat.push_back( 0 );
   status = iConfig.getUntrackedParameter< vector< int > >( "Status", defstat );

   vector< double > defDecayRadiusmin;
   defDecayRadiusmin.push_back( -1. );
   decayRadiusMin = iConfig.getUntrackedParameter< vector< double > >( "MinDecayRadius", defDecayRadiusmin );

   vector< double > defDecayRadiusmax;
   defDecayRadiusmax.push_back( 1.e5 );
   decayRadiusMax = iConfig.getUntrackedParameter< vector< double > >( "MaxDecayRadius", defDecayRadiusmax );

   vector< double > defDecayZmin;
   defDecayZmin.push_back( -1.e5 );
   decayZMin = iConfig.getUntrackedParameter< vector< double > >( "MinDecayZ", defDecayZmin );

   vector< double > defDecayZmax;
   defDecayZmax.push_back( 1.e5 );
   decayZMax = iConfig.getUntrackedParameter< vector< double > >( "MaxDecayZ", defDecayZmax );

    // check for same size
   if ( ( pMin.size() > 1   && particleID.size() != pMin.size() )
      ||  ( ptMin.size() > 1  && particleID.size() != ptMin.size() )
      ||  ( etaMin.size() > 1 && particleID.size() != etaMin.size() )
      ||  ( etaMax.size() > 1 && particleID.size() != etaMax.size() )
      ||  ( status.size() > 1 && particleID.size() != status.size() )
      ||  ( decayRadiusMin.size() > 1 && particleID.size() != decayRadiusMin.size() )
      ||  ( decayRadiusMax.size() > 1 && particleID.size() != decayRadiusMax.size() )
      ||  ( decayZMin.size() > 1 && particleID.size() != decayZMin.size() )
      ||  ( decayZMax.size() > 1 && particleID.size() != decayZMax.size() ) ) {
      cout << "WARNING: MCPROCESSFILTER : size of MinPthat and/or MaxPthat not matching with ProcessID size!!" << endl;
   }

   // if pMin size smaller than particleID , fill up further with defaults
   if( particleID.size() > pMin.size() ){
      vector< double > defpmin2;
      for( unsigned int i = 0; i < particleID.size(); i++ ){ defpmin2.push_back( 0. ); }
      pMin = defpmin2;
   }

   // if ptMin size smaller than particleID , fill up further with defaults
   if( particleID.size() > ptMin.size() ){
      vector< double > defptmin2;
      for( unsigned int i = 0; i < particleID.size(); i++ ){ defptmin2.push_back( 0. ); }
      ptMin = defptmin2;
   }

   // if etaMin size smaller than particleID , fill up further with defaults
   if( particleID.size() > etaMin.size() ){
      vector< double > defetamin2;
      for( unsigned int i = 0; i < particleID.size(); i++ ){ defetamin2.push_back( -10. ); }
      etaMin = defetamin2;
   }
   // if etaMax size smaller than particleID , fill up further with defaults
   if( particleID.size() > etaMax.size() ){
      vector< double > defetamax2;
      for( unsigned int i = 0; i < particleID.size(); i++ ){ defetamax2.push_back( 10. ); }
      etaMax = defetamax2;
   }

   // if status size smaller than particleID , fill up further with defaults
   if( particleID.size() > status.size() ){
      vector< int > defstat2;
      for( unsigned int i = 0; i < particleID.size(); i++ ){ defstat2.push_back( 0 ); }
      status = defstat2;
   }

   // if decayRadiusMin size smaller than particleID , fill up further with defaults
   if( particleID.size() > decayRadiusMin.size() ){
      vector< double > decayRadiusmin2;
      for( unsigned int i = 0; i < particleID.size(); i++ ){ decayRadiusmin2.push_back( -10. ); }
      decayRadiusMin = decayRadiusmin2;
   }

   // if decayRadiusMax size smaller than particleID , fill up further with defaults
   if( particleID.size() > decayRadiusMax.size() ){
      vector< double > decayRadiusmax2 ;
      for( unsigned int i = 0; i < particleID.size(); i++ ){ decayRadiusmax2.push_back( 1.e5 ); }
      decayRadiusMax = decayRadiusmax2;
   }

   // if decayZMin size smaller than particleID , fill up further with defaults
   if( particleID.size() > decayZMin.size() ){
      vector< double > decayZmin2 ;
      for( unsigned int i = 0; i < particleID.size(); i++ ){ decayZmin2.push_back( -1.e5 ); }
      decayZMin = decayZmin2;
   }

   // if decayZMax size smaller than particleID , fill up further with defaults
   if( particleID.size() > decayZMax.size() ){
      vector< double > decayZmax2 ;
      for( unsigned int i = 0; i < particleID.size(); i++ ){ decayZmax2.push_back( 1.e5 ); }
      decayZMax = decayZmax2;
   }
}


MCSmartSingleGenParticleFilter::~MCSmartSingleGenParticleFilter()
{
}

bool MCSmartSingleGenParticleFilter::filter( Event& iEvent, const EventSetup& iSetup )
{
   bool accepted = false;
   Handle< reco::GenParticleCollection > genParticleHandel;
   iEvent.getByLabel( label_, genParticleHandel );

   for( reco::GenParticleCollection::const_iterator pa = genParticleHandel->begin(); pa != genParticleHandel->end(); ++ pa ) {
      const reco::GenParticle* p = ( const reco::GenParticle* ) &( *pa );
      const reco::Candidate::LorentzVector p4 = p->p4();

      for( unsigned int i = 0; i < particleID.size(); i++ ) {
         if( p->pdgId() == particleID.at( i ) || particleID.at( i ) == 0 ) {

            if( p4.P() > pMin.at( i ) && p4.Pt() > ptMin.at( i )
               && p4.Eta() > etaMin.at( i ) && p4.Eta() < etaMax.at( i )
               && ( p->status() == status.at( i ) || status.at( i ) == 0 ) ) {

               //
               // WARNING: When converted from HepMCProduct to GenParticle, the
               // vertex is set to (0, 0, 0) if production_vertex returns a
               // null pointer
               //
               double decx = p->vertex().X();
               double decy = p->vertex().Y();
               double decrad = sqrt( decx*decx + decy*decy );

               if( decrad < decayRadiusMin.at( i ) ) continue;
               if( decrad > decayRadiusMax.at( i ) ) continue;

               double decz = p->vertex().Y();
               if( decz < decayZMin.at( i ) ) continue;
               if( decz > decayZMax.at( i ) ) continue;

               accepted = true;
            }
         }
      }
   }

   return accepted;
}

#include "FWCore/Framework/interface/MakerMacros.h"

//define this as a plug-in
DEFINE_FWK_MODULE(MCSmartSingleGenParticleFilter);

