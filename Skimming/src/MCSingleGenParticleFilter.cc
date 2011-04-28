#include "MUSiCProject/Skimming/interface/MCSingleGenParticleFilter.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <iostream>

using namespace edm;
using namespace std;


MCSingleGenParticleFilter::MCSingleGenParticleFilter( const ParameterSet& iConfig ) :
   label_( iConfig.getParameter< InputTag >( "genParSource" ) )
{
   //here do whatever other initialization is needed
   vector< int > defpid;
   defpid.push_back( 0 );
   particleID = iConfig.getUntrackedParameter< vector< int > >( "ParticleID", defpid );

   vector< double > defptmin;
   defptmin.push_back( 0. );
   ptMin = iConfig.getUntrackedParameter< vector< double > >( "MinPt", defptmin );

   vector< double > defetamin ;
   defetamin.push_back( -10. );
   etaMin = iConfig.getUntrackedParameter< vector< double > >( "MinEta", defetamin );

   vector< double > defetamax;
   defetamax.push_back( 10. );
   etaMax = iConfig.getUntrackedParameter< vector< double > >( "MaxEta", defetamax );

   vector< int > defstat ;
   defstat.push_back( 0 );
   status = iConfig.getUntrackedParameter< vector< int > >( "Status", defstat );


    // check for same size
   if ( ( ptMin.size() > 1 &&  particleID.size() != ptMin.size() )
      ||  ( etaMin.size() > 1 && particleID.size() != etaMin.size() )
      ||  ( etaMax.size() > 1 && particleID.size() != etaMax.size() )
      ||  ( status.size() > 1 && particleID.size() != status.size() ) ) {
      cout << "WARNING: MCPROCESSFILTER : size of MinPthat and/or MaxPthat not matching with ProcessID size!!" << endl;
   }

   // if ptMin size smaller than particleID , fill up further with defaults
   if( particleID.size() > ptMin.size() ){
      vector< double > defptmin2;
      for( unsigned int i = 0; i < particleID.size(); i++ ){ defptmin2.push_back( 0. ); }
      ptMin = defptmin2;
   }

   // if etaMin size smaller than particleID , fill up further with defaults
   if( particleID.size() > etaMin.size() ){
      vector< double > defetamin2 ;
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
}


MCSingleGenParticleFilter::~MCSingleGenParticleFilter()
{
}


// ------------ method called to skim the data  ------------
bool MCSingleGenParticleFilter::filter( Event& iEvent, const EventSetup& iSetup)
{
   bool accepted = false;
   Handle< reco::GenParticleCollection > genParticleHandel;
   iEvent.getByLabel( label_, genParticleHandel );

   for( reco::GenParticleCollection::const_iterator pa = genParticleHandel->begin(); pa != genParticleHandel->end(); ++ pa ) {
      const reco::GenParticle* p = ( const reco::GenParticle* ) &( *pa );
      const reco::Candidate::LorentzVector p4 = p->p4();

      for( unsigned int i = 0; i < particleID.size(); i++ ) {
         if( p->pdgId() == particleID.at( i ) || particleID.at( i ) == 0 ) {
            if ( p4.Pt() > ptMin.at( i ) && p4.Eta() > etaMin.at( i ) && p4.Eta() < etaMax.at( i )
               && ( p->status() == status.at( i ) || status.at( i ) == 0 ) ) {

               accepted = true;
            }
         }
      }
   }

   return accepted;
}

//define this as a plug-in
DEFINE_FWK_MODULE(MCSingleGenParticleFilter);

