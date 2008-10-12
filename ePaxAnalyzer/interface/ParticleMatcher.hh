#ifndef ParticleMatcher_hh
#define ParticleMatcher_hh

/*
Class which perform the matching between generator level particle with
reconstructed particles. The matching is based on a delta R algo. Each gen
particle points to the best matching rec particle and vice versa. If the best
matching particle has a distance large than deltaR > 0.2 ( > 0.5 for MET) the particle is
declared to has no match. For unmatched particles Match UserRecord is set to -1.
*/

#include "ePaxPxl/ePax/interface/ePax.h"
#include <iostream> 
#include "TMatrixT.h"

class ParticleMatcher {

   public: 
      // Konstruktor
      ParticleMatcher() {fDebug = 0;};
      // Destruktor
      ~ParticleMatcher() {;};
      // Match method
      void matchObjects(pxl::EventView* GenView, pxl::EventView* RecView); 
      void makeMatching(std::vector<pxl::Particle*>& gen_particles, std::vector<pxl::Particle*>& rec_particles);
      
   private:
      // Some helper methods
      int SmallestRowElement(TMatrixT<double>* matrix, unsigned int row);   
      int SmallestColumnElement(TMatrixT<double>* matrix, unsigned int col);
      int fDebug; 
      //variable to define dR which decides matching
      double DeltaRMatching;
      double DeltaRMET;
      double DeltaRParticles;
};
#endif
