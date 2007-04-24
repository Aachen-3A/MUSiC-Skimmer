#ifndef ParticleMatcher_hh
#define ParticleMatcher_hh

/*
Class which perform the matching between generator level particle with
reconstructed particles. The matching is based on a delta R algo. Each gen
particle points to the best matching rec particle and vice versa. If the best
matching particle has a distance large than deltaR > 0.2 the particle is
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
      void matchObjects(pxl::EventViewRef GenView, pxl::EventViewRef RecView); 
      void makeMatching(pxl::ParticleFilter& GenFilter, pxl::ParticleFilter& RecFilter);
      
   private:
      // Some helper methods
      int SmallestRowElement(TMatrixT<double>* matrix, int row);   
      int SmallestColumnElement(TMatrixT<double>* matrix, int col);
      int fDebug; 
};
#endif
