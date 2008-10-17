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
      ParticleMatcher(double DeltaR_Particles = 0.2, double DeltaR_MET = 0.5) : 
            _DeltaR_Particles(DeltaR_Particles), _DeltaR_MET(DeltaR_MET), _fDebug(0) {};
      // Destruktor
      ~ParticleMatcher() {;};
      // Match method
      void matchObjects(pxl::EventView* GenView, pxl::EventView* RecView); 
      void makeMatching(std::vector<pxl::Particle*>& gen_particles, std::vector<pxl::Particle*>& rec_particles);
      
   private:
      // Some helper methods
      int SmallestRowElement(TMatrixT<double>* matrix, const unsigned int& row, const double& DeltaRMatching);   
      int SmallestColumnElement(TMatrixT<double>* matrix, const unsigned int& col, const double& DeltaRMatching);
      //variable to define dR which decides matching
      double _DeltaR_Particles;
      double _DeltaR_MET;
      int _fDebug; 
};
#endif
