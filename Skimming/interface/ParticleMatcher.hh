#ifndef ParticleMatcher_hh
#define ParticleMatcher_hh

/*
Class which perform the matching between generator level particle with
reconstructed particles. The matching is based on a delta R algo. Each gen
particle points to the best matching rec particle and vice versa. If the best
matching particle has a distance large than the given limits for DeltaR, DeltaPtoPt or DeltaCharge the particle is
declared to have no match. For unmatched particles Match UserRecord is set to -1.
*/

#include "MUSiCProject/Pxl/interface/Pxl.h"
#include "MUSiCProject/Skimming/interface/jet_def.h"
#include <iostream> 
#include "TMatrixT.h"
#include <vector>

class ParticleMatcher {

   public: 
      // Konstruktor
      ParticleMatcher(double DeltaR_Particles = 0.2, double DeltaR_MET = 0.5, double DeltaPtoPt = 1000000.0, double DeltaCharge = 10.0, int fDebug = 0) : 
            _DeltaR_Particles(DeltaR_Particles), _DeltaR_MET(DeltaR_MET), _DeltaPtoPt(DeltaPtoPt), _DeltaCharge(DeltaCharge), _fDebug(fDebug) {};
      // Destruktor
      ~ParticleMatcher() {;};
      // Match method
      void matchObjects( pxl::EventView *GenView, pxl::EventView *RecView, const std::vector< jet_def > &jet_infos, const std::string &_METType); 
      void makeMatching(std::vector<pxl::Particle*>& gen_particles, std::vector<pxl::Particle*>& rec_particles, const std::string& _METType);
      
   private:
      // Some helper methods
      int SmallestRowElement(TMatrixT<double>* matrixDR, TMatrixT<double>* matrixDp, TMatrixT<double>* matrixDC, const unsigned int& row, const double& DeltaRMatching, const double& DeltaChargeMatching, const double& DeltaPtoPtMatching);   
      int SmallestColumnElement(TMatrixT<double>* matrixDR, TMatrixT<double>* matrixDp,  TMatrixT<double>* matrixDC, const unsigned int& col, const double& DeltaRMatching, const double& DeltaChargeMatching, const double& DeltaPtoPtMatching);
      //variable to define dR which decides matching
	double _DeltaR_Particles;
	double _DeltaR_MET;
      	double _DeltaPtoPt;
	double _DeltaCharge;
      int _fDebug; 
};
#endif
