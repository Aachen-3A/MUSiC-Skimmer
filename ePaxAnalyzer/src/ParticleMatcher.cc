#include "ePaxDemo/ePaxAnalyzer/interface/ParticleMatcher.hh"

using namespace std;
using namespace pxl;

// ------------ matching Method ------------

void ParticleMatcher::matchObjects(EventView* GenView, EventView* RecView, const std::string& _JetAlgo, const std::string& _METType) {
   // FIXME: Make code more generic! Generate a list of all Particle types
   std::vector<std::string> typeList;
   typeList.push_back("Muon"); 
   typeList.push_back("Ele"); 
   typeList.push_back("Gamma"); 
   typeList.push_back(_JetAlgo); 
   typeList.push_back(_METType);

   // containers to keep the filtered gen/rec particles
   vector<Particle*> gen_particles;
   vector<Particle*> rec_particles;
   for (std::vector<std::string>::const_iterator partType = typeList.begin(); partType != typeList.end(); ++partType) {
      // Choose name filter criterion
      gen_particles.clear();
      rec_particles.clear();
      ParticleNameCriterion crit(*partType);
      GenView->getObjectsOfType<Particle, PtComparator>(gen_particles, crit);
      RecView->getObjectsOfType<Particle, PtComparator>(rec_particles, crit);
      makeMatching(gen_particles, rec_particles, _METType);
   }   
}

// ------------ implementation of the matching Gen <--> Rec ------------

void ParticleMatcher::makeMatching(std::vector<Particle*>& gen_particles, std::vector<Particle*>& rec_particles, const string& _METType) {
   // First set for Gen all Matches to -1 and reset bools:
   for (std::vector<Particle*>::iterator gen_iter = gen_particles.begin(); gen_iter != gen_particles.end(); gen_iter++) {
      (*gen_iter)->setUserRecord<int>("Match", -1);
      (*gen_iter)->setUserRecord<bool>("hctaM", false);      
   }
   // same for Rec
   for (std::vector<Particle*>::iterator rec_iter = rec_particles.begin(); rec_iter != rec_particles.end(); rec_iter++) {
      (*rec_iter)->setUserRecord<int>("Match", -1);
      (*rec_iter)->setUserRecord<bool>("hctaM", false);      
   }
   unsigned int num_gen = gen_particles.size();
   unsigned int num_rec = rec_particles.size();
   
   // we need at least one Gen and one Rec to perform matching!      
   if (num_gen > 0 && num_rec > 0) {
      unsigned int col = 0;
      unsigned int row = 0;
      std::string particle;

      if (_fDebug > 1) cout << "Found " << num_gen << " Gen Objects and " << num_rec << " Rec Objects" << endl;

      TMatrixT<double> DistanzMatrix(num_gen, num_rec);
      for (std::vector<Particle*>::iterator gen_iter = gen_particles.begin(); gen_iter != gen_particles.end(); gen_iter++) {
         col = 0;
	 for (std::vector<Particle*>::iterator rec_iter = rec_particles.begin(); rec_iter != rec_particles.end(); rec_iter++) {
	    // Calculate the distance
	    if (_fDebug > 0) { 
	       cout << "Gen: " << (*gen_iter)->print(0); 
	       cout << "Rec: " << (*rec_iter)->print(0);
               cout << "Distance: " << (*gen_iter)->getVector().deltaR(&((*rec_iter)->getVector())) << endl;
	    }
	    DistanzMatrix(row,col) = (*gen_iter)->getVector().deltaR(&((*rec_iter)->getVector()));
	    col++;
	 }
	 row++;
      }
      
      
      if (_fDebug > 0) DistanzMatrix.Print(); 

      //define value in dR which defines matching
      double DeltaRMatching = _DeltaR_Particles;
      particle = (gen_particles.front())->getName();
      if (particle == _METType) DeltaRMatching = _DeltaR_MET;

      // go through every row and pushback index of Rec with smallest Distance
      for (unsigned int irow = 0; irow < num_gen; irow++) {
       	 int matched = SmallestRowElement(&DistanzMatrix, irow, DeltaRMatching);
	 gen_particles[irow]->setUserRecord<int>("Match", matched);
	 if (_fDebug > 0) cout << "GenObject " << irow << " is matched with " << matched << endl;

	 if (matched != -1){
	 	//redundant information with softlink, should replace the UserRecords after testing
	 	gen_particles[irow]->linkSoft(rec_particles[matched],"priv-gen-rec");

	 	//cout << "pt of the private matched " << rec_particles[matched]->getName() << " rec: " << rec_particles[matched]->getPt() << endl; //temporary!
	 	//cout << "pt of the private matched " << gen_particles[irow]->getName() << " gen: " << gen_particles[irow]->getPt() << endl; //temporary!

      	 	rec_particles[matched]->setUserRecord<bool>("hctaM", true);
	        if (_fDebug > 0) cout << "RecObject " << matched << " has matching Gen " << endl;      
         }
      }

      for (unsigned int icol = 0; icol < num_rec; icol++) {
         //define value in dR which defines matching
	 int matched = SmallestColumnElement(&DistanzMatrix, icol, DeltaRMatching);
         rec_particles[icol]->setUserRecord<int>("Match", matched);
	 if (_fDebug > 0) cout << "RecObject " << icol << " is matched with " << matched << endl;
	
	 if (matched != -1) {
	    //redundant information with softlink, should replace the UserRecords after testing
	    rec_particles[icol]->linkSoft(gen_particles[matched],"priv-rec-gen");
            gen_particles[matched]->setUserRecord<bool>("hctaM", true);
	    if (_fDebug > 0) cout << "GenObject " << matched << " has matching Rec " << endl;           
         }
      }
   }
}

// ---------------------- Helper Method ------------------------------

int ParticleMatcher::SmallestRowElement(TMatrixT<double>* matrix, const unsigned int& row, const double& DeltaRMatching) {

   // loop over row and return index of smallest element
   double element = (*matrix)(row, 0);
   int index = 0;
   for (int i = 1; i < matrix->GetNcols(); i++) {
      if ((*matrix)(row, i) < element) {
         element = (*matrix)(row,i);
	 index = i;
      }
   }
   if (element > DeltaRMatching) index = -1;    
   return index;
}

// ---------------------- Helper Method ------------------------------

int ParticleMatcher::SmallestColumnElement(TMatrixT<double>* matrix, const unsigned int& col, const double& DeltaRMatching) {

   // loop over row and return index of smallest element
   double element = (*matrix)(0, col);
   int index = 0;
   for (int i = 1; i < matrix->GetNrows(); i++) {
      if ((*matrix)(i, col) < element) {
         element = (*matrix)(i,col);
	 index = i;
      }
   }    
   if (element > DeltaRMatching) index = -1;
   return index;
}
