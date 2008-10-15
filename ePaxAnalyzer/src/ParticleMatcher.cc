#include "ePaxDemo/ePaxAnalyzer/interface/ParticleMatcher.hh"

using namespace std;
// ------------ matching Method ------------

void ParticleMatcher::matchObjects(pxl::EventView* GenView, pxl::EventView* RecView) {
   // FIXME: Make code more generic! Generate a list of all Particle types
   std::vector<std::string> typeList;
   typeList.push_back("Muon"); typeList.push_back("Ele"); typeList.push_back("Gamma"); typeList.push_back("KtJet"); 
   typeList.push_back("ItCone5Jet"); typeList.push_back("L2L3JESic5Jet");
   typeList.push_back("MET"); typeList.push_back("METCorr");

   //define dR for matching of different physics objects
   DeltaRMET = 0.5;
   DeltaRParticles = 0.2;
   
   pxl::ParticleFilter filter;
   // containers to keep the filtered gen/rec particles
   std::vector<pxl::Particle*> gen_particles;
   std::vector<pxl::Particle*> rec_particles;

   // getTypeList(&typeList, GenView, RecView);
   for (std::vector<std::string>::const_iterator partType = typeList.begin(); partType != typeList.end(); ++partType) {
      // Choose name filter criterion
      gen_particles.clear();
      rec_particles.clear();
      pxl::ParticleNameCriterion crit(*partType);
      filter.apply(&GenView->getObjectOwner(), gen_particles, crit);
      filter.apply(&RecView->getObjectOwner(), rec_particles, crit);
      makeMatching(gen_particles, rec_particles);
      //pxl::ParticleFilter GenFilter(GenView->getObjectOwner(), (*partType));
      //pxl::ParticleFilter RecFilter(RecView->getObjectOwner(), (*partType));
      //makeMatching(GenFilter, RecFilter);
   }   
}

// ------------ implementation of the matching Gen <--> Rec ------------

void ParticleMatcher::makeMatching(std::vector<pxl::Particle*>& gen_particles, std::vector<pxl::Particle*>& rec_particles) {
//void ParticleMatcher::makeMatching(pxl::ParticleFilter& GenFilter, pxl::ParticleFilter& RecFilter) {

   // First set for Gen all Matches to -1 and reset bools:
   for (std::vector<pxl::Particle*>::iterator gen_iter = gen_particles.begin(); gen_iter != gen_particles.end(); gen_iter++) {
      (*gen_iter)->setUserRecord<int>("Match", -1);
      (*gen_iter)->setUserRecord<bool>("hctaM", false);      
   }
   // same for Rec
   for (std::vector<pxl::Particle*>::iterator rec_iter = rec_particles.begin(); rec_iter != rec_particles.end(); rec_iter++) {
      (*rec_iter)->setUserRecord<int>("Match", -1);
      (*rec_iter)->setUserRecord<bool>("hctaM", false);      
   }
   unsigned int num_gen = gen_particles.size();
   unsigned int num_rec = rec_particles.size();
   
   // we need at least one Gen and one Rec to perform matching!      
   if (num_gen > 0 && num_rec > 0) {
   //if (GenFilter.getSize() > 0 && RecFilter.getSize() > 0) {
      unsigned int col = 0;
      unsigned int row = 0;
      std::string particle;

//if (fDebug > 1)
cout << "Found " << num_gen << " Gen Objects and " << num_rec << " Rec Objects" << endl;

      TMatrixT<double> DistanzMatrix(num_gen, num_rec);
      for (std::vector<pxl::Particle*>::iterator gen_iter = gen_particles.begin(); gen_iter != gen_particles.end(); gen_iter++) {
      //for (pxl::ParticleFilterIterator gen_iter(GenFilter); !gen_iter.isDone(); gen_iter.next()) {
         col = 0;
	 for (std::vector<pxl::Particle*>::iterator rec_iter = rec_particles.begin(); rec_iter != rec_particles.end(); rec_iter++) {
         //for (pxl::ParticleFilterIterator rec_iter(RecFilter); !rec_iter.isDone(); rec_iter.next()) {
	    // Calculate the distance
            //pxl::ParticleWkPtr gen_pa = gen_iter.wkPtr();
	    //pxl::ParticleWkPtr rec_pa = rec_iter.wkPtr();
	    if (fDebug > 0) { 
	       cout << "Gen: (" << (*gen_iter)->getPx() << ","
	            << (*gen_iter)->getPy() << ","
		    << (*gen_iter)->getPz() << ","
		    << (*gen_iter)->getE() << ")" << endl; 
	       cout << "Rec: (" << (*rec_iter)->getPx() << ","
	            << (*rec_iter)->getPy() << ","
		    << (*rec_iter)->getPz() << ","
		    << (*rec_iter)->getE() << ")" << endl;
               cout << "Distance: " << (*gen_iter)->getVector().deltaR(&((*rec_iter)->getVector())) << endl;
	    }
	    DistanzMatrix(row,col) = (*gen_iter)->getVector().deltaR(&((*rec_iter)->getVector()));
	    //DistanzMatrix(row,col) = gen_pa.get().vector(pxl::get).deltaR(rec_pa.get().vector(pxl::get));
	    col++;
	 }
	 row++;
      }
      
      
      if (fDebug > 0) DistanzMatrix.Print(); 

	//define value in dR which defines matching
	DeltaRMatching = DeltaRParticles;
	particle = (gen_particles.front())->getName();
	if (particle == "MET" || particle == "METCorr") DeltaRMatching = DeltaRMET;

      // go through every row and pushback index of Rec with smallest Distance
      for (unsigned int irow = 0; irow < num_gen; irow++) {
	 //pxl::ParticleFilterIterator tmp_gen_iter(GenFilter);
	 //pxl::ParticleWkPtr tmp_gen_pa = tmp_gen_iter.wkPtr();
       	 int matched = SmallestRowElement(&DistanzMatrix, irow);
	 gen_particles[irow]->setUserRecord<int>("Match", matched);
	 if (fDebug > 0) cout << "GenObject " << irow << " is matched with " << matched << endl;

	 if (matched != -1){
	 	//redundant information with softlink, should replace the UserRecords after testing
	 	gen_particles[irow]->linkSoft(rec_particles[matched],"priv-gen-rec");

	 	cout << "pt of the private matched rec: " << rec_particles[matched]->getPt() << endl; //temporary!
	 	cout << "pt of the private matched gen: " << gen_particles[irow]->getPt() << endl; //temporary!

      	 	rec_particles[matched]->setUserRecord<bool>("hctaM", true);
	        if (fDebug > 0) cout << "RecObject " << matched << " has matching Gen " << endl;      
         }
      }

      for (unsigned int icol = 0; icol < num_rec; icol++) {
         //define value in dR which defines matching
	 //pxl::ParticleFilterIterator tmp_gen_iter(RecFilter);
	 //pxl::ParticleWkPtr tmp_gen_pa = tmp_gen_iter.wkPtr();
	 int matched = SmallestColumnElement(&DistanzMatrix, icol);
         rec_particles[icol]->setUserRecord<int>("Match", matched);
	 if (fDebug > 0) cout << "RecObject " << icol << " is matched with " << matched << endl;
	
	 if (matched != -1) {
	 	//redundant information with softlink, should replace the UserRecords after testing
	 	rec_particles[icol]->linkSoft(gen_particles[matched],"priv-rec-gen");
	
         	gen_particles[matched]->setUserRecord<bool>("hctaM", true);
		if (fDebug > 0) cout << "GenObject " << matched << " has matching Rec " << endl;           
         }
     }
  }
}


// ------------ implementation of the matching Gen <--> Rec ------------
/*
void ParticleMatcher::makeMatching(pxl::ParticleFilter& GenFilter, pxl::ParticleFilter& RecFilter) {
   // First set all Matches to -1 and reset bools:
   for (pxl::ParticleFilterIterator gen_iter(GenFilter); !gen_iter.isDone(); gen_iter.next()) {
      pxl::ParticleWkPtr pa = gen_iter.wkPtr();
      pa.set().setUserRecord<int>("Match", -1);
      pa.set().setUserRecord<bool>("hctaM", false);
   }
   for (pxl::ParticleFilterIterator rec_iter(RecFilter); !rec_iter.isDone(); rec_iter.next()) {
      pxl::ParticleWkPtr pa = rec_iter.wkPtr();
      pa.set().setUserRecord<int>("Match", -1);
      pa.set().setUserRecord<bool>("hctaM", false);
   }
   // we need at least one Gen and one Rec to perform matching!      
   if (GenFilter.getSize() > 0 && RecFilter.getSize() > 0) {
      int col = 0;
      int row = 0;
      std::string particle;
      if (fDebug > 1) cout << "Found " << GenFilter.getSize() << " Gen Objects and " << RecFilter.getSize() << " Rec Objects" << endl;
      TMatrixT<double> DistanzMatrix(GenFilter.getSize(),RecFilter.getSize());
      for (pxl::ParticleFilterIterator gen_iter(GenFilter); !gen_iter.isDone(); gen_iter.next()) {
         col = 0;
         for (pxl::ParticleFilterIterator rec_iter(RecFilter); !rec_iter.isDone(); rec_iter.next()) {
	    // Calculate the distance
            pxl::ParticleWkPtr gen_pa = gen_iter.wkPtr();
	    pxl::ParticleWkPtr rec_pa = rec_iter.wkPtr();
	    if (fDebug > 0) { 
	       cout << "Gen: (" << gen_pa.get().vector(pxl::get).getPx() << ","
	            << gen_pa.get().vector(pxl::get).getPy() << ","
		    << gen_pa.get().vector(pxl::get).getPz() << ","
		    << gen_pa.get().vector(pxl::get).getE() << ")" << endl; 
	       cout << "Rec: (" << rec_pa.get().vector(pxl::get).getPx() << ","
	            << rec_pa.get().vector(pxl::get).getPy() << ","
		    << rec_pa.get().vector(pxl::get).getPz() << ","
		    << rec_pa.get().vector(pxl::get).getE() << ")" << endl;
               cout << "Distance: " << gen_pa.get().vector(pxl::get).deltaR(rec_pa.get().vector(pxl::get)) << endl;
	    }
	    DistanzMatrix(row,col) = gen_pa.get().vector(pxl::get).deltaR(rec_pa.get().vector(pxl::get));
	    col++;
	 }
	 row++;
      }
      if (fDebug > 0) DistanzMatrix.Print(); 
      // go through every row and pushback index of Rec with smallest Distance
      for (int irow = 0; irow < GenFilter.getSize(); irow++) {
	 //define value in dR which defines matching
	 DeltaRMatching = DeltaRParticles;
	 pxl::ParticleFilterIterator tmp_gen_iter(GenFilter);
	 pxl::ParticleWkPtr tmp_gen_pa = tmp_gen_iter.wkPtr();
	 particle = tmp_gen_pa.get().getName();
	 if(particle == "MET" || particle == "METCorr") DeltaRMatching = DeltaRMET;
         // better implementation then always iterate over all particles?????
	 int matched = SmallestRowElement(&DistanzMatrix, irow);
*/
	 //TEMPORARY FIX to always match met
	 /*pxl::ParticleFilterIterator tmp_gen_iter(GenFilter);
	 pxl::ParticleWkPtr tmp_gen_pa = tmp_gen_iter.wkPtr();
	 std::string particle = tmp_gen_pa.get().getName();
	 if(particle == "MET" || particle == "METCorr") matched = 0;*/
/*
	 int count = 0;
	 bool found = false;
         pxl::ParticleFilterIterator gen_iter(GenFilter);
	 while (!gen_iter.isDone() && !found) {
	    if (count == irow) {
	       pxl::ParticleWkPtr pa = gen_iter.wkPtr();
               pa.set().setUserRecord<int>("Match", matched);
	       found = true;
               if (fDebug > 0) cout << "GenObject " << irow << " is matched with " << matched << endl;     
	    }
	    count++;
	    gen_iter.next();
	 }
         if (matched != -1) {
            count = 0;
            found = false;
            pxl::ParticleFilterIterator rec_iter(RecFilter);
            while (!rec_iter.isDone() && !found) {
               if (count == matched) {
                  pxl::ParticleWkPtr pa = rec_iter.wkPtr();
                  pa.set().setUserRecord<bool>("hctaM", true);
                  found = true;
                  if (fDebug > 0) cout << "RecObject " << matched << " has matching Gen " << endl;      
               }
               count++;
               rec_iter.next();
            }
         } 
      }
      for (int icol = 0; icol < RecFilter.getSize(); icol++) {
         //define value in dR which defines matching
	 DeltaRMatching = DeltaRParticles;
	 pxl::ParticleFilterIterator tmp_gen_iter(RecFilter);
	 pxl::ParticleWkPtr tmp_gen_pa = tmp_gen_iter.wkPtr();
	 particle = tmp_gen_pa.get().getName();
	 if(particle == "MET" || particle == "METCorr") DeltaRMatching = DeltaRMET;
         // better implementation then always iterate over all particles?????
	 int matched = SmallestColumnElement(&DistanzMatrix, icol);
*/	 
	 //TEMPORARY FIX to always match met
	 /*pxl::ParticleFilterIterator tmp_rec_iter(RecFilter);
	 pxl::ParticleWkPtr tmp_rec_pa = tmp_rec_iter.wkPtr();
	 std::string particle = tmp_rec_pa.get().getName();
	 if(particle == "MET" || particle == "METCorr") matched = 0;*/
/*
         int count = 0;
	 bool found = false;
         pxl::ParticleFilterIterator rec_iter(RecFilter); 
	 while (!rec_iter.isDone() && !found) {
	    if (count == icol) {
	       pxl::ParticleWkPtr pa = rec_iter.wkPtr();
               pa.set().setUserRecord<int>("Match", matched);
	       found = true;
	       if (fDebug > 0) cout << "RecObject " << icol << " is matched with " << matched << endl;
	    }
	    count++;
	    rec_iter.next();
	 }
         if (matched != -1) {
            count = 0;
            found = false;
            pxl::ParticleFilterIterator gen_iter(GenFilter);
            while (!gen_iter.isDone() && !found) {
               if (count == matched) {
                  pxl::ParticleWkPtr pa = gen_iter.wkPtr();
                  pa.set().setUserRecord<bool>("hctaM", true);
                  found = true;
                  if (fDebug > 0) cout << "GenObject " << matched << " has matching Rec " << endl;           
               }
               count++;
               gen_iter.next();
            }
         }
      }
   }
}
*/

// ---------------------- Helper Method ------------------------------

int ParticleMatcher::SmallestRowElement(TMatrixT<double>* matrix, unsigned int row) {

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

int ParticleMatcher::SmallestColumnElement(TMatrixT<double>* matrix, unsigned int col) {

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
