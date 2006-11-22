#include "ePaxPxl/ePax/interface/ePax.h"
#ifndef MERGED_PXL
	#include "pcl.hh"
#endif

#include <iostream>
#include <sys/times.h>
#include <unistd.h>
//----------------------------------------------------------------------
namespace pcl {
//----------------------------------------------------------------------
double getCpuTime()
{
   // static method returning system CPU time.
   struct tms cpt;
   times(&cpt);
   return (double)(cpt.tms_utime+cpt.tms_stime)/100.;
}
//----------------------------------------------------------------------
void exception(const std::string& rout, const std::string& msg) {
  std::cerr << rout << ": " << msg << std::endl; 
  std::cerr << "pcl::exception(): throwing pcl::Exception." << std::endl; 
  throw pcl::Exception(rout, msg);
}
//----------------------------------------------------------------------
} // namespace pcl
//----------------------------------------------------------------------
#include <assert.h>

#ifndef MERGED_PXL
	#include "ptl.hh"
#endif
//
//
//----------------------------------------------------------------------
iotl__declareDataTypeExplicit(char,        "\1c",      __char, data, {pcl::BasicIoStreamer::storeBasicTypeChar(_buffer, data);}   , {pcl::BasicIoStreamer::restoreBasicTypeChar(_buffer, data);})
iotl__declareDataTypeExplicit(std::string, "\1s", std__string, data, {pcl::BasicIoStreamer::storeBasicTypeString(_buffer, data);} , {pcl::BasicIoStreamer::restoreBasicTypeString(_buffer, data);})
iotl__declareDataTypeExplicit(bool,        "\1b",      __bool, data, {pcl::BasicIoStreamer::storeBasicTypeBool(_buffer, data);}  , {pcl::BasicIoStreamer::restoreBasicTypeBool(_buffer, data);})
iotl__declareDataTypeExplicit(int,         "\1i",       __int, data, {pcl::BasicIoStreamer::storeBasicTypeInt(_buffer, data);}   , {pcl::BasicIoStreamer::restoreBasicTypeInt(_buffer, data);})
iotl__declareDataTypeExplicit(float,       "\1f",     __float, data, {pcl::BasicIoStreamer::storeBasicTypeFloat(_buffer, data);}  , {pcl::BasicIoStreamer::restoreBasicTypeFloat(_buffer, data);})
iotl__declareDataTypeExplicit(double,      "\1d",    __double, data, {pcl::BasicIoStreamer::storeBasicTypeDouble(_buffer, data);} , {pcl::BasicIoStreamer::restoreBasicTypeDouble(_buffer, data);})
//----------------------------------------------------------------------
void iotl::oStreamer::getEvent(std::ostream& cxxx, char compressionMode) {

  const std::string& strBuffer = _buffer.str();
  
  const char* cBuffer = strBuffer.c_str();
  int   lengthBuffer     = strBuffer.length();
  
  const char* cZip = cBuffer;
  int   lengthZip = lengthBuffer;

  char* cZipSpace = 0;
  unsigned long lengthZipSpace = 0;

  if (compressionMode == ' ') {;} // no compression requires no action...
  else if (compressionMode == '6') {
    // data compression a la Gero, i.e. compression level = 6:
    lengthZipSpace = int(double(lengthBuffer)*1.01+13);
    cZipSpace = new char[lengthZipSpace];
    switch (compress2((Bytef *) cZipSpace, (uLongf *) &lengthZipSpace, (const Bytef *) cBuffer, lengthBuffer, 6)) {
      case Z_MEM_ERROR:
	pcl::exception("iotl::oStreamer::getEvent()", "zlib: not enough memory");
       case Z_BUF_ERROR:
	 pcl::exception("iotl::oStreamer::getEvent()", "zlib: buffer too small");
       case Z_STREAM_ERROR:
	 pcl::exception("iotl::oStreamer::getEvent()", "zlib: level parameter invalid");
       default:
	 break;
       }
    cZip = cZipSpace;
    lengthZip = lengthZipSpace;
  }
  else {pcl::exception("iotl::oStreamer::getEvent()","Invalid compression mode.");}

  dumpMemory(cxxx, iotl__eventMarker, 4);      // Event-Marker
  dumpMemory(cxxx, &compressionMode, 1);        // Compression-Mode
  dumpMemory(cxxx, (const char*)&lengthZip, 4); // Fwd-Jumper
  
  cxxx.rdbuf()->sputn(cZip, lengthZip);          // compressed data block
  cxxx.rdbuf()->pubsync();
  
  dumpMemory(cxxx, (const char*)&lengthZip, 4); // Bwd-Jumper & consistency check
  
  if (cZipSpace) {delete [] cZipSpace;}
  
  _buffer.str("");

}
//
//
//----------------------------------------------------------------------
namespace iotl {
template<> void oStreamer::storeData<ptl::Relations>(const ptl::Relations& relations) {

  storeData<int>(relations.getSize());
  for (ptl::Relations::StlConstIterator iter = relations.begin(); 
                                      iter != relations.end(); 
				      iter++) {
       const ptl::WkPtrBase& ptr = (*iter->second);
       if (ptr.valid()) {storeId(ptr.pointer()->id());}
       else             {storeId(0);}
       }

}
//----------------------------------------------------------------------
template<> void iotl::oStreamer::storeData<ptl::Objects>(const ptl::Objects& objects) {

  // objects
  storeData<int>(objects.getSize());
  for (ptl::Objects::StlConstIterator iter = objects.begin(); iter != objects.end(); iter++) {
  	  storeAbstractObject(**iter);
  	  }

  // index
  storeData<int>(objects._index.getSize());
  for (ptl::Index::StlConstIterator iter = objects._index.begin(); iter != objects._index.end(); iter++) {
      storeData<std::string>(iter->first);
      storeId(iter->second->id());
      }

}
//----------------------------------------------------------------------
template<> void iotl::oStreamer::storeData<ptl::Layout>(const ptl::Layout& properties) {
  storeData<int>(properties._type);
  storeData<int>(properties._style);
  storeData<int>(properties._color);
  storeData<double>(properties._a);
  storeData<double>(properties._b);
  storeData<double>(properties._c);
  storeData<double>(properties._d);
}
//----------------------------------------------------------------------
template<> void iotl::oStreamer::storeData<ptl::ObjectBase>(const ptl::ObjectBase& object) {

  ptl::Id id = object.id();
  storeId(id);
  
  storeData<ptl::Relations>(object._motherRelations);
  storeData<ptl::Relations>(object._daughterRelations);
  
  if (object._ptrLayout) {storeData<bool>( true  ); storeData<ptl::Layout>(*object._ptrLayout);}
  else                   {storeData<bool>( false );}  

}
//----------------------------------------------------------------------
} // namespace iotl
//----------------------------------------------------------------------
//
//
//
//----------------------------------------------------------------------
bool iotl::iStreamer::putEvent(std::istream& cxxx, bool ignore) {

  _buffer.str("");
  
  char eventMarker[4] = {' ',' ',' ','\0'};
  char compressionMode;
  int lengthZip;
  int lengthZipVerify;
   
  redumpMemory(cxxx, eventMarker, 4);       // Event-Marker
  if (!cxxx.good()) {return false;}
  if (strcmp(iotl__eventMarker, eventMarker)) {pcl::exception("iotl::iStreamer::getEvent()","Invalid event marker while reading event header.");}
  
  redumpMemory(cxxx, &compressionMode, 1);        // Compression-Mode
  redumpMemory(cxxx, (char*)&lengthZip, 4);       // Fwd-Jumper
  
  if (ignore) {
  	// ignore the data
  	cxxx.ignore(lengthZip);
  	}
  else if (compressionMode == ' ') {
     // just read the data:
     int bytes = lengthZip;
     int have;
     while (bytes > 0) {
		if (iotl__iStreamer__lengthUnzipBuffer > bytes) {have = bytes;}
		else {have = iotl__iStreamer__lengthUnzipBuffer;}
		redumpMemory(cxxx, (char*)_outputBuffer, have);
	        _buffer.write((char*)_outputBuffer, have);
		bytes -= have;
        }
     }
  else {
     // unzip the data:
     unzipEventData(cxxx, lengthZip);
     }
  
  redumpMemory(cxxx, (char*)&lengthZipVerify, 4); // Bwd-Jumper & consistency check
  
  return (lengthZip == lengthZipVerify);
}
//----------------------------------------------------------------------
bool iotl::iStreamer::previous(std::istream& cxxx) {
  
  int pos = cxxx.tellg();   
  int lengthZipVerify;

  // seek previous block length information:
  pos -= 4;
  if (pos >= 0) {cxxx.seekg(pos);} 
  else {return false;}

  // read block length information:
  redumpMemory(cxxx, (char*)&lengthZipVerify, 4); // Bwd-Jumper & consistency check
  
  if (lengthZipVerify < 0) {pcl::exception("iotl::iStreamer::previous()","Invalid block length information.");}
  
  pos -= lengthZipVerify; // block
  pos -= 4; // block length
  pos -= 1; // compression mode
  pos -= 4; // event marker
  
  if (pos >= 0) {cxxx.seekg(pos); return true;}
   
  pcl::exception("iotl::iStreamer::previous()","Invalid block length information.");
  return false;
  
}
//----------------------------------------------------------------------
int iotl::iStreamer::unzipEventData( std::istream &in, int nBytes) {
	
  int ret, length = 0;
  unsigned int have;
  z_stream strm;

  /* allocate inflate state */
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  strm.avail_in = 0;
  strm.next_in = Z_NULL;
  
  ret = inflateInit(&strm);
  if (ret != Z_OK)
    return 0;
		
  /* decompress until deflate stream ends or end of file */
  do
  {
    strm.avail_in = in.read( (char *)_inputBuffer, iotl__iStreamer__lengthUnzipBuffer < nBytes ? iotl__iStreamer__lengthUnzipBuffer : nBytes ).gcount();
    if( in.bad() )
    {
      (void)inflateEnd(&strm);
      return 0;
    }
		
    nBytes -= strm.avail_in;
    if( in.eof() )
      nBytes = 0;
		
    if (strm.avail_in == 0)
      break;
				
    strm.next_in = _inputBuffer;
	
    /* run inflate() on input until output buffer not full */
    do
    {
      strm.avail_out = iotl__iStreamer__lengthUnzipBuffer;
      strm.next_out = _outputBuffer;
	
      ret = inflate(&strm, Z_NO_FLUSH);
      assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
				
      switch (ret)
      {
      case Z_NEED_DICT:
        ret = Z_DATA_ERROR;     /* and fall through */
      case Z_DATA_ERROR:
      case Z_MEM_ERROR:
        (void)inflateEnd(&strm);
        return 0;
      default:
        break;
      }
	
      have = iotl__iStreamer__lengthUnzipBuffer - strm.avail_out;
      _buffer.write( (char *)_outputBuffer, have);
		    
      length += have;
    } while (strm.avail_out == 0);

	/* done when inflate() says it's done */
  } while ( nBytes > 0 );

  /* clean up and return */
  inflateEnd(&strm);
    
  return length;
}
//----------------------------------------------------------------------
//
//
//----------------------------------------------------------------------
namespace iotl {
//----------------------------------------------------------------------
template<> ptl::Id iotl::iStreamer::restoreData<ptl::Relations>(ptl::Relations& relations) {

  relations.clearContainer();

  int size; restoreData<int>(size);
  for (int i=0; i<size; i++) {
       ptl::MutableId persistentId = 0; 
       restoreId(persistentId);
       // provide orphan relation
       if (persistentId != 0) relations.set(persistentId, new ptl::WkPtr<iotl::Orphan>);
       }
  return 0;       
}
//----------------------------------------------------------------------
template<> ptl::Id iotl::iStreamer::restoreData<ptl::Objects>(ptl::Objects& objects) {

  //std::cerr << "iotl::iStreamer::restoreData<ptl::Objects>(ptl::Objects& objects): START" << std::endl;

  objects.clearContainer();

  int size;   
  std::string key;
  ptl::MutableId persistentId = 0;

  // reload objects (with orphan relations)  
  restoreData<int>(size);
  for (int i=0; i<size; i++) {
      ptl::ObjectBase* item = 0;
      persistentId = restoreAbstractObject(&item);
      objects.set(*item);
      objects._copyHistory.set(persistentId, item);
      }

  // redirect relations: loop in PTL style
  for (ptl::Objects::PtlIterator iter(objects); !iter.isDone(); iter.next()) {
       ptl::ObjectBase* pNew = iter.item();
       
       // mother relations
       for (ptl::Relations::PtlTypeIterator<ptl::WkPtr<iotl::Orphan> > iterRelations(pNew->_motherRelations); !iterRelations.isDone(); iterRelations.next()) {
		   ptl::ObjectBase* pNewRel = objects._copyHistory.find(iterRelations.key(), 0);
		   if (pNewRel) {pNew->linkMother(*pNewRel);}
		   else {std::cerr << "iotl::iStreamer::restoreData<ptl::Objects>(...): WARNING: some original objects seem to have invalid relations. (persistent id is " << iterRelations.key() << ")" << std::endl;}
       }

       // daughter relations
       // have been set automatically above

       // remove orphan relations
       // mother relations
       ptl::Relations::PtlTypeIterator<ptl::WkPtr<iotl::Orphan> > iterM(pNew->_motherRelations); 
       while (!iterM.isDone()) {
		   ptl::Id id = iterM.key();
		   iterM.next();
		   pNew->_motherRelations.remove(id);
           }
       
       // daughter relations
       ptl::Relations::PtlTypeIterator<ptl::WkPtr<iotl::Orphan> > iterD(pNew->_daughterRelations); 
       while (!iterD.isDone()) {
		   ptl::Id id = iterD.key();
		   iterD.next();
		   pNew->_daughterRelations.remove(id);
           }
      }

  // reload & redirect index
  restoreData<int>(size);
  for (int i=0; i<size; i++) {
      restoreData<std::string>(key);
      restoreId(persistentId);
            
      ptl::ObjectBase* pNew = objects._copyHistory.find(persistentId, 0);
      
      if (pNew) {objects._index.set(key, pNew);}
	  else {
	  	std::cerr << "iotl::iStreamer::restoreData<ptl::Objects>(...): WARNING: some original indices could not be reconstructed." << std::endl;
	  	std::cerr << key << ": " << persistentId << std::endl;
	  	}
      }

  //std::cerr << "iotl::iStreamer::restoreData<ptl::Objects>(ptl::Objects& objects): END" << std::endl;
      
  return 0;       
}
//----------------------------------------------------------------------
template<> ptl::Id iotl::iStreamer::restoreData<ptl::Layout>(ptl::Layout& properties) {
  restoreData<int>(properties._type);
  restoreData<int>(properties._style);
  restoreData<int>(properties._color);
  restoreData<double>(properties._a);
  restoreData<double>(properties._b);
  restoreData<double>(properties._c);
  restoreData<double>(properties._d);
  return 0;
}
//----------------------------------------------------------------------
template<> ptl::Id iotl::iStreamer::restoreData<ptl::ObjectBase>(ptl::ObjectBase& object) {

  ptl::MutableId persistentId;
  bool hasLayout;
  
  restoreId(persistentId);
  
  restoreData<ptl::Relations>(object._motherRelations);
  restoreData<ptl::Relations>(object._daughterRelations);
  
  restoreData<bool>(hasLayout);
  if (hasLayout) {restoreData<ptl::Layout>( object.layout() );}
  else           {delete object._ptrLayout; object._ptrLayout = 0;}
  
  return persistentId;
}
//----------------------------------------------------------------------
ptl::Id iotl::iStreamer::restoreAbstractObject(ptl::ObjectBase** ppobj) {

  std::string objectTypeId; 
  std::string dataTypeId;
   
  pcl::BasicIoStreamer::restoreBasicTypeCStr(_buffer, objectTypeId);   
  pcl::BasicIoStreamer::restoreBasicTypeCStr(_buffer, dataTypeId); 
  
  return iotl::TypeManager::instance().restoreObject(*this, ppobj, objectTypeId, dataTypeId);
}
//----------------------------------------------------------------------
} // namespace iotl
//----------------------------------------------------------------------
void iotl::TypeManager::registerAgent(iotl::TypeAgentBase* agent) {
  
  // for debugging purposes:
//  std::cout << "iotl::TypeManager::registerAgent(): type registered: [" 
//            << agent->getObjectTypeId() << "," << agent->getDataTypeId() 
//            << "]" << std::endl;
//            //<< "] = '" << agent->getCppTypeId() << "'" << std::endl;

  // add entries to maps:
  _agentsByIotlTypeId.set(iotl::TypeIdKey(agent->getObjectTypeId(), agent->getDataTypeId()), agent);
  _agentsByCppTypeId.set(agent->getCppTypeId(), agent);
}
//----------------------------------------------------------------------
ptl::Id iotl::TypeManager::restoreObject(iStreamer& input, ptl::ObjectBase** ppobj, const std::string& objectTypeId, const std::string& dataTypeId) const {

  TypeAgentBase* agent = _agentsByIotlTypeId.find(iotl::TypeIdKey(objectTypeId, dataTypeId), 0);
  
  if (agent) {return agent->restoreObject(input, ppobj);}
  
  std::string message; 
  message += "undeclared IOTL type: [";
  message += objectTypeId;
  message += ","; 
  message += dataTypeId; 
  message += "]"; 
  pcl::exception("iotl::TypeManager::restoreObject()", message); 

  return 0;
 }
//----------------------------------------------------------------------
ptl::Id iotl::TypeManager::restoreObject(iotl::iStreamer& input,  ptl::ObjectBase& obj,  const std::string& cppTypeId) const {

  TypeAgentBase* agent = _agentsByCppTypeId.find(cppTypeId, 0);
  
  if (agent) {return agent->restoreObject(input, obj);}
  
  std::string message; 
  message += "undeclared C++ type: '";
  message += cppTypeId;
  message += "'";
  pcl::exception("iotl::TypeManager::restoreObject()", message); 

  return 0;
 }
//----------------------------------------------------------------------
void iotl::TypeManager::storeObject(oStreamer& output, const ptl::ObjectBase& obj, const std::string& cppTypeId) const {

  TypeAgentBase* agent = _agentsByCppTypeId.find(cppTypeId, 0);
  
  if (agent) {agent->storeObject(output, obj); return;}
  
  std::string message; 
  message += "undeclared C++ type: '";
  message += cppTypeId;
  message += "'";
  pcl::exception("iotl::TypeManager::storeObject()", message); 
 }
//----------------------------------------------------------------------
namespace iotl {TypeManager* TypeManager::_instance = 0;}
//----------------------------------------------------------------------
iotl::TypeManager& iotl::TypeManager::instance() {
  if (!_instance) {_instance = new TypeManager;}
  return *_instance;	
}
//----------------------------------------------------------------------
//
//
//
template class iotl::oDiskFileVx<iotl::oStreamer>;
template class iotl::iDiskFileVx<iotl::iStreamer>;

template class iotl::TypeAgent<ptl::Object<int> >;
template class iotl::TypeAgent<ptl::CowObject<int> >;

#ifndef MERGED_PXL
	#include "iotl.pol.cci"
#endif
#ifndef MERGED_PXL
	#include "iotl.hh"
	#include "pol.hh"
#endif
//
//
//----------------------------------------------------------------------
iotl__declareObjectTypeExplicit(pol::BasicObject, "\3Bo", pol__BasicObject)
iotl__declareObjectTypeExplicit(pol::BasicObjectManager, "\3Bom", pol__BasicObjectManager)
iotl__declareObjectTypeExplicit(pol::Particle, "\3Pa", pol__Particle)
iotl__declareObjectTypeExplicit(pol::Vertex, "\3Vx", pol__Vertex)
// iotl__declareObjectExplicit(pol::Collision, "\3Co", pol__Collision)
iotl__declareObjectTypeExplicit(pol::EventView, "\3Ev", pol__EventView)
iotl__declareObjectTypeExplicit(pol::AnalysisProcess, "\3Ap", pol__AnalysisProcess)
iotl__declareObjectTypeExplicit(pol::AnalysisFork, "\3Af", pol__AnalysisFork)
//----------------------------------------------------------------------
iotl__declareDataTypeExplicit(pol::BasicObjectData, "\3BoD", pol__BasicObjectData, data,
{
  storeData<bool>(data._locked);
  storeData<int>(data._monteCarloMode);
  storeData<std::string>(data._name);
  storeData<int>(data._status);
  storeData<int>(data._workflag);
  storeData(data._userRecords);
  // storeData<pol::CppPointers>(_cppPointers); // we don't store the c++ pointers!
},{
  restoreData<bool>(data._locked);
  restoreData<int>(data._monteCarloMode);
  restoreData<std::string>(data._name);
  restoreData<int>(data._status);
  restoreData<int>(data._workflag);
  restoreData(data._userRecords);
  // restoreData<pol::CppPointers>(_cppPointers); // we don't restore the c++ pointers!
})
//----------------------------------------------------------------------
iotl__declareDataTypeExplicit(pol::BasicObjectManagerData, "\3BomD", pol__BasicObjectManagerData, data,
{
  storeData<pol::Objects>(data._objects);
  storeData<pol::BasicObjectData>(data);
},{
  restoreData<pol::Objects>(data._objects);
  restoreData<pol::BasicObjectData>(data);
})
//----------------------------------------------------------------------
iotl__declareDataTypeExplicit(pol::Basic4VectorData, "\3B4vD", pol__Basic4VectorData, data,
{
  storeData<double>(data._x);
  storeData<double>(data._y);
  storeData<double>(data._z);
  storeData<double>(data._t);
},{
  restoreData<double>(data._x);
  restoreData<double>(data._y);
  restoreData<double>(data._z);
  restoreData<double>(data._t);
})
//----------------------------------------------------------------------
iotl__declareDataTypeExplicit(pol::Basic3VectorData, "\3B3vD", pol__Basic3VectorData, data,
{
  storeData<double>(data._x);
  storeData<double>(data._y);
  storeData<double>(data._z);
},{
  restoreData<double>(data._x);
  restoreData<double>(data._y);
  restoreData<double>(data._z);
})
//----------------------------------------------------------------------
iotl__declareDataTypeExplicit(pol::ParticleData, "\3PaD", pol__ParticleData, data,
{
  storeData<double>(data._charge);
  storeData<int>(data._particleId);
  storeData(data.vector());
  storeData<pol::BasicObjectData>(data);
},{
  restoreData<double>(data._charge);
  restoreData<int>(data._particleId);
  restoreData(data.vector(pol::set));
  restoreData<pol::BasicObjectData>(data);
})
//----------------------------------------------------------------------
iotl__declareDataTypeExplicit(pol::VertexData, "\3VxD", pol__VertexData, data,
{
  storeData(data.vector());
  storeData<pol::BasicObjectData>(data);
},{
  restoreData(data.vector(pol::set));
  restoreData<pol::BasicObjectData>(data);
})
//----------------------------------------------------------------------
// iotl__declareDataExplicit(pol::CollisionData, "\3CoD", pol__CollisionData, data,
// {
//   storeData<pol::BasicObjectData>(data);
// },{
//   restoreData<pol::BasicObjectData>(data);
// })
//----------------------------------------------------------------------
iotl__declareDataTypeExplicit(pol::EventViewData, "\3EvD", pol__EventViewData, data,
{
  storeData<pol::BasicObjectManagerData>(data);
},{
  restoreData<pol::BasicObjectManagerData>(data);
})
//----------------------------------------------------------------------
iotl__declareDataTypeExplicit(pol::AnalysisProcessData, "\3ApD", pol__AnalysisProcessData, data,
{
  storeData<pol::BasicObjectManagerData>(data);
},{
  restoreData<pol::BasicObjectManagerData>(data);
})
//----------------------------------------------------------------------
iotl__declareDataTypeExplicit(pol::AnalysisForkData, "\3AfD", pol__AnalysisForkData, data, 
{
  storeData<pol::BasicObjectManagerData>(data);
},{
  restoreData<pol::BasicObjectManagerData>(data);
})
//----------------------------------------------------------------------
#ifndef MERGED_PXL
	#include "ptl.hh"
#endif
//
//
//----------------------------------------------------------------------
// MEMBER FUNCTIONS
//----------------------------------------------------------------------
ptl::Objects::Objects() : ptl::Vector<ptl::ObjectBase*>(), _copyHistory(), _index() {;}
ptl::Objects::Objects(const ptl::Objects& original) : ptl::Vector<ptl::ObjectBase*>(), _copyHistory(), _index() {

  // copy objects: loop in STL style
  for (ptl::Objects::StlConstIterator iter = original.begin(); 
                                      iter != original.end(); 
				      iter++) {
       ptl::ObjectBase* pOld = *iter;
       ptl::ObjectBase* pNew = pOld->clone();

       set(*pNew);
       _copyHistory.set(pOld->id(), pNew);
      }
  
  // redirect relations: loop in PTL style
  for (ptl::Objects::PtlIterator iter(original); !iter.isDone(); iter.next()) {
       ptl::ObjectBase* pOld = iter.item();
       ptl::ObjectBase* pNew = _copyHistory.find(pOld->id(), 0);
       
       // mother relations
       for (ptl::Relations::PtlIterator iter(pOld->getMotherRelations()); !iter.isDone(); iter.next()) {
	   ptl::ObjectBase* pOldRel = iter.item()->pointer();
	   ptl::ObjectBase* pNewRel = _copyHistory.find(pOldRel->id(), 0);
	   
	   if (pOldRel) {	   
	      if (pNewRel) {pNew->linkMother(*pNewRel);}
	      else {std::cerr << "ptl::Objects::Objects(...): WARNING: some original objects had relations to objects of other holders." << std::endl;}
	      }
	   else {
	      std::cerr << "ptl::Objects::Objects(...): WARNING: some originally related objects no longer exist" << std::endl;
	      }   
       }

       // daughter relations
       // have been set automatically above
      }
  
  // redirect index:
  for (ptl::Index::PtlIterator iter(original._index); !iter.isDone(); iter.next()) {
       ptl::ObjectBase* pOld = iter.item();
       ptl::ObjectBase* pNew = _copyHistory.find(pOld->id(), 0);
       
       if (pNew) {_index.set(iter.key(),pNew);}
       else {std::cerr << "ptl::Objects::Objects(...): WARNING: some original indices pointed to objects of other holders." << std::endl;}
      }
}
//----------------------------------------------------------------------
void ptl::Objects::clearContainer() {
  for (StlConstIterator iter = _container.begin(); iter != _container.end(); iter++) {
      delete (*iter);
      }
  _container.clear();
  _copyHistory.clearContainer();
  _index.clearContainer();
}
//----------------------------------------------------------------------
void ptl::Objects::set(ptl::ObjectBase& item) {item._refObjects = this; _container.push_back(&item);}
//----------------------------------------------------------------------
void ptl::Objects::remove(ptl::ObjectBase& item) {

  // search & remove possible indices (multiple occurrences possible!)
  for (ptl::Index::StlConstIterator iter = _index.begin(); iter != _index.end(); iter++) {
       if (&item == iter->second) {_index.remove(iter->first);}
      }
  
  // search & remove possible copy history (multiple occurrences *not* possible!)
  for (ptl::CopyHistory::StlConstIterator iter = _copyHistory.begin(); iter != _copyHistory.end(); iter++) {
       if (&item == iter->second) {_copyHistory.remove(iter->first);break;}
      }

  // remove all relations:
  item.unlinkMothers();
  item.unlinkDaughters();

  // remove object
  for (StlIterator iter = _container.begin(); iter != _container.end(); iter++) {
  	   if (&item == (*iter)) {delete *iter; _container.erase(iter); break;}
  	  }
  	  
}
//----------------------------------------------------------------------
bool ptl::Objects::has(const ptl::ObjectBase& item) const {return (item._refObjects == this);}
//----------------------------------------------------------------------
ptl::Relations::Relations() : ptl::Map<ptl::Id, ptl::WkPtrBase*>() {;}
ptl::Relations::Relations(const ptl::Relations& original) : ptl::Map<ptl::Id, ptl::WkPtrBase*>() {
  for (StlConstIterator iter = original._container.begin(); iter != original._container.end(); iter++) {
      set(iter->first, iter->second->clone());
      }
}

void ptl::Relations::set(ptl::Id id, ptl::WkPtrBase* wptr) {
	ptl::Relations::remove(id); 
	_container.insert(StlPair(id, wptr));
}

void ptl::Relations::clearContainer() {
  for (StlConstIterator iter = _container.begin(); iter != _container.end(); iter++) {
      delete iter->second;
      }
  _container.clear();
}

void ptl::Relations::remove(ptl::ObjectBase& object) {delete find(object.id(), 0); _container.erase(object.id());}
void ptl::Relations::remove(ptl::Id id)              {delete find(id, 0); _container.erase(id);}

bool ptl::Relations::has(const ptl::ObjectBase& object) const {return (0 != find(object.id(), 0) );}
bool ptl::Relations::has(ptl::Id id) const                {return (0 != find(id, 0) );}
//----------------------------------------------------------------------
void ptl::WkPtrBase::notifyDeleted() {

  _objectRef = 0; 
  if (_notifyChainOut) _notifyChainOut->notifyDeleted(); 
  _notifyChainIn = 0; 
  _notifyChainOut = 0; 

}
//----------------------------------------------------------------------
void ptl::WkPtrBase::connect(ptl::ObjectBase* pointer) {

  // disconnect:
  if (_objectRef) {
     if (_objectRef->_refWkPtrSpec == this) {_objectRef->_refWkPtrSpec = _notifyChainOut;}
     if (_notifyChainIn && _notifyChainOut) {_notifyChainIn->_notifyChainOut = _notifyChainOut; _notifyChainOut->_notifyChainIn = _notifyChainIn;}
     else {
	if (_notifyChainIn)  _notifyChainIn->_notifyChainOut = 0; 
	if (_notifyChainOut) _notifyChainOut->_notifyChainIn = 0;
	}
     }
  _notifyChainOut = 0; 
  _notifyChainIn = 0;

  // connect:
  if (pointer) {
     _notifyChainIn = 0;
     _notifyChainOut = pointer->_refWkPtrSpec;
     if (_notifyChainOut) _notifyChainOut->_notifyChainIn = this;
     pointer->_refWkPtrSpec = this;
     }
  _objectRef = pointer;
  
}
//----------------------------------------------------------------------
void ptl::ObjectBase::linkMother(ptl::ObjectBase& target)  {

  ptl::WkPtrBase* pTarget = target.createSelfWkPtr();
  ptl::WkPtrBase* pThis   = this->createSelfWkPtr();
  
  if (target._refObjects != this->_refObjects) {std::cerr << "ptl::ObjectBase::linkDaughter(...): WARNING: mother and daughter have not the same object holder!" << std::endl;}

  this->_motherRelations.set(target.id(), pTarget);
  target._daughterRelations.set(this->id(), pThis);

}
//----------------------------------------------------------------------
void ptl::ObjectBase::linkDaughter(ptl::ObjectBase& target) {

  ptl::WkPtrBase* pTarget = target.createSelfWkPtr();
  ptl::WkPtrBase* pThis   = this->createSelfWkPtr();
  
  if (target._refObjects != this->_refObjects) {std::cerr << "ptl::ObjectBase::linkDaughter(...): WARNING: mother and daughter have not the same object holder!" << std::endl;}
  
  this->_daughterRelations.set(target.id(), pTarget);
  target._motherRelations.set(this->id(), pThis);

}
//----------------------------------------------------------------------
void ptl::ObjectBase::unlinkMother(ptl::ObjectBase& target)  {

  this->_motherRelations.remove(target.id());
  target._daughterRelations.remove(this->id());

}
//----------------------------------------------------------------------
void ptl::ObjectBase::unlinkDaughter(ptl::ObjectBase& target) {

  this->_daughterRelations.remove(target.id());
  target._motherRelations.remove(this->id());

}
//----------------------------------------------------------------------
void ptl::ObjectBase::unlinkMothers() {

  for (ptl::Relations::PtlIterator iter(_motherRelations); !iter.isDone(); iter.next()) {
       if (iter.item()->_objectRef) {
       	  iter.item()->_objectRef->_daughterRelations.remove(this->id());
       	  }
       }
  _motherRelations.clearContainer();
  
}  
//----------------------------------------------------------------------
void ptl::ObjectBase::unlinkDaughters() {

  for (ptl::Relations::PtlIterator iter(_daughterRelations); !iter.isDone(); iter.next()) {
       if (iter.item()->_objectRef) {
       	  iter.item()->_objectRef->_motherRelations.remove(this->id());
       	  }
       }
  _daughterRelations.clearContainer();

}
//----------------------------------------------------------------------
std::ostream& ptl::ObjectBase::printDecayTree(int level, std::ostream& os, int pan) const {

  int daughters = 0;
  
  print(level, os, pan);

  for (ptl::Relations::PtlIterator iter(_daughterRelations); !iter.isDone(); iter.next()) {
       if (iter.item()->_objectRef) {
          iter.item()->_objectRef->printDecayTree(level, os, pan+1);
	  daughters++;
	  }
      }

   if (daughters && pan > 1) {
      for (int p=0;p<pan;p++) {os << "|  ";}
      os << "*" << std::endl;
      }

  return os;

}
//----------------------------------------------------------------------
std::ostream& ptl::ObjectBase::print(int level, std::ostream& os, int pan) const {return printPan1st(os, pan) << "ptl::ObjectBase with id: " << id() << std::endl;}
//----------------------------------------------------------------------
std::ostream& ptl::ObjectBase::printPan1st(std::ostream& os, int pan) const {for (int p=0;p<pan-2;p++) {os << "|  ";} if (pan-1 > 0) os << "+--"; if (pan) os << "{ "; return os;}
std::ostream& ptl::ObjectBase::printPan(std::ostream& os, int pan) const {for (int p=0;p<pan-1;p++) {os << "|  ";}; if (pan) os << "| "; return os;}
//----------------------------------------------------------------------
namespace ptl
{
template <> std::ostream& ptl::Object<int>::print(int level, std::ostream& os, int pan) const {return printPan1st(os, pan) << "ptl::Object<int> with id: " << id() << " is " << get() << std::endl;}
template <> std::ostream& ptl::Object<unsigned int>::print(int level, std::ostream& os, int pan) const {return printPan1st(os, pan) << "ptl::Object<unsigned int> with id: " << id() << " is " << get() << std::endl;}
template <> std::ostream& ptl::Object<bool>::print(int level, std::ostream& os, int pan) const {return printPan1st(os, pan) << "ptl::Object<bool> with id: " << id() << " is " << get() << std::endl;}
template <> std::ostream& ptl::Object<double>::print(int level, std::ostream& os, int pan) const {return printPan1st(os, pan) << "ptl::Object<double> with id: " << id() << " is " << get() << std::endl;}
template <> std::ostream& ptl::Object<float>::print(int level, std::ostream& os, int pan) const {return printPan1st(os, pan) << "ptl::Object<float> with id: " << id() << " is " << get() << std::endl;}
//----------------------------------------------------------------------
template <> std::ostream& ptl::CowObject<int>::print(int level, std::ostream& os, int pan) const {return printPan1st(os, pan) << "ptl::CowObject<int> with id: " << id() << " is " << get() << std::endl;}
template <> std::ostream& ptl::CowObject<unsigned int>::print(int level, std::ostream& os, int pan) const {return printPan1st(os, pan) << "ptl::CowObject<unsigned int> with id: " << id() << " is " << get() << std::endl;}
template <> std::ostream& ptl::CowObject<bool>::print(int level, std::ostream& os, int pan) const {return printPan1st(os, pan) << "ptl::CowObject<bool> with id: " << id() << " is " << get() << std::endl;}
template <> std::ostream& ptl::CowObject<double>::print(int level, std::ostream& os, int pan) const {return printPan1st(os, pan) << "ptl::CowObject<double> with id: " << id() << " is " << get() << std::endl;}
template <> std::ostream& ptl::CowObject<float>::print(int level, std::ostream& os, int pan) const {return printPan1st(os, pan) << "ptl::CowObject<float> with id: " << id() << " is " << get() << std::endl;}
}
//----------------------------------------------------------------------
std::ostream& operator << (std::ostream& cxxx, const ptl::ObjectBase& obj) {return obj.print(0, cxxx, 0);}
//----------------------------------------------------------------------
// EXPLICIT INSTANCIATION
template class ptl::Vector<int>;
template class ptl::Map<int,int>;
template class ptl::Objects::PtlTypeIterator<ptl::Object<int> >;
template class ptl::Relations::PtlTypeIterator<ptl::Object<int> >;

template class ptl::Ptr<int>;
template class ptl::WkPtrSpec<int, ptl::Object<int> >;
template class ptl::WkPtr<int>;
template class ptl::CowWkPtr<int>;
template class ptl::SpyWkPtr<int>;

template class ptl::Object<int>;
template class ptl::CowObject<int>;
template class ptl::SpyObject<int>;
//----------------------------------------------------------------------
#ifndef MERGED_PXL
	#include "pol.hh"
#endif
//
//
//----------------------------------------------------------------------
namespace pol {
//----------------------------------------------------------------------
bool const operator==(const pol::Basic4VectorData& obj1, const pol::Basic4VectorData& obj2)
{return (obj1.getX() == obj2.getX() && obj1.getY() == obj2.getY() && obj1.getZ() == obj2.getZ() && obj1.getE() == obj2.getE());}

bool const operator!=(const pol::Basic4VectorData& obj1, const pol::Basic4VectorData& obj2)
{return (obj1.getX() != obj2.getX() || obj1.getY() != obj2.getY() || obj1.getZ() != obj2.getZ() || obj1.getE() != obj2.getE());}
//----------------------------------------------------------------------
bool const operator==(const pol::Basic3VectorData& obj1, const pol::Basic3VectorData& obj2) 
{return (obj1.getX() == obj2.getX() && obj1.getY() == obj2.getY() && obj1.getZ() == obj2.getZ());}

bool const operator!=(const pol::Basic3VectorData& obj1, const pol::Basic3VectorData& obj2) 
{return (obj1.getX() != obj2.getX() || obj1.getY() != obj2.getY() || obj1.getZ() != obj2.getZ());}
//----------------------------------------------------------------------
bool const operator==(const pol::VertexData& obj1, const pol::VertexData& obj2)
{return (obj1.vector() == obj2.vector());}

bool const operator!=(const pol::VertexData& obj1, const pol::VertexData& obj2) 
{return (obj1.vector() != obj2.vector());}
//----------------------------------------------------------------------
bool const operator==(const pol::ParticleData& obj1, const pol::ParticleData& obj2) 
{return (obj1.vector() == obj2.vector() && obj1.getCharge() != obj2.getCharge());}

bool const operator!=(const pol::ParticleData& obj1, const pol::ParticleData& obj2) 
{return (obj1.vector() != obj2.vector() || obj1.getCharge() != obj2.getCharge());}
//----------------------------------------------------------------------
} // namespace pol
//----------------------------------------------------------------------
ptl::ObjectBase* pol::AnalysisProcess::clone() const {return new AnalysisProcess(*this);}
ptl::WkPtrBase* pol::AnalysisProcess::createSelfWkPtr() {return new AnalysisProcessWkPtr(*this);}
void pol::AnalysisProcess::storeYourSelf(iotl::oStreamer& output) const {output.storeObject(*this);}
//----------------------------------------------------------------------
void pol::AnalysisFork::beginJob(const pol::Objects* input) {

  for (pol::Objects::PtlTypeIterator<pol::AnalysisFork> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().beginJob(input);
      }
      
  for (pol::Objects::PtlTypeIterator<pol::AnalysisProcess> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().beginJob(input);
      }
      
}
//----------------------------------------------------------------------
void pol::AnalysisFork::buildTemplate(int mode) {

  for (pol::Objects::PtlTypeIterator<pol::AnalysisFork> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().buildTemplate(mode);
      }
      
  for (pol::Objects::PtlTypeIterator<pol::AnalysisProcess> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().buildTemplate(mode);
      }
      
}
//----------------------------------------------------------------------
void pol::AnalysisFork::beginRun(const pol::Objects* input) {

  for (pol::Objects::PtlTypeIterator<pol::AnalysisFork> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().beginRun(input);
      }
      
  for (pol::Objects::PtlTypeIterator<pol::AnalysisProcess> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().beginRun(input);
      }
      
}
//----------------------------------------------------------------------
void pol::AnalysisFork::analyseEvent(const pol::Objects* input) {

  for (pol::Objects::PtlTypeIterator<pol::AnalysisFork> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().analyseEvent(input);
      }
      
  for (pol::Objects::PtlTypeIterator<pol::AnalysisProcess> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().analyseEvent(input);
      }
      
}
//----------------------------------------------------------------------
void pol::AnalysisFork::finishEvent(const pol::Objects* input) {

  for (pol::Objects::PtlTypeIterator<pol::AnalysisFork> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().finishEvent(input);
      }
      
  for (pol::Objects::PtlTypeIterator<pol::AnalysisProcess> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().finishEvent(input);
      }
      
}
//----------------------------------------------------------------------
void pol::AnalysisFork::endRun(const pol::Objects* input) {

  for (pol::Objects::PtlTypeIterator<pol::AnalysisFork> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().endRun(input);
      }
      
  for (pol::Objects::PtlTypeIterator<pol::AnalysisProcess> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().endRun(input);
      }
      
}
//----------------------------------------------------------------------
void pol::AnalysisFork::endJob(const pol::Objects* input) {

  for (pol::Objects::PtlTypeIterator<pol::AnalysisFork> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().endJob(input);
      }
      
  for (pol::Objects::PtlTypeIterator<pol::AnalysisProcess> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().endJob(input);
      }
      
}
//----------------------------------------------------------------------
ptl::ObjectBase* pol::AnalysisFork::clone() const {return new AnalysisFork(*this);}
ptl::WkPtrBase* pol::AnalysisFork::createSelfWkPtr() {return new AnalysisForkWkPtr(*this);}
void pol::AnalysisFork::storeYourSelf(iotl::oStreamer& output) const {output.storeObject(*this);}
//----------------------------------------------------------------------
//
//
//
//----------------------------------------------------------------------
namespace ptl {
//----------------------------------------------------------------------
template <>  std::ostream& ptl::CowObject<pol::ParticleData>::print(int level, std::ostream& os, int pan) const {

//  printPan1st(os, pan) << "ptl::CowObject<pol::ParticleData> with id: " << id() << " (data socket currently located at " << _dataSocket << ")" << std::endl;
//  printPan(os, pan)    << "     name: " << get().getName() << std::endl;
  printPan1st(os, pan) << "Particle: '" << get().getName() << "', p = (" 
                             << get().vector().getPt() << ", " 
                             << get().vector().getPz() << ") m = " 
                             << get().vector().getMass() << std::endl;
  return os;

}
//----------------------------------------------------------------------
template <>  std::ostream& ptl::CowObject<pol::VertexData>::print(int level, std::ostream& os, int pan) const {

//  printPan1st(os, pan) << "ptl::CowObject<pol::VertexData> with id: " << id() << " (data socket currently located at " << _dataSocket << ")" << std::endl;
//  printPan(os, pan)    << "     name: " << get().getName() << std::endl;
  printPan1st(os, pan) << "Vertex: '" << get().getName() << "', x = (" 
                             << get().vector().getX() << ", " 
                             << get().vector().getY() << ", " 
                             << get().vector().getZ() << ")" << std::endl;
  return os;

}
//----------------------------------------------------------------------
template <>  std::ostream& ptl::CowObject<pol::CollisionData>::print(int level, std::ostream& os, int pan) const {

//  printPan1st(os, pan) << "ptl::CowObject<pol::CollisionData> with id: " << id() << " (data socket currently located at " << _dataSocket << ")" << std::endl;
//  printPan(os, pan)    << "     name: " << get().getName() << std::endl;
  printPan1st(os, pan) << "Collision: " << get().getName() << std::endl;
  return os;

}
//----------------------------------------------------------------------
template <>  std::ostream& ptl::Object<pol::EventViewData>::print(int level, std::ostream& os, int pan) const {

//  printPan1st(os, pan) << "ptl::CowObject<pol::EventViewData> with id: " << id() << std::endl;
//  printPan(os, pan)    << "     name: " << get().getName() << std::endl;
  printPan1st(os, pan) << "EventView: " << get().getName() << std::endl;
  for (pol::Objects::PtlIterator iter(get().getObjects()); !iter.isDone(); iter.next()) {
      if (iter.object().getMotherRelations().getSize() == 0) {iter.object().printDecayTree(level, os, pan);}
      }
  return os;

}
//----------------------------------------------------------------------
template <>  std::ostream& ptl::Object<pol::AnalysisProcessData>::print(int level, std::ostream& os, int pan) const {

//  printPan1st(os, pan) << "ptl::CowObject<pol::AnalysisProcessData> with id: " << id() << std::endl;
//  printPan(os, pan)    << "     name: " << get().getName() << std::endl;
  printPan1st(os, pan) << "AnalysisProcess: " << get().getName() << std::endl;
  for (pol::Objects::PtlIterator iter(get().getObjects()); !iter.isDone(); iter.next()) {
      if (iter.object().getMotherRelations().getSize() == 0) {iter.object().printDecayTree(level, os, pan);}
      }
  return os;

}
//----------------------------------------------------------------------
template <>  std::ostream& ptl::Object<pol::AnalysisForkData>::print(int level, std::ostream& os, int pan) const {

//  printPan1st(os, pan) << "ptl::CowObject<pol::AnalysisForkData> with id: " << id() << std::endl;
//  printPan(os, pan)    << "     name: " << get().getName() << std::endl;
  printPan1st(os, pan) << "AnalysisFork: " << get().getName() << std::endl;
  for (pol::Objects::PtlIterator iter(get().getObjects()); !iter.isDone(); iter.next()) {
      if (iter.object().getMotherRelations().getSize() == 0) {iter.object().printDecayTree(level, os, pan);}
      }
  return os;

}
//----------------------------------------------------------------------
}
//----------------------------------------------------------------------
#ifndef MERGED_PXL
	#include "ePax.hh"
#endif	
//
//
//----------------------------------------------------------------------
// ePax
//----------------------------------------------------------------------
iotl__declareObjectTypeExplicit(ePax::ePaxParticle,  "\4ePa", ePax__ePaxParticle)
iotl__declareObjectTypeExplicit(ePax::ePaxVertex,    "\4eVe", ePax__ePaxVertex)
iotl__declareObjectTypeExplicit(ePax::ePaxCollision, "\4eCo", ePax__ePaxCollision)
iotl__declareObjectTypeExplicit(ePax::ePaxEventView, "\4eEv", ePax__ePaxEventView)
iotl__declareObjectTypeExplicit(ePax::ePaxAnalysisProcess, "\4eAp", ePax__ePaxAnalysisProcess)
iotl__declareObjectTypeExplicit(ePax::ePaxAnalysisFork,    "\4eAf", ePax__ePaxAnalysisFork)
//----------------------------------------------------------------------
ptl::ObjectBase* ePax::ePaxParticle::clone() const {return new ePaxParticle(*this);}
std::ostream& ePax::ePaxParticle::print(int level, std::ostream& os, int pan) const {return ptl::CowObject<pol::ParticleData>::print(level, os, pan);}
ptl::WkPtrBase* ePax::ePaxParticle::createSelfWkPtr() {return new ePaxParticleWkPtr(*this);}
void ePax::ePaxParticle::storeYourSelf(iotl::oStreamer& output) const {output.storeObject(*this);}
//----------------------------------------------------------------------
ptl::ObjectBase* ePax::ePaxVertex::clone() const {return new ePaxVertex(*this);}
std::ostream& ePax::ePaxVertex::print(int level, std::ostream& os, int pan) const {return ptl::CowObject<pol::VertexData>::print(level, os, pan);}
ptl::WkPtrBase* ePax::ePaxVertex::createSelfWkPtr() {return new ePaxVertexWkPtr(*this);}
void ePax::ePaxVertex::storeYourSelf(iotl::oStreamer& output) const {output.storeObject(*this);}
//----------------------------------------------------------------------
ptl::ObjectBase* ePax::ePaxCollision::clone() const {return new ePaxCollision(*this);}
std::ostream& ePax::ePaxCollision::print(int level, std::ostream& os, int pan) const {return ptl::CowObject<pol::CollisionData>::print(level, os, pan);}
ptl::WkPtrBase* ePax::ePaxCollision::createSelfWkPtr() {return new ePaxCollisionWkPtr(*this);}
void ePax::ePaxCollision::storeYourSelf(iotl::oStreamer& output) const {output.storeObject(*this);}
//----------------------------------------------------------------------
ptl::ObjectBase* ePax::ePaxEventView::clone() const {return new ePaxEventView(*this);}
std::ostream& ePax::ePaxEventView::print(int level, std::ostream& os, int pan) const {return ptl::Object<pol::EventViewData>::print(level, os, pan);}
ptl::WkPtrBase* ePax::ePaxEventView::createSelfWkPtr() {return new ePaxEventViewWkPtr(*this);}
void ePax::ePaxEventView::storeYourSelf(iotl::oStreamer& output) const {output.storeObject(*this);}
//----------------------------------------------------------------------
ptl::ObjectBase* ePax::ePaxAnalysisProcess::clone() const {return new ePaxAnalysisProcess(*this);}
std::ostream& ePax::ePaxAnalysisProcess::print(int level, std::ostream& os, int pan) const {return ptl::Object<pol::AnalysisProcessData>::print(level, os, pan);}
ptl::WkPtrBase* ePax::ePaxAnalysisProcess::createSelfWkPtr() {return new ePaxAnalysisProcessWkPtr(*this);}
void ePax::ePaxAnalysisProcess::storeYourSelf(iotl::oStreamer& output) const {output.storeObject(*this);}
//----------------------------------------------------------------------
ptl::ObjectBase* ePax::ePaxAnalysisFork::clone() const {return new ePaxAnalysisFork(*this);}
std::ostream& ePax::ePaxAnalysisFork::print(int level, std::ostream& os, int pan) const {return ptl::Object<pol::AnalysisForkData>::print(level, os, pan);}
ptl::WkPtrBase* ePax::ePaxAnalysisFork::createSelfWkPtr() {return new ePaxAnalysisForkWkPtr(*this);}
void ePax::ePaxAnalysisFork::storeYourSelf(iotl::oStreamer& output) const {output.storeObject(*this);}
//----------------------------------------------------------------------
