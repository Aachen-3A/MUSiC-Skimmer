#include "ePaxPxl/ePax/interface/ePax.h"
#ifndef MERGED_PXL
	#include "pcl.hh"
#endif

#include <iostream>
#include <sys/times.h>
#include <unistd.h>
//----------------------------------------------------------------------
namespace pxl {
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
  std::cerr << "pxl::exception(): throwing pxl::Exception." << std::endl; 
  throw pxl::Exception(rout, msg);
}
//----------------------------------------------------------------------
} // namespace pxl
//----------------------------------------------------------------------
#include <assert.h>

#ifndef MERGED_PXL
	#include "ptl.hh"
#endif
//
//
//----------------------------------------------------------------------
iotl__declareDataTypeExplicit(char,        "\1c",      __char, data, {pxl::BasicIoStreamer::storeBasicTypeChar(_buffer, data);}   , {pxl::BasicIoStreamer::restoreBasicTypeChar(_buffer, data);})
iotl__declareDataTypeExplicit(std::string, "\1s", std__string, data, {pxl::BasicIoStreamer::storeBasicTypeString(_buffer, data);} , {pxl::BasicIoStreamer::restoreBasicTypeString(_buffer, data);})
iotl__declareDataTypeExplicit(bool,        "\1b",      __bool, data, {pxl::BasicIoStreamer::storeBasicTypeBool(_buffer, data);}  , {pxl::BasicIoStreamer::restoreBasicTypeBool(_buffer, data);})
iotl__declareDataTypeExplicit(int,         "\1i",       __int, data, {pxl::BasicIoStreamer::storeBasicTypeInt(_buffer, data);}   , {pxl::BasicIoStreamer::restoreBasicTypeInt(_buffer, data);})
iotl__declareDataTypeExplicit(float,       "\1f",     __float, data, {pxl::BasicIoStreamer::storeBasicTypeFloat(_buffer, data);}  , {pxl::BasicIoStreamer::restoreBasicTypeFloat(_buffer, data);})
iotl__declareDataTypeExplicit(double,      "\1d",    __double, data, {pxl::BasicIoStreamer::storeBasicTypeDouble(_buffer, data);} , {pxl::BasicIoStreamer::restoreBasicTypeDouble(_buffer, data);})
//----------------------------------------------------------------------
void pxl::oStreamer::getEvent(std::ostream& cxxx, const std::string& info, char compressionMode) {

  // info block:
  const char* cInfo = info.c_str(); 
  int lengthInfo = info.length();

  // zip block:
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
	pxl::exception("pxl::oStreamer::getEvent()", "zlib: not enough memory");
       case Z_BUF_ERROR:
	 pxl::exception("pxl::oStreamer::getEvent()", "zlib: buffer too small");
       case Z_STREAM_ERROR:
	 pxl::exception("pxl::oStreamer::getEvent()", "zlib: level parameter invalid");
       default:
	 break;
       }
    cZip = cZipSpace;
    lengthZip = lengthZipSpace;
  }
  else {pxl::exception("pxl::oStreamer::getEvent()","Invalid compression mode.");}
  
  // write event
  dumpMemory(cxxx, iotl__eventMarker, 4);          // Event-Marker
  storeBasicTypeChar(cxxx, compressionMode);       // Compression-Mode
  
  storeBasicTypeInt(cxxx, lengthInfo);             // Fwd-Jumper
  cxxx.rdbuf()->sputn(cInfo, lengthInfo);          // info block
  cxxx.rdbuf()->pubsync();
  storeBasicTypeInt(cxxx, lengthInfo);             // Bwd-Jumper & consistency check
  
  storeBasicTypeInt(cxxx, lengthZip);              // Fwd-Jumper
  cxxx.rdbuf()->sputn(cZip, lengthZip);            // compressed zip block
  cxxx.rdbuf()->pubsync();
  storeBasicTypeInt(cxxx, lengthZip);               // Bwd-Jumper & consistency check
  
  if (cZipSpace) {delete [] cZipSpace;}
  
  _buffer.str("");
  _buffer.clear();

}
//
//
//----------------------------------------------------------------------
namespace pxl {
template<> void oStreamer::storeData<pxl::Relations>(const pxl::Relations& relations) {

  storeData<int>(relations.getSize());
  for (pxl::Relations::StlConstIterator iter = relations.begin(); 
                                      iter != relations.end(); 
				      iter++) {
       const pxl::WkPtrBase& ptr = (*iter->second);
       if (ptr.valid()) {storeId(ptr.pointer()->id());}
       else             {storeId(0);}
       }

}
//----------------------------------------------------------------------
template<> void pxl::oStreamer::storeData<pxl::Objects>(const pxl::Objects& objects) {

  // objects
  storeData<int>(objects.getSize());
  for (pxl::Objects::StlConstIterator iter = objects.begin(); iter != objects.end(); iter++) {
  	  storeAbstractObject(**iter);
  	  }

  // index
  storeData<int>(objects._index.getSize());
  for (pxl::Index::StlConstIterator iter = objects._index.begin(); iter != objects._index.end(); iter++) {
      storeData<std::string>(iter->first);
      storeId(iter->second->id());
      }

}
//----------------------------------------------------------------------
template<> void pxl::oStreamer::storeData<pxl::Layout>(const pxl::Layout& properties) {
  storeData<int>(properties._type);
  storeData<int>(properties._style);
  storeData<int>(properties._color);
  storeData<double>(properties._a);
  storeData<double>(properties._b);
  storeData<double>(properties._c);
  storeData<double>(properties._d);
}
//----------------------------------------------------------------------
template<> void pxl::oStreamer::storeData<pxl::ObjectBase>(const pxl::ObjectBase& object) {

  pxl::Id id = object.id();
  storeId(id);
  
  storeData<pxl::Relations>(object._motherRelations);
  storeData<pxl::Relations>(object._daughterRelations);
  
  if (object._ptrLayout) {storeData<bool>( true  ); storeData<pxl::Layout>(*object._ptrLayout);}
  else                   {storeData<bool>( false );}  

}
//----------------------------------------------------------------------
} // namespace pxl
//----------------------------------------------------------------------
//
//
//
//----------------------------------------------------------------------
bool pxl::iStreamer::putEvent(std::istream& cxxx, int mode, const std::string& infoCondition) {

  // clear the buffer stringstream
  _buffer.str("");
  _buffer.clear();

  char eventMarker[4] = {' ',' ',' ','\0'};
  std::string info = "";
  char compressionMode = ' ';
  
  int lengthInfo;
  int lengthInfoVerify;

  int lengthZip;
  int lengthZipVerify;
   
  redumpMemory(cxxx, eventMarker, 4);            // Event-Marker
  if (!cxxx.good()) {return false;}
  if (strcmp(iotl__eventMarker, eventMarker)) {pxl::exception("pxl::iStreamer::getEvent()","Invalid event marker while reading event header.");}
  restoreBasicTypeChar(cxxx, compressionMode);   // Compression-Mode
  
  restoreBasicTypeInt(cxxx, lengthInfo);         // Fwd-Jumper
  if (lengthInfo) {                              // info block
	    char* buffer = new char[lengthInfo];
	    cxxx.read(buffer, lengthInfo);
	    info.assign(buffer, lengthInfo);
	    delete [] buffer;
	  }
  restoreBasicTypeInt(cxxx, lengthInfoVerify);   // Bwd-Jumper & consistency check
  
  if (mode == -1 && infoCondition != info) {mode = 0;}
  
  restoreBasicTypeInt(cxxx, lengthZip);          // Fwd-Jumper
  if (mode == 0) {                               // zip block
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
  restoreBasicTypeInt(cxxx, lengthZipVerify);        // Bwd-Jumper & consistency check

  // seek the beginning:
  _buffer.seekg(0);
  
  return ((lengthInfo == lengthInfoVerify) && (lengthZip == lengthZipVerify));
}
//----------------------------------------------------------------------
bool pxl::iStreamer::previous(std::istream& cxxx) {
  
  int pos = cxxx.tellg();   
  int lengthInfoVerify;
  int lengthZipVerify;

  // read zip length information:
  pos -= 4;
  if (pos >= 0) {cxxx.seekg(pos);} 
  else {return false;}
  restoreBasicTypeInt(cxxx, lengthZipVerify); // Bwd-Jumper & consistency check
  if (lengthZipVerify < 0) {pxl::exception("pxl::iStreamer::previous()","Invalid zip block length information.");}
  
  pos -= lengthZipVerify; // zip block
  pos -= 4; // zip block length
  
  // read info length information:
  pos -= 4;
  if (pos >= 0) {cxxx.seekg(pos);} 
  else {return false;}
  restoreBasicTypeInt(cxxx, lengthInfoVerify); // Bwd-Jumper & consistency check
  if (lengthInfoVerify < 0) {pxl::exception("pxl::iStreamer::previous()","Invalid info block length information.");}
  
  pos -= lengthInfoVerify; // info block
  pos -= 4; // info block length
  
  pos -= 1; // compression mode
  pos -= 4; // event marker
  
  if (pos >= 0) {cxxx.seekg(pos); return true;}
   
  pxl::exception("pxl::iStreamer::previous()","Invalid block length information.");
  return false;
  
}
//----------------------------------------------------------------------
int pxl::iStreamer::unzipEventData( std::istream &in, int nBytes) {
	
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
namespace pxl {
//----------------------------------------------------------------------
template<> pxl::Id pxl::iStreamer::restoreData<pxl::Relations>(pxl::Relations& relations) {

  relations.clearContainer();

  int size; restoreData<int>(size);
  for (int i=0; i<size; i++) {
       pxl::MutableId persistentId = 0; 
       restoreId(persistentId);
       // provide orphan relation
       if (persistentId != 0) relations.set(persistentId, new pxl::WkPtr<pxl::Orphan>);
       }
  return 0;       
}
//----------------------------------------------------------------------
template<> pxl::Id pxl::iStreamer::restoreData<pxl::Objects>(pxl::Objects& objects) {

  //std::cerr << "pxl::iStreamer::restoreData<pxl::Objects>(pxl::Objects& objects): START" << std::endl;

  objects.clearContainer();

  int size;   
  std::string key;
  pxl::MutableId persistentId = 0;

  // reload objects (with orphan relations)  
  restoreData<int>(size);
  for (int i=0; i<size; i++) {
      pxl::ObjectBase* item = 0;
      persistentId = restoreAbstractObject(&item);
      objects.set(*item);
      objects._copyHistory.set(persistentId, item);
      }

  // redirect relations: loop in PTL style
  for (pxl::Objects::Iterator iter(objects); !iter.isDone(); iter.next()) {
       pxl::ObjectBase* pNew = iter.item();
       
       // mother relations
       for (pxl::Relations::TypeIterator<pxl::WkPtr<pxl::Orphan> > iterRelations(pNew->_motherRelations); !iterRelations.isDone(); iterRelations.next()) {
		   pxl::ObjectBase* pNewRel = objects._copyHistory.find(iterRelations.key(), 0);
		   if (pNewRel) {pNew->linkMother(*pNewRel);}
		   else {std::cerr << "pxl::iStreamer::restoreData<pxl::Objects>(...): WARNING: some original objects seem to have invalid relations. (persistent id is " << iterRelations.key() << ")" << std::endl;}
       }

       // daughter relations
       // have been set automatically above

       // remove orphan relations
       // mother relations
       pxl::Relations::TypeIterator<pxl::WkPtr<pxl::Orphan> > iterM(pNew->_motherRelations); 
       while (!iterM.isDone()) {
		   pxl::Id id = iterM.key();
		   iterM.next();
		   pNew->_motherRelations.remove(id);
           }
       
       // daughter relations
       pxl::Relations::TypeIterator<pxl::WkPtr<pxl::Orphan> > iterD(pNew->_daughterRelations); 
       while (!iterD.isDone()) {
		   pxl::Id id = iterD.key();
		   iterD.next();
		   pNew->_daughterRelations.remove(id);
           }
      }

  // reload & redirect index
  restoreData<int>(size);
  for (int i=0; i<size; i++) {
      restoreData<std::string>(key);
      restoreId(persistentId);
            
      pxl::ObjectBase* pNew = objects._copyHistory.find(persistentId, 0);
      
      if (pNew) {objects._index.set(key, pNew);}
	  else {
	  	std::cerr << "pxl::iStreamer::restoreData<pxl::Objects>(...): WARNING: some original indices could not be reconstructed." << std::endl;
	  	std::cerr << key << ": " << persistentId << std::endl;
	  	}
      }

  //std::cerr << "pxl::iStreamer::restoreData<pxl::Objects>(pxl::Objects& objects): END" << std::endl;
      
  return 0;       
}
//----------------------------------------------------------------------
template<> pxl::Id pxl::iStreamer::restoreData<pxl::Layout>(pxl::Layout& properties) {
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
template<> pxl::Id pxl::iStreamer::restoreData<pxl::ObjectBase>(pxl::ObjectBase& object) {

  pxl::MutableId persistentId;
  bool hasLayout;
  
  restoreId(persistentId);
  
  restoreData<pxl::Relations>(object._motherRelations);
  restoreData<pxl::Relations>(object._daughterRelations);
  
  restoreData<bool>(hasLayout);
  if (hasLayout) {restoreData<pxl::Layout>( object.layout() );}
  else           {delete object._ptrLayout; object._ptrLayout = 0;}
  
  return persistentId;
}
//----------------------------------------------------------------------
pxl::Id pxl::iStreamer::restoreAbstractObject(pxl::ObjectBase** ppobj) {

  std::string objectTypeId; 
  std::string dataTypeId;
   
  pxl::BasicIoStreamer::restoreBasicTypeCStr(_buffer, objectTypeId);   
  pxl::BasicIoStreamer::restoreBasicTypeCStr(_buffer, dataTypeId); 
  
  return pxl::TypeManager::instance().restoreObject(*this, ppobj, objectTypeId, dataTypeId);
}
//----------------------------------------------------------------------
} // namespace pxl
//----------------------------------------------------------------------
void pxl::TypeManager::registerAgent(pxl::TypeAgentBase* agent) {
  
  // for debugging purposes:
//  std::cout << "pxl::TypeManager::registerAgent(): type registered: [" 
//            << agent->getObjectTypeId() << "," << agent->getDataTypeId() 
//            << "]" << std::endl;
//            //<< "] = '" << agent->getCppTypeId() << "'" << std::endl;

  // add entries to maps:
  _agentsByIotlTypeId.set(pxl::TypeIdKey(agent->getObjectTypeId(), agent->getDataTypeId()), agent);
  _agentsByCppTypeId.set(agent->getCppTypeId(), agent);
}
//----------------------------------------------------------------------
pxl::Id pxl::TypeManager::restoreObject(iStreamer& input, pxl::ObjectBase** ppobj, const std::string& objectTypeId, const std::string& dataTypeId) const {

  TypeAgentBase* agent = _agentsByIotlTypeId.find(pxl::TypeIdKey(objectTypeId, dataTypeId), 0);
  
  if (agent) {return agent->restoreObject(input, ppobj);}
  
  std::string message; 
  message += "undeclared IOTL type: [";
  message += objectTypeId;
  message += ","; 
  message += dataTypeId; 
  message += "]"; 
  pxl::exception("pxl::TypeManager::restoreObject()", message); 

  return 0;
 }
//----------------------------------------------------------------------
pxl::Id pxl::TypeManager::restoreObject(pxl::iStreamer& input,  pxl::ObjectBase& obj,  const std::string& cppTypeId) const {

  TypeAgentBase* agent = _agentsByCppTypeId.find(cppTypeId, 0);
  
  if (agent) {return agent->restoreObject(input, obj);}
  
  std::string message; 
  message += "undeclared C++ type: '";
  message += cppTypeId;
  message += "'";
  pxl::exception("pxl::TypeManager::restoreObject()", message); 

  return 0;
 }
//----------------------------------------------------------------------
void pxl::TypeManager::storeObject(oStreamer& output, const pxl::ObjectBase& obj, const std::string& cppTypeId) const {

  TypeAgentBase* agent = _agentsByCppTypeId.find(cppTypeId, 0);
  
  if (agent) {agent->storeObject(output, obj); return;}
  
  std::string message; 
  message += "undeclared C++ type: '";
  message += cppTypeId;
  message += "'";
  pxl::exception("pxl::TypeManager::storeObject()", message); 
 }
//----------------------------------------------------------------------
namespace pxl {TypeManager* TypeManager::_instance = 0;}
//----------------------------------------------------------------------
pxl::TypeManager& pxl::TypeManager::instance() {
  if (!_instance) {_instance = new TypeManager;}
  return *_instance;	
}
//----------------------------------------------------------------------
//
//
//
template class pxl::oDiskFileVx<pxl::oStreamer>;
template class pxl::iDiskFileVx<pxl::iStreamer>;

template class pxl::TypeAgent<pxl::Object<int> >;
template class pxl::TypeAgent<pxl::CowObject<int> >;

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
iotl__declareObjectTypeExplicit(pxl::BasicObject, "\3Bo", pol__BasicObject)
iotl__declareObjectTypeExplicit(pxl::BasicObjectManager, "\3Bom", pol__BasicObjectManager)
iotl__declareObjectTypeExplicit(pxl::Particle, "\3Pa", pol__Particle)
iotl__declareObjectTypeExplicit(pxl::Vertex, "\3Vx", pol__Vertex)
// iotl__declareObjectExplicit(pxl::Collision, "\3Co", pol__Collision)
iotl__declareObjectTypeExplicit(pxl::EventView, "\3Ev", pol__EventView)
iotl__declareObjectTypeExplicit(pxl::AnalysisProcess, "\3Ap", pol__AnalysisProcess)
iotl__declareObjectTypeExplicit(pxl::AnalysisFork, "\3Af", pol__AnalysisFork)
//----------------------------------------------------------------------
iotl__declareDataTypeExplicit(pxl::BasicObjectData, "\3BoD", pol__BasicObjectData, data,
{
  storeData<bool>(data._locked);
  storeData<int>(data._monteCarloMode);
  storeData<std::string>(data._name);
  storeData<int>(data._status);
  storeData<int>(data._workflag);
  storeData(data._userRecords);
  // storeData<pxl::CppPointers>(_cppPointers); // we don't store the c++ pointers!
},{
  restoreData<bool>(data._locked);
  restoreData<int>(data._monteCarloMode);
  restoreData<std::string>(data._name);
  restoreData<int>(data._status);
  restoreData<int>(data._workflag);
  restoreData(data._userRecords);
  // restoreData<pxl::CppPointers>(_cppPointers); // we don't restore the c++ pointers!
})
//----------------------------------------------------------------------
iotl__declareDataTypeExplicit(pxl::BasicObjectManagerData, "\3BomD", pol__BasicObjectManagerData, data,
{
  storeData<pxl::Objects>(data._objects);
  storeData<pxl::BasicObjectData>(data);
},{
  restoreData<pxl::Objects>(data._objects);
  restoreData<pxl::BasicObjectData>(data);
})
//----------------------------------------------------------------------
iotl__declareDataTypeExplicit(pxl::Basic4VectorData, "\3B4vD", pol__Basic4VectorData, data,
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
iotl__declareDataTypeExplicit(pxl::Basic3VectorData, "\3B3vD", pol__Basic3VectorData, data,
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
iotl__declareDataTypeExplicit(pxl::ParticleData, "\3PaD", pol__ParticleData, data,
{
  storeData<double>(data._charge);
  storeData<int>(data._particleId);
  storeData(data.vector());
  storeData<pxl::BasicObjectData>(data);
},{
  restoreData<double>(data._charge);
  restoreData<int>(data._particleId);
  restoreData(data.vector(pxl::set));
  restoreData<pxl::BasicObjectData>(data);
})
//----------------------------------------------------------------------
iotl__declareDataTypeExplicit(pxl::VertexData, "\3VxD", pol__VertexData, data,
{
  storeData(data.vector());
  storeData<pxl::BasicObjectData>(data);
},{
  restoreData(data.vector(pxl::set));
  restoreData<pxl::BasicObjectData>(data);
})
//----------------------------------------------------------------------
// iotl__declareDataExplicit(pxl::CollisionData, "\3CoD", pol__CollisionData, data,
// {
//   storeData<pxl::BasicObjectData>(data);
// },{
//   restoreData<pxl::BasicObjectData>(data);
// })
//----------------------------------------------------------------------
iotl__declareDataTypeExplicit(pxl::EventViewData, "\3EvD", pol__EventViewData, data,
{
  storeData<pxl::BasicObjectManagerData>(data);
},{
  restoreData<pxl::BasicObjectManagerData>(data);
})
//----------------------------------------------------------------------
iotl__declareDataTypeExplicit(pxl::AnalysisProcessData, "\3ApD", pol__AnalysisProcessData, data,
{
  storeData<pxl::BasicObjectManagerData>(data);
},{
  restoreData<pxl::BasicObjectManagerData>(data);
})
//----------------------------------------------------------------------
iotl__declareDataTypeExplicit(pxl::AnalysisForkData, "\3AfD", pol__AnalysisForkData, data, 
{
  storeData<pxl::BasicObjectManagerData>(data);
},{
  restoreData<pxl::BasicObjectManagerData>(data);
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
pxl::ObjectOwner::ObjectOwner() : pxl::Vector<pxl::ObjectBase*>(), _copyHistory(), _index() {;}
pxl::ObjectOwner::ObjectOwner(const pxl::ObjectOwner& original) : pxl::Vector<pxl::ObjectBase*>(), _copyHistory(), _index() {

  // copy objects: loop in STL style
  for (pxl::ObjectOwner::StlConstIterator iter = original.begin(); 
                                      iter != original.end(); 
				      iter++) {
       pxl::ObjectBase* pOld = *iter;
       pxl::ObjectBase* pNew = pOld->clone();

       set(*pNew);
       _copyHistory.set(pOld->id(), pNew);
      }
  
  // redirect relations: loop in PTL style
  for (pxl::ObjectOwner::Iterator iter(original); !iter.isDone(); iter.next()) {
       pxl::ObjectBase* pOld = iter.item();
       pxl::ObjectBase* pNew = _copyHistory.find(pOld->id(), 0);
       
       // mother relations
       for (pxl::Relations::Iterator iter(pOld->getMotherRelations()); !iter.isDone(); iter.next()) {
	   pxl::ObjectBase* pOldRel = iter.item()->pointer();
	   pxl::ObjectBase* pNewRel = _copyHistory.find(pOldRel->id(), 0);
	   
	   if (pOldRel) {	   
	      if (pNewRel) {pNew->linkMother(*pNewRel);}
	      else {std::cerr << "pxl::ObjectOwner::ObjectOwner(...): WARNING: some original objects had relations to objects of other owners." << std::endl;}
	      }
	   else {
	      std::cerr << "pxl::ObjectOwner::ObjectOwner(...): WARNING: some originally related objects no longer exist" << std::endl;
	      }   
       }

       // daughter relations
       // have been set automatically above
      }
  
  // redirect index:
  for (pxl::Index::Iterator iter(original._index); !iter.isDone(); iter.next()) {
       pxl::ObjectBase* pOld = iter.item();
       pxl::ObjectBase* pNew = _copyHistory.find(pOld->id(), 0);
       
       if (pNew) {_index.set(iter.key(),pNew);}
       else {std::cerr << "pxl::ObjectOwner::ObjectOwner(...): WARNING: some original indices pointed to objects of other owners." << std::endl;}
      }
}
//----------------------------------------------------------------------
void pxl::Objects::clearContainer() {
  for (StlConstIterator iter = _container.begin(); iter != _container.end(); iter++) {
      delete (*iter);
      }
  _container.clear();
  _copyHistory.clearContainer();
  _index.clearContainer();
}
//----------------------------------------------------------------------
void pxl::Objects::set(pxl::ObjectBase& item) {item._refObjectOwner = this; _container.push_back(&item);}
//----------------------------------------------------------------------
void pxl::Objects::remove(pxl::ObjectBase& item) {

  // search & remove possible indices (multiple occurrences possible!)
  for (pxl::Index::StlConstIterator iter = _index.begin(); iter != _index.end(); iter++) {
       if (&item == iter->second) {_index.remove(iter->first);}
      }
  
  // search & remove possible copy history (multiple occurrences *not* possible!)
  for (pxl::CopyHistory::StlConstIterator iter = _copyHistory.begin(); iter != _copyHistory.end(); iter++) {
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
bool pxl::Objects::has(const pxl::ObjectBase& item) const {return (item._refObjectOwner == this);}
//----------------------------------------------------------------------
void pxl::WkPtrBase::notifyDeleted() {

  _objectRef = 0; 
  if (_notifyChainOut) _notifyChainOut->notifyDeleted(); 
  _notifyChainIn = 0; 
  _notifyChainOut = 0; 

}
//----------------------------------------------------------------------
void pxl::WkPtrBase::connect(pxl::ObjectBase* pointer) {

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
void pxl::ObjectBase::linkMother(pxl::ObjectBase& target)  {

  pxl::WkPtrBase* pTarget = target.createSelfWkPtr();
  pxl::WkPtrBase* pThis   = this->createSelfWkPtr();
  
  if (target._refObjectOwner != this->_refObjectOwner) {std::cerr << "pxl::ObjectBase::linkDaughter(...): WARNING: mother and daughter have not the same object owner!" << std::endl;}

  this->_motherRelations.set(target.id(), pTarget);
  target._daughterRelations.set(this->id(), pThis);

}
//----------------------------------------------------------------------
void pxl::ObjectBase::linkDaughter(pxl::ObjectBase& target) {

  pxl::WkPtrBase* pTarget = target.createSelfWkPtr();
  pxl::WkPtrBase* pThis   = this->createSelfWkPtr();
  
  if (target._refObjectOwner != this->_refObjectOwner) {std::cerr << "pxl::ObjectBase::linkDaughter(...): WARNING: mother and daughter have not the same object owner!" << std::endl;}
  
  this->_daughterRelations.set(target.id(), pTarget);
  target._motherRelations.set(this->id(), pThis);

}
//----------------------------------------------------------------------
void pxl::ObjectBase::unlinkMother(pxl::ObjectBase& target)  {

  this->_motherRelations.remove(target.id());
  target._daughterRelations.remove(this->id());

}
//----------------------------------------------------------------------
void pxl::ObjectBase::unlinkDaughter(pxl::ObjectBase& target) {

  this->_daughterRelations.remove(target.id());
  target._motherRelations.remove(this->id());

}
//----------------------------------------------------------------------
void pxl::ObjectBase::unlinkMothers() {

  for (pxl::Relations::Iterator iter(_motherRelations); !iter.isDone(); iter.next()) {
       if (iter.item()->_objectRef) {
       	  iter.item()->_objectRef->_daughterRelations.remove(this->id());
       	  }
       }
  _motherRelations.clearContainer();
  
}  
//----------------------------------------------------------------------
void pxl::ObjectBase::unlinkDaughters() {

  for (pxl::Relations::Iterator iter(_daughterRelations); !iter.isDone(); iter.next()) {
       if (iter.item()->_objectRef) {
       	  iter.item()->_objectRef->_motherRelations.remove(this->id());
       	  }
       }
  _daughterRelations.clearContainer();

}
//----------------------------------------------------------------------
std::ostream& pxl::ObjectBase::printDecayTree(int level, std::ostream& os, int pan) const {

  int daughters = 0;
  
  print(level, os, pan);

  for (pxl::Relations::Iterator iter(_daughterRelations); !iter.isDone(); iter.next()) {
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
std::ostream& pxl::ObjectBase::print(int level, std::ostream& os, int pan) const {return printPan1st(os, pan) << "pxl::ObjectBase with id: " << id() << std::endl;}
//----------------------------------------------------------------------
std::ostream& pxl::ObjectBase::printPan1st(std::ostream& os, int pan) const {for (int p=0;p<pan-2;p++) {os << "|  ";} if (pan-1 > 0) os << "+--"; if (pan) os << "{ "; return os;}
std::ostream& pxl::ObjectBase::printPan(std::ostream& os, int pan) const {for (int p=0;p<pan-1;p++) {os << "|  ";}; if (pan) os << "| "; return os;}
//----------------------------------------------------------------------
namespace pxl
{
template <> std::ostream& pxl::Object<int>::print(int level, std::ostream& os, int pan) const {return printPan1st(os, pan) << "pxl::Object<int> with id: " << id() << " is " << get() << std::endl;}
template <> std::ostream& pxl::Object<unsigned int>::print(int level, std::ostream& os, int pan) const {return printPan1st(os, pan) << "pxl::Object<unsigned int> with id: " << id() << " is " << get() << std::endl;}
template <> std::ostream& pxl::Object<bool>::print(int level, std::ostream& os, int pan) const {return printPan1st(os, pan) << "pxl::Object<bool> with id: " << id() << " is " << get() << std::endl;}
template <> std::ostream& pxl::Object<double>::print(int level, std::ostream& os, int pan) const {return printPan1st(os, pan) << "pxl::Object<double> with id: " << id() << " is " << get() << std::endl;}
template <> std::ostream& pxl::Object<float>::print(int level, std::ostream& os, int pan) const {return printPan1st(os, pan) << "pxl::Object<float> with id: " << id() << " is " << get() << std::endl;}
//----------------------------------------------------------------------
template <> std::ostream& pxl::CowObject<int>::print(int level, std::ostream& os, int pan) const {return printPan1st(os, pan) << "pxl::CowObject<int> with id: " << id() << " is " << get() << std::endl;}
template <> std::ostream& pxl::CowObject<unsigned int>::print(int level, std::ostream& os, int pan) const {return printPan1st(os, pan) << "pxl::CowObject<unsigned int> with id: " << id() << " is " << get() << std::endl;}
template <> std::ostream& pxl::CowObject<bool>::print(int level, std::ostream& os, int pan) const {return printPan1st(os, pan) << "pxl::CowObject<bool> with id: " << id() << " is " << get() << std::endl;}
template <> std::ostream& pxl::CowObject<double>::print(int level, std::ostream& os, int pan) const {return printPan1st(os, pan) << "pxl::CowObject<double> with id: " << id() << " is " << get() << std::endl;}
template <> std::ostream& pxl::CowObject<float>::print(int level, std::ostream& os, int pan) const {return printPan1st(os, pan) << "pxl::CowObject<float> with id: " << id() << " is " << get() << std::endl;}
}
//----------------------------------------------------------------------
std::ostream& operator << (std::ostream& cxxx, const pxl::ObjectBase& obj) {return obj.print(0, cxxx, 0);}
//----------------------------------------------------------------------
// EXPLICIT INSTANCIATION
template class pxl::Vector<int>;
template class pxl::Map<int,int>;
template class pxl::Objects::TypeIterator<pxl::WkPtr<int> >;
template class pxl::Relations::TypeIterator<pxl::WkPtr<int> >;

template class pxl::Ptr<int>;
template class pxl::WkPtrSpec<int, pxl::Object<int> >;
template class pxl::WkPtr<int>;
template class pxl::CowWkPtr<int>;
template class pxl::SpyWkPtr<int>;

template class pxl::Object<int>;
template class pxl::CowObject<int>;
template class pxl::SpyObject<int>;
//----------------------------------------------------------------------
#ifndef MERGED_PXL
	#include "pol.hh"
#endif
//
//
//----------------------------------------------------------------------
namespace pxl {
//----------------------------------------------------------------------
bool const operator==(const pxl::Basic4VectorData& obj1, const pxl::Basic4VectorData& obj2)
{return (obj1.getX() == obj2.getX() && obj1.getY() == obj2.getY() && obj1.getZ() == obj2.getZ() && obj1.getE() == obj2.getE());}

bool const operator!=(const pxl::Basic4VectorData& obj1, const pxl::Basic4VectorData& obj2)
{return (obj1.getX() != obj2.getX() || obj1.getY() != obj2.getY() || obj1.getZ() != obj2.getZ() || obj1.getE() != obj2.getE());}
//----------------------------------------------------------------------
bool const operator==(const pxl::Basic3VectorData& obj1, const pxl::Basic3VectorData& obj2) 
{return (obj1.getX() == obj2.getX() && obj1.getY() == obj2.getY() && obj1.getZ() == obj2.getZ());}

bool const operator!=(const pxl::Basic3VectorData& obj1, const pxl::Basic3VectorData& obj2) 
{return (obj1.getX() != obj2.getX() || obj1.getY() != obj2.getY() || obj1.getZ() != obj2.getZ());}
//----------------------------------------------------------------------
bool const operator==(const pxl::VertexData& obj1, const pxl::VertexData& obj2)
{return (obj1.vector() == obj2.vector());}

bool const operator!=(const pxl::VertexData& obj1, const pxl::VertexData& obj2) 
{return (obj1.vector() != obj2.vector());}
//----------------------------------------------------------------------
bool const operator==(const pxl::ParticleData& obj1, const pxl::ParticleData& obj2) 
{return (obj1.vector() == obj2.vector() && obj1.getCharge() == obj2.getCharge());}

bool const operator!=(const pxl::ParticleData& obj1, const pxl::ParticleData& obj2) 
{return (obj1.vector() != obj2.vector() || obj1.getCharge() != obj2.getCharge());}
//----------------------------------------------------------------------
} // namespace pxl
//----------------------------------------------------------------------
pxl::ObjectBase* pxl::AnalysisProcess::clone() const {return new AnalysisProcess(*this);}
pxl::WkPtrBase* pxl::AnalysisProcess::createSelfWkPtr() {return new AnalysisProcessWkPtr(*this);}
void pxl::AnalysisProcess::storeYourSelf(pxl::oStreamer& output) const {output.storeObject(*this);}
//----------------------------------------------------------------------
void pxl::AnalysisFork::beginJob(const pxl::Objects* input) {

  for (pxl::Objects::TypeIterator<pxl::AnalysisFork> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().beginJob(input);
      }
      
  for (pxl::Objects::TypeIterator<pxl::AnalysisProcess> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().beginJob(input);
      }
      
}
//----------------------------------------------------------------------
void pxl::AnalysisFork::buildTemplate(int mode) {

  for (pxl::Objects::TypeIterator<pxl::AnalysisFork> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().buildTemplate(mode);
      }
      
  for (pxl::Objects::TypeIterator<pxl::AnalysisProcess> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().buildTemplate(mode);
      }
      
}
//----------------------------------------------------------------------
void pxl::AnalysisFork::beginRun(const pxl::Objects* input) {

  for (pxl::Objects::TypeIterator<pxl::AnalysisFork> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().beginRun(input);
      }
      
  for (pxl::Objects::TypeIterator<pxl::AnalysisProcess> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().beginRun(input);
      }
      
}
//----------------------------------------------------------------------
void pxl::AnalysisFork::analyseEvent(const pxl::Objects* input) {

  for (pxl::Objects::TypeIterator<pxl::AnalysisFork> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().analyseEvent(input);
      }
      
  for (pxl::Objects::TypeIterator<pxl::AnalysisProcess> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().analyseEvent(input);
      }
      
}
//----------------------------------------------------------------------
void pxl::AnalysisFork::finishEvent(const pxl::Objects* input) {

  for (pxl::Objects::TypeIterator<pxl::AnalysisFork> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().finishEvent(input);
      }
      
  for (pxl::Objects::TypeIterator<pxl::AnalysisProcess> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().finishEvent(input);
      }
      
}
//----------------------------------------------------------------------
void pxl::AnalysisFork::endRun(const pxl::Objects* input) {

  for (pxl::Objects::TypeIterator<pxl::AnalysisFork> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().endRun(input);
      }
      
  for (pxl::Objects::TypeIterator<pxl::AnalysisProcess> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().endRun(input);
      }
      
}
//----------------------------------------------------------------------
void pxl::AnalysisFork::endJob(const pxl::Objects* input) {

  for (pxl::Objects::TypeIterator<pxl::AnalysisFork> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().endJob(input);
      }
      
  for (pxl::Objects::TypeIterator<pxl::AnalysisProcess> iter(get().getObjects()); !iter.isDone(); iter.next()) {
      iter.object().endJob(input);
      }
      
}
//----------------------------------------------------------------------
pxl::ObjectBase* pxl::AnalysisFork::clone() const {return new AnalysisFork(*this);}
pxl::WkPtrBase* pxl::AnalysisFork::createSelfWkPtr() {return new AnalysisForkWkPtr(*this);}
void pxl::AnalysisFork::storeYourSelf(pxl::oStreamer& output) const {output.storeObject(*this);}
//----------------------------------------------------------------------
//
//
//
//----------------------------------------------------------------------
namespace pxl {
//----------------------------------------------------------------------
template <>  std::ostream& pxl::CowObject<pxl::ParticleData>::print(int level, std::ostream& os, int pan) const {

//  printPan1st(os, pan) << "pxl::CowObject<pxl::ParticleData> with id: " << id() << " (data socket currently located at " << _dataSocket << ")" << std::endl;
//  printPan(os, pan)    << "     name: " << get().getName() << std::endl;
  printPan1st(os, pan) << "Particle: '" << get().getName() << "', p = (" 
                             << get().vector().getPt() << ", " 
                             << get().vector().getPz() << ") m = " 
                             << get().vector().getMass() << std::endl;
  return os;

}
//----------------------------------------------------------------------
template <>  std::ostream& pxl::CowObject<pxl::VertexData>::print(int level, std::ostream& os, int pan) const {

//  printPan1st(os, pan) << "pxl::CowObject<pxl::VertexData> with id: " << id() << " (data socket currently located at " << _dataSocket << ")" << std::endl;
//  printPan(os, pan)    << "     name: " << get().getName() << std::endl;
  printPan1st(os, pan) << "Vertex: '" << get().getName() << "', x = (" 
                             << get().vector().getX() << ", " 
                             << get().vector().getY() << ", " 
                             << get().vector().getZ() << ")" << std::endl;
  return os;

}
//----------------------------------------------------------------------
template <>  std::ostream& pxl::CowObject<pxl::CollisionData>::print(int level, std::ostream& os, int pan) const {

//  printPan1st(os, pan) << "pxl::CowObject<pxl::CollisionData> with id: " << id() << " (data socket currently located at " << _dataSocket << ")" << std::endl;
//  printPan(os, pan)    << "     name: " << get().getName() << std::endl;
  printPan1st(os, pan) << "Collision: " << get().getName() << std::endl;
  return os;

}
//----------------------------------------------------------------------
template <>  std::ostream& pxl::Object<pxl::EventViewData>::print(int level, std::ostream& os, int pan) const {

//  printPan1st(os, pan) << "pxl::CowObject<pxl::EventViewData> with id: " << id() << std::endl;
//  printPan(os, pan)    << "     name: " << get().getName() << std::endl;
  printPan1st(os, pan) << "EventView: " << get().getName() << std::endl;
  for (pxl::Objects::Iterator iter(get().getObjects()); !iter.isDone(); iter.next()) {
      if (iter.object().getMotherRelations().getSize() == 0) {iter.object().printDecayTree(level, os, pan);}
      }
  return os;

}
//----------------------------------------------------------------------
template <>  std::ostream& pxl::Object<pxl::AnalysisProcessData>::print(int level, std::ostream& os, int pan) const {

//  printPan1st(os, pan) << "pxl::CowObject<pxl::AnalysisProcessData> with id: " << id() << std::endl;
//  printPan(os, pan)    << "     name: " << get().getName() << std::endl;
  printPan1st(os, pan) << "AnalysisProcess: " << get().getName() << std::endl;
  for (pxl::Objects::Iterator iter(get().getObjects()); !iter.isDone(); iter.next()) {
      if (iter.object().getMotherRelations().getSize() == 0) {iter.object().printDecayTree(level, os, pan);}
      }
  return os;

}
//----------------------------------------------------------------------
template <>  std::ostream& pxl::Object<pxl::AnalysisForkData>::print(int level, std::ostream& os, int pan) const {

//  printPan1st(os, pan) << "pxl::CowObject<pxl::AnalysisForkData> with id: " << id() << std::endl;
//  printPan(os, pan)    << "     name: " << get().getName() << std::endl;
  printPan1st(os, pan) << "AnalysisFork: " << get().getName() << std::endl;
  for (pxl::Objects::Iterator iter(get().getObjects()); !iter.isDone(); iter.next()) {
      if (iter.object().getMotherRelations().getSize() == 0) {iter.object().printDecayTree(level, os, pan);}
      }
  return os;

}
//----------------------------------------------------------------------
}
//----------------------------------------------------------------------
#ifndef MERGED_PXL
	#include "ePaxPxl/ePax/interface/ePax.h"
#endif	
//
//
//----------------------------------------------------------------------
// ePax
//----------------------------------------------------------------------
iotl__declareObjectTypeExplicit(ePaxParticle,  "\4ePa", ePaxParticle)
iotl__declareObjectTypeExplicit(ePaxVertex,    "\4eVe", ePaxVertex)
iotl__declareObjectTypeExplicit(ePaxCollision, "\4eCo", ePaxCollision)
iotl__declareObjectTypeExplicit(ePaxEventView, "\4eEv", ePaxEventView)
iotl__declareObjectTypeExplicit(ePaxAnalysisProcess, "\4eAp", ePaxAnalysisProcess)
iotl__declareObjectTypeExplicit(ePaxAnalysisFork,    "\4eAf", ePaxAnalysisFork)
//----------------------------------------------------------------------
pxl::ObjectBase* ePaxParticle::clone() const {return new ePaxParticle(*this);}
std::ostream& ePaxParticle::print(int level, std::ostream& os, int pan) const {return pxl::CowObject<pxl::ParticleData>::print(level, os, pan);}
pxl::WkPtrBase* ePaxParticle::createSelfWkPtr() {return new ePaxParticleWkPtr(*this);}
void ePaxParticle::storeYourSelf(pxl::oStreamer& output) const {output.storeObject(*this);}
//----------------------------------------------------------------------
pxl::ObjectBase* ePaxVertex::clone() const {return new ePaxVertex(*this);}
std::ostream& ePaxVertex::print(int level, std::ostream& os, int pan) const {return pxl::CowObject<pxl::VertexData>::print(level, os, pan);}
pxl::WkPtrBase* ePaxVertex::createSelfWkPtr() {return new ePaxVertexWkPtr(*this);}
void ePaxVertex::storeYourSelf(pxl::oStreamer& output) const {output.storeObject(*this);}
//----------------------------------------------------------------------
pxl::ObjectBase* ePaxCollision::clone() const {return new ePaxCollision(*this);}
std::ostream& ePaxCollision::print(int level, std::ostream& os, int pan) const {return pxl::CowObject<pxl::CollisionData>::print(level, os, pan);}
pxl::WkPtrBase* ePaxCollision::createSelfWkPtr() {return new ePaxCollisionWkPtr(*this);}
void ePaxCollision::storeYourSelf(pxl::oStreamer& output) const {output.storeObject(*this);}
//----------------------------------------------------------------------
pxl::ObjectBase* ePaxEventView::clone() const {return new ePaxEventView(*this);}
std::ostream& ePaxEventView::print(int level, std::ostream& os, int pan) const {return pxl::Object<pxl::EventViewData>::print(level, os, pan);}
pxl::WkPtrBase* ePaxEventView::createSelfWkPtr() {return new ePaxEventViewWkPtr(*this);}
void ePaxEventView::storeYourSelf(pxl::oStreamer& output) const {output.storeObject(*this);}
//----------------------------------------------------------------------
pxl::ObjectBase* ePaxAnalysisProcess::clone() const {return new ePaxAnalysisProcess(*this);}
std::ostream& ePaxAnalysisProcess::print(int level, std::ostream& os, int pan) const {return pxl::Object<pxl::AnalysisProcessData>::print(level, os, pan);}
pxl::WkPtrBase* ePaxAnalysisProcess::createSelfWkPtr() {return new ePaxAnalysisProcessWkPtr(*this);}
void ePaxAnalysisProcess::storeYourSelf(pxl::oStreamer& output) const {output.storeObject(*this);}
//----------------------------------------------------------------------
pxl::ObjectBase* ePaxAnalysisFork::clone() const {return new ePaxAnalysisFork(*this);}
std::ostream& ePaxAnalysisFork::print(int level, std::ostream& os, int pan) const {return pxl::Object<pxl::AnalysisForkData>::print(level, os, pan);}
pxl::WkPtrBase* ePaxAnalysisFork::createSelfWkPtr() {return new ePaxAnalysisForkWkPtr(*this);}
void ePaxAnalysisFork::storeYourSelf(pxl::oStreamer& output) const {output.storeObject(*this);}
//----------------------------------------------------------------------
