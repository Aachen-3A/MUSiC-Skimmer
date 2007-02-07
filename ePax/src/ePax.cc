#include <iostream>
#include <sys/times.h>
#include <unistd.h>

#include "ePaxPxl/ePax/interface/ePax.h"

namespace pxl {

double getCpuTime()
{
    // static method returning system CPU time.
    struct tms cpt;
    times(&cpt);
    return (double)(cpt.tms_utime + cpt.tms_stime) / 100.0;
}

void exception(const std::string& routine, const std::string& message)
{
    std::cerr << routine << ": " << message << std::endl;
    std::cerr << "pxl::exception(): throwing pxl::Exception." << std::endl;
    throw pxl::Exception(routine, message);
}

} // namespace pxl


namespace pxl {

#define PRINT_NATIVE(type) \
template<> \
std::ostream& CowObject<type>::print(int level, std::ostream& os, int pan) const \
{ \
    return printPan1st(os, pan) << "pxl::CowObject<" #type "> with id: " \
                                << id() << " is " << get() << std::endl; \
}

PRINT_NATIVE(int)
PRINT_NATIVE(unsigned int)
PRINT_NATIVE(bool)
PRINT_NATIVE(double)
PRINT_NATIVE(float)

#undef PRINT_NATIVE

// EXPLICIT INSTANCIATION
template class CowObject<int>;
template class CowObject<unsigned int>;
template class CowObject<bool>;
template class CowObject<double>;
template class CowObject<float>;

} // namespace pxl

namespace pxl {

// EXPLICIT INSTANCIATION
template class CowWkPtr<int>;

} // namespace pxl

namespace pxl {

const pxl::Get get;
const pxl::Set set;

} // namespace pxl

namespace pxl {

// EXPLICIT INSTANCIATION
template class Map<int, int>;

} // namespace pxl


namespace pxl {

void ObjectBase::linkMother(pxl::ObjectBase& target)
{
    pxl::WkPtrBase* pTarget = target.createSelfWkPtr();
    pxl::WkPtrBase* pThis   = this->createSelfWkPtr(); 

    if (target._refObjectOwner != this->_refObjectOwner)
        // FIXME: exception?
        std::cerr << "pxl::ObjectBase::linkDaughter(...): WARNING: mother and daughter have not the same object holder!" << std::endl;

    this->_motherRelations.set(target.id(), pTarget);
    target._daughterRelations.set(this->id(), pThis);
}

void ObjectBase::linkDaughter(pxl::ObjectBase& target)
{
    pxl::WkPtrBase* pTarget = target.createSelfWkPtr();
    pxl::WkPtrBase* pThis   = this->createSelfWkPtr(); 

   if (target._refObjectOwner != this->_refObjectOwner)
       // FIXME: exception?
       std::cerr << "pxl::ObjectBase::linkDaughter(...): WARNING: mother and daughter have not the same object holder!" << std::endl;

   this->_daughterRelations.set(target.id(), pTarget);
   target._motherRelations.set(this->id(), pThis);
}

void ObjectBase::unlinkMother(pxl::ObjectBase& target)
{
    this->_motherRelations.remove(target.id());
    target._daughterRelations.remove(this->id());
}

void ObjectBase::unlinkDaughter(pxl::ObjectBase& target)
{
    this->_daughterRelations.remove(target.id());
    target._motherRelations.remove(this->id());  
}

void ObjectBase::unlinkMothers() {
    for(pxl::Relations::Iterator iter(_motherRelations);
        !iter.isDone(); iter.next()) {

        if (iter.item()->_objectRef)
            iter.item()->_objectRef->_daughterRelations.remove(this->id());
    }

    _motherRelations.clearContainer();
}

void ObjectBase::unlinkDaughters() {
    for(pxl::Relations::Iterator iter(_daughterRelations);
       !iter.isDone(); iter.next()) {

       if (iter.item()->_objectRef)
           iter.item()->_objectRef->_motherRelations.remove(this->id());
    }

    _daughterRelations.clearContainer();
}

std::ostream& ObjectBase::printDecayTree(int level, std::ostream& os, int pan) const
{
    int daughters = 0;

    print(level, os, pan);

    for(pxl::Relations::Iterator iter(_daughterRelations);
        !iter.isDone(); iter.next()) {

        if (iter.item()->_objectRef) {
            iter.item()->_objectRef->printDecayTree(level, os, pan + 1);
            daughters++;
        }
    }

   if (daughters && pan > 1) {
       for(int p = 0; p < pan; p++)
           os << "|  ";
       os << "*" << std::endl;
   }

   return os;
}

std::ostream& ObjectBase::print(int level, std::ostream& os, int pan) const
{
    return printPan1st(os, pan) << "pxl::ObjectBase with id: " << id() << std::endl;
}

std::ostream& ObjectBase::printPan1st(std::ostream& os, int pan) const
{
    for(int p = 0; p < pan - 2; p++)
        os << "|  ";
    if (pan - 1 > 0)
        os << "+--";
    if (pan)
        os << "{ ";
    return os;
}

std::ostream& ObjectBase::printPan(std::ostream& os, int pan) const
{
    for(int p = 0; p < pan - 1; p++)
        os << "|  ";
    if (pan)
        os << "| ";
    return os;
}

std::ostream& operator<<(std::ostream& cxxx, const pxl::ObjectBase& obj)
{
    return obj.print(0, cxxx, 0);
}

} // namespace pxl


namespace pxl {

#define PRINT_NATIVE(type) \
template<> \
std::ostream& Object<type>::print(int level, std::ostream& os, int pan) const \
{ \
    return printPan1st(os, pan) << "pxl::Object<" #type "> with id: " \
                                << id() << " is " << get() << std::endl; \
}

PRINT_NATIVE(int)
PRINT_NATIVE(unsigned int)
PRINT_NATIVE(bool)
PRINT_NATIVE(double)
PRINT_NATIVE(float)

#undef PRINT_NATIVE

// EXPLICIT INSTANCIATION
template class Object<int>;
template class Object<unsigned int>;
template class Object<bool>;
template class Object<double>;
template class Object<float>;

} // namespace pxl


namespace pxl {

ObjectOwner::ObjectOwner() :
    pxl::Vector<pxl::ObjectBase*>(),
    _copyHistory(),
    _index()
{
}

ObjectOwner::ObjectOwner(const pxl::ObjectOwner& original) :
    pxl::Vector<pxl::ObjectBase*>(),
    _copyHistory(),
    _index()
{
    // copy objects: loop in STL style
    for(StlConstIterator iter = original._container.begin();
        iter != original._container.end(); iter++) {

        pxl::ObjectBase* pOld = *iter;
        pxl::ObjectBase* pNew = pOld->clone();

        set(*pNew);
        _copyHistory.set(pOld->id(), pNew);
    }

    // FIXME: possibly, inefficient, might be done all in one loop
    // redirect relations: loop in PTL style
    for(Iterator iter(original); !iter.isDone(); iter.next()) {
        pxl::ObjectBase* pOld = iter.item();
        pxl::ObjectBase* pNew = _copyHistory.find(pOld->id(), 0);

        // mother relations
        for(pxl::Relations::Iterator iter(pOld->getMotherRelations());
            !iter.isDone(); iter.next()) {

            pxl::ObjectBase* pOldRel = iter.item()->pointer();
            pxl::ObjectBase* pNewRel = _copyHistory.find(pOldRel->id(), 0);

            if (pOldRel) {
                if (pNewRel)
                    pNew->linkMother(*pNewRel);
                else
                    // FIXME: cerr again?
                    std::cerr << "pxl::ObjectOwner::ObjectOwner(...): WARNING: some original objects had relations to objects of other owners." << std::endl;
            } else
                std::cerr << "pxl::ObjectOwner::ObjectOwner(...): WARNING: some originally related objects no longer exist." << std::endl;
        }

        // daughter relations
        // have been set automatically above
    }

    // redirect index:
    for(pxl::Index::Iterator iter(original._index);
        !iter.isDone(); iter.next()) {

        pxl::ObjectBase* pOld = iter.item();
        pxl::ObjectBase* pNew = _copyHistory.find(pOld->id(), 0);

        if (pNew)
            _index.set(iter.key(), pNew);
        else
            std::cerr << "pxl::ObjectOwner::ObjectOwner(...): WARNING: some original indices pointed to objects of other owners." << std::endl;
    }
}

void ObjectOwner::clearContainer() {
    for(StlConstIterator iter = _container.begin();
        iter != _container.end(); iter++)
        delete (*iter);
    _container.clear();
    _copyHistory.clearContainer();
    _index.clearContainer();
}

void ObjectOwner::set(pxl::ObjectBase& item)
{
    item._refObjectOwner = this;
    _container.push_back(&item);
}

void ObjectOwner::remove(pxl::ObjectBase& item)
{
    // search & remove possible indices (multiple occurrences possible!)
    for(pxl::Index::StlConstIterator iter = _index.begin();
        iter != _index.end(); iter++) {

        // FIXME: inefficient
        if (&item == iter->second)
            _index.remove(iter->first);
    }

    // search & remove possible copy history
    for(pxl::CopyHistory::StlConstIterator iter = _copyHistory.begin();
        iter != _copyHistory.end(); iter++) {

        // FIXME: inefficient
        if (&item == iter->second) {
            _copyHistory.remove(iter->first);
            break; // multiple occurrences *not* possible!
        }
    }

    // remove all relations:
    item.unlinkMothers();
    item.unlinkDaughters();

    for(StlIterator iter = _container.begin();
        iter != _container.end(); iter++) {

        if (&item == (*iter)) {
            delete *iter;
            _container.erase(iter);
            break;
        }
    }
}

bool ObjectOwner::has(const pxl::ObjectBase& item) const
{
    return item._refObjectOwner == this;
}

// EXPLICIT INSTANCIATION
template class ObjectOwner::TypeIterator<pxl::Object<int> >;

} // namespace pxl

namespace pxl {

// EXPLICIT INSTANCIATION
template class Ptr<int>;

} // namespace pxl

namespace pxl {

// EXPLICIT INSTANCIATION
template class pxl::Relations::TypeIterator<pxl::WkPtr<int> >;

} // namespace pxl


namespace pxl {

// EXPLICIT INSTANCIATION
template class SpyObject<int>;   

} // namespace pxl

namespace pxl {

// EXPLICIT INSTANCIATION
template class Vector<int>;

} // namespace pxl

namespace pxl {

void WkPtrBase::notifyDeleted()
{
    _objectRef = 0;
    if (_notifyChainOut)
        _notifyChainOut->notifyDeleted();
    _notifyChainIn = 0; 
    _notifyChainOut = 0;
}

void WkPtrBase::connect(pxl::ObjectBase* pointer)
{
    // disconnect:
    if (_objectRef) {
        if (_objectRef->_refWkPtrSpec == this)
            _objectRef->_refWkPtrSpec = _notifyChainOut;
        if (_notifyChainIn && _notifyChainOut) {
            _notifyChainIn->_notifyChainOut = _notifyChainOut;
            _notifyChainOut->_notifyChainIn = _notifyChainIn;
        } else {
            if (_notifyChainIn)
                _notifyChainIn->_notifyChainOut = 0;
            if (_notifyChainOut)
                _notifyChainOut->_notifyChainIn = 0;
        }
    }   
    _notifyChainOut = 0;
    _notifyChainIn = 0; 

    // connect
    if (pointer) {
        _notifyChainIn = 0;
        _notifyChainOut = pointer->_refWkPtrSpec;
        if (_notifyChainOut)
            _notifyChainOut->_notifyChainIn = this;
        pointer->_refWkPtrSpec = this;
    }

    _objectRef = pointer;
}

} // namespace pxl

namespace pxl {

// EXPLICIT INSTANCIATION
template class WkPtr<int>;

} // namespace pxl

namespace pxl {

// EXPLICIT INSTANCIATION
template class WkPtrSpec<int, pxl::Object<int> >;

} // namespace pxl

namespace pxl {

// EXPLICIT INSTANCIATION
template class pxl::iDiskFileVx<pxl::iStreamer>;

} // namespace pxl
#include <string>
#include <cstring>

#include <zlib.h>




#define iotl__iStreamer__lengthUnzipBuffer	65536

namespace pxl {

// dummy class
class Orphan {};

iStreamer::iStreamer() :
    pxl::BasicIoStreamer(), _buffer()  
{
    _inputBuffer = new unsigned char[iotl__iStreamer__lengthUnzipBuffer];
    _outputBuffer = new unsigned char[iotl__iStreamer__lengthUnzipBuffer];
}

iStreamer::~iStreamer()
{
    delete[] _inputBuffer;
    delete[] _outputBuffer;
}

int iStreamer::unzipEventData(std::istream &in, int nBytes)
{
    int ret, length = 0;
    unsigned int have;

    z_stream strm;
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL; 
    strm.opaque = Z_NULL;
    strm.avail_in = 0;   
    strm.next_in = Z_NULL;

    ret = inflateInit(&strm);
    if (ret != Z_OK)
        return 0;

    // decompress until deflate stream ends or end of file
    do {
        int size = nBytes;
        if (size > iotl__iStreamer__lengthUnzipBuffer)
            size = iotl__iStreamer__lengthUnzipBuffer;

        strm.avail_in = in.read((char*)_inputBuffer, size).gcount();
        if (in.bad()) {
            inflateEnd(&strm);
            return 0;
        }

        nBytes -= strm.avail_in;
        if (in.eof())
            nBytes = 0;

        if (strm.avail_in == 0)
            break;

        strm.next_in = _inputBuffer;

        // run inflate() on input until output buffer not full
        do {
            strm.avail_out = iotl__iStreamer__lengthUnzipBuffer;
            strm.next_out = _outputBuffer;

            ret = inflate(&strm, Z_NO_FLUSH);
            switch(ret) {
                case Z_STREAM_ERROR:
                    pxl::exception("pxl::unzipEventData()", "Internal inflate stream error.");
                case Z_NEED_DICT:
                    ret = Z_DATA_ERROR;	// fall through
                case Z_DATA_ERROR:
                case Z_MEM_ERROR:
                    inflateEnd(&strm);
                    return 0;
                default:
                    break;
            }

            have = iotl__iStreamer__lengthUnzipBuffer - strm.avail_out;
            _buffer.write((char *)_outputBuffer, have);

            length += have;
        } while(strm.avail_out == 0);
    } while(nBytes > 0); // done when inflate() says it's done

    inflateEnd(&strm);

    return length;
}

bool iStreamer::putEvent(std::istream& cxxx, int mode, const std::string& infoCondition)
{
    _buffer.clear();
    _buffer.str("");

    char eventMarker[4] = "   ";
    std::string info;
    char compressionMode = 0;

    int lengthInfo = 0;
    int lengthInfoVerify = 0;

    int lengthZip = 0;
    int lengthZipVerify = 0;

    restoreMemory(cxxx, eventMarker, 4);	// event marker
    if (!cxxx.good())
        return false;
    if (std::strcmp(iotl__eventMarker, eventMarker) != 0)
        pxl::exception("pxl::iStreamer::getEvent()",
                       "Invalid event marker while reading event header.");

    restoreBasicTypeChar(cxxx, compressionMode); // compression mode

    restoreBasicTypeInt(cxxx, lengthInfo);       // forward jumper
    if (lengthInfo) {
        char* buffer = new char[lengthInfo];
        cxxx.read(buffer, lengthInfo);
        info.assign(buffer, lengthInfo);
        delete[] buffer;
    }
    restoreBasicTypeInt(cxxx, lengthInfoVerify); // backward jumper & consistency check

    if (lengthInfo != lengthInfoVerify)
        return false;

    if (mode == -1 && infoCondition != info)
        mode = 0;

    restoreBasicTypeInt(cxxx, lengthZip);        // forward jumper
    if (mode == 0)
        cxxx.ignore(lengthZip); // ignore data
    else if (compressionMode == ' ') {
        // just read the data:
        int size;
        for(int bytes = lengthZip; bytes > 0; bytes -= size) {
            size = bytes;
            if (size > iotl__iStreamer__lengthUnzipBuffer)
                size = iotl__iStreamer__lengthUnzipBuffer;

            restoreMemory(cxxx, (char*)_outputBuffer, size);  
            _buffer.write((char*)_outputBuffer, size);
        }
    } else // unzip the data:
        unzipEventData(cxxx, lengthZip);

    restoreBasicTypeInt(cxxx, lengthZipVerify); // backward jumper & consistency check

    _buffer.seekg(0);

    return lengthZip == lengthZipVerify;
}

bool iStreamer::previous(std::istream& cxxx)
{
    int pos = cxxx.tellg();	// FIXME: Yikes! offset_t or something
    int lengthInfoVerify;
    int lengthZipVerify;

    // read zip length information:
    pos -= 4;
    if (pos < 0)
        return false;
    cxxx.seekg(pos);

    restoreBasicTypeInt(cxxx, lengthZipVerify); // backward jumper & consistency check
    if (lengthZipVerify < 0)
        pxl::exception("pxl::iStreamer::previous()", "Invalid zip block length information.");

    pos -= lengthZipVerify;
    pos -= 4;

    // read info length information:
    pos -= 4;
    if (pos < 0)
        return false;
    cxxx.seekg(pos);

    restoreBasicTypeInt(cxxx, lengthInfoVerify); // backward jumper & consistency check
    if (lengthInfoVerify < 0)
        pxl::exception("pxl::iStreamer::previous()", "Invalid info block length information.");

    pos -= lengthInfoVerify;
    pos -= 4;

    pos -= 1; // compression mode;
    pos -= 4; // event marker

    if (pos < 0) {
        pxl::exception("pxl::iStreamer::previous()", "Invalid block length information.");
        return false;
    }
    cxxx.seekg(pos);
    return true;
}

template<>
pxl::Id iStreamer::restoreData<pxl::Relations>(pxl::Relations& relations)
{
    relations.clearContainer();

    int size = 0;
    restoreData<int>(size);
    for(int i = 0; i < size; i++) {
        pxl::MutableId persistentId = 0; 

        restoreId(persistentId);
        // provide orphan relation
        if (persistentId != 0)
            relations.set(persistentId, new pxl::WkPtr<pxl::Orphan>);
    }

    return 0;
}

template<>
pxl::Id iStreamer::restoreData<pxl::ObjectOwner>(pxl::ObjectOwner& objects)
{
    objects.clearContainer();

    int size;
    std::string key;
    pxl::MutableId persistentId = 0;

    // reload objects (with orphan relations)
    size = 0;
    restoreData<int>(size);
    for(int i = 0; i < size; i++) {
        pxl::ObjectBase* item = 0;
        persistentId = restoreAbstractObject(&item);
        objects.set(*item);
        objects._copyHistory.set(persistentId, item);
    }

    // redirect relations: loop in PTL style
    for(pxl::ObjectOwner::Iterator iter(objects);
        !iter.isDone(); iter.next()) {
        pxl::ObjectBase* pNew = iter.item();

        // mother relations
        for(pxl::Relations::TypeIterator<pxl::WkPtr<pxl::Orphan> > iterRelations(pNew->_motherRelations);
            !iterRelations.isDone(); iterRelations.next()) {

            pxl::ObjectBase* pNewRel = objects._copyHistory.find(iterRelations.key(), 0);
            if (pNewRel)
                pNew->linkMother(*pNewRel);
            else  // FIXME:: std::cerr again?
                std::cerr << "pxl::iStreamer::restoreData<pxl::ObjectOwner>(...): WARNING: some original objects seem have invalid relations."
                          << " (persistent id is " << iterRelations.key() << ")" << std::endl;
        }

        // daughter relations have been set automatically above
        // remove orphan relations

        // mother relations
        pxl::Relations::TypeIterator<pxl::WkPtr<pxl::Orphan> > iterM(pNew->_motherRelations); 
        while(!iterM.isDone()) {
            pxl::Id id = iterM.key();
            iterM.next();
            pNew->_motherRelations.remove(id);
        }

        // daughter relations
        pxl::Relations::TypeIterator<pxl::WkPtr<pxl::Orphan> > iterD(pNew->_daughterRelations);
        while(!iterD.isDone()) {
            pxl::Id id = iterD.key();
            iterD.next();
            pNew->_daughterRelations.remove(id);
        }
    }

    // reload & redirect index
    size = 0;
    restoreData<int>(size);
    for(int i = 0; i < size; i++) {
        restoreData<std::string>(key);
        restoreId(persistentId);

        pxl::ObjectBase* pNew = objects._copyHistory.find(persistentId, 0);

        if (pNew) {
            objects._index.set(key, pNew);
        } else {
            // FIXME: my standard cerr complaint
            std::cerr << "pxl::iStreamer::restoreData<pxl::Objects>(...): WARNING: some original indices could not be reconstructed." << std::endl;
            std::cerr << key << ": " << persistentId << std::endl;
        }
    }

    return 0;
}

template<>
pxl::Id iStreamer::restoreData<pxl::Layout>(pxl::Layout& layout)
{
    restoreData<int>(layout._type);
    restoreData<int>(layout._style);
    restoreData<int>(layout._color);
    restoreData<double>(layout._a);
    restoreData<double>(layout._b);
    restoreData<double>(layout._c);
    restoreData<double>(layout._d);
    return 0;
}

template<>
pxl::Id iStreamer::restoreData<pxl::ObjectBase>(pxl::ObjectBase& object)
{
    pxl::MutableId persistentId ;
    bool hasLayout;

    restoreId(persistentId);

    restoreData<pxl::Relations>(object._motherRelations);
    restoreData<pxl::Relations>(object._daughterRelations);

    restoreData<bool>(hasLayout);
    if (hasLayout) {
        restoreData<pxl::Layout>(object.layout());
    } else {
        delete object._ptrLayout;
        object._ptrLayout = 0;
    }

    return persistentId;
}

pxl::Id iStreamer::restoreAbstractObject(pxl::ObjectBase** ppobj)
{
    std::string objectTypeId; 
    std::string dataTypeId;

    pxl::BasicIoStreamer::restoreBasicTypeCStr(_buffer, objectTypeId);   
    pxl::BasicIoStreamer::restoreBasicTypeCStr(_buffer, dataTypeId); 

    return pxl::TypeManager::instance().restoreObject(*this, ppobj, objectTypeId, dataTypeId);
}

} // namespace pxl

namespace pxl {

// EXPLICIT INSTANCIATION
template class pxl::oDiskFileVx<pxl::oStreamer>;

} // namespace pxl





namespace pxl {

void oStreamer::getEvent(std::ostream& cxxx, const std::string& info, char compressionMode)
{
    // info block:
    const char* cInfo = info.c_str();
    int lengthInfo = info.length();

    // zip block:
    const std::string& strBuffer = _buffer.str();

    const char* cBuffer = strBuffer.c_str();
    int lengthBuffer = strBuffer.length();

    const char* cZip = cBuffer;
    int lengthZip = lengthBuffer;

    char* cZipSpace = 0;
    unsigned long lengthZipSpace = 0;

    if (compressionMode == ' ') {
        // no compression requires no action...
    } else if (compressionMode >= '0' && compressionMode <= '9') {
        // data compression a la Gero, i.e. compression level = 6:
        lengthZipSpace = int(double(lengthBuffer) * 1.05 + 16);
        cZipSpace = new char[lengthZipSpace];

        int status = compress2((Bytef*)cZipSpace,
                               (uLongf*)&lengthZipSpace,
                               (const Bytef*)cBuffer,
                               lengthBuffer,
                               compressionMode - '0');
        switch(status) {
          case Z_MEM_ERROR:
            pxl::exception("pxl::oStreamer::getEvent()", "zlib: not enough memory");
            break;
          case Z_BUF_ERROR:
            pxl::exception("pxl::oStreamer::getEvent()", "zlib: buffer too small");
            break;
          case Z_STREAM_ERROR:
            pxl::exception("pxl::oStreamer::getEvent()", "zlib: level parameter invalid");
            break;
          default:
            break;
        }

        cZip = cZipSpace;
        lengthZip = lengthZipSpace;
    } else
        pxl::exception("pxl::oStreamer::getEvent()","Invalid compression mode.");

    // write event
    dumpMemory(cxxx, iotl__eventMarker, 4);       // event marker
    storeBasicTypeChar(cxxx, compressionMode);    // compression mode

    storeBasicTypeInt(cxxx, lengthInfo);          // forward jumper
    cxxx.rdbuf()->sputn(cInfo, lengthInfo);       // info block
    cxxx.rdbuf()->pubsync();
    storeBasicTypeInt(cxxx, lengthInfo);          // backward jumper & consistency check

    storeBasicTypeInt(cxxx, lengthZip);           // forward jumper
    cxxx.rdbuf()->sputn(cZip, lengthZip);         // compressed data block
    cxxx.rdbuf()->pubsync();
    storeBasicTypeInt(cxxx, lengthZip);           // backward jumper & consistency check

    if (cZipSpace)
        delete[] cZipSpace;

    _buffer.str("");
    _buffer.clear();
}

template<>
void oStreamer::storeData<pxl::Relations>(const pxl::Relations& relations)
{
    storeData<int>(relations.getSize());
    for(pxl::Relations::StlConstIterator iter = relations.begin();
        iter != relations.end(); iter++) {

        const pxl::WkPtrBase& ptr = (*iter->second);
        pxl::Id id = ptr.valid() ? ptr.pointer()->id() : 0;

        storeId(id);
    }
}

template<>
void oStreamer::storeData<pxl::ObjectOwner>(const pxl::ObjectOwner& objects)
{
    // objects
    storeData<int>(objects.getSize());
    for(pxl::ObjectOwner::StlConstIterator iter = objects.begin(); 
        iter != objects.end(); iter++)
        storeAbstractObject(**iter);

    // index
    storeData<int>(objects._index.getSize());
    for(pxl::Index::StlConstIterator iter = objects._index.begin(); 
        iter != objects._index.end(); iter++) {

        storeData<std::string>(iter->first);
        storeId(iter->second->id());
    }
}

template<>
void oStreamer::storeData<pxl::Layout>(const pxl::Layout& layout)
{
    storeData<int>(layout._type);
    storeData<int>(layout._style);
    storeData<int>(layout._color);
    storeData<double>(layout._a);
    storeData<double>(layout._b);
    storeData<double>(layout._c);
    storeData<double>(layout._d);
}

template<>
void oStreamer::storeData<pxl::ObjectBase>(const pxl::ObjectBase& object)
{
    pxl::Id id = object.id();
    storeId(id);

    storeData<pxl::Relations>(object._motherRelations);
    storeData<pxl::Relations>(object._daughterRelations);

    if (object._ptrLayout) {
        storeData<bool>(true);
        storeData<pxl::Layout>(*object._ptrLayout);
    } else
        storeData<bool>(false);
}

} // namespace pxl


namespace pxl {

// EXPLICIT INSTANCIATION
template class TypeAgent<pxl::Object<int> >;
template class TypeAgent<pxl::CowObject<int> >;

} // namespace pxl




namespace pxl {

// The singleton instance
TypeManager* TypeManager::_instance = 0;

TypeManager& TypeManager::instance()
{
    if (!_instance)
        _instance = new TypeManager;

    return *_instance;    
}

void pxl::TypeManager::registerAgent(pxl::TypeAgentBase* agent) {
// for debugging purposes:
//    std::cout << "pxl::TypeManager::registerAgent(): type registered: [" 
//              << agent->getObjectTypeId() << "," << agent->getDataTypeId() 
//              << "]" << std::endl;
//              //<< "] = '" << agent->getCppTypeId() << "'" << std::endl;

    // add entries to maps:
    _agentsByIotlTypeId.set(pxl::TypeIdKey(agent->getObjectTypeId(), agent->getDataTypeId()), agent);
    _agentsByCppTypeId.set(agent->getCppTypeId(), agent);
}

pxl::Id TypeManager::restoreObject(iStreamer& input, pxl::ObjectBase** ppobj, const std::string& objectTypeId, const std::string& dataTypeId) const
{
    TypeAgentBase* agent = _agentsByIotlTypeId.find(pxl::TypeIdKey(objectTypeId, dataTypeId), 0);
    if (!agent) {
        std::string message; 
        message += "undeclared IOTL type: [";
        message += objectTypeId;
        message += ", "; 
        message += dataTypeId; 
        message += "]"; 
        pxl::exception("pxl::TypeManager::restoreObject()", message);

        return 0;
    }

    return agent->restoreObject(input, ppobj);
}

pxl::Id TypeManager::restoreObject(pxl::iStreamer& input, pxl::ObjectBase& obj, const std::string& cppTypeId) const
{
    TypeAgentBase* agent = _agentsByCppTypeId.find(cppTypeId, 0);
    if (!agent) {
        std::string message; 
        message += "undeclared C++ type: '";
        message += cppTypeId;
        message += "'";
        pxl::exception("pxl::TypeManager::restoreObject()", message); 

        return 0;
    }

    return agent->restoreObject(input, obj);
}

void TypeManager::storeObject(pxl::oStreamer& output, const pxl::ObjectBase& obj, const std::string& cppTypeId) const
{
    TypeAgentBase* agent = _agentsByCppTypeId.find(cppTypeId, 0);
    if (!agent) {
        std::string message; 
        message += "undeclared C++ type: '";
        message += cppTypeId;
        message += "'";
        pxl::exception("pxl::TypeManager::storeObject()", message);
 
        return;
    }

    agent->storeObject(output, obj);
}

} // namespace pxl


iotl__declareDataTypeExplicit(char,        "\1c",      __char, data, { pxl::BasicIoStreamer::storeBasicTypeChar(_buffer, data); },   { pxl::BasicIoStreamer::restoreBasicTypeChar(_buffer, data); })
iotl__declareDataTypeExplicit(std::string, "\1s", std__string, data, { pxl::BasicIoStreamer::storeBasicTypeString(_buffer, data); }, { pxl::BasicIoStreamer::restoreBasicTypeString(_buffer, data); })
iotl__declareDataTypeExplicit(bool,        "\1b",      __bool, data, { pxl::BasicIoStreamer::storeBasicTypeBool(_buffer, data); },   { pxl::BasicIoStreamer::restoreBasicTypeBool(_buffer, data); })
iotl__declareDataTypeExplicit(int,         "\1i",       __int, data, { pxl::BasicIoStreamer::storeBasicTypeInt(_buffer, data); },    { pxl::BasicIoStreamer::restoreBasicTypeInt(_buffer, data); })
iotl__declareDataTypeExplicit(float,       "\1f",     __float, data, { pxl::BasicIoStreamer::storeBasicTypeFloat(_buffer, data); },  { pxl::BasicIoStreamer::restoreBasicTypeFloat(_buffer, data); })
iotl__declareDataTypeExplicit(double,      "\1d",    __double, data, { pxl::BasicIoStreamer::storeBasicTypeDouble(_buffer, data); }, { pxl::BasicIoStreamer::restoreBasicTypeDouble(_buffer, data); })




namespace pxl {

void AnalysisFork::beginJob(const pxl::Objects* input)
{
    for(pxl::Objects::TypeIterator<pxl::AnalysisFork> iter(get().getObjects());
       !iter.isDone(); iter.next())
        iter.object().beginJob(input);

    for(pxl::Objects::TypeIterator<pxl::AnalysisProcess> iter(get().getObjects());
        !iter.isDone(); iter.next())
        iter.object().beginJob(input);
}

void AnalysisFork::buildTemplate(int mode)
{
    for(pxl::Objects::TypeIterator<pxl::AnalysisFork> iter(get().getObjects());
        !iter.isDone(); iter.next())
        iter.object().buildTemplate(mode);

    for(pxl::Objects::TypeIterator<pxl::AnalysisProcess> iter(get().getObjects());
        !iter.isDone(); iter.next())
        iter.object().buildTemplate(mode);
}

void pxl::AnalysisFork::beginRun(const pxl::Objects* input)
{
    for(pxl::Objects::TypeIterator<pxl::AnalysisFork> iter(get().getObjects());
        !iter.isDone(); iter.next())
        iter.object().beginRun(input);

    for(pxl::Objects::TypeIterator<pxl::AnalysisProcess> iter(get().getObjects());
        !iter.isDone(); iter.next())
        iter.object().beginRun(input);
}

void pxl::AnalysisFork::analyseEvent(const pxl::Objects* input)
{
    for(pxl::Objects::TypeIterator<pxl::AnalysisFork> iter(get().getObjects());
        !iter.isDone(); iter.next())
        iter.object().analyseEvent(input);
      
    for(pxl::Objects::TypeIterator<pxl::AnalysisProcess> iter(get().getObjects());
        !iter.isDone(); iter.next())
        iter.object().analyseEvent(input);
}

void pxl::AnalysisFork::finishEvent(const pxl::Objects* input)
{
    for(pxl::Objects::TypeIterator<pxl::AnalysisFork> iter(get().getObjects());
        !iter.isDone(); iter.next())
        iter.object().finishEvent(input);

    for(pxl::Objects::TypeIterator<pxl::AnalysisProcess> iter(get().getObjects());
        !iter.isDone(); iter.next())
        iter.object().finishEvent(input);
}

void pxl::AnalysisFork::endRun(const pxl::Objects* input)
{
    for(pxl::Objects::TypeIterator<pxl::AnalysisFork> iter(get().getObjects());
        !iter.isDone(); iter.next())
        iter.object().endRun(input);

    for(pxl::Objects::TypeIterator<pxl::AnalysisProcess> iter(get().getObjects());
        !iter.isDone(); iter.next())
        iter.object().endRun(input);
}

void pxl::AnalysisFork::endJob(const pxl::Objects* input)
{
    for(pxl::Objects::TypeIterator<pxl::AnalysisFork> iter(get().getObjects());
        !iter.isDone(); iter.next())
        iter.object().endJob(input);

    for(pxl::Objects::TypeIterator<pxl::AnalysisProcess> iter(get().getObjects());
        !iter.isDone(); iter.next())
        iter.object().endJob(input);
}

pxl::ObjectBase* AnalysisFork::clone() const
{
    return new AnalysisFork(*this);
}

pxl::WkPtrBase* AnalysisFork::createSelfWkPtr()
{
    return new AnalysisForkWkPtr(*this);
}

void AnalysisFork::storeYourSelf(pxl::oStreamer& output) const
{
    output.storeObject(*this);
}

template<>
std::ostream& pxl::Object<pxl::AnalysisForkData>::print(int level, std::ostream& os, int pan) const
{
//    printPan1st(os, pan) << "pxl::CowObject<pxl::AnalysisForkData> with id: " << id() << std::endl;
//    printPan(os, pan)    << "     name: " << get().getName() << std::endl;
    printPan1st(os, pan) << "AnalysisFork: " << get().getName() << std::endl;
    for(pxl::Objects::Iterator iter(get().getObjects());
        !iter.isDone(); iter.next()) {

        if (iter.object().getMotherRelations().getSize() == 0)
            iter.object().printDecayTree(level, os, pan);
    }
    return os;
}

} // namespace pxl

iotl__declareObjectTypeExplicit(pxl::AnalysisFork, "\3Af", pol__AnalysisFork)

iotl__declareDataTypeExplicit(pxl::AnalysisForkData, "\3AfD", pol__AnalysisForkData, data, {
    storeData<pxl::BasicObjectManagerData>(data);
}, {
    restoreData<pxl::BasicObjectManagerData>(data);
})




namespace pxl {

pxl::ObjectBase* AnalysisProcess::clone() const
{
    return new AnalysisProcess(*this);
}

pxl::WkPtrBase* AnalysisProcess::createSelfWkPtr()
{
    return new AnalysisProcessWkPtr(*this);
}

void AnalysisProcess::storeYourSelf(pxl::oStreamer& output) const
{
    output.storeObject(*this);
}

template<>
std::ostream& pxl::Object<pxl::AnalysisProcessData>::print(int level, std::ostream& os, int pan) const
{
//    printPan1st(os, pan) << "pxl::CowObject<pxl::AnalysisProcessData> with id: " << id() << std::endl;
//    printPan(os, pan)    << "     name: " << get().getName() << std::endl;
    printPan1st(os, pan) << "AnalysisProcess: " << get().getName() << std::endl;
    for(pxl::Objects::Iterator iter(get().getObjects());
        !iter.isDone(); iter.next()) {

        if (iter.object().getMotherRelations().getSize() == 0)
            iter.object().printDecayTree(level, os, pan);
    }
    return os;
}

} // namespace pxl

iotl__declareObjectTypeExplicit(pxl::AnalysisProcess, "\3Ap", pol__AnalysisProcess)

iotl__declareDataTypeExplicit(pxl::AnalysisProcessData, "\3ApD", pol__AnalysisProcessData, data, {
    storeData<pxl::BasicObjectManagerData>(data);
}, {
    restoreData<pxl::BasicObjectManagerData>(data);
})



namespace pxl {

bool const operator==(const pxl::Basic3VectorData& obj1, const pxl::Basic3VectorData& obj2) 
{
    return obj1.getX() == obj2.getX() &&
           obj1.getY() == obj2.getY() &&
           obj1.getZ() == obj2.getZ();
}

bool const operator!=(const pxl::Basic3VectorData& obj1, const pxl::Basic3VectorData& obj2) 
{
    return obj1.getX() != obj2.getX() ||
           obj1.getY() != obj2.getY() ||
           obj1.getZ() != obj2.getZ();
}

} // namespace pxl

iotl__declareDataTypeExplicit(pxl::Basic3VectorData, "\3B3vD", pol__Basic3VectorData, data, {
    storeData<double>(data._x);
    storeData<double>(data._y);
    storeData<double>(data._z);
}, {
    restoreData<double>(data._x);
    restoreData<double>(data._y);
    restoreData<double>(data._z);
})



namespace pxl {

bool const operator==(const pxl::Basic4VectorData& obj1, const pxl::Basic4VectorData& obj2)
{
    return obj1.getX() == obj2.getX() &&
           obj1.getY() == obj2.getY() &&
           obj1.getZ() == obj2.getZ() &&
           obj1.getE() == obj2.getE();
}

bool const operator!=(const pxl::Basic4VectorData& obj1, const pxl::Basic4VectorData& obj2)
{
    return obj1.getX() != obj2.getX() ||
           obj1.getY() != obj2.getY() ||
           obj1.getZ() != obj2.getZ() ||
           obj1.getE() != obj2.getE();
}

} // namespace pxl

iotl__declareDataTypeExplicit(pxl::Basic4VectorData, "\3B4vD", pol__Basic4VectorData, data, {
    storeData<double>(data._x);
    storeData<double>(data._y);
    storeData<double>(data._z);
    storeData<double>(data._t);
}, {
    restoreData<double>(data._x);
    restoreData<double>(data._y);
    restoreData<double>(data._z);
    restoreData<double>(data._t);
})



iotl__declareObjectTypeExplicit(pxl::BasicObject, "\3Bo", pol__BasicObject)

iotl__declareDataTypeExplicit(pxl::BasicObjectData, "\3BoD", pol__BasicObjectData, data, {
    storeData<bool>(data._locked);
    storeData<int>(data._monteCarloMode);
    storeData<std::string>(data._name);
    storeData<int>(data._status);
    storeData<int>(data._workflag);
    storeData(data._userRecords);
    // storeData<pxl::CppPointers>(_cppPointers); // we don't store the c++ pointers!
}, {
    restoreData<bool>(data._locked);
    restoreData<int>(data._monteCarloMode);
    restoreData<std::string>(data._name);
    restoreData<int>(data._status);
    restoreData<int>(data._workflag);
    restoreData(data._userRecords);
    // restoreData<pxl::CppPointers>(_cppPointers); // we don't restore the c++ pointers!
})



iotl__declareObjectTypeExplicit(pxl::BasicObjectManager, "\3Bom", pol__BasicObjectManager)

iotl__declareDataTypeExplicit(pxl::BasicObjectManagerData, "\3BomD", pol__BasicObjectManagerData, data, {
    storeData<pxl::Objects>(data._objects);
    storeData<pxl::BasicObjectData>(data);
}, {
    restoreData<pxl::Objects>(data._objects);
    restoreData<pxl::BasicObjectData>(data);
})



namespace pxl {

template<>
std::ostream& pxl::CowObject<pxl::CollisionData>::print(int level, std::ostream& os, int pan) const
{
//    printPan1st(os, pan) << "pxl::CowObject<pxl::CollisionData> with id: " << id() << " (data socket currently located at " << _dataSocket << ")" << std::endl;
//    printPan(os, pan)    << "     name: " << get().getName() << std::endl;
    printPan1st(os, pan) << "Collision: " << get().getName() << std::endl;
    return os;
}

} // namespace pxl




namespace pxl {

template<>
std::ostream& pxl::Object<pxl::EventViewData>::print(int level, std::ostream& os, int pan) const
{
//    printPan1st(os, pan) << "pxl::CowObject<pxl::EventViewData> with id: " << id() << std::endl;
//    printPan(os, pan)    << "     name: " << get().getName() << std::endl;
    printPan1st(os, pan) << "EventView: " << get().getName() << std::endl;
    for(pxl::Objects::Iterator iter(get().getObjects());
        !iter.isDone(); iter.next()) {

        if (iter.object().getMotherRelations().getSize() == 0)
            iter.object().printDecayTree(level, os, pan);
    }
    return os;
}

} // namespace pxl

iotl__declareObjectTypeExplicit(pxl::EventView, "\3Ev", pol__EventView)

iotl__declareDataTypeExplicit(pxl::EventViewData, "\3EvD", pol__EventViewData, data, {
    storeData<pxl::BasicObjectManagerData>(data);
}, {
    restoreData<pxl::BasicObjectManagerData>(data);
})




namespace pxl {

bool const operator==(const pxl::ParticleData& obj1, const pxl::ParticleData& obj2) 
{
    return obj1.vector() == obj2.vector() &&
           obj1.getCharge() == obj2.getCharge();
}

bool const operator!=(const pxl::ParticleData& obj1, const pxl::ParticleData& obj2) 
{
    return obj1.vector() != obj2.vector() ||
           obj1.getCharge() != obj2.getCharge();
}

template<>
std::ostream& pxl::CowObject<pxl::ParticleData>::print(int level, std::ostream& os, int pan) const
{
//    printPan1st(os, pan) << "pxl::CowObject<pxl::ParticleData> with id: " << id() << " (data socket currently located at " << _dataSocket << ")" << std::endl;
//    printPan(os, pan)    << "     name: " << get().getName() << std::endl;
    printPan1st(os, pan) << "Particle: '" << get().getName() << "', p = ("
                         << get().vector().getPt() << ", " 
                         << get().vector().getPz() << ") m = " 
                         << get().vector().getMass() << std::endl;
    return os;
}

} // namespace pxl

iotl__declareObjectTypeExplicit(pxl::Particle, "\3Pa", pol__Particle)

iotl__declareDataTypeExplicit(pxl::ParticleData, "\3PaD", pol__ParticleData, data, {
    storeData<double>(data._charge);
    storeData<int>(data._particleId);
    storeData(data.vector());
    storeData<pxl::BasicObjectData>(data);
}, {
    restoreData<double>(data._charge);
    restoreData<int>(data._particleId);
    restoreData(data.vector(pxl::set));
    restoreData<pxl::BasicObjectData>(data);
})




namespace pxl {

bool const operator==(const pxl::VertexData& obj1, const pxl::VertexData& obj2)
{
    return obj1.vector() == obj2.vector();
}

bool const operator!=(const pxl::VertexData& obj1, const pxl::VertexData& obj2) 
{
    return obj1.vector() != obj2.vector();
}

template<>
std::ostream& pxl::CowObject<pxl::VertexData>::print(int level, std::ostream& os, int pan) const
{
//    printPan1st(os, pan) << "pxl::CowObject<pxl::VertexData> with id: " << id() << " (data socket currently located at " << _dataSocket << ")" << std::endl;
//    printPan(os, pan)    << "     name: " << get().getName() << std::endl;
    printPan1st(os, pan) << "Vertex: '" << get().getName() << "', x = ("
                         << get().vector().getX() << ", "
                         << get().vector().getY() << ", "
                         << get().vector().getZ() << ")" << std::endl;
    return os;
}

} // namespace pxl

iotl__declareObjectTypeExplicit(pxl::Vertex, "\3Vx", pol__Vertex)

iotl__declareDataTypeExplicit(pxl::VertexData, "\3VxD", pol__VertexData, data, {
    storeData(data.vector());
    storeData<pxl::BasicObjectData>(data);
}, {
    restoreData(data.vector(pxl::set));
    restoreData<pxl::BasicObjectData>(data);
})
//
//
//----------------------------------------------------------------------
// ePax
//----------------------------------------------------------------------
const pxl::Get ePaxGet;
const pxl::Set ePaxSet;
//----------------------------------------------------------------------
iotl__declareObjectTypeExplicit(ePaxParticle,  "\4ePa", ePaxParticle)
iotl__declareObjectTypeExplicit(ePaxVertex,    "\4eVe", ePaxVertex)
iotl__declareObjectTypeExplicit(ePaxCollision, "\4eCo", ePaxCollision)
iotl__declareObjectTypeExplicit(ePaxEventView, "\4eEv", ePaxEventView)
iotl__declareObjectTypeExplicit(ePaxAnalysisProcess, "\4eAp", ePaxAnalysisProcess)
iotl__declareObjectTypeExplicit(ePaxAnalysisFork,    "\4eAf", ePaxAnalysisFork)
//----------------------------------------------------------------------
pxl::ObjectBase* ePaxParticle::clone() const { return new ePaxParticle(*this); }
std::ostream& ePaxParticle::print(int level, std::ostream& os, int pan) const { return pxl::CowObject<pxl::ParticleData>::print(level, os, pan); }
pxl::WkPtrBase* ePaxParticle::createSelfWkPtr() { return new ePaxParticleWkPtr(*this); }
void ePaxParticle::storeYourSelf(pxl::oStreamer& output) const { output.storeObject(*this); }
//----------------------------------------------------------------------
pxl::ObjectBase* ePaxVertex::clone() const { return new ePaxVertex(*this); }
std::ostream& ePaxVertex::print(int level, std::ostream& os, int pan) const { return pxl::CowObject<pxl::VertexData>::print(level, os, pan); }
pxl::WkPtrBase* ePaxVertex::createSelfWkPtr() { return new ePaxVertexWkPtr(*this); }
void ePaxVertex::storeYourSelf(pxl::oStreamer& output) const { output.storeObject(*this); }
//----------------------------------------------------------------------
pxl::ObjectBase* ePaxCollision::clone() const { return new ePaxCollision(*this); }
std::ostream& ePaxCollision::print(int level, std::ostream& os, int pan) const { return pxl::CowObject<pxl::CollisionData>::print(level, os, pan); }
pxl::WkPtrBase* ePaxCollision::createSelfWkPtr() { return new ePaxCollisionWkPtr(*this); }
void ePaxCollision::storeYourSelf(pxl::oStreamer& output) const { output.storeObject(*this); }
//----------------------------------------------------------------------
pxl::ObjectBase* ePaxEventView::clone() const { return new ePaxEventView(*this); }
std::ostream& ePaxEventView::print(int level, std::ostream& os, int pan) const { return pxl::Object<pxl::EventViewData>::print(level, os, pan); }
pxl::WkPtrBase* ePaxEventView::createSelfWkPtr() { return new ePaxEventViewWkPtr(*this); }
void ePaxEventView::storeYourSelf(pxl::oStreamer& output) const { output.storeObject(*this); }
//----------------------------------------------------------------------
pxl::ObjectBase* ePaxAnalysisProcess::clone() const { return new ePaxAnalysisProcess(*this); }
std::ostream& ePaxAnalysisProcess::print(int level, std::ostream& os, int pan) const { return pxl::Object<pxl::AnalysisProcessData>::print(level, os, pan); }
pxl::WkPtrBase* ePaxAnalysisProcess::createSelfWkPtr() { return new ePaxAnalysisProcessWkPtr(*this); }
void ePaxAnalysisProcess::storeYourSelf(pxl::oStreamer& output) const { output.storeObject(*this); }
//----------------------------------------------------------------------
pxl::ObjectBase* ePaxAnalysisFork::clone() const { return new ePaxAnalysisFork(*this); }
std::ostream& ePaxAnalysisFork::print(int level, std::ostream& os, int pan) const { return pxl::Object<pxl::AnalysisForkData>::print(level, os, pan); }
pxl::WkPtrBase* ePaxAnalysisFork::createSelfWkPtr() { return new ePaxAnalysisForkWkPtr(*this); }
void ePaxAnalysisFork::storeYourSelf(pxl::oStreamer& output) const { output.storeObject(*this); }
//----------------------------------------------------------------------

