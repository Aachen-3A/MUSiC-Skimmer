#ifndef epax_hh
#define epax_hh

#ifndef pxl_hh
#define pxl_hh

#ifndef pxl_pcl_hh
#define pxl_pcl_hh

#ifndef pxl_pcl_Exception_hh
#define pxl_pcl_Exception_hh

#include <string>

namespace pxl {

class Exception {
  public: 
    Exception() :
        routine("unspecified routine"),
        message("unspecified error") {}

    Exception(const std::string& routine, const std::string& message) :
        routine(routine),
        message(message) {}

    ~Exception() {} 

    const std::string &getRoutine() const { return routine; }
    const std::string &getMessage() const { return message; }

  private:
    std::string routine;
    std::string message;
};

} // namespace pxl

#endif // pxl_pcl_Exception_hh
#ifndef pxl_pcl_functions_hh
#define pxl_pcl_functions_hh


namespace pxl {

double getCpuTime();
void exception(const std::string& routine, const std::string& message);

} // namespace pxl

#endif // pxl_pcl_functions_hh
#ifndef pxl_pcl_BasicIoStreamer_hh
#define pxl_pcl_BasicIoStreamer_hh

#include <iostream>

#ifndef pxl_pcl_BasicLinuxIoStreamer_hh
#define pxl_pcl_BasicLinuxIoStreamer_hh


/*
 * FIXME / General comment:
 * BasicLinuxIoStreamer seems generic enough, not at all Linux specific?
 */

namespace pxl {

class BasicLinuxIoStreamer {
  protected:
    // writing data
    inline void dumpMemory(std::ostream& cxxx, const char* address, int bytes)
    {
        cxxx.rdbuf()->sputn(address, bytes);
        cxxx.rdbuf()->pubsync();
    }

    inline void storeBasicTypeChar(std::ostream& cxxx, char data)
    {
        dumpMemory(cxxx, (const char*)&data, 1); 
    }

    inline void storeBasicTypeBool(std::ostream& cxxx, bool data)
    {
        if (data) dumpMemory(cxxx, "Y", 1);
        else dumpMemory(cxxx, "N", 1);
    }

    inline void storeBasicTypeInt(std::ostream& cxxx, int data)
    {
        dumpMemory(cxxx, (const char*)&data, 4); // FIXME: endian/64bit
    }

    inline void storeBasicTypeFloat(std::ostream& cxxx, float data)
    {
        dumpMemory(cxxx, (const char*)&data, 4); // FIXME: endian/64bit
    }

    inline void storeBasicTypeDouble(std::ostream& cxxx, double data)
    {
        dumpMemory(cxxx, (const char*)&data, 8); // FIXME: endian/64bit
    }

    inline void storeBasicTypeCStr(std::ostream& cxxx, const char* address)
    {
        for(; (*address) != '\0'; address++) // FIXME: performance!!!
      	    cxxx.rdbuf()->sputc(*address);   // better: len, data ? (as below)
        cxxx.rdbuf()->sputc('\0');
    }

    inline void storeBasicTypeString(std::ostream& cxxx, const std::string& data)
    {
        int length = data.length();
        storeBasicTypeInt(cxxx, length);
        dumpMemory(cxxx, data.c_str(), length);
    }
  
    // reading data
    inline bool restoreMemory(std::istream& cxxx, char* address, int bytes)
    {
        cxxx.read(address, bytes);
        return bytes == cxxx.gcount();
    }

    inline bool restoreBasicTypeChar(std::istream& cxxx, char& data)
    {
        return restoreMemory(cxxx, (char*)&data, 1);
    }

    inline bool restoreBasicTypeBool(std::istream& cxxx, bool& data)
    {
      char cYesNo = ' ';
      bool success = restoreMemory(cxxx, &cYesNo, 1);
      data = (cYesNo == 'Y');
      return success;
    }

    inline bool restoreBasicTypeInt(std::istream& cxxx, int& data)
    {
        return restoreMemory(cxxx, (char*)&data, 4); // FIXME: endian/64bit
    }

    inline bool restoreBasicTypeFloat(std::istream& cxxx, float& data)
    {
        return restoreMemory(cxxx, (char*)&data, 4); // FIXME: endian/64bit
    }

    inline bool restoreBasicTypeDouble(std::istream& cxxx, double& data)
    {
        return restoreMemory(cxxx, (char*)&data, 8); // FIXME: endian/64bit
    }

    inline bool restoreBasicTypeCStr(std::istream& cxxx, std::string& data, int bytes = -1)
    {
        int read; // see comment above about persistent format
        data.clear();
        cxxx.get(_operationsBuffer, BUF_SIZE - 1, '\0');
        read = cxxx.gcount();
        _operationsBuffer[read + 1] = '\0';
        while (read == BUF_SIZE - 1) {
             data += _operationsBuffer;
             cxxx.get(_operationsBuffer, BUF_SIZE - 1, '\0');
             read = cxxx.gcount();
             _operationsBuffer[read+1] = '\0';
        }
        data += _operationsBuffer;
        cxxx.ignore(1);
        return true;
    }

    inline bool restoreBasicTypeString(std::istream& cxxx, std::string& data)
    {
        int length;
        restoreBasicTypeInt(cxxx, length); 
        if (length) {
            char *buffer = new char[length];	// FIXME: performance
            cxxx.read(buffer, length);
            data.assign(buffer, length);
            delete[] buffer;
         } else
            data = "";
        return 1;
    }

  private:
    static const unsigned int BUF_SIZE = 1024;

    char _operationsBuffer[BUF_SIZE];
};

} // namespace pxl

#endif // pxl_pcl_BasicLinuxIoStreamer_hh

namespace pxl {

typedef BasicLinuxIoStreamer BasicIoStreamer;

} // namespace pxl

#endif // pxl_pcl_BasicIoStreamer_hh

#endif // pxl_pcl_hh
#ifndef pxl_ptl_hh
#define pxl_ptl_hh

#ifndef pxl_ptl_MutableId_hh
#define pxl_ptl_MutableId_hh

namespace pxl {

typedef int MutableId;

} // namespace pxl

#endif // pxl_ptl_MutableId_hh
#ifndef pxl_ptl_Id_hh
#define pxl_ptl_Id_hh


namespace pxl {

typedef const pxl::MutableId Id;

} // namespace pxl

#endif // pxl_ptl_Id_hh
#ifndef pxl_ptl_Ptr_hh
#define pxl_ptl_Ptr_hh


namespace pxl {

// ptl

template<class datatype>
class Ptr {
  public:
    Ptr() : _cppPtr(0) {}
    Ptr(const pxl::Ptr<datatype>& original) : _cppPtr(original._cppPtr) {}
    Ptr(datatype* original) : _cppPtr(original) {}
    ~Ptr() {}

    inline datatype* pointer() const  { return _cppPtr; }
    inline bool valid() const         { return _cppPtr != 0; }

    inline datatype& object() const   { return *access(); }

    inline void operator=(const pxl::Ptr<datatype>& pptr) { _cppPtr = pptr._cppPtr; }
    inline void operator=(datatype& data)                 { _cppPtr = &data; }
    inline void operator=(datatype* dataptr)              { _cppPtr = dataptr; }

    inline datatype* operator->() const { return _cppPtr; }

  protected:
    // safe access to object
    inline datatype* access() const
    {
        if (_cppPtr)
            return _cppPtr; 
        std::cerr << "pxl::Ptr::access(): FATAL: The object you intend to access does not exist!" << std::endl;
        return 0;
    }

    datatype* _cppPtr;
};


template<class datatype>
datatype& operator*(pxl::Ptr<datatype>& ptr)
{ return ptr.object(); }

template<class datatype>
const datatype& operator*(const pxl::Ptr<datatype>& ptr)
{ return ptr.object(); }

} // namespace pxl

#endif // pxl_ptl_Ptr_hh
#ifndef pxl_ptl_Vector_hh
#define pxl_ptl_Vector_hh

#include <vector>

namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// ptl

template<class itemtype>
class Vector {
  public:
    virtual ~Vector() {}

    typedef          std::vector<itemtype>                 StlContainer;
    typedef typename std::vector<itemtype>::const_iterator StlConstIterator;
    typedef typename std::vector<itemtype>::iterator       StlIterator;

    // navigation
    const StlConstIterator begin() const { return _container.begin(); }
    const StlConstIterator end()   const { return _container.end(); }

    inline StlContainer& getContainer() { return _container; }
    virtual void clearContainer() { return _container.clear(); }

    inline int getSize() const { return _container.size(); }

    class Iterator {
      public:
        Iterator(const pxl::Vector<itemtype>::Iterator& original) :
            _iter(original._iter), _containerRef(original._containerRef) {}

        Iterator(const pxl::Vector<itemtype>& vector) :
            _containerRef(&vector) { first(); }

        inline void first()  { _iter = _containerRef->begin(); }
        inline void next()   { _iter++; }
        inline bool isDone() { return _iter == _containerRef->end(); }

	inline itemtype item() { return *_iter; }

      private:
        Iterator() : _containerRef(0) {}

        StlConstIterator             _iter;
        const pxl::Vector<itemtype>* _containerRef;
    };

  protected:
    StlContainer _container;

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;
};

} // namespace pxl

#endif // pxl_ptl_Vector_hh
#ifndef pxl_ptl_Map_hh
#define pxl_ptl_Map_hh

#include <map>


namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// pcl

template<class keytype, class itemtype>
class Map {
  public:
    virtual ~Map() {}

    typedef          std::pair<keytype, itemtype>                StlPair;
    typedef          std::map<keytype, itemtype>                 StlContainer;
    typedef typename std::map<keytype, itemtype>::const_iterator StlConstIterator;
    typedef typename std::map<keytype, itemtype>::iterator       StlIterator;

    void set(const keytype& key, const itemtype& item)
    {
        _container.erase(key);
        _container.insert(StlPair(key, item));
    }

    void remove(const keytype& key)
    { _container.erase(key); }

    itemtype find(const keytype& key, itemtype defaultitem) const
    {
        StlConstIterator found = _container.find(key);
        if (found != _container.end())
            return found->second;
        return defaultitem;
    }

    itemtype find(const keytype& key) const
    {
        StlConstIterator found = _container.find(key);
        if (found != _container.end())
            return found->second;
        pxl::exception("ptl::Map::find(...)", "key not found and no default item provided");
        throw;
    }

    // navigation
    const StlConstIterator begin() const { return _container.begin(); }
    const StlConstIterator end()   const { return _container.end(); }

    inline StlContainer& getContainer() { return _container; }
    virtual void clearContainer() { return _container.clear(); }

    inline int getSize() const { return _container.size(); }

    class Iterator {
      public:
        Iterator(const pxl::Map<keytype, itemtype>::Iterator& original) :
            _iter(original._iter), _containerRef(original._containerRef) {}

        Iterator(const pxl::Map<keytype, itemtype>& map) :
            _containerRef(&map) { first(); }

        inline void first()  { _iter = _containerRef->begin(); }
        inline void next()   { _iter++; }
        inline bool isDone() { return _iter == _containerRef->end(); }

        inline keytype  key()  { return _iter->first; }
	inline itemtype item() { return _iter->second; }

      private:
        Iterator() : _containerRef(0) {}

        StlConstIterator                   _iter;
        const pxl::Map<keytype, itemtype>* _containerRef;
    };

  protected:
    StlContainer _container;

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;
};

} // namespace pxl

#endif // pxl_ptl_Map_hh
#ifndef pxl_ptl_CopyHistory_hh
#define pxl_ptl_CopyHistory_hh


namespace pxl {

// ptl

class ObjectBase;

typedef pxl::Map<pxl::Id, pxl::ObjectBase*> CopyHistory;

} // namespace pxl

#endif // pxl_ptl_CopyHistory_hh
#ifndef pxl_ptl_Index_hh
#define pxl_ptl_Index_hh


namespace pxl {

// ptl

class ObjectBase;

typedef pxl::Map<std::string, pxl::ObjectBase*> Index;

} // namespace pxl

#endif // pxl_ptl_Index_hh
#ifndef pxl_ptl_WkPtrBase_hh
#define pxl_ptl_WkPtrBase_hh

namespace pxl {

// ptl

class ObjectBase;

class WkPtrBase {
  public:
    virtual ~WkPtrBase() { connect(0); }

    virtual pxl::WkPtrBase* clone() const
    { return new pxl::WkPtrBase(*this); }

    inline pxl::ObjectBase* pointer() const { return _objectRef; }
    inline bool valid() const { return _objectRef != 0; }
    inline pxl::ObjectBase* operator->() const { return _objectRef; }

  protected:
    WkPtrBase() :
        _notifyChainIn(0),
        _notifyChainOut(0),
        _objectRef(0) {}

    void notifyDeleted();

    void connect(pxl::ObjectBase* pointer);

    pxl::WkPtrBase* _notifyChainIn;
    pxl::WkPtrBase* _notifyChainOut;

    pxl::ObjectBase* _objectRef;

  friend class ObjectBase;
};

} // namespace pxl

#endif // pxl_ptl_WkPtrBase_hh
#ifndef pxl_ptl_WkPtrSpec_hh
#define pxl_ptl_WkPtrSpec_hh


#ifndef pxl_ptl_GetSet_hh
#define pxl_ptl_GetSet_hh

namespace pxl {

// ptl

class Get { public: Get() {} };
class Set { public: Set() {} };

extern const pxl::Get get;
extern const pxl::Set set;

} // namespace pxl

#endif // pxl_ptl_GetSet_hh

namespace pxl {

// iotl

class iStreamer;
class oStreamer;

// ptl

template<class datatype, class objecttype>
class WkPtrSpec : public pxl::WkPtrBase {
  public:
    WkPtrSpec() : pxl::WkPtrBase() {}
    WkPtrSpec(objecttype* ptr) : pxl::WkPtrBase()
    { pxl::WkPtrBase::connect(ptr); }
    WkPtrSpec(objecttype& object) : pxl::WkPtrBase()
    { pxl::WkPtrBase::connect(&object); }
    WkPtrSpec(const pxl::WkPtrSpec<datatype, objecttype>& original) : pxl::WkPtrBase()
    { pxl::WkPtrBase::connect((objecttype*) original._objectRef); }
    virtual ~WkPtrSpec()
    { pxl::WkPtrBase::connect(0); }

    virtual pxl::WkPtrBase* clone() const
    { return new pxl::WkPtrSpec<datatype, objecttype>(*this); }

    // for deep copies:
    inline void operator=(const pxl::WkPtrSpec<datatype, objecttype>& pptr)
    { connect(pptr._objectRef); }
    inline void operator=(objecttype& object)
    { connect(&object); }
    inline void operator=(objecttype* objectptr)
    { connect(objectptr); }

    // methods to grant object & data access
    inline objecttype& object() const  { return *access(); }
    inline const datatype& get() const { return access()->get(); }
    inline       datatype& set()       { return access()->set(); }

    inline const datatype& operator()()                const { return get(); }
    inline const datatype& operator()(const pxl::Get&) const { return get(); }
    inline       datatype& operator()(const pxl::Set&)       { return set(); }

    inline objecttype* operator->() const { return access(); }

  protected:
    // safe access to object
    inline objecttype* access() const
    {
        if (_objectRef) return (objecttype*)_objectRef;
        std::cerr << "pxl::WkPtrSpec::access(): FATAL: The object you intend to access does not exist!" << std::endl;
        return 0; // FIXME: Exception?
    }

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;
};

template<class datatype, class objecttype>
objecttype& operator*(pxl::WkPtrSpec<datatype, objecttype>& wkPtr)
{ return wkPtr.object(); }

template<class datatype, class objecttype>
const objecttype& operator*(const pxl::WkPtrSpec<datatype, objecttype>& wkPtr)
{ return wkPtr.object(); }

} // namespace pxl

#endif // pxl_ptl_WkPtrSpec_hh
#ifndef pxl_ptl_WkPtrOwner_hh
#define pxl_ptl_WkPtrOwner_hh



namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// ptl

class WkPtrBase;
class ObjectBase;

template<class keytype>
class WkPtrOwner : public pxl::Map<keytype, pxl::WkPtrBase*> {
  public:
    typedef          pxl::Map<keytype, pxl::WkPtrBase*>                   PtlMap;
    typedef typename pxl::Map<keytype, pxl::WkPtrBase*>::StlConstIterator StlConstIterator;
    typedef typename pxl::Map<keytype, pxl::WkPtrBase*>::StlIterator      StlIterator;

    WkPtrOwner() : pxl::Map<keytype, pxl::WkPtrBase*>() {}
    virtual ~WkPtrOwner() { pxl::WkPtrOwner<keytype>::clearContainer(); }

    void set(const keytype &key, pxl::WkPtrBase* wptr)
    {
        // FIXME: efficiency, remove & set in one?
        pxl::WkPtrOwner<keytype>::remove(key);
        PtlMap::set(key, wptr);
    }
    void set(const keytype& key, pxl::ObjectBase& obj);

    void remove(const keytype& key)
    {
        delete find(key, 0);
        PtlMap::_container.erase(key);
    }

    bool has(const keytype& key) const { return find(key, 0) != 0; }

    virtual void clearContainer();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - -
    class Iterator {
      public:
        Iterator(const pxl::WkPtrOwner<keytype>::Iterator& original) :
            _iter(original._iter), _containerRef(original._containerRef) {}
        Iterator(const pxl::WkPtrOwner<keytype>& map) :
            _containerRef(&map) { first(); }

        inline void first()  { _iter = _containerRef->begin(); }
        inline void next()   { _iter++; }
        inline bool isDone() { return _iter == _containerRef->end(); }

        inline keytype key()           { return _iter->first; }
        inline pxl::WkPtrBase* item()  { return _iter->second; }
        inline pxl::WkPtrBase& wkPtr() { return *item(); }
        inline pxl::ObjectBase& object()
        {
            if (!wkPtr().valid()) // FIXME: exception?
                std::cerr << "pxl::WkPtrOwner::object(): FATAL: The object you intend to access does not exist!" << std::endl;
            return *wkPtr().pointer();
        }

      private:
        Iterator() : _containerRef(0) {}

        StlConstIterator                _iter;
        const pxl::WkPtrOwner<keytype>* _containerRef;
    };
    // - - - - - - - - - - - - - - - - - - - - - - - - - - -
    template<class wkptrtype>
    class TypeIterator {
      public:
        TypeIterator(const pxl::WkPtrOwner<keytype>::TypeIterator<wkptrtype>& original) :
            _iter(original._iter), _containerRef(original._containerRef) {}

        TypeIterator(const pxl::WkPtrOwner<keytype>& oo) :
            _containerRef(&oo) { first(); }

        inline void first()  { for(_iter = _containerRef->begin(); !tryItem(); _iter++); }
        inline void next()   { do _iter++; while(!tryItem()); }
        inline bool isDone() { return _iter == _containerRef->end(); }

        inline keytype key()            { return _iter->first; }
        inline wkptrtype* item()        { return dynamic_cast<wkptrtype*>(_iter->second); }
        inline const wkptrtype& wkPtr() { return *item(); }
        inline pxl::ObjectBase& object()
        {
            if (!wkPtr().valid()) // FIXME: exception?
                std::cerr << "pxl::WkPtrOwner::object(): FATAL: The object you intend to access does not exist!" << std::endl;
            return *wkPtr().pointer();
        }

      private:
        TypeIterator() : _containerRef(0) {}

        inline bool tryItem()
        {
            if (isDone()) return true;
            return dynamic_cast<wkptrtype*>(_iter->second) != 0;
        }

        StlConstIterator                _iter;
        const pxl::WkPtrOwner<keytype>* _containerRef;
    };
    // - - - - - - - - - - - - - - - - - - - - - - - - - - -

  private:
    WkPtrOwner(const pxl::WkPtrOwner<keytype>& original);

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;
};

} // namespace pxl


#endif // pxl_ptl_WkPtrOwner_hh
#ifndef pxl_ptl_Relations_hh
#define pxl_ptl_Relations_hh


namespace pxl {

typedef pxl::WkPtrOwner<pxl::Id> Relations;

} // namespace pxl

#endif // pxl_ptl_Relations_hh
#ifndef pxl_ptl_ObjectBase_hh
#define pxl_ptl_ObjectBase_hh


#ifndef pxl_ptl_Layout_hh
#define pxl_ptl_Layout_hh

namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// ptl

class Layout {
  public:
    Layout() : _type(0), _style(0), _color(0), _a(0.), _b(0.), _c(0.), _d(0.) {}
    ~Layout() {}

    inline int getType() { return _type; }
    inline int getStyle() { return _style; }
    inline int getColor() { return _color; }

    inline double getA() { return _a; }
    inline double getB() { return _b; }
    inline double getC() { return _c; }
    inline double getD() { return _d; }

    inline void setType(int v){ _type = v; }
    inline void setStyle(int v) { _style = v; }
    inline void setColor(int v) { _color = v; }

    inline void setA(double v) { _a = v; }
    inline void setB(double v) { _b = v; }
    inline void setC(double v) { _c = v; }
    inline void setD(double v) { _d = v; }

  protected:
    int _type;
    int _style;
    int _color;
 
    double _a;
    double _b;
    double _c;
    double _d;

  friend class iStreamer;
  friend class oStreamer;
};

} // namespace pxl

#endif // pxl_ptl_Layout_hh

namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// ptl

class ObjectOwner;

class ObjectBase {
  public:
    virtual ~ObjectBase()
    {
        if (_refWkPtrSpec)
            _refWkPtrSpec->notifyDeleted();
        delete _ptrLayout;
    }

    // FIXME: not good!
    inline pxl::Id id() const { return (int)(void*)this;  }

    inline pxl::ObjectOwner* holder() const { return _refObjectOwner; }

    virtual pxl::ObjectBase* clone() const { return new pxl::ObjectBase(*this); }

    inline pxl::Relations& getMotherRelations()   { return _motherRelations; }
    inline pxl::Relations& getDaughterRelations() { return _daughterRelations; }

    void linkMother(pxl::ObjectBase& target);
    void linkDaughter(pxl::ObjectBase& target);

    inline void linkMother(pxl::WkPtrBase& target)   { if (target._objectRef) linkMother(*(target._objectRef)); }
    inline void linkDaughter(pxl::WkPtrBase& target) { if (target._objectRef) linkDaughter(*(target._objectRef)); }

    void unlinkMother(pxl::ObjectBase& target);
    void unlinkDaughter(pxl::ObjectBase& target);

    inline void unlinkMother(pxl::WkPtrBase& target)   { if (target._objectRef) unlinkMother(*(target._objectRef)); }
    inline void unlinkDaughter(pxl::WkPtrBase& target) { if (target._objectRef) unlinkDaughter(*(target._objectRef)); }

    void unlinkMothers();
    void unlinkDaughters();

    inline pxl::Layout& layout()
    {
        if (!_ptrLayout)
            _ptrLayout = new pxl::Layout;
        return *_ptrLayout;
    }

    std::ostream& printDecayTree(int level = 0, std::ostream& os = std::cout, int pan = 1) const;
    virtual std::ostream&  print(int level = 0, std::ostream& os = std::cout, int pan = 0) const;

  protected:
    ObjectBase() :
        _refWkPtrSpec(0),
        _refObjectOwner(0),
        _motherRelations(),
        _daughterRelations(),
        _ptrLayout(0) {}

    ObjectBase(const ObjectBase& original) :
        _refWkPtrSpec(0),
        _refObjectOwner(0),
        _motherRelations(),
        _daughterRelations(),
        _ptrLayout(0)
    {
        if (original._ptrLayout)
            _ptrLayout = new pxl::Layout;
    }

    virtual pxl::WkPtrBase* createSelfWkPtr()
    { pxl::exception("pxl::ObjectBase::createSelfWkPtr()", "ATTENTION! Inheriting class must reimplement this virtual method."); return 0; }
    virtual void storeYourSelf(pxl::oStreamer& output) const
    { pxl::exception("pxl::ObjectBase::storeYourSelf()", "ATTENTION! Inheriting class must reimplement this virtual method."); }

    std::ostream& printPan1st(std::ostream& os, int pan) const;
    std::ostream& printPan(std::ostream& os, int pan) const;

    pxl::WkPtrBase*   _refWkPtrSpec;
    pxl::ObjectOwner* _refObjectOwner;

    pxl::Relations _motherRelations;
    pxl::Relations _daughterRelations;

    pxl::Layout* _ptrLayout;

  friend class WkPtrBase;
  template<class keytype>
  friend class WkPtrOwner;
  friend class ObjectOwner; 
      
  friend class iStreamer;
  friend class oStreamer;
};

} // namespace pxl

// ptl operators

std::ostream& operator<<(std::ostream& cxxx, const pxl::ObjectBase& obj);

// FIXME: template methods for WkPtrOwner
namespace pxl {

template<class keytype>
void WkPtrOwner<keytype>::clearContainer()
{
    for(StlConstIterator iter = PtlMap::_container.begin();
        iter != PtlMap::_container.end(); iter++)
        delete iter->second;
    PtlMap::_container.clear();
}

template<class keytype>
void WkPtrOwner<keytype>::set(const keytype& key, pxl::ObjectBase& obj)
{
    WkPtrOwner<keytype>::remove(key);
    PtlMap::set(key, obj.createSelfWkPtr());
}

} // namespace pxl

#endif // pxl_ptl_ObjectBase_hh
#ifndef pxl_ptl_Object_hh
#define pxl_ptl_Object_hh


#ifndef pxl_ptl_WkPtr_hh
#define pxl_ptl_WkPtr_hh


namespace pxl {

// ptl

template<class datatype>
class Object;

template<class datatype>
class WkPtr : public pxl::WkPtrSpec<datatype, pxl::Object<datatype> > {
  public:
    WkPtr() :
    	pxl::WkPtrSpec<datatype, pxl::Object<datatype> >() {}
    WkPtr(pxl::Object<datatype>* ptr) :
        pxl::WkPtrSpec<datatype, pxl::Object<datatype> >(ptr) {}
    WkPtr(pxl::Object<datatype>& object) :
        pxl::WkPtrSpec<datatype, pxl::Object<datatype> >(object) {}
    WkPtr(const pxl::WkPtrSpec<datatype, pxl::Object<datatype> >& original) :
        pxl::WkPtrSpec<datatype, pxl::Object<datatype> >(original) {}
};

} // namespace pxl


#endif // pxl_ptl_WkPtr_hh

namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// ptl

template<class datatype>
class Object : public pxl::ObjectBase {
  public:
    Object() : pxl::ObjectBase(), _data() {}
    Object(const datatype& original) : pxl::ObjectBase(), _data(original) {}
    Object(const pxl::Object<datatype>& original) :
        pxl::ObjectBase(original), _data(original._data) {}
    virtual ~Object() {}

    // for deep copies:
    virtual pxl::ObjectBase* clone() const { return new pxl::Object<datatype>(*this); }

    inline const datatype& get() const { return _data; }
    inline       datatype& set()       { return _data; }

    inline const datatype& operator()()                const { return get(); }
    inline const datatype& operator()(const pxl::Get&) const { return get(); }
    inline       datatype& operator()(const pxl::Set&)       { return set(); }

    inline pxl::Object<datatype>& operator=(const datatype& original)
    { _data = original; return *this; }
    inline pxl::Object<datatype>& operator=(const pxl::Object<datatype>& original)
    { _data = original._data; return *this; }

    virtual std::ostream& print(int level = 0, std::ostream& os = std::cout, int pan = 0) const
    {
        os << "called by pxl::Object<...>: ";
        return pxl::ObjectBase::print(level, os, pan);
    }

  protected:
    virtual pxl::WkPtrBase* createSelfWkPtr()
    { return new pxl::WkPtr<datatype>(*this); }

    virtual void storeYourSelf(pxl::oStreamer& output) const;

    datatype _data;

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;
};

} // namespace pxl

#ifndef pxl_iotl_oStreamer_hh
#define pxl_iotl_oStreamer_hh

#include <typeinfo>
#include <sstream>



#ifndef pxl_iotl_internal_hh
#define pxl_iotl_internal_hh

// defaults

#define iotl__default__compressionMode '6'

// file format

#define iotl__eventMarker     "<E>"

#endif // pxl_iotl_internal_hh
#ifndef pxl_iotl_TypeManager_hh
#define pxl_iotl_TypeManager_hh



#ifndef pxl_iotl_TypeIdKey_hh
#define pxl_iotl_TypeIdKey_hh


namespace pxl {

class TypeIdKey {
  public:
    TypeIdKey(const std::string& objectTypeId, const std::string& dataTypeId) : 
        _objectTypeId(objectTypeId), _dataTypeId(dataTypeId) {}

    inline bool operator<(const pxl::TypeIdKey& other) const
    {
        if (_dataTypeId != other._dataTypeId)
            return _dataTypeId < other._dataTypeId;
        else
            return _objectTypeId < other._objectTypeId;
    }

    inline bool operator==(const pxl::TypeIdKey& other) const
    {
        return _dataTypeId   == other._dataTypeId &&
               _objectTypeId == other._objectTypeId;
    }

    // FIXME: why public?
    const std::string& _objectTypeId;
    const std::string& _dataTypeId;

  private:
    TypeIdKey();
};

} // namespace pxl

#endif // pxl_iotl_TypeIdKey_hh

namespace pxl {

class TypeAgentBase;

class TypeManager {
  public:
    static TypeManager& instance();

    void registerAgent(pxl::TypeAgentBase* interface);

    void storeObject(pxl::oStreamer& output, const pxl::ObjectBase& obj, const std::string& cppTypeId) const;
    pxl::Id restoreObject(pxl::iStreamer& input, pxl::ObjectBase& obj, const std::string& cppTypeId) const;
    pxl::Id restoreObject(pxl::iStreamer& input, pxl::ObjectBase** ppobj, const std::string& cppTypeId, const std::string& dataTypeId) const;

  protected:
    TypeManager() {} // locked for singletons
    ~TypeManager() {} // locked for singletons

    pxl::Map<TypeIdKey, TypeAgentBase*>   _agentsByIotlTypeId;
    pxl::Map<std::string, TypeAgentBase*> _agentsByCppTypeId;

  private:
    TypeManager(const TypeManager&); // locked for singletons
    TypeManager& operator=(const TypeManager&); // locked for singletons

    static TypeManager* _instance;
};

} // namespace pxl

#endif // pxl_iotl_TypeManager_hh

namespace pxl {

// iotl

class oStreamer : public pxl::BasicIoStreamer {
  public:
    oStreamer() : pxl::BasicIoStreamer(), _buffer() {}
    ~oStreamer() {}

    void getEvent(std::ostream& cxxx, const std::string& info, char compressionMode = iotl__default__compressionMode);

    template<class objecttype>
    void storeObject(const objecttype& obj)
    { pxl::TypeManager::instance().storeObject(*this, obj, typeid(obj).name()); }

    inline void storeAbstractObject(const pxl::ObjectBase& data)
    { data.storeYourSelf(*this); }

    template<class datatype>
    void storeData(const datatype& data);
    template<class itemtype>
    void storeData(const pxl::Vector<itemtype>& vector);
    template<class keytype, class itemtype>
    void storeData(const pxl::Map<keytype, itemtype>& map);

    inline void storeTypeId(const char* typeId)
    { pxl::BasicIoStreamer::storeBasicTypeCStr(_buffer, typeId); }
    inline void storeTypeId(std::ostream& cxxx, const char* typeId)
    { pxl::BasicIoStreamer::storeBasicTypeCStr(cxxx, typeId); }

  protected:
    inline void storeId(pxl::Id& id)
    { storeMemory((const char*)&id, 4); } // FIXME: endianess/64bit
    inline void storeMemory(const char* address, int bytes)
    { pxl::BasicIoStreamer::dumpMemory(_buffer, address, bytes); }

    std::stringstream _buffer;
};

// template methods

template<class itemtype>
void oStreamer::storeData(const pxl::Vector<itemtype>& vector)
{
    typedef typename pxl::Vector<itemtype>::StlConstIterator PtlVectorConstIterator;

    storeData<int>(vector.getSize());
    for(PtlVectorConstIterator iter = vector.begin();
        iter != vector.end(); iter++)
        storeData<itemtype>(*iter);
}

template<class keytype, class itemtype>
void oStreamer::storeData(const pxl::Map<keytype, itemtype>& map)
{
    typedef typename pxl::Map<keytype, itemtype>::StlConstIterator PtlMapConstIterator;

    storeData<int>(map.getSize());
    for(PtlMapConstIterator iter = map.begin();
        iter != map.end(); iter++) {

        storeData<keytype>(iter->first);
        storeData<itemtype>(iter->second);
    }
}

template<> void oStreamer::storeData<pxl::ObjectOwner>(const pxl::ObjectOwner& objects);

} // namespace pxl

#endif // pxl_iotl_oStreamer_hh

namespace pxl {

// template methods

template<class datatype>
void pxl::Object<datatype>::storeYourSelf(pxl::oStreamer& output) const
{ output.storeObject(*this); }

#define PCL_PRINT_NATIVE(type) \
template<> \
std::ostream& Object<type>::print(int level, std::ostream& os, int pan) const;

PCL_PRINT_NATIVE(int)
PCL_PRINT_NATIVE(unsigned int)
PCL_PRINT_NATIVE(bool)
PCL_PRINT_NATIVE(double)
PCL_PRINT_NATIVE(float)

#undef PCL_PRINT_NATIVE

} // namespace pxl

#endif // pxl_ptl_Object_hh
#ifndef pxl_ptl_ObjectOwner_hh
#define pxl_ptl_ObjectOwner_hh


namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// ptl

class ObjectBase;

class ObjectOwner : public pxl::Vector<pxl::ObjectBase*> {
  public:
    ObjectOwner();
    ObjectOwner(const pxl::ObjectOwner& original);
    virtual ~ObjectOwner() { pxl::ObjectOwner::clearContainer(); }

    template<class objecttype>
    inline objecttype& create()
    {
        objecttype* pitem = new objecttype;
        pitem->_refObjectOwner = this;
        _container.push_back(static_cast<pxl::ObjectBase*>(pitem));
        return *pitem;
    }

    template<class objecttype>
    inline objecttype& create(const objecttype& original)
    {
        objecttype* pitem = new objecttype(original);
        pitem->_refObjectOwner = this;
        _container.push_back(static_cast<pxl::ObjectBase*>(pitem));
        return *pitem;
    }

    template<class objecttype, class ctrtype>
    inline objecttype& create(const ctrtype& original)
    {
        objecttype* pitem = new objecttype(original);
        pitem->_refObjectOwner = this;
        _container.push_back(static_cast<pxl::ObjectBase*>(pitem));
        return *pitem;
    }

    void set(pxl::ObjectBase& item);
    void remove(pxl::ObjectBase& item);
    bool has(const pxl::ObjectBase& item) const;

    virtual void clearContainer();

    template<class objecttype>
    inline objecttype* findObject(const std::string& idx) const	// goes via Index & casts
    { return dynamic_cast<objecttype*>(_index.find(idx, 0)); }

    template<class objecttype>
    inline objecttype* findCopyOf(const pxl::ObjectBase& original) const // goes via CopyHistory & casts
    { return dynamic_cast<objecttype*>(_copyHistory.find(original.id(), 0)); }    

    inline const pxl::CopyHistory& getCopyHistory() const { return _copyHistory; }
    inline void clearCopyHistory() { _copyHistory.clearContainer(); }

    inline bool setIndex(const std::string& idx, pxl::ObjectBase& obj)
    {
        if (!idx.length() || !has(obj))
            return false;
        _index.set(idx, &obj);
        return true;
    }

    inline const pxl::Index& getIndex() const { return _index; }
    inline void removeIndex(const std::string& idx) { _index.remove(idx); }
    inline void clearIndex() { _index.clearContainer(); }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - -
    class Iterator {
      public:
        Iterator(const pxl::ObjectOwner::Iterator& original) :
            _iter(original._iter), _containerRef(original._containerRef) {}
        Iterator(const pxl::ObjectOwner& vector) :
            _containerRef(&vector) { first(); }

        inline void first()  { _iter = _containerRef->begin(); }
        inline void next()   { _iter++; }
        inline bool isDone() { return _iter == _containerRef->end(); }

        inline pxl::ObjectBase* item() { return *_iter; }
        inline pxl::ObjectBase& object() { return *item(); }

      private:
        Iterator() : _containerRef(0) {;}

        StlConstIterator    _iter;
        const pxl::ObjectOwner* _containerRef;
    };
    // - - - - - - - - - - - - - - - - - - - - - - - - - - -
    template<class objecttype>
    class TypeIterator {
      public:
        TypeIterator(const pxl::ObjectOwner::TypeIterator<objecttype>& original) :
            _iter(original._iter), _containerRef(original._containerRef) {}

        TypeIterator(const pxl::ObjectOwner& oo) :
            _containerRef(&oo) { first(); }

        inline void first()  { for(_iter = _containerRef->begin(); !tryItem(); _iter++); }
        inline void next()   { do _iter++; while(!tryItem()); }
        inline bool isDone() { return _iter == _containerRef->end(); }

        inline objecttype* item() { return dynamic_cast<objecttype*>(*_iter); }
        inline objecttype& object() { return *item(); }

      private:
        TypeIterator() : _containerRef(0) {}

        inline bool tryItem()
        {
            if (isDone()) return true;
            return dynamic_cast<objecttype*>(*_iter) != 0;
        }

        StlConstIterator    _iter;
        const pxl::ObjectOwner* _containerRef;
    };
    // - - - - - - - - - - - - - - - - - - - - - - - - - - -

  protected:
    pxl::CopyHistory _copyHistory;
    pxl::Index       _index;

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;
};

} // namespace pxl

#endif // pxl_ptl_ObjectOwner_hh
#ifndef pxl_ptl_Objects_hh
#define pxl_ptl_Objects_hh


namespace pxl {

typedef ObjectOwner Objects;

} // namespace pxl

#endif // pxl_ptl_Objects_hh
#ifndef pxl_ptl_CowObject_hh
#define pxl_ptl_CowObject_hh


#ifndef pxl_ptl_CowWkPtr_hh
#define pxl_ptl_CowWkPtr_hh


namespace pxl {

// ptl

template<class datatype>
class CowObject;

template<class datatype>
class CowWkPtr : public pxl::WkPtrSpec<datatype, pxl::CowObject<datatype> > {
  public:
    CowWkPtr() :
    	pxl::WkPtrSpec<datatype, pxl::CowObject<datatype> >() {}
    CowWkPtr(pxl::CowObject<datatype>* ptr) :
        pxl::WkPtrSpec<datatype, pxl::CowObject<datatype> >(ptr) {}
    CowWkPtr(pxl::CowObject<datatype>& object) :
        pxl::WkPtrSpec<datatype, pxl::CowObject<datatype> >(object) {}
    CowWkPtr(const pxl::WkPtrSpec<datatype, pxl::CowObject<datatype> >& original) :
        pxl::WkPtrSpec<datatype, pxl::CowObject<datatype> >(original) {}
};

} // namespace pxl


#endif // pxl_ptl_CowWkPtr_hh

namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// ptl

template<class datatype>
class CowObject : public pxl::ObjectBase {
  private:
    class DataSocket {
      public:
        DataSocket() : _references(1) {}
        DataSocket(const datatype& original) :
            _references(1), _data(original) {}
        DataSocket(const pxl::CowObject<datatype>::DataSocket& original) :
            _references(1), _data(original._data) {}
        virtual ~DataSocket() {}

        // for deep copies
        virtual DataSocket* clone() const { return new DataSocket(*this); }

        // methods to grant data access
        inline datatype& getData() { return _data; }
        inline void setData(const datatype& object) { _data = object; }

        unsigned int _references;
        datatype     _data;

      friend class pxl::iStreamer;
      friend class pxl::oStreamer;
    };

  public:
    CowObject() : pxl::ObjectBase() { _dataSocket = new DataSocket; }
    CowObject(const datatype& original) : pxl::ObjectBase() { _dataSocket = new DataSocket; }
    CowObject(const pxl::CowObject<datatype>& original) : pxl::ObjectBase(original)
    { _dataSocket = original._dataSocket; _dataSocket->_references++; }
    virtual ~CowObject() { dropDataSocket(); }

    // for deep copies:
    virtual pxl::ObjectBase* clone() const { return new pxl::CowObject<datatype>(*this); }

    inline const datatype& get() const { return _dataSocket->getData(); }
    inline datatype& set()
    {
        if (_dataSocket->_references > 1) {
            _dataSocket->_references--;
            _dataSocket = new DataSocket(*_dataSocket);
        }
        return _dataSocket->getData(); 
    }

    inline const datatype& operator()()                const { return get(); }
    inline const datatype& operator()(const pxl::Get&) const { return get(); }
    inline       datatype& operator()(const pxl::Set&)       { return set(); }

    inline pxl::CowObject<datatype>& operator=(const datatype& original)
    {
        dropDataSocket();
        _dataSocket = new DataSocket(original);
        return *this;
    }
    inline pxl::CowObject<datatype>& operator=(const pxl::CowObject<datatype>& original)
    {
        dropDataSocket();
        _dataSocket = original._dataSocket;
        _dataSocket->_references++;
        return *this;
    }

    virtual std::ostream& print(int level = 0, std::ostream& os = std::cout, int pan = 0) const
    {
        os << "called by pxl::CowObject<...>: ";
        return pxl::ObjectBase::print(level, os, pan);
    }

  protected:
    CowObject(DataSocket& original) : pxl::ObjectBase()
    { _dataSocket = &original; _dataSocket->_references++; }

    inline void dropDataSocket()
    {
        if (_dataSocket->_references-- == 1)
            delete _dataSocket;
    }

    virtual pxl::WkPtrBase* createSelfWkPtr()
    { return new pxl::CowWkPtr<datatype>(*this); }

    virtual void storeYourSelf(pxl::oStreamer& output) const;

    DataSocket* _dataSocket;

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;
};

} // namespace pxl


namespace pxl {

// template methods

template<class datatype>
void pxl::CowObject<datatype>::storeYourSelf(pxl::oStreamer& output) const
{ output.storeObject(*this); }

#define PCL_PRINT_NATIVE(type) \
template<> \
std::ostream& CowObject<type>::print(int level, std::ostream& os, int pan) const;

PCL_PRINT_NATIVE(int)
PCL_PRINT_NATIVE(unsigned int)
PCL_PRINT_NATIVE(bool)
PCL_PRINT_NATIVE(double)
PCL_PRINT_NATIVE(float)

#undef PCL_PRINT_NATIVE

} // namespace pxl

#endif // pxl_ptl_CowObject_hh
#ifndef pxl_ptl_SpyObject_hh
#define pxl_ptl_SpyObject_hh


#ifndef pxl_ptl_SpyWkPtr_hh
#define pxl_ptl_SpyWkPtr_hh


namespace pxl {

// ptl

template<class datatype>
class SpyObject;

template<class datatype>
class SpyWkPtr : public pxl::WkPtrSpec<datatype, pxl::SpyObject<datatype> > {
  public:
    SpyWkPtr() : pxl::WkPtrSpec<datatype, pxl::SpyObject<datatype> >() {}
    SpyWkPtr(pxl::SpyObject<datatype>* ptr) :
        pxl::WkPtrSpec<datatype, pxl::SpyObject<datatype> >(ptr) {}
    SpyWkPtr(pxl::SpyObject<datatype>& object) :
        pxl::WkPtrSpec<datatype, pxl::SpyObject<datatype> >(object) {}
    SpyWkPtr(const pxl::WkPtrSpec<datatype, pxl::SpyObject<datatype> >& original) :
        pxl::WkPtrSpec<datatype, pxl::SpyObject<datatype> >(original) {}
};

} // namespace pxl


#endif // pxl_ptl_SpyWkPtr_hh

namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// ptl

template<class datatype>
class SpyObject : public pxl::Object<pxl::Ptr<datatype> > {
  public:
    SpyObject() :
        pxl::Object<pxl::Ptr<datatype> >() {}
    SpyObject(datatype* original) :
        pxl::Object<pxl::Ptr<datatype> >(pxl::Ptr<datatype>(original)) {}
    SpyObject(const pxl::Ptr<datatype>& original) :
        pxl::Object<pxl::Ptr<datatype> >(original) {}
    SpyObject(const pxl::SpyObject<datatype>& original) :
        pxl::Object<pxl::Ptr<datatype> >(original.get()) {}
    virtual ~SpyObject() {}

    // for deep copies:
    virtual pxl::ObjectBase* clone() const { return new pxl::SpyObject<datatype>(*this); }

    inline pxl::SpyObject<datatype>& operator=(datatype *original)
    {
        pxl::Object<pxl::Ptr<datatype> >::operator=(pxl::Ptr<datatype>(original));
        return *this;
    }
    inline pxl::SpyObject<datatype>& operator=(const pxl::Ptr<datatype>& original)
    {
        pxl::Object<pxl::Ptr<datatype> >::operator=(original);
        return *this;
    }
    inline pxl::SpyObject<datatype>& operator=(const pxl::SpyObject<datatype>& original)
    {
        pxl::Object<pxl::Ptr<datatype> >::operator=(pxl::Ptr<datatype>(original.get()));
        return *this;
    }

    virtual std::ostream& print(int level = 0, std::ostream& os = std::cout, int pan = 0) const
    {
        os << "called by pxl::SpyObject<...>: ";
        return pxl::ObjectBase::print(level, os, pan);
    }

  protected:
    virtual pxl::WkPtrBase* createSelfWkPtr()
    { return new pxl::SpyWkPtr<datatype>(*this); }

    virtual void storeYourSelf(pxl::oStreamer& output) const;

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;
};

} // namespace pxl


namespace pxl {

// template methods

template<class datatype>
void SpyObject<datatype>::storeYourSelf(pxl::oStreamer& output) const
{ output.storeObject(*this); }

} // namespace pxl

#endif // pxl_ptl_SpyObject_hh

#ifndef pxl_iotl_hh
#define pxl_iotl_hh

#ifndef pxl_iotl_iStreamer_hh
#define pxl_iotl_iStreamer_hh





namespace pxl {

// iotl

class iStreamer : public pxl::BasicIoStreamer {
  public:
    iStreamer();
    virtual ~iStreamer();

    inline bool putEventIf(std::istream& cxxx, const std::string& infoCondition)
    { return putEvent(cxxx, -1, infoCondition); }
    inline bool putEvent(std::istream& cxxx)
    { return putEvent(cxxx, +1, ""); }

    inline bool next(std::istream& cxxx) { return putEvent(cxxx,  0, ""); }
    bool previous(std::istream& cxxx);
    inline bool endOfEvent() { return _buffer.peek() == EOF; }

    template<class objecttype>
    pxl::Id restoreObject(objecttype& obj)
    { return pxl::TypeManager::instance().restoreObject(*this, obj, typeid(obj).name()); }

    inline pxl::ObjectBase& restoreAbstractObject()
    {
        pxl::ObjectBase* pobj;
        restoreAbstractObject(&pobj);
        return *pobj;
    }

    template<class datatype>
    pxl::Id restoreData(datatype& data);
    template<class itemtype>
    pxl::Id restoreData(pxl::Vector<itemtype>& vector);
    template<class keytype, class itemtype>
    pxl::Id restoreData(pxl::Map<keytype, itemtype>& map);

    inline void restoreTypeId(const char* expectedTypeId)
    {
        std::string read;
        pxl::BasicIoStreamer::restoreBasicTypeCStr(_buffer, read);
        if (read != expectedTypeId)
            pxl::exception("pxl::iStreamer::restoreTypeId()",
                           std::string("Unexpected object type: ") + read);
    }

    inline void restoreTypeId(std::istream& cxxx, const char* expectedTypeId)
    {
        std::string read;
        pxl::BasicIoStreamer::restoreBasicTypeCStr(cxxx, read);
        if (read != expectedTypeId)
            pxl::exception("pxl::iStreamer::restoreTypeId()",
                           std::string("Unexpected object type: ") + read);
    }

  protected:
    bool putEvent(std::istream& cxxx, int mode, const std::string& infoCondition);
    pxl::Id restoreAbstractObject(pxl::ObjectBase** ppobj);

    inline void restoreId(pxl::MutableId& persistentId)
    {
        restoreMemory((char*)&persistentId, 4); // FIXME: endianess/64bit
        persistentId *= -1;	// FIXME: What's this for?
    }

    using BasicIoStreamer::restoreMemory;
    inline void restoreMemory(char* address, int bytes)
    { pxl::BasicIoStreamer::restoreMemory(_buffer, address, bytes); }

    std::stringstream _buffer;

  private:
    int unzipEventData(std::istream &in, int nBytes);

    unsigned char *_inputBuffer;
    unsigned char *_outputBuffer;
};

// template methods

template<class itemtype>
pxl::Id iStreamer::restoreData(pxl::Vector<itemtype>& vector)
{
    vector.clearContainer();
    int size = 0;
    pxl::iStreamer::restoreData<int>(size);
    for(int i = 0; i < size; i++) {
        itemtype item;
        pxl::iStreamer::restoreData<itemtype>(item);
        vector.set(item);
    }
    return 0;
}

template<class keytype, class itemtype>
pxl::Id iStreamer::restoreData(pxl::Map<keytype, itemtype>& map)
{
    map.clearContainer();
    int size = 0;
    pxl::iStreamer::restoreData<int>(size);
    for(int i = 0; i < size; i++) {
        keytype key;
        itemtype item;
        pxl::iStreamer::restoreData<keytype>(key);
        pxl::iStreamer::restoreData<itemtype>(item);
        map.set(key, item);
    }
    return 0;
}

class ObjectOwner;
class Layout;

template<> pxl::Id iStreamer::restoreData<pxl::Relations>(pxl::Relations& relations);
template<> pxl::Id iStreamer::restoreData<pxl::ObjectOwner>(pxl::ObjectOwner& objects);
template<> pxl::Id iStreamer::restoreData<pxl::Layout>(pxl::Layout& layout);
template<> pxl::Id iStreamer::restoreData<pxl::ObjectBase>(pxl::ObjectBase& object);

} // namespace pxl

#endif // pxl_iotl_iStreamer_hh
#ifndef pxl_iotl_iFile_hh
#define pxl_iotl_iFile_hh

namespace pxl {

class iFile {
  public:
    virtual ~iFile() {}

    virtual bool open(const std::string& filename)
    { return false; }

    virtual void close() {}
    virtual bool readEvent() { return false; }
    virtual bool readEventIf(const std::string& info) { return false; }
    inline  int  skipEvent() { return skipEvents(1); }
    virtual int  skipEvents(int n) { return 0; }
    virtual bool endOfFile() { return true; }

  protected:
    iFile() {}

  private:
    iFile(const pxl::iFile&) {}
};

} // namespace pxl

#endif // pxl_iotl_iFile_hh
#ifndef pxl_iotl_iDiskFileVx_hh
#define pxl_iotl_iDiskFileVx_hh

#include <fstream>
#include <ios>



namespace pxl {

template<class istreamertype>
class iDiskFileVx : public pxl::iFile, public istreamertype {
  public:
    iDiskFileVx() : pxl::iFile(), istreamertype() {}
    virtual ~iDiskFileVx() { close(); }

    virtual bool open(const std::string& filename)
    {
        close();
        _file.open(filename.c_str(), std::fstream::in | std::fstream::binary);
        if (!_file.is_open()) {
            // FIXME: cerr...
            std::cerr << "pxl::oDiskFileVx::open(): error while opening file "
                      << filename << std::endl;
            return false;
        }
        return true;
    }

    virtual void close()
    {
        if (!_file.is_open()) return;
        _file.close();
    }

    virtual bool readEvent()
    {
        if (_file.is_open()) {
            return istreamertype::putEvent(_file);
        } else {
            pxl::exception("pxl::iDiskFileVx<>::readEvent()", "No file open.");
            return false;
        }
    }

    virtual bool readEventIf(const std::string& info)
    {
        if (_file.is_open()) {
            return istreamertype::putEventIf(_file, info);
        } else {
            pxl::exception("pxl::iDiskFileVx<>::readEventIf()", "No file open.");
            return false;
        }
    }

    virtual int skipEvents(int n)
    {
        int s = 0;
        for(; n < 0 && istreamertype::previous(_file); n++) s--;
        for(; n > 0 && istreamertype::next(_file);     n--) s++;
        return s;
    }

    virtual bool endOfFile()
    { return _file.peek() == EOF; }

  protected:
    std::fstream _file;

  private:
    iDiskFileVx(const pxl::iDiskFileVx<istreamertype>&) {}
};

} // namespace pxl

#endif // pxl_iotl_iDiskFileVx_hh
#ifndef pxl_iotl_iDiskFile_hh
#define pxl_iotl_iDiskFile_hh


namespace pxl {

typedef pxl::iDiskFileVx<pxl::iStreamer> iDiskFile;

} // namespace pxl

#endif // pxl_iotl_iDiskFile_hh
#ifndef pxl_iotl_oFile_hh
#define pxl_iotl_oFile_hh


namespace pxl {

class oFile {
  public:
    virtual ~oFile() {}

    virtual bool open(const std::string& filename, bool append = false)
    { return false; }

    virtual void close() {}

    virtual void writeEvent(const std::string& info = "", char compressionMode = iotl__default__compressionMode) {}

  protected:
    oFile() {}

  private:
    oFile(const pxl::oFile&) {}
};

} // namespace pxl

#endif // pxl_iotl_oFile_hh
#ifndef pxl_iotl_oDiskFileVx_hh
#define pxl_iotl_oDiskFileVx_hh




namespace pxl {

template<class ostreamertype>
class oDiskFileVx : public pxl::oFile, public ostreamertype {
  public:
    oDiskFileVx() : pxl::oFile(), ostreamertype() {}
    virtual ~oDiskFileVx() { close(); }

    virtual bool open(const std::string& filename, bool append = false)
    {
        close();
        if (append)
            _file.open(filename.c_str(), std::fstream::app | std::fstream::out | std::fstream::binary);
        else
            _file.open(filename.c_str(), std::fstream::trunc | std::fstream::out | std::fstream::binary);
        if (!_file.good()) {
            // FIXME: cerr...
            std::cerr << "pxl::oDiskFileVx::open(): error while opening file "
                      << filename << std::endl;
            return false;
        }
        return true;
    }

    virtual void close()
    {
        if (!_file.is_open()) return;
        _file.rdbuf()->pubsync();
        _file.close();
    }

    virtual void writeEvent(const std::string& info = "", char compressionMode = iotl__default__compressionMode)
    {
        if (_file.is_open())
            ostreamertype::getEvent(_file, info, compressionMode);
        else
            pxl::exception("pxl::oDiskFileVx<>::writeEvent()", "No file open.");
    }

  protected:
    std::fstream _file;

  private:
    oDiskFileVx(const pxl::oDiskFileVx<ostreamertype>&) {}
};

} // namespace pxl

#endif // pxl_iotl_oDiskFileVx_hh
#ifndef pxl_iotl_oDiskFile_hh
#define pxl_iotl_oDiskFile_hh


namespace pxl {

typedef pxl::oDiskFileVx<pxl::oStreamer> oDiskFile;

} // namespace pxl

#endif // pxl_iotl_oDiskFile_hh
#ifndef pxl_iotl_TypeAgentBase_hh
#define pxl_iotl_TypeAgentBase_hh



namespace pxl {

class TypeAgentBase {
  public:
    virtual ~TypeAgentBase() {}
    virtual void storeObject(pxl::oStreamer& output, const pxl::ObjectBase& objectbase) {}
    virtual pxl::Id restoreObject(pxl::iStreamer& input, pxl::ObjectBase& obj) { return 0; }
    virtual pxl::Id restoreObject(pxl::iStreamer& input, pxl::ObjectBase** ppobj) { return 0; }

    inline const std::string& getObjectTypeId() { return _objectTypeId; }
    inline const std::string& getDataTypeId() { return _dataTypeId; }
    inline const std::string& getCppTypeId() { return _cppTypeId; }

  protected:
    std::string _objectTypeId;
    std::string _dataTypeId;
    std::string _cppTypeId;
};

} // namespace pxl

#endif // pxl_iotl_TypeAgentBase_hh
#ifndef pxl_iotl_TypeAgent_hh
#define pxl_iotl_TypeAgent_hh




namespace pxl {

template<class objecttype>
class TypeAgent : public TypeAgentBase {
  public:
    TypeAgent(const std::string& objectTypeId);
    virtual ~TypeAgent() {}

    virtual void storeObject(pxl::oStreamer &output, const pxl::ObjectBase& obj);
    virtual pxl::Id restoreObject(pxl::iStreamer& input, pxl::ObjectBase& obj);
    virtual pxl::Id restoreObject(pxl::iStreamer& input, pxl::ObjectBase** ppobj);
};

template<class type>
const char* getIotlTypeId(const type* = 0);

// template methods

template<class objecttype>
TypeAgent<objecttype>::TypeAgent(const std::string& objectTypeId)
{
    objecttype obj;

    _objectTypeId = objectTypeId;
    _dataTypeId   = pxl::getIotlTypeId(&obj.get());
    _cppTypeId    = std::string(typeid(obj).name());

    pxl::TypeManager::instance().registerAgent(this);
}

template<class objecttype>
void TypeAgent<objecttype>::storeObject(pxl::oStreamer& output, const pxl::ObjectBase& obj)
{
    const objecttype* ptr = dynamic_cast<const objecttype*>(&obj);

    output.storeTypeId(_objectTypeId.c_str());
    output.storeTypeId(_dataTypeId.c_str());
    output.storeData(ptr->get());
    output.storeData(obj);
}

template<class objecttype>
pxl::Id TypeAgent<objecttype>::restoreObject(pxl::iStreamer& input, pxl::ObjectBase& obj)
{
    objecttype* ptr = dynamic_cast<objecttype*>(&obj); // FIXME static_cast?

    input.restoreTypeId(_objectTypeId.c_str());
    input.restoreTypeId(_dataTypeId.c_str());
    input.restoreData(ptr->set());
    return input.restoreData(obj);
}

template<class objecttype>
pxl::Id TypeAgent<objecttype>::restoreObject(pxl::iStreamer& input, pxl::ObjectBase** ppobj)
{
    objecttype* ptr = new objecttype;
    // please notice: object & data type ids have been read already! 
    input.restoreData(ptr->set());
    (*ppobj) = static_cast<pxl::ObjectBase*>(ptr);
    return input.restoreData(**ppobj);
}

} // namespace pxl

#endif // pxl_iotl_TypeAgent_hh
#ifndef pxl_iotl_types_hh
#define pxl_iotl_types_hh




// complex data types

#define iotl__declareDataTypeExplicit(type, id, interfacename, data, storecode, restorecode) \
namespace pxl { \
\
template<> \
const char* getIotlTypeId<type>(const type*) \
{ return id; } \
\
template<> \
void oStreamer::storeData<type>(const type& data) \
{ storecode; } \
\
template<> \
pxl::Id iStreamer::restoreData<type>(type& data) \
{ do { restorecode; } while(false); return 0; } \
\
TypeAgent<pxl::Object<type> > _typeInterface_ptl__Object_ ## interfacename("\1OO"); \
TypeAgent<pxl::CowObject<type> > _typeInterface_ptl__CowObject_ ## interfacename("\1CO"); \
\
} // namespace pxl

#define iotl__declareDataType(type, data, storecode, restorecode) \
iotl__declareDataTypeExplicit(type, #type, type, data, storecode, restorecode)

#define iotl__declareDataType(type, data, storecode, restorecode) \
iotl__declareDataTypeExplicit(type, #type, type, data, storecode, restorecode)

#define iotl__declareDataTypeProto(type) \
namespace pxl { \
\
template<> \
const char* getIotlTypeId<type>(const type*); \
\
template<> \
void oStreamer::storeData<type>(const type& data); \
\
template<> \
pxl::Id iStreamer::restoreData<type>(type& data); \
\
} // namespace pxl

#define iotl__declareSpyTypeExplicit(type, id, interfacename, ptr, storecode, restorecode) \
namespace pxl { \
\
template<> \
const char* getIotlTypeId<pxl::Ptr<type> >(const pxl::Ptr<type>*) \
{ return id; } \
\
template<> \
void oStreamer::storeData<pxl::Ptr<type> >(const pxl::Ptr<type>& ptr) \
{ storecode; } \
\
template<> \
pxl::Id iStreamer::restoredata<pxl::Ptr<type> >(pxl::Ptr<type>& ptr) \
{ do { restoredata; } while(false); return 0; } \
\
TypeAgent<pxl::SpyObject<type> > _typeInterface_ptl__SpyObject_ ## interfacename("\1SO"); \
\
} // namespace pxl

#define iotl__declareSpyType(type, ptr, storecode, restorecode) \
iotl__declareSpyTypeExplicit(type, #type, type, ptr, storecode, restorecode)

#define iotl__declareSpyTypeProto(type) \
namespace pxl { \
\
template<> \
const char* getIotlTypeId<pxl::Ptr<type> >(const pxl::Ptr<type>*); \
\
template<> \
void oStreamer::storeData<pxl::Ptr<type> >(const pxl::Ptr<type>& ptr); \
\
template<> \
pxl::Id iStreamer::restoredata<pxl::Ptr<type> >(pxl::Ptr<type>& ptr); \
\
} // namespace pxl

#define iotl__declareObjectTypeExplicit(type, id, interfacename) \
namespace pxl { \
\
TypeAgent<type> _typeInterface_ ## interfacename(id); \
\
}

#define iotl__declareObjectType(type) \
iotl__declareObjectTypeExplicit(type, #type, type)

#define iotl__declareObjectTypeProto(type)


iotl__declareDataTypeProto(char)
iotl__declareDataTypeProto(std::string)
iotl__declareDataTypeProto(bool)
iotl__declareDataTypeProto(int)
iotl__declareDataTypeProto(float)
iotl__declareDataTypeProto(double)

#endif // pxl_iotl_types_hh

#endif // pxl_iotl_hh

#endif // pxl_ptl_hh
#ifndef pxl_pol_hh
#define pxl_pol_hh

#ifndef pxl_pol_Basic3Vector_hh
#define pxl_pol_Basic3Vector_hh

#include <cmath>


namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// pol

#define EPSILON 1.0e-9

class Basic3VectorData {
  public:
    Basic3VectorData() :
        _x(0), _y(0), _z(0) {}
    Basic3VectorData(const Basic3VectorData& orig) :
        _x(orig._x), _y(orig._y), _z(orig._z) {}
    Basic3VectorData(double x, double y, double z) :
        _x(x), _y(y), _z(z) {}

    inline void setX(double x) { _x = x; }
    inline void setY(double y) { _y = y; }
    inline void setZ(double z) { _z = z; }

    inline void setRhoPhi(double perp, double phi)
    { _x = std::cos(phi) * perp; _y = std::sin(phi) * perp; }
    inline void setRhoPhiZ(double perp, double phi, double z)
    { setRhoPhi(perp, phi); _z = z; }
    inline void setRThetaPhi(double r, double theta, double phi)
    { setRhoPhiZ(std::cos(theta) * r, phi, std::sin(theta) * r); }

    inline double getX() const { return _x; }
    inline double getY() const { return _y; }
    inline double getZ() const { return _z; }

    inline bool isNullPerp() const
    { return _x > -EPSILON && _x < EPSILON && _y > -EPSILON && _y < EPSILON; }
    inline bool isNull() const
    { return isNullPerp() && _z > -EPSILON && _z < EPSILON; }

    inline double getPerp2() const
    { return _x*_x + _y*_y; }
    inline double getPerp() const
    { return std::sqrt(getPerp2()); }
    inline double getPhi() const
    { return isNullPerp() ? 0.0 : std::atan2(_y, _x); }

    inline double getMag2() const
    { return _x*_x + _y*_y + _z*_z; }
    inline double getMag() const
    { return std::sqrt(getMag2()); }

    inline double getCosTheta() const
    { double mag = getMag(); return mag < EPSILON ? 1.0 : _z / mag; }
    inline double getCos2Theta() const
    { double mag2 = getMag2(); return mag2 < EPSILON ? 1.0 : _z*_z / mag2; }

    inline double getTheta() const
    { return isNull() ? 0.0 : std::atan2(getPerp(), _z); }

    inline double deltaRho(const Basic3VectorData& fv) const
    {
        double dDtheta = deltaTheta(fv);
        double dDphi = deltaPhi(fv);
        return std::sqrt(dDtheta*dDtheta + dDphi*dDphi);
    }

    inline double deltaPhi(const Basic3VectorData& fv) const
    { 
        double dDphi = getPhi() - fv.getPhi();
        while(dDphi > M_PI)
            dDphi -= 2 * M_PI;
        while(dDphi < -M_PI)
            dDphi += 2 * M_PI;
        return dDphi;
    }

    inline double deltaTheta(const Basic3VectorData& fv) const
    { 
        double dDtheta = getTheta() - fv.getTheta();
        while(dDtheta > M_PI)
            dDtheta -= 2 * M_PI;
        while(dDtheta < -M_PI)
            dDtheta += 2 * M_PI;
        return dDtheta;
    }

    inline const pxl::Basic3VectorData& operator=(const pxl::Basic3VectorData& vec)
    { _x = vec._x; _y = vec._y; _z = vec._z; return *this; }

    inline const pxl::Basic3VectorData& operator+=(const pxl::Basic3VectorData& vec)
    { _x += vec._x; _y += vec._y; _z += vec._z; return *this; }
    inline const pxl::Basic3VectorData& operator-=(const pxl::Basic3VectorData& vec)
    { _x -= vec._x; _y -= vec._y; _z -= vec._z; return *this; }

  private:
    double _x;
    double _y;
    double _z;

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;      
};

#undef EPSILON

// non-member operators
bool const operator==(const pxl::Basic3VectorData& obj1, const pxl::Basic3VectorData& obj2);
bool const operator!=(const pxl::Basic3VectorData& obj1, const pxl::Basic3VectorData& obj2);

} // namespace pxl

iotl__declareDataTypeProto(pxl::Basic3VectorData)

#endif // pxl_pol_Basic3Vector_hh
#ifndef pxl_pol_Basic4Vector_hh
#define pxl_pol_Basic4Vector_hh




namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// pol

class Basic4VectorData : public Basic3VectorData {
  public:
    Basic4VectorData() :
        Basic3VectorData(0, 0, 0), _t(0) {}
    Basic4VectorData(const Basic3VectorData &orig, double t = 0) :
        Basic3VectorData(orig), _t(t) {}
    Basic4VectorData(const Basic4VectorData &orig) :
        Basic3VectorData(orig), _t(orig._t) {}
    Basic4VectorData(double x, double y, double z, double t = 0) :
        Basic3VectorData(x, y, z), _t(t) {}

    // setX inherited
    // setY inherited
    // setZ inherited
    inline void setT (double t) { _t = t; }

    inline void setPx(double px) { setX(px); }
    inline void setPy(double py) { setY(py); }
    inline void setPz(double pz) { setZ(pz); }

    inline void setE (double e)  { _t = e; }
    inline void setMass(double m)
    { _t = std::sqrt(m*m + getMag2()); }

    // getX inherited
    // getY inherited
    // getZ inherited
    inline double getT() const { return _t; }

    inline double getPx() const { return getX(); }
    inline double getPy() const { return getY(); }
    inline double getPz() const { return getZ(); }

    inline double getE()  const { return _t; }

    inline double getMass2() const
    { return _t*_t - getMag2(); }

    inline double getMass() const
    { double m2 = getMass2(); return m2 < 0.0 ? 0.0 : std::sqrt(m2); }
    // getPerp inherited
    inline double getPt() const
    { return getPerp(); }

    // getPhi inherited
    // getTheta inherited
    // deltaRho inherited
    // deltaPhi inherited
    // deltaTheta inherited

    inline double getEta() const
    { return -std::log(std::tan(getTheta()*0.5)); }

    inline double getEt2() const
    { double pt2 = getPerp2(); return pt2 == 0.0 ? 0.0 : _t*_t * pt2 / getMag2(); }
    inline double getEt() const
    { return std::sqrt(getEt2()); }

// FIXME: Carsten Hof transverse mass
inline double getMT() const { return sqrt(getEt()*getEt() - getPerp2()); }

    inline double deltaR(const Basic4VectorData& fv) const
    { 
        double dDeta = deltaEta(fv);
        double dDphi = deltaPhi(fv);
        return std::sqrt(dDeta*dDeta + dDphi*dDphi);
    }

    inline double deltaEta(const Basic4VectorData& fv) const
    { return getEta() - fv.getEta(); }

    inline const pxl::Basic4VectorData& operator=(const pxl::Basic3VectorData& vec)
    { Basic3VectorData::operator=(vec); return *this; }
    inline const pxl::Basic4VectorData& operator=(const pxl::Basic4VectorData& vec)
    { Basic3VectorData::operator=(vec); _t = vec._t; return *this; }

    inline const pxl::Basic4VectorData& operator+=(const pxl::Basic4VectorData& vec)
    { Basic3VectorData::operator+=(vec); _t += vec._t; return *this; }
    inline const pxl::Basic4VectorData& operator-=(const pxl::Basic4VectorData& vec)
    { Basic3VectorData::operator-=(vec); _t -= vec._t; return *this; }

  private:
    double _t;

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;      
};

// non-member operators
bool const operator==(const pxl::Basic4VectorData& obj1, const pxl::Basic4VectorData& obj2);
bool const operator!=(const pxl::Basic4VectorData& obj1, const pxl::Basic4VectorData& obj2);

} // namespace pxl

iotl__declareDataTypeProto(pxl::Basic4VectorData)

#endif // pxl_pol_Basic4Vector_hh
#ifndef pxl_pol_BasicObject_hh
#define pxl_pol_BasicObject_hh




#ifndef pxl_pol_types_hh
#define pxl_pol_types_hh


namespace pxl {

typedef CopyHistory& CopyHistoryRef;
typedef const CopyHistory& CopyHistoryConstRef;

typedef Index& IndexRef;
typedef const Index& IndexConstRef;

typedef Objects& ObjectsRef;
typedef const Objects& ObjectsConstRef;

typedef Relations& RelationsRef;
typedef const Relations& RelationsConstRef;

typedef pxl::Map<std::string, double> UserRecords;
typedef UserRecords& UserRecordsRef;
typedef const UserRecords& UserRecordsConstRef;

typedef pxl::Map<std::string, void*> CppPointers;
typedef CppPointers& CppPointersRef;
typedef const CppPointers& CppPointersConstRef;

} // namespace pxl

#endif // pxl_pol_types_hh

namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// pol

class BasicObjectData {
  public: 
    BasicObjectData() :
        _locked(0),
        _monteCarloMode(0),
        _name("default"),
        _status(0),
        _workflag(0),
        _userRecords(),
        _cppPointers()
    {}

    inline bool        getLocked()   const { return _locked; }
    inline int         getMonteCarloMode() const { return _monteCarloMode; }
    inline std::string getName()     const { return _name; }
    inline int         getStatus()   const { return _status; }
    inline int         getWorkflag() const { return _workflag; }

    inline void setLocked(bool v)      { _locked = v; }
    inline void setMonteCarloMode(int v) { _monteCarloMode = v; }
    inline void setName(std::string v) { _name = v; }
    inline void setStatus(int v)       { _status = v; }
    inline void setWorkflag(int v)     { _workflag = v; }

    inline const pxl::UserRecords& getUserRecords() const { return _userRecords; }
    inline const pxl::CppPointers& getCppPointers() const { return _cppPointers; }

    inline void setUserRecord(const std::string& key, double item) { _userRecords.set(key, item); }
    inline void setCppPointer(const std::string& key, void* cptr)  { _cppPointers.set(key, cptr); }

    inline void removeUserRecord(const std::string& key) { _userRecords.remove(key); }
    inline void removeCppPointer(const std::string& key) { _cppPointers.remove(key); }

    inline double findUserRecord(const std::string& key, double defaultitem) const
    { return _userRecords.find(key, defaultitem); }
    inline double findUserRecord(const std::string& key)                     const
    { return _userRecords.find(key); }

    template<class datatype>
    inline datatype* findCppPointer(const std::string& key) const
    { return dynamic_cast<datatype*>(_cppPointers.find(key, 0)); }

  protected:
    bool        _locked;
    int         _monteCarloMode;
    std::string _name;
    int         _status;
    int         _workflag;

    pxl::UserRecords _userRecords;
    pxl::CppPointers _cppPointers;

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;  
};

typedef pxl::CowObject<pxl::BasicObjectData> BasicObject;
typedef pxl::WkPtrSpec<pxl::BasicObjectData, pxl::BasicObject> BasicObjectWkPtr;
typedef BasicObject& BasicObjectRef;
typedef const BasicObject& BasicObjectConstRef;

} // namespace pxl

iotl__declareObjectTypeProto(pxl::BasicObject)
iotl__declareDataTypeProto(pxl::BasicObjectData)

#endif // pxl_pol_BasicObject_hh
#ifndef pxl_pol_BasicManager_hh
#define pxl_pol_BasicManager_hh





namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// pol

class BasicObjectManagerData : public pxl::BasicObjectData {
  public: 
    BasicObjectManagerData() :
        BasicObjectData(), _objects() {}
    BasicObjectManagerData(const pxl::BasicObjectManagerData& original) :
        BasicObjectData(original), _objects(original._objects) {}

    // create
    template<class datatype>
    datatype& create()
    { return _objects.create<datatype>(); }

    template<class datatype, class ctrdatatype>
    datatype& create(const ctrdatatype& ori)
    { return _objects.create<datatype,ctrdatatype>(ori); }

    // crateIndexed
    template<class datatype>
    datatype& createIndexed(const std::string& idx)
    { datatype& obj = _objects.create<datatype>(); setIndex(idx, obj); return obj; }

    template<class datatype, class ctrdatatype>
    datatype& createIndexed(const ctrdatatype& ori, const std::string& idx)
    { datatype& obj = _objects.create<datatype,ctrdatatype>(ori); setIndex(idx, obj); return obj; }


    inline void setObject(pxl::ObjectBase& obj, const std::string& idx)
    {_objects.set(obj); setIndex(idx, obj);}

    inline bool setIndex(const std::string& idx, pxl::ObjectBase& obj)
    { return _objects.setIndex(idx, obj); }

    inline const pxl::Objects& getObjects() const
    { return _objects; }

    inline void removeObject(pxl::ObjectBase& obj)
    { _objects.remove(obj); }

    inline void clearObjects()
    { _objects.clearContainer(); }


    template<class objecttype>
    inline objecttype* findObject(const std::string idx) const
    { return _objects.findObject<objecttype>(idx); }

    template<class objecttype>
    inline objecttype* findCopyOf(const pxl::ObjectBase& original) const
    { return _objects.findCopyOf<objecttype>(original); }


    inline const pxl::CopyHistory& getCopyHistory() const
    { return _objects.getCopyHistory(); }

    inline void clearCopyHistory()
    { _objects.clearCopyHistory(); }


    inline const pxl::Index& getIndex() const
    { return _objects.getIndex(); }

    inline void removeIndex(const std::string& idx)
    { _objects.removeIndex(idx); }

    inline void clearIndex()
    { _objects.clearIndex(); }

  protected:
    pxl::Objects _objects;

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;      
};

typedef pxl::Object<pxl::BasicObjectManagerData> BasicObjectManager; // NOTICE: BasicObjectManagers cannot be managed by Copy On Write!
typedef pxl::WkPtrSpec<pxl::BasicObjectManagerData, pxl::BasicObjectManager> BasicObjectManagerWkPtr;
typedef BasicObjectManager& BasicObjectManagerRef;
typedef const BasicObjectManager& BasicObjectManagerConstRef;

} // namespace pxl

iotl__declareObjectTypeProto(pxl::BasicObjectManager)
iotl__declareDataTypeProto(pxl::BasicObjectManagerData)

#endif // pxl_pol_BasicManager_hh
#ifndef pxl_pol_Vertex_hh
#define pxl_pol_Vertex_hh





namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// pol

class VertexData : public pxl::BasicObjectData {
  public: 
    inline const pxl::Basic3VectorData& vector()                const { return _vector; }
    inline const pxl::Basic3VectorData& vector(const pxl::Get&) const { return _vector; }
    inline       pxl::Basic3VectorData& vector(const pxl::Set&)       { return _vector; }

    inline const pxl::VertexData& operator+=(const pxl::VertexData& vec) { this->vector(pxl::set) += vec.vector(); return *this; }
    inline const pxl::VertexData& operator-=(const pxl::VertexData& vec) { this->vector(pxl::set) -= vec.vector(); return *this; }

  protected:
    pxl::Basic3VectorData _vector;

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;
};

// non-member operators
bool const operator==(const pxl::VertexData& obj1, const pxl::VertexData& obj2);
bool const operator!=(const pxl::VertexData& obj1, const pxl::VertexData& obj2);

// typedefs
/// data = >pxl::VertexData
typedef pxl::CowObject<pxl::VertexData> Vertex;
typedef pxl::WkPtrSpec<pxl::VertexData, pxl::Vertex> VertexWkPtr;
typedef Vertex& VertexRef;
typedef const Vertex& VertexConstRef;

template<> std::ostream& CowObject<pxl::VertexData>::print(int level, std::ostream& os, int pan) const;

} // namespace pxl

iotl__declareObjectTypeProto(pxl::Vertex)
iotl__declareDataTypeProto(pxl::VertexData)

#endif // pxl_pol_Vertex_hh
#ifndef pxl_pol_Particle_hh
#define pxl_pol_Particle_hh




namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// pol

class ParticleData : public pxl::BasicObjectData {
  public: 
    ParticleData() :
        BasicObjectData(), _vector(), _charge(0), _particleId(0) {}

    inline const pxl::Basic4VectorData& vector()                const { return _vector; }
    inline const pxl::Basic4VectorData& vector(const pxl::Get&) const { return _vector; }
    inline       pxl::Basic4VectorData& vector(const pxl::Set&)       { return _vector; }

    inline double getCharge() const { return _charge; }
    inline void setCharge(double v) { _charge = v; }

    inline int getParticleId() const { return _particleId; }
    inline void setParticleId(int v) { _particleId = v; }

    inline const pxl::ParticleData& operator+=(const pxl::ParticleData& vec)
    { this->vector(pxl::set) += vec.vector(); _charge += vec._charge; return *this; }

    inline const pxl::ParticleData& operator-=(const pxl::ParticleData& vec)
    { this->vector(pxl::set) -= vec.vector(); _charge += vec._charge; return *this; }

  protected:
    pxl::Basic4VectorData _vector;
    double                _charge;
    int                   _particleId;

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;      
};    

// non-member operators
bool const operator==(const pxl::ParticleData& obj1, const pxl::ParticleData& obj2);
bool const operator!=(const pxl::ParticleData& obj1, const pxl::ParticleData& obj2);

// typedefs
/// data => pxl::ParticleData
typedef pxl::CowObject<pxl::ParticleData> Particle;
typedef pxl::WkPtrSpec<pxl::ParticleData, pxl::Particle> ParticleWkPtr;
typedef Particle& ParticleRef;
typedef const Particle& ParticleConstRef;

template<> std::ostream& CowObject<pxl::ParticleData>::print(int level, std::ostream& os, int pan) const;

} // namespace pxl

iotl__declareObjectTypeProto(pxl::Particle)
iotl__declareDataTypeProto(pxl::ParticleData)

#endif // pxl_pol_Particle_hh
#ifndef pxl_pol_Collision_hh
#define pxl_pol_Collision_hh



namespace pxl {

typedef pxl::BasicObjectData CollisionData;

// typedefs
/// data => pxl::CollisionData
typedef pxl::CowObject<pxl::CollisionData> Collision;
typedef pxl::WkPtrSpec<pxl::CollisionData, pxl::Collision> CollisionWkPtr;
typedef Collision& CollisionRef;
typedef const Collision& CollisionConstRef;

template<> std::ostream& CowObject<pxl::CollisionData>::print(int level, std::ostream& os, int pan) const;

} // namespace pxl

#endif // pxl_pol_Collision_hh
#ifndef pxl_pol_EventView_hh
#define pxl_pol_EventView_hh




namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// pol

class EventViewData : public pxl::BasicObjectManagerData {
  public: 
    EventViewData()
        : BasicObjectManagerData() {}
    EventViewData(const pxl::EventViewData& original)
        : BasicObjectManagerData(original) {}

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;      
};

/// data => pxl::EventViewData
typedef pxl::Object<pxl::EventViewData> EventView; // NOTICE: EventViews cannot be managed by Copy On Write!
typedef pxl::WkPtrSpec<pxl::EventViewData, pxl::EventView> EventViewWkPtr;
typedef EventView& EventViewRef;
typedef const EventView& EventViewConstRef;

template<> std::ostream& Object<pxl::EventViewData>::print(int level, std::ostream& os, int pan) const;

} // namespace pxl

iotl__declareObjectTypeProto(pxl::EventView)
iotl__declareDataTypeProto(pxl::EventViewData)

#endif // pxl_pol_EventView_hh
#ifndef pxl_pol_AnalysisProcess_hh
#define pxl_pol_AnalysisProcess_hh




namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// pol

class AnalysisProcessData : public pxl::BasicObjectManagerData {
  public:
    AnalysisProcessData() :
        BasicObjectManagerData() {}
    AnalysisProcessData(const pxl::AnalysisProcessData& original) :
        BasicObjectManagerData(original) {}

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;  
};

/// data=>pxl::AnalysisProcessData
class AnalysisProcess : public pxl::Object<pxl::AnalysisProcessData> {
	// NOTICE: AnalysisProcesses cannot be managed by Copy On Write!
  public: 
    AnalysisProcess() :
        pxl::Object<pxl::AnalysisProcessData>() {}
    AnalysisProcess(const pxl::AnalysisProcess& original) :
        pxl::Object<pxl::AnalysisProcessData>(original) {}
    virtual ~AnalysisProcess() {}

    virtual void buildTemplate(int mode = 0) {}

    virtual void beginJob(const pxl::Objects* input = 0) {}
    virtual void beginRun(const pxl::Objects* input = 0) {}
    virtual void analyseEvent(const pxl::Objects* input = 0) {}
    virtual void finishEvent(const pxl::Objects* input = 0) { set().clearObjects(); }
    virtual void endRun(const pxl::Objects* input = 0) {}
    virtual void endJob(const pxl::Objects* input = 0) {}

    template<class objecttype>
    const objecttype& castInput(const pxl::Objects* input)
    {
        const objecttype* obj = dynamic_cast<const objecttype*>(input);
        if (!obj)
            // FIXME: cerr vs exception vs possibly other rant
            std::cerr << "pxl::AnalysisProcess::castInput(): FATAL: The pointer you intend to cast does not exist!" << std::endl;
        return *obj;
    }

    virtual pxl::ObjectBase* clone() const;

  protected:
    virtual pxl::WkPtrBase* createSelfWkPtr();
    virtual void storeYourSelf(pxl::oStreamer& output) const;
};

typedef pxl::WkPtrSpec<pxl::AnalysisProcessData, pxl::AnalysisProcess> AnalysisProcessWkPtr;
typedef AnalysisProcess& AnalysisProcessRef;
typedef const AnalysisProcess& AnalysisProcessConstRef;

template<> std::ostream& Object<pxl::AnalysisProcessData>::print(int level, std::ostream& os, int pan) const;

} // namespace pxl

iotl__declareObjectTypeProto(pxl::AnalysisProcess)
iotl__declareDataTypeProto(pxl::AnalysisProcessData)

#endif // pxl_pol_AnalysisProcess_hh
#ifndef pxl_pol_AnalysisFork_hh
#define pxl_pol_AnalysisFork_hh





namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// pol

class AnalysisForkData : public pxl::BasicObjectManagerData {
  public: 
    AnalysisForkData() :
        BasicObjectManagerData() {}
    AnalysisForkData(const pxl::AnalysisForkData& original) :
        BasicObjectManagerData(original) {}

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;
};

/// data = >pxl::AnalysisForkData
class AnalysisFork : public pxl::Object<pxl::AnalysisForkData> {
	// NOTICE: AnalysisForks cannot be managed by Copy On Write!
  public: 
    AnalysisFork() :
        pxl::Object<pxl::AnalysisForkData>() {}
    AnalysisFork(const pxl::AnalysisFork& original) :
        pxl::Object<pxl::AnalysisForkData>(original) {}
    virtual ~AnalysisFork() {}

    virtual void buildTemplate(int mode = 0);

    virtual void beginJob(const pxl::Objects*input = 0);
    virtual void beginRun(const pxl::Objects*input = 0);
    virtual void analyseEvent(const pxl::Objects*input = 0);
    virtual void finishEvent(const pxl::Objects*input = 0);
    virtual void endRun(const pxl::Objects*input = 0);
    virtual void endJob(const pxl::Objects*input = 0);

    virtual pxl::ObjectBase* clone() const;

  protected:
    virtual pxl::WkPtrBase* createSelfWkPtr();
    virtual void storeYourSelf(pxl::oStreamer& output) const;
};

typedef pxl::WkPtrSpec<pxl::AnalysisForkData, pxl::AnalysisFork> AnalysisForkWkPtr;
typedef AnalysisFork& AnalysisForkRef;
typedef const AnalysisFork& AnalysisForkConstRef;

template<> std::ostream& Object<pxl::AnalysisForkData>::print(int level, std::ostream& os, int pan) const;

} // namespace pxl

iotl__declareObjectTypeProto(pxl::AnalysisFork)
iotl__declareDataTypeProto(pxl::AnalysisForkData)

#endif // pxl_pol_AnalysisFork_hh

#endif // pxl_pol_hh

#endif // pxl_hh
//
//----------------------------------------------------------------------
extern const pxl::Get ePaxGet;
extern const pxl::Set ePaxSet;
//----------------------------------------------------------------------
typedef pxl::Objects ePaxObjects;
typedef ePaxObjects& ePaxObjectsRef;
typedef const ePaxObjects& ePaxObjectsConstRef;
//----------------------------------------------------------------------
typedef pxl::Relations ePaxRelations;
typedef ePaxRelations& ePaxRelationsRef;
typedef const ePaxRelations& ePaxRelationsConstRef;
//----------------------------------------------------------------------
typedef pxl::CopyHistory ePaxCopyHistory;
typedef ePaxCopyHistory& ePaxCopyHistoryRef;
typedef const ePaxCopyHistory& ePaxCopyHistoryConstRef;
//----------------------------------------------------------------------
typedef pxl::Index ePaxIndex;
typedef ePaxIndex& ePaxIndexRef;
typedef const ePaxIndex& ePaxIndexConstRef;
//----------------------------------------------------------------------
typedef pxl::UserRecords ePaxUserRecords;
typedef ePaxUserRecords& ePaxUserRecordsRef;
typedef const ePaxUserRecords& ePaxUserRecordsConstRef;
//----------------------------------------------------------------------
typedef pxl::CppPointers ePaxCppPointers;
typedef ePaxCppPointers& ePaxCppPointersRef;
typedef const ePaxCppPointers& ePaxCppPointersConstRef;
//----------------------------------------------------------------------
// ePAXify pxl::AnalysisFork
/// object=>pxl::AnalysisFork data=>pxl::AnalysisForkData
class ePaxAnalysisFork : public pxl::AnalysisFork {
  public: 
    // Gero-style service methods go here...
//     template<class datatype> datatype& create(const std::string idx)                                                 { return set().create<datatype>(idx); }
//     template<class datatype> datatype& create(const datatype& original, const std::string idx)                       { return set().create<datatype>(original, idx); }
//     template<class datatype, class ctrdatatype> datatype& create(const ctrdatatype& original, const std::string idx) { return set().create<datatype,ctrdatatype>(original, idx); }
//     
//     inline void setUserRecord(const std::string& key, double item) { set().setUserRecord(key, item); }
//     inline void setCppPointer(const std::string& key, void* cptr)  { set().setCppPointer(key, cptr); }
//     ...


    virtual pxl::ObjectBase* clone() const;
    virtual std::ostream&  print(int level = 0, std::ostream& os = std::cout, int pan = 0) const;
  protected:    
    virtual pxl::WkPtrBase* createSelfWkPtr();
    virtual void storeYourSelf(pxl::oStreamer& output) const;
};
typedef pxl::WkPtrSpec<pxl::AnalysisForkData, ePaxAnalysisFork> ePaxAnalysisForkWkPtr;
typedef ePaxAnalysisFork& ePaxAnalysisForkRef;
typedef const ePaxAnalysisFork& ePaxAnalysisForkConstRef;
//----------------------------------------------------------------------
// ePAXify pxl::AnalysisProcess
/// object=>pxl::AnalysisProcess data=>pxl::AnalysisProcessData
class ePaxAnalysisProcess : public pxl::AnalysisProcess {
  public: 
    // Gero-style service methods go here...
    // ...

    virtual pxl::ObjectBase* clone() const;
    virtual std::ostream&  print(int level = 0, std::ostream& os = std::cout, int pan = 0) const;
  protected:    
    virtual pxl::WkPtrBase* createSelfWkPtr();
    virtual void storeYourSelf(pxl::oStreamer& output) const;
};
typedef pxl::WkPtrSpec<pxl::AnalysisProcessData, ePaxAnalysisProcess> ePaxAnalysisProcessWkPtr;
typedef ePaxAnalysisProcess& ePaxAnalysisProcessRef;
typedef const ePaxAnalysisProcess& ePaxAnalysisProcessConstRef;
//----------------------------------------------------------------------
// ePAXify pxl::EventView
/// object=>pxl::EventView data=>pxl::EventViewData
class ePaxEventView : public pxl::EventView {
  public: 
    // Gero-style service methods go here...
    // ...

    virtual pxl::ObjectBase* clone() const;
    virtual std::ostream&  print(int level = 0, std::ostream& os = std::cout, int pan = 0) const;
  protected:    
    virtual pxl::WkPtrBase* createSelfWkPtr();
    virtual void storeYourSelf(pxl::oStreamer& output) const;
};
typedef pxl::WkPtrSpec<pxl::EventViewData, ePaxEventView> ePaxEventViewWkPtr;
typedef ePaxEventView& ePaxEventViewRef;
typedef const ePaxEventView& ePaxEventViewConstRef;
//----------------------------------------------------------------------
// ePAXify pxl::Particle
/// object=>pxl::Particle data=>pxl::ParticleData
class ePaxParticle : public pxl::Particle {
  public: 
    // Gero-style service methods go here...
    // ...

    virtual pxl::ObjectBase* clone() const;
    virtual std::ostream&  print(int level = 0, std::ostream& os = std::cout, int pan = 0) const;
  protected:    
    virtual pxl::WkPtrBase* createSelfWkPtr();
    virtual void storeYourSelf(pxl::oStreamer& output) const;
};
typedef pxl::WkPtrSpec<pxl::ParticleData, ePaxParticle> ePaxParticleWkPtr;
typedef ePaxParticle& ePaxParticleRef;
typedef const ePaxParticle& ePaxParticleConstRef;
//----------------------------------------------------------------------
// ePAXify pxl::Vertex
/// object=>pxl::Vertex data=>pxl::VertexData
class ePaxVertex : public pxl::Vertex {
  public: 
    // Gero-style service methods go here...
    // ...

    virtual pxl::ObjectBase* clone() const;
    virtual std::ostream&  print(int level = 0, std::ostream& os = std::cout, int pan = 0) const;
  protected:    
    virtual pxl::WkPtrBase* createSelfWkPtr();
    virtual void storeYourSelf(pxl::oStreamer& output) const;
};
typedef pxl::WkPtrSpec<pxl::VertexData, ePaxVertex> ePaxVertexWkPtr;
typedef ePaxVertex& ePaxVertexRef;
typedef const ePaxVertex& ePaxVertexConstRef;
//----------------------------------------------------------------------
// ePAXify pxl::Collision
/// object=>pxl::Collision data=>pxl::CollisionData
class ePaxCollision : public pxl::Collision {
  public: 
    // Gero-style service methods go here...
    // ...

    virtual pxl::ObjectBase* clone() const;
    virtual std::ostream&  print(int level = 0, std::ostream& os = std::cout, int pan = 0) const;
  protected:    
    virtual pxl::WkPtrBase* createSelfWkPtr();
    virtual void storeYourSelf(pxl::oStreamer& output) const;
};
typedef pxl::WkPtrSpec<pxl::CollisionData, ePaxCollision> ePaxCollisionWkPtr;
typedef ePaxCollision& ePaxCollisionRef;
typedef const ePaxCollision& ePaxCollisionConstRef;
//----------------------------------------------------------------------
iotl__declareObjectTypeProto(ePaxParticle)
iotl__declareObjectTypeProto(ePaxVertex)
iotl__declareObjectTypeProto(ePaxCollision)
iotl__declareObjectTypeProto(ePaxEventView)
iotl__declareObjectTypeProto(ePaxAnalysisProcess)
iotl__declareObjectTypeProto(ePaxAnalysisFork)
//----------------------------------------------------------------------
#endif // epax_hh

