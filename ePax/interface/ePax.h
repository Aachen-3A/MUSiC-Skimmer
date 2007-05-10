//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_hh
#define pxl_hh

//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_pcl_hh
#define pxl_pcl_hh

//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_pcl_macros_hh
#define pxl_pcl_macros_hh

#include <cstddef>

#undef PXL_LIKELY
#undef PXL_UNLIKELY
#if __GNUC__ > 2 || (__GNUC__ == 2 && __GNUC_MINOR__ > 4)
#	define PXL_NORETURN		__attribute__((__noreturn__))
#	define PXL_LIKELY(expr)		(__builtin_expect((bool)(expr), true))
#	define PXL_UNLIKELY(expr)	(__builtin_expect((bool)(expr), false))
#else
#	define PXL_NORETURN
#	define PXL_LIKELY(expr)		(expr)
#	define PXL_UNLIKELY(expr)	(expr)
#endif

#ifdef offsetof
#	define PXL_OFFSETOF(t, f)	((std::ptrdiff_t)offsetof(t, f))
#else
#	define PXL_OFFSETOF(t, f)	((std::ptrdiff_t)((char *) &((t*)0)->f))
#endif

#define PXL_BASE(t, f, v)		(reinterpret_cast<t*>(reinterpret_cast<char*>(v) - PXL_OFFSETOF(t, f)))

#endif // pxl_pcl_macros_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_pcl_Exception_hh
#define pxl_pcl_Exception_hh

#include <string>

namespace pxl {

enum Particles { Electron = 1, Muon = 2, Gamma = 3, MET = 4, KtJet = 5, ItCone5Jet = 7, MidCone5Jet = 8, MidCone7Jet = 9 };

/**
This class is provided for PXL exception and error handling.
*/
class Exception {
  public: 
    /// Constructor indicating unspecified routine and message.
    Exception() :
        routine("unspecified routine"),
        message("unspecified error") {}

    /// Standard constructor, routine and message can be specified as arguments.
    Exception(const std::string& routine, const std::string& message) :
        routine(routine),
        message(message) {}

    ~Exception() {} 

    /// This method returns the routine the exception appeared in.
    const std::string &getRoutine() const { return routine; }
    /// This method returns a clear text message describing the exception.
    const std::string &getMessage() const { return message; }

  private:
    std::string routine;
    std::string message;
};

} // namespace pxl

#endif // pxl_pcl_Exception_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_pcl_BasicIoStreamer_hh
#define pxl_pcl_BasicIoStreamer_hh

#include <iostream>

//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_pcl_BasicLinuxIoStreamer_hh
#define pxl_pcl_BasicLinuxIoStreamer_hh


/*
 * FIXME / General comment:
 * BasicLinuxIoStreamer seems generic enough, not at all Linux specific?
 */

namespace pxl {

/**
This internal class provides binary streaming of basic data types on a Linux platform.  
*/ 
class BasicLinuxIoStreamer {
  protected:
    // writing data
    /// This methods stores a number of bytes \p bytes at address \p address to the \p cxxx stream.
    inline void dumpMemory(std::ostream& cxxx, const char* address, int bytes)
    {
        cxxx.rdbuf()->sputn(address, bytes);
        cxxx.rdbuf()->pubsync();
    }

    /// This methods stores a char \p data in a device-independent representation to the \p cxxx stream.
    inline void storeBasicTypeChar(std::ostream& cxxx, char data)
    {
        dumpMemory(cxxx, (const char*)&data, 1); 
    }

    /// This methods stores a bool \p data in a device-independent representation to the \p cxxx stream.
    inline void storeBasicTypeBool(std::ostream& cxxx, bool data)
    {
        if (data) dumpMemory(cxxx, "Y", 1);
        else dumpMemory(cxxx, "N", 1);
    }

    /// This methods stores an int \p data in a device-independent representation to the \p cxxx stream.
    inline void storeBasicTypeInt(std::ostream& cxxx, int data)
    {
        dumpMemory(cxxx, (const char*)&data, 4); // FIXME: endian/64bit
    }

    /// This methods stores a float \p data in a device-independent representation to the \p cxxx stream.
    inline void storeBasicTypeFloat(std::ostream& cxxx, float data)
    {
        dumpMemory(cxxx, (const char*)&data, 4); // FIXME: endian/64bit
    }

    /// This methods stores a double \p data in a device-independent representation to the \p cxxx stream.
    inline void storeBasicTypeDouble(std::ostream& cxxx, double data)
    {
        dumpMemory(cxxx, (const char*)&data, 8); // FIXME: endian/64bit
    }

    /// This methods stores a C string \p data in a device-independent representation to the \p cxxx stream.
    inline void storeBasicTypeCStr(std::ostream& cxxx, const char* address)
    {
        for(; (*address) != '\0'; address++) // FIXME: performance!!!
      	    cxxx.rdbuf()->sputc(*address);   // better: len, data ? (as below)
        cxxx.rdbuf()->sputc('\0');
    }

    /// This methods stores a C++ std::string \p data in a device-independent representation to the \p cxxx stream.
    inline void storeBasicTypeString(std::ostream& cxxx, const std::string& data)
    {
        int length = data.length();
        storeBasicTypeInt(cxxx, length);
        dumpMemory(cxxx, data.c_str(), length);
    }
  
    // reading data
    /// This methods restores a number of \p bytes bytes from the \p cxxx stream to the address \p address.
    inline bool restoreMemory(std::istream& cxxx, char* address, int bytes)
    {
        cxxx.read(address, bytes);
        return bytes == cxxx.gcount();
    }

    /// This methods restores a char \p data in a device-independent representation from the \p cxxx stream.
    inline bool restoreBasicTypeChar(std::istream& cxxx, char& data)
    {
        return restoreMemory(cxxx, (char*)&data, 1);
    }

    /// This methods restores a bool \p data in a device-independent representation from the \p cxxx stream.
    inline bool restoreBasicTypeBool(std::istream& cxxx, bool& data)
    {
      char cYesNo = ' ';
      bool success = restoreMemory(cxxx, &cYesNo, 1);
      data = (cYesNo == 'Y');
      return success;
    }

    /// This methods restores a int \p data in a device-independent representation from the \p cxxx stream.
    inline bool restoreBasicTypeInt(std::istream& cxxx, int& data)
    {
        return restoreMemory(cxxx, (char*)&data, 4); // FIXME: endian/64bit
    }

    /// This methods restores a char \p data in a device-independent representation from the \p cxxx stream.
    inline bool restoreBasicTypeFloat(std::istream& cxxx, float& data)
    {
        return restoreMemory(cxxx, (char*)&data, 4); // FIXME: endian/64bit
    }

    /// This methods restores a double \p data in a device-independent representation from the \p cxxx stream.
    inline bool restoreBasicTypeDouble(std::istream& cxxx, double& data)
    {
        return restoreMemory(cxxx, (char*)&data, 8); // FIXME: endian/64bit
    }

    /// This methods restores a C string \p data in a device-independent representation from the \p cxxx stream.
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

    /// This methods restores a C++ std::string \p data in a device-independent representation from the \p cxxx stream.
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

/// This typedef is for platform selection (currently: pxl::BasicLinuxIoStreamer)
typedef BasicLinuxIoStreamer BasicIoStreamer;

} // namespace pxl

#endif // pxl_pcl_BasicIoStreamer_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_pcl_VariantBase_hh
#define pxl_pcl_VariantBase_hh

#include <cstdlib>
#include <vector>


namespace pxl {

/** 
This base class serves the PXL variant data type. 
It does not contain the information about the type store, this task
is delegated to the caller. This way memory can be saved, if, for instance,
table-like structures are to be stored, where the type information is stored
in a separate header structure.
*/
// Important notice: pxl::VariantBase is not a fully opaque type, you have to
// take care of calling clear() and dup() after copying or freeing yourself!
class VariantBase {
  public:
    /// This enum represents the possible value types that can be stored by the PXL variant data type.
    enum Type {
        TYPE_NULL = 0,
        TYPE_BOOL,
        TYPE_CHAR,
        TYPE_UCHAR,
        TYPE_SHORT,
        TYPE_USHORT,
        TYPE_INT,
        TYPE_UINT,
        TYPE_LONG,
        TYPE_ULONG,
        TYPE_FLOAT,
        TYPE_DOUBLE,
        TYPE_STRING,
        TYPE_PTR,
        TYPE_USER
    };

    inline VariantBase() { std::memset(&v, 0, sizeof v); }
    inline VariantBase(const pxl::VariantBase& orig) : v(orig.v) {}

    /// This method template returns the variant content and tests for compatibility of the template expansion with the type given by \p t.
    template<typename T>
    inline T get(Type t) const;

    /// This method template sets the variant content and tests for compatibility of the template expansion with the type given by \p t.
    template<typename T>
    inline void set(Type t, T arg);

    /// This method clears the variant content and releases possible allocations, depending on the type given by \p t.
    inline void clear(Type t)
    { 
        if (t == TYPE_STRING)
            delete[] (char*)v.p;
        else if (PXL_UNLIKELY(t >= TYPE_USER)) {
            const TypeInfo& tInfo = getTypeInfo(t);
            if (tInfo.clear)
                tInfo.clear(&v);
        }
        std::memset(&v, 0, sizeof v);
    }

    /// This method ensures correct duplication of possible internal memory allocations after a copy of the variant, depending on the type given by \p t.
    inline void dup(Type t)
    { 
        if (t == TYPE_STRING) {
            std::size_t len = std::strlen((char*)v.p) + 1;
            char* copy = new char[len];
            std::memcpy(copy, v.p, len);
            v.p = copy;
        } else if (PXL_UNLIKELY(t >= TYPE_USER)) {
            const TypeInfo& tInfo = getTypeInfo(t);
            if (tInfo.dup)
                tInfo.dup(&v);
        }
    }

    /// This method template returns a suitable type id for the used template instantiation
    template<typename T>
    static Type findType();

  protected:
    /// This union serves as actual value storage for the PXL variant data type.
    union Value {
        bool           b;
        char           c;
        unsigned char  uc;
        short          s;
        unsigned short us;
        int            i;
        unsigned int   ui;
        long           l;
        unsigned long  ul;
        float          f;
        double         d;
        void           *p;
    } v;

    /// This subclass serves as internal type information the PXL variant data type.
    class TypeInfo {
      public:
        inline TypeInfo(const char* name) :
            clear(0), dup(0), name(name) {}

        void        (*clear)(Value* v);
        void        (*dup)(Value* v);
        const char* name;
    };

    /// This method returns internal type information for the type given by \p t .
    static inline const TypeInfo& getTypeInfo(Type t)
    {
        if (PXL_UNLIKELY((std::size_t)t >= types.size()))
            return fallbackGetTypeInfo(t);
        return types[t];
    }

    /// This method throws a detailed exception about a type mismatch between \p tShould and \p tIs .
    void wrongType(Type tShould, Type tIs) const;

  private:
    static const TypeInfo& fallbackGetTypeInfo(Type t);

    static std::vector<TypeInfo> types;
};

#define PCL_ANYBASE_CHECK(t_, t) \
if (TYPE_##t_ != t) \
    wrongType(TYPE_##t_, t);

#define PCL_ANYBASE_SIMPLE(t_, type, idx) \
template<> \
inline type VariantBase::get<type>(Type t) const \
{ \
    PCL_ANYBASE_CHECK(t_, t) \
    return v.idx; \
} \
\
template<> \
inline void VariantBase::set<type>(Type t, type arg) \
{ \
    PCL_ANYBASE_CHECK(t_, t) \
        v.idx = arg; \
} \
\
template<> \
inline VariantBase::Type VariantBase::findType<type>() \
{ return TYPE_##t_; }

PCL_ANYBASE_SIMPLE(BOOL, bool, b)
PCL_ANYBASE_SIMPLE(CHAR, char, c)
PCL_ANYBASE_SIMPLE(UCHAR, unsigned char, uc)
PCL_ANYBASE_SIMPLE(SHORT, short, s)
PCL_ANYBASE_SIMPLE(USHORT, unsigned short, us)
PCL_ANYBASE_SIMPLE(INT, int, i)
PCL_ANYBASE_SIMPLE(UINT, unsigned int, ui)
PCL_ANYBASE_SIMPLE(LONG, long, l)
PCL_ANYBASE_SIMPLE(ULONG, unsigned long, ul)
PCL_ANYBASE_SIMPLE(FLOAT, float, f)
PCL_ANYBASE_SIMPLE(DOUBLE, double, f)
PCL_ANYBASE_SIMPLE(PTR, void*, p)

#undef PCL_ANYBASE_SIMPLE

// specializations for complex data types

// TYPE_STRING

template<>
inline std::string VariantBase::get(Type t) const
{
    PCL_ANYBASE_CHECK(STRING, t)
    return std::string((const char*)v.p);
}

template<>
inline void VariantBase::set(Type t, std::string arg)
{
    PCL_ANYBASE_CHECK(STRING, t)
    if (v.p)
        delete[] (char*)v.p;
    std::string::size_type size = arg.size();
    v.p = new char[size + 1];
    std::memcpy(v.p, arg.c_str(), size);
    ((char*)v.p)[size] = 0;
}

template<>
inline VariantBase::Type VariantBase::findType<std::string>()
{ return TYPE_STRING; }

#undef PCL_ANYBASE_CHECK

} // namespace pxl

#endif // pxl_pcl_VariantBase_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_pcl_Variant_hh
#define pxl_pcl_Variant_hh


namespace pxl {

/** 
This class represents the PXL variant data type.
It carries a generic data holder along with type information about
the type of the currently stored value.
*/
class Variant : public VariantBase {
  public:
    inline Variant() : type(TYPE_NULL) {}
    inline Variant(Type t) : type(t) {}
    inline Variant(const pxl::Variant& orig) :
    	VariantBase(orig), type(orig.type) { dup(type); }
    inline ~Variant() { if (type != TYPE_NULL) VariantBase::clear(type); }

    inline pxl::Variant& operator=(const pxl::Variant& orig)
    {
        if (type != TYPE_NULL) {
            VariantBase::clear(type);
            type = orig.type;
        } else if (PXL_UNLIKELY(type != orig.type))
            wrongType(orig.type, type);
        v = orig.v;
        dup(type);
        return *this;
    }

    /// This method returns the type id of the currently stored value.
    inline Type getType() const { return type; }
    /// This method returns the type name of the currently stored value.
    inline const char *getTypeName() const { return getTypeInfo(type).name; }

    /// This method sets the type of the value to \p t and throws a pxl::Exception if a type is already assigned.
    inline void setType(Type t)
    {
        if (PXL_UNLIKELY(type != TYPE_NULL))
            wrongType(TYPE_NULL, type);
        type = t;
    }

    /// This method clears the contents and type of the currently stored value.
    inline void clear()
    { VariantBase::clear(type); type = TYPE_NULL; }

    /// This method initialises the variant with an empty value of the type that matches the template instantiation.
    template<typename datatype>
    inline void init()
    { setType(findType<datatype>()); }

    /// This method returns the value of the variant or throws a pxl::Exception if the currently stored type doesn't match the template instantiation.
    template<typename datatype>
    inline datatype get() const
    { return VariantBase::get<datatype>(type); }

    /// This method set the value of the variant or throws a pxl::Exception if the currently assigned type doesn't match the template instantiation.
    template<typename datatype>
    inline void set(datatype arg)
    { VariantBase::set<datatype>(type, arg); }

  protected:
    Type type;
};

} // namespace pxl

#endif // pxl_pcl_Variant_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_pcl_functions_hh
#define pxl_pcl_functions_hh


namespace pxl {

/// @internal This function returns a platform-specific CPU timestamp; internally used for performance tests.
double getCpuTime();

// cheaper to call with <const char*> arguments, no inline string construction
/// This function can be used to throw a pxl::Exception with routine and message information.
void exception(const char* routine, const char* message);

// <std::string>.c_str() is inlined as a simple pointer dereference, cheaper
/// This function can be used to throw a pxl::Exception with routine and message information.
inline void exception(const std::string& routine, const char* message)
{ pxl::exception(routine.c_str(), message); }
/// This function can be used to throw a pxl::Exception with routine and message information.
inline void exception(const char *routine, const std::string& message)
{ pxl::exception(routine, message.c_str()); }
/// This function can be used to throw a pxl::Exception with routine and message information.
inline void exception(const std::string& routine, const std::string& message)
{ pxl::exception(routine.c_str(), message.c_str()); }

} // namespace pxl

#endif // pxl_pcl_functions_hh

#endif // pxl_pcl_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_hh
#define pxl_ptl_hh

//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_MutableId_hh
#define pxl_ptl_MutableId_hh

namespace pxl {

/// This typedef is intended to provide a data type for the PXL unique object-id
typedef int MutableId;

} // namespace pxl

#endif // pxl_ptl_MutableId_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_Id_hh
#define pxl_ptl_Id_hh


namespace pxl {

/// This typedef is intended to provide a read-only data type for the PXL unique object-id
typedef const pxl::MutableId Id;

} // namespace pxl

#endif // pxl_ptl_Id_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_Ptr_hh
#define pxl_ptl_Ptr_hh


namespace pxl {

// ptl

/** 
This class template implements a PXL style dumb C++ pointer to \p datatype; 
it is intended to parallelize pointer-related analysis code for PXL objects and 
user-defined C++ objects.   
*/
template<class datatype>
class Ptr {
  public:
    Ptr() : _cppPtr(0) {}
    Ptr(const pxl::Ptr<datatype>& original) : _cppPtr(original._cppPtr) {}
    Ptr(datatype* original) : _cppPtr(original) {}
    ~Ptr() {}

    /// This method returns the C++ pointer to the referenced object  
    inline datatype* pointer() const  { return _cppPtr; }
    /// This method returns true, if the pointer is not a null pointer; 
    /// CAUTION: valid() cannot check, if the referenced object (still) exists.   
    inline bool valid() const         { return _cppPtr != 0; }

    /// This method provides direct access to the referenced object;
    /// CAUTION: object() will not check, if the referenced object (still) exists.   
    inline datatype& object() const   { return *access(); }

    /// This assignment operator causes the pointer to reference the C++ object referenced by \p pptr.
    inline void operator=(const pxl::Ptr<datatype>& pptr) { _cppPtr = pptr._cppPtr; }
    /// This assignment operator causes the pointer to reference the C++ object \p data.
    inline void operator=(datatype& data)                 { _cppPtr = &data; }
    /// This assignment operator causes the pointer to reference the C++ object pointed to by \p dataptr.
    inline void operator=(datatype* dataptr)              { _cppPtr = dataptr; }

    /// This arrow operator de-references the pointer.   
    /// CAUTION: it will not be checked, if the referenced object (still) exists.   
    inline datatype* operator->() const { return _cppPtr; }

  protected:
    // access to object
    inline datatype* access() const
    {
        if (_cppPtr)
            return _cppPtr; 
        std::cerr << "pxl::Ptr::access(): FATAL: The object you intend to access does not exist!" << std::endl;
        return 0; // FIXME: excexption?
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
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_Vector_hh
#define pxl_ptl_Vector_hh


namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// ptl

/** 
This class template aggregates a STL vector and extends its functionalites. 
*/
template<class itemtype>
class Vector {
  public:
    virtual ~Vector() {}

    /// For STL-style access: std::vector<> data type of the aggregated STL vector
    typedef          std::vector<itemtype>                 StlContainer;
    /// For STL-style iteration: const_iterator data type of the aggregated STL vector 
    typedef typename std::vector<itemtype>::const_iterator StlConstIterator;
    /// For STL-style iteration: iterator data type of the aggregated STL vector 
    typedef typename std::vector<itemtype>::iterator       StlIterator;

    /// This method inserts \p item.
    void pushBack(const itemtype& item) {_container.push_back(item);}

    // navigation
    /// For STL-style iteration: returns the STL begin iterator of the aggregated STL vector. 
    const StlConstIterator begin() const { return _container.begin(); }
    /// For STL-style iteration: returns the STL end iterator of the aggregated STL vector. 
    const StlConstIterator end()   const { return _container.end(); }

    /// This method grants access to the aggregated STL vector. 
    inline StlContainer& getContainer() { return _container; }
    /// This method clears the aggregated STL vector. 
    virtual void clearContainer() { return _container.clear(); }

    /// This method returns the number of items contained in the aggregated STL vector. 
    inline int getSize() const { return _container.size(); }

    /// For PXL-style iteration: PXL iterator class 
    class Iterator {
      public:
        /// Go-to constructor moves iterator to the same positon as \p original. 
        Iterator(const Iterator& original) :
            _iter(original._iter), _containerRef(original._containerRef) {}

        /// Standard constructor moves iterator to the beginning of \p vector. 
        Iterator(const pxl::Vector<itemtype>& vector) :
            _containerRef(&vector) { first(); }

        /// This method moves the iterator to the beginning of the vector. 
        inline void first()  { _iter = _containerRef->begin(); }
        /// This method moves the iterator to the next pair in the vector. 
        inline void next()   { _iter++; }
         /// This method returns true, if the iterator hit upper or lower boundary of the vector. 
        inline bool isDone() { return _iter == _containerRef->end(); }

        /// This method returns the current item. 
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
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_Map_hh
#define pxl_ptl_Map_hh

#include <map>


namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// pcl

/** 
This class template aggregates a STL map and extends its functionalites; 
a map is a list of \p keytype-itemtype pairs sorted by a unique key of \p keytype. 
*/
template<class keytype, class itemtype>
class Map {
  public:
    virtual ~Map() {}

    /// For STL-style access: std::pair<> data type of the aggregated STL map
    typedef          std::pair<keytype, itemtype>                StlPair;
    /// For STL-style access: std::map<> data type of the aggregated STL map
    typedef          std::map<keytype, itemtype>                 StlContainer;
    /// For STL-style iteration: const_iterator data type of the aggregated STL map 
    typedef typename std::map<keytype, itemtype>::const_iterator StlConstIterator;
    /// For STL-style iteration: iterator data type of the aggregated STL map 
    typedef typename std::map<keytype, itemtype>::iterator       StlIterator;

    /// This method inserts (or replaces) the \p key-item pair indetified by \p key.
    void set(const keytype& key, const itemtype& item)
    {
        StlIterator insertPos = _container.lower_bound(key);
        if (insertPos == _container.end() || insertPos->first != key)
            _container.insert(insertPos, StlPair(key, item));
        else
            insertPos->second = item;
    }

    /// This method removes (if existing) the \p key-item pair indetified by \p key. 
    void remove(const keytype& key)
    { _container.erase(key); }

    /// This method searches and returns the item indetified by \p key; \p defaultitem is returned in case the key is not found. 
    itemtype find(const keytype& key, itemtype defaultitem) const
    {
        StlConstIterator found = _container.find(key);
        if (found != _container.end())
            return found->second;
        return defaultitem;
    }

    /// This method searches and returns the item indetified by \p key; a pxl::Exception is thrown in case the key is not found. 
    itemtype find(const keytype& key) const
    {
        StlConstIterator found = _container.find(key);
        if (found != _container.end())
            return found->second;
        pxl::exception("pxl::Map::find(...)", "key not found and no default item provided");
        throw;
    }

  protected:
    itemtype& findOrAlloc(const keytype &key)
    {
        StlIterator insertPos = _container.lower_bound(key);
        if (insertPos == _container.end() || insertPos->first != key)
            return _container.insert(insertPos, StlPair(key, itemtype()))->second;
        else
            return insertPos->second;
    }

    const itemtype* findOrReturn(const keytype &key) const
    {
        StlConstIterator found = _container.find(key);
        if (found == _container.end())
           return 0;
        return &found->second;
    }

  public:
    // navigation
    /// For STL-style iteration: returns the STL begin iterator of the aggregated STL map. 
    const StlConstIterator begin() const { return _container.begin(); }
    /// For STL-style iteration: returns the STL end iterator of the aggregated STL map. 
    const StlConstIterator end()   const { return _container.end(); }

    /// This method grants access to the aggregated STL map. 
    inline StlContainer& getContainer() { return _container; }
    /// This method clears the aggregated STL map. 
    virtual void clearContainer() { return _container.clear(); }

    /// This method returns the number of pairs contained in the aggregated STL map. 
    inline int getSize() const { return _container.size(); }

    /// For PXL-style iteration: PXL iterator class 
    class Iterator {
      public:
        /// Go-to constructor moves iterator to the same positon as the original. 
        Iterator(const Iterator& original) :
            _iter(original._iter), _containerRef(original._containerRef) {}

        /// Standard constructor moves iterator to the beginning of the map. 
        Iterator(const pxl::Map<keytype, itemtype>& map) :
            _containerRef(&map) { first(); }

        /// This method moves the iterator to the beginning of the map. 
        inline void first()  { _iter = _containerRef->begin(); }
        /// This method moves the iterator to the next pair in the map. 
        inline void next()   { _iter++; }
        /// This method returns true, if the iterator hit upper or lower boundary of the map. 
        inline bool isDone() { return _iter == _containerRef->end(); }

        /// This method returns the key of the current pair. 
        inline keytype  key()  { return _iter->first; }
        /// This method returns the item of the current pair. 
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
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_CopyHistory_hh
#define pxl_ptl_CopyHistory_hh


namespace pxl {

// ptl

class ObjectBase;

/// This typedef is intended to hold the origin information of copied objects in pxl::ObjectOwner instances. 
typedef pxl::Map<pxl::Id, pxl::ObjectBase*> CopyHistory;
/// This typedef defines a reference for pxl::CopyHistory
typedef CopyHistory& CopyHistoryRef;
/// This typedef defines a const reference for pxl::CopyHistory
typedef const CopyHistory& CopyHistoryConstRef;

} // namespace pxl

#endif // pxl_ptl_CopyHistory_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_Index_hh
#define pxl_ptl_Index_hh


namespace pxl {

// ptl

class ObjectBase;

/// This typedef is intended to hold the index information (string-object associations) in pxl::ObjectOwner instances. 
typedef pxl::Map<std::string, pxl::ObjectBase*> Index;
/// This typedef defines a reference for pxl::Index
typedef Index& IndexRef;
/// This typedef defines a const reference for pxl::Index
typedef const Index& IndexConstRef;


} // namespace pxl

#endif // pxl_ptl_Index_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_WkPtrBase_hh
#define pxl_ptl_WkPtrBase_hh

namespace pxl {

// ptl

class ObjectBase;

/** 
This base class provides common functionalities for all derived PXL weak pointers. 
*/
class WkPtrBase {
  public:
    virtual ~WkPtrBase() { connect(0); }

    /// This virtual method creates a deep copy and returns a C++ pointer to the newly created weak pointer instance.  
    virtual pxl::WkPtrBase* clone() const
    { return new pxl::WkPtrBase(*this); }

    /// This method returns a C++ pointer of type pxl::ObjectBase to the referenced object  
    inline pxl::ObjectBase* pointer() const { return _objectRef; }
    /// This method returns true, if the referenced object exists.  
    inline bool valid() const { return _objectRef != 0; }
    /// This arrow operator de-references the weak pointer.   
    inline pxl::ObjectBase* operator->() const { return _objectRef; }

    /// This method attempts a dynamic cast on the referenced object
    static inline WkPtrBase* cast_dynamic(WkPtrBase* orig)
    { return orig; }

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
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_WkPtrSpec_hh
#define pxl_ptl_WkPtrSpec_hh



//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_GetSet_hh
#define pxl_ptl_GetSet_hh

namespace pxl {

// ptl

/// This class is used for overloading the bracket operator of PXL objects and weak pointers with read access. 
class Get { public: Get() {} };
/// This class is used for overloading the bracket operator of PXL objects and weak pointers with write access. 
class Set { public: Set() {} };

/// This constant can be used to select the bracket operator of PXL objects and weak pointers for read access. 
extern const pxl::Get get;
/// This constant can be used to select the bracket operator of PXL objects and weak pointers for write access. 
extern const pxl::Set set;

} // namespace pxl

#endif // pxl_ptl_GetSet_hh

namespace pxl {

// iotl

class iStreamer;
class oStreamer;

// ptl

/** 
This class template represents a weak pointer to PXL objects of \p objecttype that aggregate data of \p datatype.   
*/
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

    /// This virtual method creates a deep copy and returns a C++ pointer to the newly created weak pointer instance.  
    virtual pxl::WkPtrBase* clone() const
    { return new pxl::WkPtrSpec<datatype, objecttype>(*this); }

    /// This assignment operator causes the weak pointer to reference the object referenced by \p pptr.
    inline void operator=(const pxl::WkPtrSpec<datatype, objecttype>& pptr)
    { connect(pptr._objectRef); }
    /// This assignment operator causes the weak pointer to reference the object.
    inline void operator=(objecttype& object)
    { connect(&object); }
    /// This assignment operator causes the weak pointer to reference the object pointed to by \p objectptr.
    inline void operator=(objecttype* objectptr)
    { connect(objectptr); }

    // methods to grant object & data access
    /// This method provides direct access to the referenced object.
    inline objecttype& object() const  { return *access(); }
    /// This method grants read access to the data aggregated by the referenced object.
    inline const datatype& get() const { return access()->get(); }
    /// This method grants write access to the data aggregated by the referenced object.
    inline       datatype& set()       { return access()->set(); }

    /// This bracket operator grants read access to the data aggregated by the referenced object.
    inline const datatype& operator()()                const { return get(); }
    /// This bracket operator grants read access to the data aggregated by the referenced object.
    inline const datatype& operator()(const pxl::Get&) const { return get(); }
    /// This bracket operator grants write access to the data aggregated by the referenced object.
    inline       datatype& operator()(const pxl::Set&)       { return set(); }

    /// This arrow operator de-references the weak pointer.   
    inline objecttype* operator->() const { return access(); }

    /// This method attempts a dynamic cast on the referenced object
    static inline WkPtrSpec<datatype, objecttype>* cast_dynamic(WkPtrBase* orig)
    {
        objecttype *object = dynamic_cast<objecttype*>(orig->pointer());
        if (!object)
            return 0;

        // FIXME: This is crude but required:
        if (PXL_UNLIKELY(reinterpret_cast<void*>(object) !=
                         reinterpret_cast<void*>(orig->pointer())))
                pxl::exception("pxl::WkPtrSpec::cast_dynamic()",
                               "Unsupportede multiple inheritance configuration.");

        return reinterpret_cast<pxl::WkPtrSpec<datatype, objecttype>*>(orig);
    }

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
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourmap Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_WkPtrOwner_hh
#define pxl_ptl_WkPtrOwner_hh



namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// ptl

class ObjectBase;

/** 
This class template is a container for weak pointers sorted (uniquely) by \p keytype. 
It contains and manages weak pointers to PXL objects by using the inherited 
functionality of the pxl::Map. Using the object-id as unique key, this class 
is used for relation management and a C++ typedef is provided to pxl::Relations.
*/
template<class keytype>
class WkPtrOwner : public pxl::Map<keytype, pxl::WkPtrBase*> {
  public:
    typedef          pxl::Map<keytype, pxl::WkPtrBase*>                   PtlMap;
    typedef typename pxl::Map<keytype, pxl::WkPtrBase*>::StlConstIterator StlConstIterator;
    typedef typename pxl::Map<keytype, pxl::WkPtrBase*>::StlIterator      StlIterator;

    WkPtrOwner() : pxl::Map<keytype, pxl::WkPtrBase*>() {}
    virtual ~WkPtrOwner() { pxl::WkPtrOwner<keytype>::clearContainer(); }

    /// This method inserts (or replaces) the weak pointer \p wptr indetified by \p key.
    void set(const keytype &key, pxl::WkPtrBase* wptr)
    {
        // FIXME: efficiency, remove & set in one?
        pxl::WkPtrOwner<keytype>::remove(key);
        PtlMap::set(key, wptr);
    }

     /// This method inserts (or replaces) the weak pointer \p wptr indetified by \p key.
    void set(const keytype &key, const pxl::WkPtrBase& wptr)
    {
        // FIXME: efficiency, remove & set in one?
        pxl::WkPtrOwner<keytype>::remove(key);
        PtlMap::set(key, wptr.clone());
    }

    /// This method inserts (or replaces) a weak pointer to the object \p obj indetified by \p key.
    void set(const keytype& key, pxl::ObjectBase& obj);

    /// This method removes the weak pointer indetified by \p key.
    void remove(const keytype& key)
    {
        delete find(key, 0);
        PtlMap::_container.erase(key);
    }

    /// This method returns true, if a weak pointer associated with \p key is found.
    bool has(const keytype& key) const { return find(key, 0) != 0; }

    /// This method removes all  weak pointers contained.
    virtual void clearContainer();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /// For PXL-style iteration: PXL iterator class 
    class Iterator {
      public:
        /// Go-to constructor moves iterator to the same positon as \p original. 
        Iterator(const Iterator& original) :
            _iter(original._iter), _containerRef(original._containerRef) {}
         /// Standard constructor moves iterator to the beginning of \p wpo. 
        Iterator(const pxl::WkPtrOwner<keytype>& wpo) :
            _containerRef(&wpo) { first(); }

        /// This method moves the iterator to the beginning of the map. 
        inline void first()  { _iter = _containerRef->begin(); }
         /// This method moves the iterator to the next pair in the map. 
        inline void next()   { _iter++; }
         /// This method returns true, if the iterator hit upper or lower boundary of the map. 
        inline bool isDone() { return _iter == _containerRef->end(); }

        /// This method returns the current key. 
        inline keytype key()           { return _iter->first; }
        /// This method returns the current item, thus a C++ pointer of type pxl::WkPtrBase* to the weak pointer; see also the wkPtr() and object() methods.  
        inline pxl::WkPtrBase* item()  { return _iter->second; }
        /// This method returns a pxl::WkPtrBase reference to the current weak pointer. 
        inline pxl::WkPtrBase& wkPtr() { return *item(); }
        /// This method returns a pxl::ObjectBase reference to the object pointed at by the current weak pointer. 
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
    /// For PXL-style iteration: PXL iterator class template; 
    /// this iterator behaves like a normal PXL iterator but ignores all weak pointers 
    /// that cannot be interpreted as type wkptrtype (tested by dynamic casts). 
    template<class wkptrtype>
    class TypeIterator {
      public:
         /// Go-to constructor moves iterator to the same positon as the original. 
        TypeIterator(const TypeIterator<wkptrtype>& original) :
            _iter(original._iter), _containerRef(original._containerRef) {}

         /// Standard constructor moves iterator to the beginning of \p wpo. 
        TypeIterator(const pxl::WkPtrOwner<keytype>& wpo) :
            _containerRef(&wpo) { first(); }

        /// This method moves the iterator to the beginning of the map. 
        inline void first()
        { for(_iter = _containerRef->begin(); !isDone() && !item(); _iter++); }

         /// This method moves the iterator to the next pair in the map. 
        inline void next()
        { do _iter++; while(!isDone() && !item()); }

         /// This method returns true, if the iterator hits upper boundary of the map.
        inline bool isDone() const
        { return _iter == _containerRef->end(); }

        /// This method returns the current key. 
        inline keytype key() const
        { return _iter->first; }

        /// This method returns the current item, thus a C++ pointer of type \p wkptrtype* to the weak pointer; see also the wkPtr() and object() methods.  
        inline wkptrtype* item() const
        { return isDone() ? 0 : wkptrtype::cast_dynamic(_iter->second); }

        /// This method returns a wkptrtype reference to the current weak pointer. 
        inline const wkptrtype& wkPtr() const
        { return *item(); }

        /// This method returns a pxl::ObjectBase reference to the object pointed at by the current weak pointer. 
        inline pxl::ObjectBase& object() const
        {
            if (!wkPtr().valid()) // FIXME: exception?
                std::cerr << "pxl::WkPtrOwner::object(): FATAL: The object you intend to access does not exist!" << std::endl;
            return *wkPtr().pointer();
        }

      private:
        TypeIterator() : _containerRef(0) {}

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
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_Relations_hh
#define pxl_ptl_Relations_hh


namespace pxl {

/**
For the convenience of an easy-to-read syntax, we provide this C++ typedef from pxl::WkPtrOwner with pxl::Id.
The class pxl::ObjectBase owns two instances of pxl::Relations 
for manageing mother and daughter relations to other pxl::ObjectBase derivatives. 
For decay tree integrity reasons, the class pxl::ObjectBase allows to 
establish relations only to objects contained in the same object owner 
(and in the case of both objects not being contained in owners). 
Please notice, that by building the relation management on the basis of weak 
pointers, it becomes fail-safe and easy to debug even in the unforeseen case 
of non-adequate object deletion.
*/ 
typedef pxl::WkPtrOwner<pxl::Id> Relations;
/// This typedef defines a reference for pxl::Relations
typedef Relations& RelationsRef;
/// This typedef defines a const reference for pxl::Relations
typedef const Relations& RelationsConstRef;


} // namespace pxl

#endif // pxl_ptl_Relations_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_ObjectBase_hh
#define pxl_ptl_ObjectBase_hh


//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_Layout_hh
#define pxl_ptl_Layout_hh

namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// ptl

/** 
This class holds layout information of PXL objects when visualized
in the Graphical User Interface VisualPxl. For internal use only, 
methods and data members are self-explanatory.  
*/
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

/** 
This base class provides common functionalities for all derived PXL objects,  
such as mother-daughter relations, weak pointer concept and related service routines. 
In addition, it aggregates the PXL unique object-id and layout information for visualization in VisualPxl.
It has a C++ pointer to the pxl::ObjectOwner it is aggregated in in order to avoid
mother/daughter relations to outside objects to be established.
*/
class ObjectBase {
  public:
    virtual ~ObjectBase()
    {
        if (_refWkPtrSpec)
            _refWkPtrSpec->notifyDeleted();
        delete _ptrLayout;
    }

    // FIXME: replace by uuid
    /// This method returns the PXL unique object-id.
    inline pxl::Id id() const { return (int)(void*)this;  }

    /// This method returns a C++ pointer to the pxl::ObjectOwner it is owned by.  
    inline pxl::ObjectOwner* owner() const { return _refObjectOwner; }

    /// This virtual method creates a deep copy and returns a C++ pointer to the newly created object.  
    virtual pxl::ObjectBase* clone() const { return new pxl::ObjectBase(*this); }

    /// This method grants access to the pxl::Relations instance manageing mother relations  
    inline pxl::Relations& getMotherRelations()   { return _motherRelations; }
    /// This method grants access to the pxl::Relations instance manageing daughter relations  
    inline pxl::Relations& getDaughterRelations() { return _daughterRelations; }

    /// This method establishes a mother relation to the \p target object; please notice, that
    /// only relations between objects owned by the same object owner will be established.   
    void linkMother(pxl::ObjectBase& target);
    /// This method establishes a daughter relation to the \p target object; please notice, that
    /// only relations between objects owned by the same object owner will be established.   
    void linkDaughter(pxl::ObjectBase& target);

    /// This method establishes a mother relation to the object referred to by \p target; please notice, that
    /// only relations between objects owned by the same object owner will be established.   
    inline void linkMother(pxl::WkPtrBase& target)   { if (target._objectRef) linkMother(*(target._objectRef)); }
    /// This method establishes a daughter relation to the object referred to by \p target; please notice, that
    /// only relations between objects owned by the same object owner will be established.   
    inline void linkDaughter(pxl::WkPtrBase& target) { if (target._objectRef) linkDaughter(*(target._objectRef)); }

    /// This method removes an existing daughter relation to the \p target object.
    void unlinkMother(pxl::ObjectBase& target);
    /// This method removes an existing daughter relation to the \p target object.
    void unlinkDaughter(pxl::ObjectBase& target);

    /// This method removes an existing mother relation to the object referred to by \p target.
    inline void unlinkMother(pxl::WkPtrBase& target)   { if (target._objectRef) unlinkMother(*(target._objectRef)); }
    /// This method removes an existing daughter relation to the object referred to by \p target.
    inline void unlinkDaughter(pxl::WkPtrBase& target) { if (target._objectRef) unlinkDaughter(*(target._objectRef)); }

    /// This method removes all existing mother relations.
    void unlinkMothers();
    /// This method removes all existing daughter relations.
    void unlinkDaughters();

    /// This method grants access to the layout information provided for visualization in VisualPxl.
    inline pxl::Layout& layout()
    {
        if (!_ptrLayout)
            _ptrLayout = new pxl::Layout;
        return *_ptrLayout;
    }

    /// This virtual method recursively invokes its own and the print() methods of all daughter objects.
    /// @param level verbosity level
    /// @param os output stream, default is std::cout
    /// @param pan print indention
    /// @return output stream
    std::ostream& printDecayTree(int level = 0, std::ostream& os = std::cout, int pan = 1) const;
    
    /// This virtual method is intended to print out object state information on various verbosity levels.
    /// @param level verbosity level
    /// @param os output stream, default is std::cout
    /// @param pan print indention
    /// @return output stream
    virtual std::ostream& print(int level = 0, std::ostream& os = std::cout, int pan = 0) const;

    /// This method creates a weak pointer to itself and returns a pxl::WkPtrBase* to the newly created weak pointer instance. 
    virtual pxl::WkPtrBase* createSelfWkPtr()
    { pxl::exception("pxl::ObjectBase::createSelfWkPtr()", "ATTENTION! Inheriting class must reimplement this virtual method."); return 0; }

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
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_Object_hh
#define pxl_ptl_Object_hh


//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_WkPtr_hh
#define pxl_ptl_WkPtr_hh


namespace pxl {

// ptl

template<class datatype>
class Object;

/** 
This class template is a weak pointer to pxl::Object instances; 
it inherits methods and operators to directly access the data
aggregated in the referenced pxl::Object objects.   
Read/write access is granted by the methods get() and set() as well as by 
the bracket operator.   
Using the assignment operator one acts directly on the aggregated data. 
This class is used to define the pxl::EventViewWkPtr, for instance.
*/
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

    /// This method attempts a dynamic cast on the referenced object
    static inline WkPtr<datatype>* cast_dynamic(WkPtrBase* orig)
    { return reinterpret_cast<WkPtr<datatype>*>(pxl::WkPtrSpec<datatype, pxl::Object<datatype> >::cast_dynamic(orig)); }
};

} // namespace pxl


#endif // pxl_ptl_WkPtr_hh

namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// ptl

/** 
This class template aggregates arbitrary data (type \p datatype); it is able to manage mother-daughter 
relations with other pxl::ObjectBase derivatives. 
Read/write access to the aggregated data is granted by the methods get() and set() as well as by 
the bracket operator. 
Using the assignment operator one acts directly on the aggregated data. 
The copy constructor of this class template copies the data to a 
new instance with no relations.
This class is used to define the pxl::EventView, for instance.
*/
template<class datatype>
class Object : public pxl::ObjectBase {
  public:
    Object() : pxl::ObjectBase(), _data() {}
    Object(const datatype& original) : pxl::ObjectBase(), _data(original) {}
    Object(const pxl::Object<datatype>& original) :
        pxl::ObjectBase(original), _data(original._data) {}
    virtual ~Object() {}

    /// This virtual method creates a deep copy and returns a C++ pointer to it.
    virtual pxl::ObjectBase* clone() const { return new pxl::Object<datatype>(*this); }

    /// This method grants read access to the aggregated data.
    inline const datatype& get() const { return _data; }
    /// This method grants write access to the aggregated data.
    inline       datatype& set()       { return _data; }

     /// This bracket operator grants read access to the aggregated data.
    inline const datatype& operator()()                const { return get(); }
     /// This bracket operator grants read access to the aggregated data.
    inline const datatype& operator()(const pxl::Get&) const { return get(); }
    /// This bracket operator grants write access to the aggregated data.
    inline       datatype& operator()(const pxl::Set&)       { return set(); }

    /// This assignment operator acts directly on the aggregated data.
    inline pxl::Object<datatype>& operator=(const datatype& original)
    { _data = original; return *this; }

    /// This assignment operator acts directly on the aggregated data.
    inline pxl::Object<datatype>& operator=(const pxl::Object<datatype>& original)
    { _data = original._data; return *this; }

    /// This virtual method is intended to print out object state information on various verbosity levels.
    /// @param level verbosity level
    /// @param os output stream, default is std::cout
    /// @param pan print indention
    /// @return output stream
    virtual std::ostream& print(int level = 0, std::ostream& os = std::cout, int pan = 0) const
    {
        os << "called by pxl::Object<...>: ";
        return pxl::ObjectBase::print(level, os, pan);
    }

    virtual pxl::WkPtrBase* createSelfWkPtr()
    { return new pxl::WkPtr<datatype>(*this); }

  protected:
    virtual void storeYourSelf(pxl::oStreamer& output) const;

    datatype _data;

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;
};

} // namespace pxl

//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_iotl_oStreamer_hh
#define pxl_iotl_oStreamer_hh

#include <typeinfo>
#include <sstream>



//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_iotl_internal_hh
#define pxl_iotl_internal_hh

// defaults

#define iotl__default__compressionMode '6'

// file format

#define iotl__eventMarker     "<E>"

#endif // pxl_iotl_internal_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_iotl_TypeManager_hh
#define pxl_iotl_TypeManager_hh



//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_iotl_TypeIdKey_hh
#define pxl_iotl_TypeIdKey_hh


namespace pxl {

/** 
This PXL-internal class is used by the I/O type manager for fast string searching; 
it carries the PXL invariant object- and data-id for a objecttype<datatype> combination 
(such as pxl::Object<double>, for instance)
*/ 
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

    inline bool operator!=(const pxl::TypeIdKey& other) const
    { return !(*this == other); }

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

/** 
This PXL-internal class represents the I/O type manager singleton for  
storing and restoring PXL and user-defined objects/data.
It holds two maps (one sorted by the pxl::TypeIdKey and one by the C++ class-id) 
with pointers to pxl::TypeAgent specializations for all supported object/data types.  
*/ 
class TypeManager {
  public:
    /// This static method provides the singleton access.
    static TypeManager& instance();

    /// This method registers a type agent pointed to by \p agent. It is called by the constructor of the class pxl::TypeAgent.
    void registerAgent(pxl::TypeAgentBase* agent);

    /// This method stores the pxl::ObjectBase derivative \p obj to the given output stream.
    /// For a unique identification of the correct type agent, the C++ class-id of obj must be specified in 
    /// \p cppTypeId.   
    void storeObject(pxl::oStreamer& output, const pxl::ObjectBase& obj, const std::string& cppTypeId) const;
    /// This method restores the pxl::ObjectBase derivative \p obj from the given input stream.
    /// and returns the persistent (=old) object-id.    
    /// For a unique identification of the correct type agent, the C++ class-id of obj must be specified in 
    /// \p cppTypeId.   
    pxl::Id restoreObject(pxl::iStreamer& input, pxl::ObjectBase& obj, const std::string& cppTypeId) const;
    /// This method restores the pxl::ObjectBase derivative obj from the given input stream.
    /// For a unique identification of the correct type agent, the PXL invariant object- and data-id of \p obj 
    /// must be specified in the parameters \p objectTypeId and \p dataTypeId.   
    /// It allocates memory for the new instance, 
    /// makes the pxl::ObjectBase pointer \p (*ppobj) point to this instance,  
    /// restores the object from the given input stream
    /// and returns the persistent (=old) object-id.    
    pxl::Id restoreObject(pxl::iStreamer& input, pxl::ObjectBase** ppobj, const std::string& objectTypeId, const std::string& dataTypeId) const;

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

/**
This class template provides event-by-event data storage to STL streams; 
it inherits from pxl::BasicIoStreamer ( = pxl::BasicLinuxIoStreamer).
This class provides an event buffer to which PXL objects or data can 
be stored using storeObject() and storeData(). The method getEvent()
can be called at the end of each event to have the buffer compressed 
(by ZLIB), streamed into a standard STL stream and finally cleared. 
*/
class oStreamer : public pxl::BasicIoStreamer {
  public:
    oStreamer() : pxl::BasicIoStreamer(), _buffer() {}
    ~oStreamer() {}

    /// This method can be called at the end of each event to have the 
    /// event buffer compressed, streamed to \p cxxx and finally cleared.
    /// An optional information string \p info can be provided for filtered restoring of the data with the pxl::iStreamer::putEventIf() method.
    /// The optional \p compressionMode identifier can be '6' for ZLIB level 6 compression or ' ' for no compression.  
    void getEvent(std::ostream& cxxx, const std::string& info="", char compressionMode = iotl__default__compressionMode);

    /// This method template stores the pxl::ObjectBase derivative \p obj of the known type \p objecttype to the event buffer.
    template<class objecttype>
    void storeObject(const objecttype& obj)
    { pxl::TypeManager::instance().storeObject(*this, obj, typeid(obj).name()); }

    /// This method stores an unknown (abstract) pxl::ObjectBase derivative \p obj to the event buffer.
    inline void storeAbstractObject(const pxl::ObjectBase& data)
    { data.storeYourSelf(*this); }

    /// This method template stores the \p data of \p datatype to the event buffer.
    template<class datatype>
    void storeData(const datatype& data);
    /// This method stores the contents of a pxl::Vector to the event buffer.
    template<class itemtype>
    void storeData(const pxl::Vector<itemtype>& vector);
    /// This method stores the contents of a pxl::Map to the event buffer.
    template<class keytype, class itemtype>
    void storeData(const pxl::Map<keytype, itemtype>& map);

    /// This internal method stores a type-id to the event buffer.
    inline void storeTypeId(const char* typeId)
    { pxl::BasicIoStreamer::storeBasicTypeCStr(_buffer, typeId); }
    /// This internal method stores a type-id to the stream \p cxxx.
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
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_ObjectOwner_hh
#define pxl_ptl_ObjectOwner_hh



namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// ptl

/** 
This class is an active container for pxl::ObjectBase derivatives, such as 
pxl::Object, pxl::CowObject, pxl::Particle or pxl::EventView; 
it has the ownership and deletion responsability for the contained objects. 
The method template create() can be used to create derivatives 
of pxl::ObjectBase within object owners. The method set() can be 
used to manually add objects, has() tests object ownerships, and 
remove() explicitely removes objects from the owner and deletes them.
The copy constructor of the class pxl::ObjectOwner also produces copies of 
the contained objects, and re-establishes corresponding mother-daughter relations 
amongst the copied objects. For the convenience of a quick and targeted object access, the
newly created owner carries a so-called copy history for mapping original and 
copied objects. This information is used by the findCopyOf() method:
provided a reference to the original object, this method returns a pointer to the copied object.
A further, powerful tool for targeted object access is the so-called index, which
allows to map objects to unique string identifiers, the index-id. The method findObject() 
can be used, to directly access objects by their index-ids or object-ids.
The pxl::ObjectOwner extends the functionality of the specialized class 
template pxl::Vector. It provides a selective iterator, the class template 
pxl::ObjectOwner::TypeIterator, that ignores all objects other than the 
specialized data type. 
For the convenience of an easy-to-read syntax, we provide a C++ typedef to pxl::Objects. 
*/
class ObjectOwner : public pxl::Vector<pxl::ObjectBase*> {
  public:
    ObjectOwner();
    /// This copy constructor performs a deep copy of object 
    /// owner \p original and all contained objects with (redirected) relations.
    /// A copy history keeps track of originals and copies
    /// and the findCopyOf() method allows quick access to the copies. 
    ObjectOwner(const pxl::ObjectOwner& original);
    /// This destructor deletes all contained objects.
    virtual ~ObjectOwner() { pxl::ObjectOwner::clearContainer(); }

    /// This method template creates a new instance of \p objecttype;
    /// objecttype must be a class inheriting from pxl::ObjectBase;
    /// the newly created instance is owned and will be deleted by this object owner. 
    template<class objecttype>
    inline objecttype& create()
    {
        objecttype* pitem = new objecttype;
        pitem->_refObjectOwner = this;
        _container.push_back(static_cast<pxl::ObjectBase*>(pitem));
        return *pitem;
    }

    /// This method template creates a copy of \p original by invoking the copy constructor of \p objecttype; 
    /// \p objecttype must be a class inheriting from pxl::ObjectBase;
    /// the newly created instance is owned and will be deleted by this object owner. 
    template<class objecttype>
    inline objecttype& create(const objecttype& original)
    {
        objecttype* pitem = new objecttype(original);
        pitem->_refObjectOwner = this;
        _container.push_back(static_cast<pxl::ObjectBase*>(pitem));
        return *pitem;
    }

    /// This method template creates a new \p objecttype instance by invoking a \p ctrtype overloaded constructor; 
    /// \p objecttype must be a class inheriting from pxl::ObjectBase;
    /// the newly created instance is owned and will be deleted by this object owner. 
    template<class objecttype, class ctrtype>
    inline objecttype& create(const ctrtype& original)
    {
        objecttype* pitem = new objecttype(original);
        pitem->_refObjectOwner = this;
        _container.push_back(static_cast<pxl::ObjectBase*>(pitem));
        return *pitem;
    }

    /// This method inserts \p item in the container of this object owner and takes deletion responsability.
    void set(pxl::ObjectBase& item);
    /// This method deletes \p item.
    void remove(pxl::ObjectBase& item);
    /// This method returns true, if \p item is owned by this object owner.
    bool has(const pxl::ObjectBase& item) const;

    /// This method clears the object owner and deletes all owned objects. 
    virtual void clearContainer();

    /// This method searches the index for the index-id \p idx and returns a dynamically casted 
    /// C++ pointer of type \p objecttype* to the corresponding object; 
    /// in case idx is not found a null pointer is returned.
    template<class objecttype>
    inline objecttype* findObject(const std::string& idx) const	// goes via Index & casts
    { return dynamic_cast<objecttype*>(_index.find(idx, 0)); }

    /// This method searches the copy history to locate the copy of \p original and 
    /// returns a dynamically casted C++ pointer of type \p objecttype* to the corresponding copy; 
    /// in case no copy can be traced a null pointer is returned.
    template<class objecttype>
    inline objecttype* findCopyOf(const pxl::ObjectBase& original) const // goes via CopyHistory & casts
    { return dynamic_cast<objecttype*>(_copyHistory.find(original.id(), 0)); }    

    /// This method provides direct access to the copy history (created by the copy constructor). 
    inline const pxl::CopyHistory& getCopyHistory() const { return _copyHistory; }
    /// This method clears the copy history  (created by the copy constructor). 
    inline void clearCopyHistory() { _copyHistory.clearContainer(); }

    /// This method registers the object \p obj with the index-id \p idx in the index and returns true in case of success;
    /// please notice, that \p obj must be owned by this object owner and \p idx must not be a zero length string.  
    inline bool setIndex(const std::string& idx, pxl::ObjectBase& obj)
    {
        if (!idx.length() || !has(obj))
            return false;
        _index.set(idx, &obj);
        return true;
    }

    /// This method provides direct access to the index. 
    inline const pxl::Index& getIndex() const { return _index; }
    /// This method removes the index entry with index-id \p idx; please notice: it does not remove the object itself. 
    inline void removeIndex(const std::string& idx) { _index.remove(idx); }
    /// This method clears the index; please notice: it does not remove the objects themselves.
    inline void clearIndex() { _index.clearContainer(); }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /// For PXL-style iteration: PXL iterator class 
    class Iterator {
      public:
        /// Go-to constructor moves iterator to the same positon as \p original. 
        Iterator(const pxl::ObjectOwner::Iterator& original) :
            _iter(original._iter), _containerRef(original._containerRef) {}
 
        /// Standard constructor moves iterator to the beginning of \p vector. 
        Iterator(const pxl::ObjectOwner& vector) :
            _containerRef(&vector) { first(); }

        /// This method moves the iterator to the beginning of the vector. 
        inline void first()  { _iter = _containerRef->begin(); }
        /// This method moves the iterator to the next pair in the vector. 
        inline void next()   { _iter++; }
         /// This method returns true, if the iterator hit upper or lower boundary of the vector. 
        inline bool isDone() { return _iter == _containerRef->end(); }

        /// This method returns the current item, thus a C++ pointer of type pxl::ObjectBase* to the object; see also the object() method. 
        inline pxl::ObjectBase* item() { return *_iter; }
        /// This method returns a pxl::ObjectBase reference to the current object. 
        inline pxl::ObjectBase& object() { return *item(); }

      private:
        Iterator() : _containerRef(0) {;}

        StlConstIterator    _iter;
        const pxl::ObjectOwner* _containerRef;
    };
    // - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /// For PXL-style iteration: PXL iterator class template; 
    /// this iterator behaves like a normal PXL iterator but ignores all objects 
    /// that cannot be interpreted as type objecttype (tested by dynamic casts). 
    template<class objecttype>
    class TypeIterator {
      public:
        /// Go-to constructor moves iterator to the same positon as the original. 
        TypeIterator(const pxl::ObjectOwner::TypeIterator<objecttype>& original) :
            _iter(original._iter), _containerRef(original._containerRef) {}

        /// Standard constructor moves iterator to the beginning of the vector of \p oo. 
        TypeIterator(const pxl::ObjectOwner& oo) :
            _containerRef(&oo) { first(); }

        /// This method moves the iterator to the beginning of the vector. 
        inline void first()  { for(_iter = _containerRef->begin(); !tryItem(); _iter++); }
        /// This method moves the iterator to the next pair in the vector. 
        inline void next()   { do _iter++; while(!tryItem()); }
          /// This method returns true, if the iterator hit upper or lower boundary of the vector. 
        inline bool isDone() { return _iter == _containerRef->end(); }

        /// This method returns the current item, thus a C++ pointer of type \p objecttype* to the object; see also the object() method. 
        inline objecttype* item() { return dynamic_cast<objecttype*>(*_iter); }
        /// This method returns a objecttype reference to the current object. 
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
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_Objects_hh
#define pxl_ptl_Objects_hh


namespace pxl {

/// For the convenience of an easy-to-read syntax, we provide this C++ typedef from pxl::ObjectOwner. 
typedef ObjectOwner Objects;
/// This typedef defines a reference for pxl::Objects
typedef Objects& ObjectsRef;
/// This typedef defines a const reference for pxl::Objects
typedef const Objects& ObjectsConstRef;


} // namespace pxl

#endif // pxl_ptl_Objects_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_CowObject_hh
#define pxl_ptl_CowObject_hh


//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_CowWkPtr_hh
#define pxl_ptl_CowWkPtr_hh


namespace pxl {

// ptl

template<class datatype>
class CowObject;

/** 
This class template is a weak pointer to pxl::CowObject instances; 
it inherits methods and operators to directly access the data
aggregated in the referenced pxl::CowObject objects.   
Read/write access is granted by the methods get() and set() as well as by 
the bracket operator.  
Using the assignment operator one acts directly on the aggregated data. 
This class is used to define the pxl::ParticleWkPtr, for instance.
*/
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

    /// This method attempts a dynamic cast on the referenced object
    static inline CowWkPtr<datatype>* cast_dynamic(WkPtrBase* orig)
    { return reinterpret_cast<CowWkPtr<datatype>*>(pxl::WkPtrSpec<datatype, pxl::CowObject<datatype> >::cast_dynamic(orig)); }
};

} // namespace pxl


#endif // pxl_ptl_CowWkPtr_hh

namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// ptl

/** 
This class template deploys a copy-on-write mechanism to aggregate arbitrary data (type \p datatype); 
it is able to manage mother-daughter relations with other pxl::ObjectBase derivatives.  
Read/write access to the aggregated data is granted by the methods get() and set() as well as by   
the bracket operator. 
Using the assignment operator one acts directly on the aggregated data. 
The copy constructor of this class template creates 
a new instance with no relations. It does *not* copy the aggregated data, that is placed in 
a pxl::CowObject::DataSocket, but increases its reference counter by one. 
Only upon write access, and if more than one object refers to this pxl::CowObject::DataSocket  
instance, the data socket including data are acutally duplicated. 
This class is used to define the pxl::Particle, for instance. 
*/
template<class datatype>
class CowObject : public pxl::ObjectBase {
  private:
    class DataSocket {
      public:
        DataSocket() : _references(1) {}
        DataSocket(const datatype& original) :
            _references(1), _data(original) {}
        DataSocket(const DataSocket& original) :
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
    CowObject(const datatype& original) : pxl::ObjectBase() { _dataSocket = new DataSocket(original); }
    CowObject(const pxl::CowObject<datatype>& original) : pxl::ObjectBase(original)
    { _dataSocket = original._dataSocket; _dataSocket->_references++; }
    virtual ~CowObject() { dropDataSocket(); }

    /// This virtual method creates a deep copy and returns a C++ pointer to it.
    virtual pxl::ObjectBase* clone() const { return new pxl::CowObject<datatype>(*this); }

    /// This method grants read access to the aggregated data.
    inline const datatype& get() const { return _dataSocket->getData(); }
    
    /// This method grants write access to the aggregated data; 
    /// if necessary, the copy-on-write mechanism performs a deep copy of the aggregated data first. 
    inline datatype& set()
    {
        if (_dataSocket->_references > 1) {
            _dataSocket->_references--;
            _dataSocket = new DataSocket(*_dataSocket);
        }
        return _dataSocket->getData(); 
    }

    /// This bracket operator grants read access to the aggregated data.
    inline const datatype& operator()()                const { return get(); }
    /// This bracket operator grants read access to the aggregated data.
    inline const datatype& operator()(const pxl::Get&) const { return get(); }
    /// This bracket operator grants write access to the aggregated data; 
    /// if necessary, the copy-on-write mechanism performs a deep copy of the aggregated data first. 
    inline       datatype& operator()(const pxl::Set&)       { return set(); }

    /// This assignment operator acts directly on the aggregated data.
    inline pxl::CowObject<datatype>& operator=(const datatype& original)
    {
        dropDataSocket();
        _dataSocket = new DataSocket(original);
        return *this;
    }

    /// This assignment operator acts directly on the aggregated data.
    inline pxl::CowObject<datatype>& operator=(const pxl::CowObject<datatype>& original)
    {
        dropDataSocket();
        _dataSocket = original._dataSocket;
        _dataSocket->_references++;
        return *this;
    }

    /// This virtual method is intended to print out object state information on various verbosity levels.
    /// @param level verbosity level
    /// @param os output stream, default is std::cout
    /// @param pan print indention
    /// @return output stream
    virtual std::ostream& print(int level = 0, std::ostream& os = std::cout, int pan = 0) const
    {
        os << "called by pxl::CowObject<...>: ";
        return pxl::ObjectBase::print(level, os, pan);
    }

    virtual pxl::WkPtrBase* createSelfWkPtr()
    { return new pxl::CowWkPtr<datatype>(*this); }

  protected:
    CowObject(DataSocket& original) : pxl::ObjectBase()
    { _dataSocket = &original; _dataSocket->_references++; }

    inline void dropDataSocket()
    {
        if (_dataSocket->_references-- == 1)
            delete _dataSocket;
    }

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
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_SpyObject_hh
#define pxl_ptl_SpyObject_hh


//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_SpyWkPtr_hh
#define pxl_ptl_SpyWkPtr_hh


namespace pxl {

// ptl

template<class datatype>
class SpyObject;

/// This is an undocumented class under development. Usage deprecated. 
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

    /// This method attempts a dynamic cast on the referenced object
    static inline SpyWkPtr<datatype>* cast_dynamic(WkPtrBase* orig)
    { return reinterpret_cast<SpyWkPtr<datatype>*>(pxl::WkPtrSpec<datatype, pxl::SpyObject<datatype> >::cast_dynamic(orig)); }
};

} // namespace pxl


#endif // pxl_ptl_SpyWkPtr_hh

namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// ptl

/// This is an undocumented class under development. Usage deprecated. 
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

    virtual pxl::WkPtrBase* createSelfWkPtr()
    { return new pxl::SpyWkPtr<datatype>(*this); }

  protected:

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
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_ptl_Filter_hh
#define pxl_ptl_Filter_hh



namespace pxl {

/** 
This class template provides a sorted filter for PXL physics objects;
it inherits from pxl::WkPtrOwner and can hold a list of weak pointers to 
objects of type \p objecttype sorted by \p sorttype. The user can write
his own filters by public inheritance from this class and reimplementation
of the methods pass() and sort().  
*/ 
template<class objecttype, class sorttype>
class Filter : public pxl::WkPtrOwner<int> {
  public:
    virtual ~Filter() {}

  protected:
    Filter() {}

    /// This method returns true in case the object \p obj passes the filter criterion. 
    virtual bool     pass(const objecttype& obj) const = 0;
    /// This method returns the sort criterion of type \p sorttype. 
    virtual sorttype sort(const objecttype& obj) const = 0;

    /// This method applies the filter by running over the \p objects container and fills
    /// itself with weak pointers to the objects passing the filter criteria.  
    int apply(const pxl::ObjectOwner& objects)
    {
        typedef std::pair<sorttype,objecttype*>     LocalStlPair;
        typedef std::multimap<sorttype,objecttype*> LocalStlMultiMap;

        LocalStlMultiMap map;
        // fill map:
        for(pxl::Objects::TypeIterator<objecttype> iter(objects);
            !iter.isDone(); iter.next()) {

            if (!pass(iter.object()))
                 continue;

            map.insert(LocalStlPair(sort(iter.object()), iter.item()));
        }

        int position = 0;
        clearContainer();
        for(typename LocalStlMultiMap::const_iterator iterMap = map.begin();
            iterMap != map.end(); ++iterMap)
            set(position++, iterMap->second->createSelfWkPtr());

        return position;   
    }
};

/// This typedef provides a PXL-style Iterator for pxl::Filter
typedef WkPtrOwner<int>::Iterator FilterIterator;

} // namespace pxl

#endif // pxl_ptl_Filter_hh

//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_iotl_hh
#define pxl_iotl_hh

//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_iotl_iStreamer_hh
#define pxl_iotl_iStreamer_hh





namespace pxl {

// iotl

/**
This class template provides event-by-event data reading from STL streams; 
it inherits from pxl::BasicIoStreamer ( = pxl::BasicLinuxIoStreamer).
This class provides an event buffer from which PXL objects or data can 
be restored using restoreObject() and restoreData(). The method putEvent()
must be called to have the buffer filled with the decompressed event data.
*/
class iStreamer : public pxl::BasicIoStreamer {
  public:
    iStreamer();
    virtual ~iStreamer();

    /// This method decompresses the next event from the input stream \p cxxx 
    /// to the event buffer the information string of which is identical to \p info; 
    /// if the end of the stream was reached, false is returned.
    inline bool putEventIf(std::istream& cxxx, const std::string& info)
    { return putEvent(cxxx, -1, info); }
    
    /// This method decompresses the next event from the input stream \p cxxx 
    /// to the event buffer; if the end of the stream was reached, false is returned.
    inline bool putEvent(std::istream& cxxx)
    { return putEvent(cxxx, +1, ""); }

    /// This method proceeds to the beginning of the next event in the stream \p cxxx; true is returned in case of success.
    inline bool next(std::istream& cxxx) { return putEvent(cxxx,  0, ""); }
    /// This method steps back to the beginning of the previous event in the stream \p cxxx; true is returned in case of success.
    bool previous(std::istream& cxxx);
    
    /// This method returns true, if the end of the stream was reached.
    inline bool endOfEvent() { return _buffer.peek() == EOF; }

    /// This method template restores the subsequent pxl::ObjectBase derivative \p obj of the known type \p objecttype from the event buffer.
    template<class objecttype>
    pxl::Id restoreObject(objecttype& obj)
    { return pxl::TypeManager::instance().restoreObject(*this, obj, typeid(obj).name()); }

    /// This method creates and restores the subsequent pxl::ObjectBase derivative from the event buffer
    /// and returns a pxl::ObjectBase reference to it.
    inline pxl::ObjectBase& restoreAbstractObject()
    {
        pxl::ObjectBase* pobj;
        restoreAbstractObject(&pobj);
        return *pobj;
    }

    /// This method template restores the subsequent \p data of \p datatype from the event buffer.
    template<class datatype>
    pxl::Id restoreData(datatype& data);
    /// This method restores the contents of the subsequent pxl::Vector from the event buffer.
    template<class itemtype>
    pxl::Id restoreData(pxl::Vector<itemtype>& vector);
    /// This method restores the contents of the subsequent pxl::Map from the event buffer.
    template<class keytype, class itemtype>
    pxl::Id restoreData(pxl::Map<keytype, itemtype>& map);

    /// This method restores the subsequent type id from the stream \p cxxx and returns it in form of a C++ std::string.
    inline std::string restoreTypeId(std::istream& cxxx)
    {
        std::string read;
        pxl::BasicIoStreamer::restoreBasicTypeCStr(cxxx, read);
        return read;
    }

    /// This method restores the subsequent type id from the event buffer and returns it in form of a C++ std::string.
    inline std::string restoreTypeId()
    { return restoreTypeId(_buffer); }

    /// This method restores the subsequent type id from the stream \p cxxx and ensures that it is identical to \p expectedTypeId.
    inline void restoreTypeId(std::istream& cxxx, const char* expectedTypeId)
    {
        std::string read;
        pxl::BasicIoStreamer::restoreBasicTypeCStr(cxxx, read);
        if (read != expectedTypeId)
            pxl::exception("pxl::iStreamer::restoreTypeId()",
                           std::string("Unexpected object type: ") + read);
    }

    /// This method restores the subsequent type id from the event buffer and ensures that it is identical to \p expectedTypeId.
    inline void restoreTypeId(const char* expectedTypeId)
    { restoreTypeId(_buffer, expectedTypeId); }

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
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_iotl_iFile_hh
#define pxl_iotl_iFile_hh

namespace pxl {

/**
This base class provides the user interface for file-based restoring of PXL objects; see pxl::iDiskFile (= pxl::iDiskFileVx) for a concrete implementation.   
pxl::iFile provides methods for opening and closing files and for reading events.
*/
class iFile {
  public:
    virtual ~iFile() {}

    /// This virtual method is intended to open a file named \p filename.
    virtual bool open(const std::string& filename)
    { return false; }

    /// This virtual method is intended to close an open file.
    virtual void close() {}
    /// This virtual method is intended to read the following event in the open file;
    virtual bool readEvent() { return false; }
    /// This virtual method is intended to read the following event in the open file the information string of which is identical to \p info;
    virtual bool readEventIf(const std::string& info) { return false; }
    /// This virtual method skips the next event and returns 1 in case of success.
    inline  int  skipEvent() { return skipEvents(1); }
    /// This virtual method is intended to try to skip the next \p n events and to return the number of events actually skipped.
    virtual int  skipEvents(int n) { return 0; }
    /// This virtual method is intended to return true if the end of the file is reached.
    virtual bool endOfFile() { return true; }

  protected:
    iFile() {}

  private:
    iFile(const pxl::iFile&) {}
};

} // namespace pxl

#endif // pxl_iotl_iFile_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_iotl_iDiskFileVx_hh
#define pxl_iotl_iDiskFileVx_hh

#include <fstream>
#include <ios>



namespace pxl {

/**
This class template allows to restore objects from PXL I/O disk files; it provides methods 
for opening and closing files, for reading/skipping events and all methods inherited by \p istreamertype
(pxl::iStreamer for the current PXL I/O file version).
*/
template<class istreamertype>
class iDiskFileVx : public pxl::iFile, public istreamertype {
  public:
    iDiskFileVx() : pxl::iFile(), istreamertype() {}
    virtual ~iDiskFileVx() { close(); }

    /// This virtual method opens a disk file with name \p filename.
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

    /// This virtual method closes an open disk file.
    virtual void close()
    {
        if (!_file.is_open()) return;
        _file.close();
    }

    /// This virtual method reads and decompresses the following event in the open file;
    /// if the end of the file was reached, false is returned. 
    virtual bool readEvent()
    {
        if (_file.is_open()) {
            return istreamertype::putEvent(_file);
        } else {
            pxl::exception("pxl::iDiskFileVx<>::readEvent()", "No file open.");
            return false;
        }
    }

    /// This virtual method reads and decompresses the next event in the open file
    /// the information string of which is identical to \p info; all other events are skipped;
    /// if the end of the file was reached, false is returned. 
    virtual bool readEventIf(const std::string& info)
    {
        if (_file.is_open()) {
            return istreamertype::putEventIf(_file, info);
        } else {
            pxl::exception("pxl::iDiskFileVx<>::readEventIf()", "No file open.");
            return false;
        }
    }

    /// This virtual method tries to skip the next \p n events and returns the number of events actually skipped.
    virtual int skipEvents(int n)
    {
        int s = 0;
        for(; n < 0 && istreamertype::previous(_file); n++) s--;
        for(; n > 0 && istreamertype::next(_file);     n--) s++;
        return s;
    }

    /// This virtual method returns true if the end of the file was reached.
    virtual bool endOfFile()
    { return _file.peek() == EOF; }

  protected:
    std::fstream _file;

  private:
    iDiskFileVx(const pxl::iDiskFileVx<istreamertype>&) {}
};

} // namespace pxl

#endif // pxl_iotl_iDiskFileVx_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_iotl_iDiskFile_hh
#define pxl_iotl_iDiskFile_hh


namespace pxl {

/// This typedef provides the current PXL I/O file version by specializing  pxl::iDiskFileVx with pxl::iStreamer.
typedef pxl::iDiskFileVx<pxl::iStreamer> iDiskFile;

} // namespace pxl

#endif // pxl_iotl_iDiskFile_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_iotl_oFile_hh
#define pxl_iotl_oFile_hh


namespace pxl {

/**
This base class provides the user interface for file-based storing of PXL objects; see pxl::oDiskFile (= pxl::oDiskFileVx) for a concrete implementation.   
pxl::oFile provides methods for opening and closing files and for writing events.
*/
class oFile {
  public:
    virtual ~oFile() {}

    /// This virtual method is intended to open a file named \p filename in append or overwrite mode.
    virtual bool open(const std::string& filename, bool append = false)
    { return false; }

    /// This virtual method is intended to close an open file.
    virtual void close() {}

    /// This virtual method is intended to write the current event buffer to the open file;
    /// an information string as well as the compression mode can be provided.  
    virtual void writeEvent(const std::string& info = "", char compressionMode = iotl__default__compressionMode) {}

  protected:
    oFile() {}

  private:
    oFile(const pxl::oFile&) {}
};

} // namespace pxl

#endif // pxl_iotl_oFile_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_iotl_oDiskFileVx_hh
#define pxl_iotl_oDiskFileVx_hh




namespace pxl {

/**
This class template provides data storage to PXL I/O disk files; it provides methods 
for opening and closing files, for writing events and all methods inherited by \p ostreamertype
(pxl::oStreamer for the current PXL I/O file version).
*/
template<class ostreamertype>
class oDiskFileVx : public pxl::oFile, public ostreamertype {
  public:
    oDiskFileVx() : pxl::oFile(), ostreamertype() {}
    virtual ~oDiskFileVx() { close(); }

    /// This virtual method opens a disk file with name \p filename in append or overwrite mode.
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

    /// This virtual method closes an open disk file.
    virtual void close()
    {
        if (!_file.is_open()) return;
        _file.rdbuf()->pubsync();
        _file.close();
    }

    /// This virtual method writes the current event buffer to the open file;
    /// an information string as well as the compression mode can be provided  
    /// (see pxl::oStreamer::getEvent() for details of the current PXL I/O file version).
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
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_iotl_oDiskFile_hh
#define pxl_iotl_oDiskFile_hh


namespace pxl {

/// This typedef provides the current PXL I/O file version by specializing  pxl::oDiskFileVx with pxl::oStreamer.
typedef pxl::oDiskFileVx<pxl::oStreamer> oDiskFile;

} // namespace pxl

#endif // pxl_iotl_oDiskFile_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_iotl_TypeAgentBase_hh
#define pxl_iotl_TypeAgentBase_hh



namespace pxl {

/** 
This PXL-internal base class is used with the I/O type manager singleton pxl::TypeManager 
for storing and restoring PXL objects with PXL or user-defined data types; see also pxl::TypeAgent.
*/ 
class TypeAgentBase {
  public:
    virtual ~TypeAgentBase() {}
    /// This virtual method is (re-)implemented in pxl::TypeAgent.  
    virtual void storeObject(pxl::oStreamer& output, const pxl::ObjectBase& objectbase) {}
    /// This virtual method is (re-)implemented in pxl::TypeAgent.  
    virtual pxl::Id restoreObject(pxl::iStreamer& input, pxl::ObjectBase& obj) { return 0; }
    /// This virtual method is (re-)implemented in pxl::TypeAgent.  
    virtual pxl::Id restoreObject(pxl::iStreamer& input, pxl::ObjectBase** ppobj) { return 0; }

    /// This method returns the PXL invariant object-id that can be handled by this instance.
    inline const std::string& getObjectTypeId() { return _objectTypeId; }
    /// This method returns the PXL invariant data-id that can be handled by this instance.
    inline const std::string& getDataTypeId() { return _dataTypeId; }
    /// This method returns the (transient!) C++ class-id that can be handled by this instance.
    inline const std::string& getCppTypeId() { return _cppTypeId; }

  protected:
    std::string _objectTypeId;
    std::string _dataTypeId;
    std::string _cppTypeId;
};

} // namespace pxl

#endif // pxl_iotl_TypeAgentBase_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_iotl_TypeAgent_hh
#define pxl_iotl_TypeAgent_hh




namespace pxl {

/** 
This PXL-internal class template is used with the I/O type manager singleton pxl::TypeManager for 
storing and restoring PXL objects of type \p objecttype.
(Please notice: pxl::TypeAgent specializations must be globally declard and register themselves 
automatically with the pxl::TypeManager.) 
*/ 
template<class objecttype>
class TypeAgent : public TypeAgentBase {
  public:
    /// This constructor accepts the PXL invariant object-id \p objectTypeId and determines 
    /// data-id as well as the C++ class-id of \p objecttype; it registers this instance with 
    /// the pxl::TypeManager.
    TypeAgent(const std::string& objectTypeId);
    virtual ~TypeAgent() {}

    /// This virtual method stores \p obj (that must be of \p objecttype!) to the given output stream.  
    virtual void storeObject(pxl::oStreamer &output, const pxl::ObjectBase& obj);
    /// This virtual method restores \p obj (that must be of \p objecttype!) from the given input stream
    /// and returns the persistent (=old) object-id.    
    virtual pxl::Id restoreObject(pxl::iStreamer& input, pxl::ObjectBase& obj);
    /// This virtual method allocates memory for a new objecttype instance, 
    /// makes the pxl::ObjectBase pointer \p (*ppobj) point to this instance,  
    /// restores the object from the given input stream
    /// and returns the persistent (=old) object-id.    
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
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
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

iotl__declareDataTypeProto(pxl::Variant)

#endif // pxl_iotl_types_hh

#endif // pxl_iotl_hh

#endif // pxl_ptl_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_pol_hh
#define pxl_pol_hh

//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_pol_Basic3Vector_hh
#define pxl_pol_Basic3Vector_hh

#include <cmath>


namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// pol

#define EPSILON 1.0e-9
#ifndef M_PI
	#define M_PI        3.14159265358979323846
#endif

/** 
This class provides a simple threevector with basic algebra. The methods provided are self-explanatory.
*/
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
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_pol_Basic4Vector_hh
#define pxl_pol_Basic4Vector_hh



namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// pol

/** 
This class provides a simple Lorentz-fourvector with basic algebra. The methods provided are self-explanatory.
*/
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
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_pol_UserRecord_hh
#define pxl_pol_UserRecord_hh





namespace pxl {
/**
This class is intented to aggregate information complementary to data members in form
of string-variant pairs. 
All PXL physics objects such as own user records and provide methods for quick 
access to individual user record entries.
*/ 
class UserRecord : public pxl::Map<std::string, pxl::Variant> {
  public:
    typedef pxl::Map<std::string, pxl::Variant> PtlMap;

    /// This method template inserts (or replaces) the user record indetified by \p key.
    template<typename datatype>
    inline void set(const std::string& key, const datatype& item)
    {
        pxl::Variant& value = PtlMap::findOrAlloc(key);
        if (PXL_UNLIKELY(value.getType() != Variant::TYPE_NULL))
            value.clear();

        value.template init<datatype>();
        value.template set<datatype>(item);
    }

    /// This method template searches and returns the user record item identified by \p key; \p defaultitem is returned in case the key is not found. 
    template<typename datatype>
    inline datatype find(const std::string& key, const datatype& defaultitem) const
    {
        const pxl::Variant* value = PtlMap::findOrReturn(key);
        if (!value)
            return defaultitem;
        return value->template get<datatype>();
    }

    /// This method template searches and returns the user record item indetified by \p key; a pxl::Exception is thrown in case the key is not found. 
    template<typename datatype>
    inline datatype find(const std::string& key) const
    {
        const pxl::Variant* value = PtlMap::findOrReturn(key);
        if (!value)
            pxl::exception("pxl::UserRecord::find(...)",
                           "key not found and no default item provided");
        return value->template get<datatype>();
    }
};

/// This typedef defines a reference for pxl::UserRecord
typedef UserRecord& UserRecordRef;
/// This typedef defines a const reference for pxl::UserRecord
typedef const UserRecord& UserRecordConstRef;

// map stream implementations to pxl::Map<...> base class

template<>
inline pxl::Id iStreamer::restoreData(pxl::UserRecord& map)
{ return restoreData(static_cast<pxl::UserRecord::PtlMap&>(map)); }

template<>
inline void oStreamer::storeData(const pxl::UserRecord& map)
{ storeData(static_cast<const pxl::UserRecord::PtlMap&>(map)); }

} // namespace pxl

#endif // pxl_pol_UserRecord_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_pol_BasicObject_hh
#define pxl_pol_BasicObject_hh






namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// pol

/** 
This class provides common functionalities of PXL physics objects like 
data members for storing an object name and flags for status, Monte-Carlo mode and 
object locking; more specific information, such as b-tags, jet cone sizes or energy 
corrections, for instance, can be stored in the so-called user records (see pxl::UserRecord). 
An integer workflag facilitates tagging of individual objects. 
*/
class BasicObjectData {
  public: 
    BasicObjectData() :
        _locked(0),
        _monteCarloMode(0),
        _name("default"),
        _status(0),
        _workflag(0),
        _userRecords()
    {}

    /// This method returns the value of the lock flag.
    inline bool        getLocked()   const { return _locked; }
    /// This method returns the value of the Monte-Carlo flag.
    inline int         getMonteCarloMode() const { return _monteCarloMode; }
    /// This method returns the name.
    inline std::string getName()     const { return _name; }
    /// This method returns the value of the status flag.
    inline int         getStatus()   const { return _status; }
    /// This method returns the value of the workflag.
    inline int         getWorkflag() const { return _workflag; }

    /// This method sets the value of the lock flag to \p v.
    inline void setLocked(bool v)      { _locked = v; }
    /// This method sets the value of the Monte-Carlo flag to \p v.
    inline void setMonteCarloMode(int v) { _monteCarloMode = v; }
    /// This method sets the name to the contents of \p v.
    inline void setName(std::string v) { _name = v; }
    /// This method sets the value of the status flag to \p v.
    inline void setStatus(int v)       { _status = v; }
    /// This method sets the value of the workflag to \p v.
    inline void setWorkflag(int v)     { _workflag = v; }

    /// This method provides access to the user records.
    inline const pxl::UserRecord& getUserRecord() const { return _userRecords; }

    /// This method sets the user record entry identified by \p key to \p item.
    template<typename datatype>
    inline void setUserRecord(const std::string& key, const datatype& item)
    { _userRecords.template set<datatype>(key, item); }

    /// This method removes the user record entry identified by \p key.
    inline void removeUserRecord(const std::string& key)
    { _userRecords.remove(key); }

    /// This method searches the user record entry identified by \p key; \p defaultitem is returned in case key is not found.
    template<typename datatype>
    inline datatype findUserRecord(const std::string& key, const datatype& defaultitem) const
    { return _userRecords.template find<datatype>(key, defaultitem); }

    /// This method searches the user record entry identified by \p key; an exception is thrown in case key is not found.
    template<typename datatype>
    inline datatype findUserRecord(const std::string& key) const
    { return _userRecords.template find<datatype>(key); }

  protected:
    bool        _locked;
    int         _monteCarloMode;
    std::string _name;
    int         _status;
    int         _workflag;

    pxl::UserRecord  _userRecords;

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;  
};

/// This typedef defines a basic PXL physics object with pxl::BasicObjectData as data.
typedef pxl::CowObject<pxl::BasicObjectData> BasicObject;
/// This typedef defines a weak pointer for pxl::BasicObject
typedef pxl::WkPtrSpec<pxl::BasicObjectData, pxl::BasicObject> BasicObjectWkPtr;
/// This typedef defines a reference for pxl::BasicObject
typedef BasicObject& BasicObjectRef;
/// This typedef defines a const reference for pxl::BasicObject
typedef const BasicObject& BasicObjectConstRef;

} // namespace pxl

iotl__declareObjectTypeProto(pxl::BasicObject)
iotl__declareDataTypeProto(pxl::BasicObjectData)

#endif // pxl_pol_BasicObject_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_pol_BasicManager_hh
#define pxl_pol_BasicManager_hh





namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// pol

/** 
This class the functionality of the pxl::BasicObjectData class by providing an object owner (see pxl::ObjectOwner) and 
corresponding service methods. This way, physics objects (like instances of the classes 
pxl::Particle, pxl::Vertex or pxl::Collision as well as other arbitrary pxl::ObjectBase derivatives can be 
aggregated and managed. 
*/
class BasicObjectManagerData : public pxl::BasicObjectData {
  public: 
    BasicObjectManagerData() :
        BasicObjectData(), _objects() {}
    /// This copy constructor performs a deep copy of \p original 
    /// with all contained objects and their (redirected) relations.
    BasicObjectManagerData(const pxl::BasicObjectManagerData& original) :
        BasicObjectData(original), _objects(original._objects) {}

    // create
    /// This method template creates a new instance of \p objecttype;
    /// objecttype must be a class inheriting from pxl::ObjectBase;
    /// the newly created instance is owned and will be deleted by the object owner. 
    template<class datatype>
    datatype& create()
    { return _objects.create<datatype>(); }

    /// This method template creates a new \p objecttype instance by invoking a \p ctrtype overloaded constructor; 
    /// \p objecttype must be a class inheriting from pxl::ObjectBase;
    /// the newly created instance is owned and will be deleted by the object owner. 
    template<class datatype, class ctrdatatype>
    datatype& create(const ctrdatatype& ori)
    { return _objects.create<datatype,ctrdatatype>(ori); }

    // crateIndexed
    /// This method template acts like create() and registers the newly created instance under \p idx in the index.
    template<class datatype>
    datatype& createIndexed(const std::string& idx)
    { datatype& obj = _objects.create<datatype>(); setIndex(idx, obj); return obj; }

    /// This method template acts like create() and registers the newly created instance under \p idx in the index.
    template<class datatype, class ctrdatatype>
    datatype& createIndexed(const ctrdatatype& ori, const std::string& idx)
    { datatype& obj = _objects.create<datatype,ctrdatatype>(ori); setIndex(idx, obj); return obj; }


    /// This method inserts \p obj in the container of the object owner and takes deletion responsability.
    inline void setObject(pxl::ObjectBase& obj, const std::string& idx)
    {_objects.set(obj); setIndex(idx, obj);}

    /// This method registers the object \p obj with the index-id \p idx in the index and returns true in case of success;
    /// please notice, that obj must be owned by this object owner and \p idx must not be a zero length string.  
    inline bool setIndex(const std::string& idx, pxl::ObjectBase& obj)
    { return _objects.setIndex(idx, obj); }

    /// This method provides access to the object owner.
    inline const pxl::Objects& getObjects() const
    { return _objects; }

    /// This method deletes the object \p obj.
    inline void removeObject(pxl::ObjectBase& obj)
    { _objects.remove(obj); }

    /// This method clears the object owner and deletes all owned objects. 
    inline void clearObjects()
    { _objects.clearContainer(); }


    /// This method searches the index for the index-id \p idx and returns a dynamically casted 
    /// C++ pointer of type \p objecttype* to the corresponding object; 
    /// in case idx is not found a null pointer is returned.
    template<class objecttype>
    inline objecttype* findObject(const std::string idx) const
    { return _objects.findObject<objecttype>(idx); }

    /// This method searches the copy history to locate the copy of \p original and 
    /// returns a dynamically casted C++ pointer of type \p objecttype* to the corresponding copy; 
    /// in case no copy can be traced a null pointer is returned.
    template<class objecttype>
    inline objecttype* findCopyOf(const pxl::ObjectBase& original) const
    { return _objects.findCopyOf<objecttype>(original); }


    /// This method provides direct access to the copy history (created by the copy constructor). 
    inline const pxl::CopyHistory& getCopyHistory() const
    { return _objects.getCopyHistory(); }

    /// This method clears the copy history  (created by the copy constructor). 
    inline void clearCopyHistory()
    { _objects.clearCopyHistory(); }


    /// This method provides direct access to the index. 
    inline const pxl::Index& getIndex() const
    { return _objects.getIndex(); }

    /// This method removes the index entry with index-id \p idx; please notice: it does not remove the object itself. 
    inline void removeIndex(const std::string& idx)
    { _objects.removeIndex(idx); }

    /// This method clears the index; please notice: it does not remove the objects themself.
    inline void clearIndex()
    { _objects.clearIndex(); }

  protected:
    pxl::Objects _objects;

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;      
};

/// This typedef defines a basic PXL physics object with pxl::BasicObjectManagerData as data.
typedef pxl::Object<pxl::BasicObjectManagerData> BasicObjectManager; // NOTICE: BasicObjectManagers cannot be managed by Copy On Write!
/// This typedef defines a weak pointer for pxl::BasicObjectManager
typedef pxl::WkPtrSpec<pxl::BasicObjectManagerData, pxl::BasicObjectManager> BasicObjectManagerWkPtr;
/// This typedef defines a reference for pxl::BasicObjectManager
typedef BasicObjectManager& BasicObjectManagerRef;
/// This typedef defines a const reference for pxl::BasicObjectManager
typedef const BasicObjectManager& BasicObjectManagerConstRef;

} // namespace pxl

iotl__declareObjectTypeProto(pxl::BasicObjectManager)
iotl__declareDataTypeProto(pxl::BasicObjectManagerData)

#endif // pxl_pol_BasicManager_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_pol_Vertex_hh
#define pxl_pol_Vertex_hh




namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// pol
/**
This class allows to store threevector and further properties of the decay vertex; see also pxl::BasicObjectData.
*/
class VertexData : public pxl::BasicObjectData {
  public:
    /// This method grants read access to the vector. 
    inline const pxl::Basic3VectorData& vector()                const { return _vector; }
    /// This method grants read access to the vector. 
    inline const pxl::Basic3VectorData& vector(const pxl::Get&) const { return _vector; }
    /// This method grants write access to the vector. 
    inline       pxl::Basic3VectorData& vector(const pxl::Set&)       { return _vector; }

    /// This method adds the vector of \p vxd. 
    inline const pxl::VertexData& operator+=(const pxl::VertexData& vxd) { this->vector(pxl::set) += vxd.vector(); return *this; }
    /// This method subtracts the vector of \p vxd. 
    inline const pxl::VertexData& operator-=(const pxl::VertexData& vxd) { this->vector(pxl::set) -= vxd.vector(); return *this; }

  protected:
    pxl::Basic3VectorData _vector;

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;
};

// non-member operators
bool const operator==(const pxl::VertexData& obj1, const pxl::VertexData& obj2);
bool const operator!=(const pxl::VertexData& obj1, const pxl::VertexData& obj2);

// typedefs
/**
This typedef represents decay vertices, thus spatial points of particle decays;
data is aggregated in pxl::VertexData.
It is intended to store three-vector and further properties of the decay vertex
and to establish relations to mother and daughter pxl::Particle objects, for instance.
*/
typedef pxl::CowObject<pxl::VertexData> Vertex;
/// This typedef defines a weak pointer for pxl::Vertex
typedef pxl::WkPtrSpec<pxl::VertexData, pxl::Vertex> VertexWkPtr;
/// This typedef defines a reference for pxl::Vertex
typedef Vertex& VertexRef;
/// This typedef defines a const reference for pxl::Vertex
typedef const Vertex& VertexConstRef;

template<> std::ostream& CowObject<pxl::VertexData>::print(int level, std::ostream& os, int pan) const;

} // namespace pxl

iotl__declareObjectTypeProto(pxl::Vertex)
iotl__declareDataTypeProto(pxl::VertexData)

#endif // pxl_pol_Vertex_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_pol_Particle_hh
#define pxl_pol_Particle_hh



namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// pol
/**
This class allows to store Lorentz-fourvector and further properties of particles or reconstructed objects
such as charge, particle-id plus the inherited properties of pxl::BasicObjectData.
*/
class ParticleData : public pxl::BasicObjectData {
  public: 
    ParticleData() :
        BasicObjectData(), _vector(), _charge(0), _particleId(0) {}

    /// This method grants read access to the vector. 
    inline const pxl::Basic4VectorData& vector()                const { return _vector; }
     /// This method grants read access to the vector. 
    inline const pxl::Basic4VectorData& vector(const pxl::Get&) const { return _vector; }
    /// This method grants write access to the vector. 
    inline       pxl::Basic4VectorData& vector(const pxl::Set&)       { return _vector; }

    /// This method returns the particle charge.
    inline double getCharge() const { return _charge; }
    /// This method sets the particle charge to v.
    inline void setCharge(double v) { _charge = v; }

    /// This method returns the particle-id.
    inline int getParticleId() const { return _particleId; }
    /// This method sets the particle-id to \p v.
    inline void setParticleId(int v) { _particleId = v; }

    /// This method adds vector and charge of \p pad. 
    inline const pxl::ParticleData& operator+=(const pxl::ParticleData& pad)
    { this->vector(pxl::set) += pad.vector(); _charge += pad._charge; return *this; }

    /// This method subtracts vector and charge of of \p pad. 
    inline const pxl::ParticleData& operator-=(const pxl::ParticleData& pad)
    { this->vector(pxl::set) -= pad.vector(); _charge += pad._charge; return *this; }

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
/**
This typedef represents particles and reconstructed 
objects such as muons, electrons, photons, jets; data is aggregated in pxl::ParticleData.
*/
typedef pxl::CowObject<pxl::ParticleData> Particle;
/// This typedef defines a weak pointer for pxl::Particle
typedef pxl::WkPtrSpec<pxl::ParticleData, pxl::Particle> ParticleWkPtr;
/// This typedef defines a reference for pxl::Particle
typedef Particle& ParticleRef;
/// This typedef defines a const reference for pxl::Particle
typedef const Particle& ParticleConstRef;

template<> std::ostream& CowObject<pxl::ParticleData>::print(int level, std::ostream& os, int pan) const;

} // namespace pxl

iotl__declareObjectTypeProto(pxl::Particle)
iotl__declareDataTypeProto(pxl::ParticleData)

#endif // pxl_pol_Particle_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_pol_Collision_hh
#define pxl_pol_Collision_hh



namespace pxl {

/**
The data aggregated in pxl::CollisionData is identical to pxl::CollisionData;
it allows the separation of different collisions in multicollision events 
(as they occur at high-rate hadron colliders) by providing the relation management 
necessary to associate pxl::Vertex or pxl::Particle objects, for instance.
*/
typedef pxl::BasicObjectData CollisionData;

/**
This typedef represents individual interactions in  multicollision events; 
data is aggregated in pxl::CollisionData (= pxl::BasicObjectData).
It allows the separation of different collisions as they occur 
at high-rate hadron colliders by providing the relation management 
necessary to associate pxl::Vertex or pxl::Particle objects, for instance. 
*/
typedef pxl::CowObject<pxl::CollisionData> Collision;
/// This typedef defines a weak pointer for pxl::Collision
typedef pxl::WkPtrSpec<pxl::CollisionData, pxl::Collision> CollisionWkPtr;
/// This typedef defines a reference for pxl::Collision
typedef Collision& CollisionRef;
/// This typedef defines a const reference for pxl::Collision
typedef const Collision& CollisionConstRef;

template<> std::ostream& CowObject<pxl::CollisionData>::print(int level, std::ostream& os, int pan) const;

} // namespace pxl

#endif // pxl_pol_Collision_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_pol_EventView_hh
#define pxl_pol_EventView_hh



namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// pol
/**
By inheritance from pxl::BasicObjectManagerData and pxl::BasicObjectData, 
this class is capable of holding the complete information of one 
multicollision event with decay trees, spatial vertex information, 
four-momenta as well as additional event-related reconstruction data 
in the user records. Physics objects (i.e. instances of the classes pxl::Particle, 
pxl::Vertex or pxl::Collision) as well as other arbitrary pxl::ObjectBase 
derivatives can be aggregated and managed.
The name 'event view'  arises from the fact, that it is 
intended to represent a distinct view of an event (e.g. connecting 
particles to the decay tree according to one out of a
number of hypotheses, applying different jet energy corrections, etc.). 
To facilitate the development of numerous 
parallel or subsequent event views, as needed for hypothesis evolution, 
for instance, this class features a copy constructor, 
which provides a deep copy of the event container with all data members, 
physics objects, and their (redirected) relations. 
Please remember, that PXL physics objects are based on the copy-on-write 
mechanism. 
Although behaving like deep copies, it is only upon write access, 
that the individual CPU and memory intense copy processes are carried out.
This way, the PXL provides a flexible generalized event container 
with comfortable high-performance duplication functionality, 
meeting the needs of HEP analyses in channels with ambiguous 
event topologies.
*/
class EventViewData : public pxl::BasicObjectManagerData {
  public:
    EventViewData()
        : BasicObjectManagerData() {}
    /// This copy constructor provides a deep copy of the event container \p original with all data members, 
    /// physics objects, and their (redirected) relations. 
    EventViewData(const pxl::EventViewData& original)
        : BasicObjectManagerData(original) {}

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;      
};

/**
This typedef represents the generalized event container provided by PXL; 
data is aggregated in pxl::EventViewData. 
The event view is capable of holding the complete information of one 
multicollision event with decay trees, spatial vertex information, 
four-momenta as well as additional event-related reconstruction data 
in the user records. Physics objects (i.e. instances of the classes pxl::Particle, 
pxl::Vertex or pxl::Collision) as well as other arbitrary pxl::ObjectBase 
derivatives can be aggregated and managed.
The name 'event view'  arises from the fact, that it is 
intended to represent a distinct view of an event (e.g. connecting 
particles to the decay tree according to one out of a
number of hypotheses, applying different jet energy corrections, etc.). 
To facilitate the development of numerous 
parallel or subsequent event views, as needed for hypothesis evolution, 
for instance, the pxl::EventViewData class features a copy constructor, 
which provides a deep copy of the event container with all data members, 
physics objects, and their (redirected) relations. 
Please remember, that PXL physics objects are based on the copy-on-write 
mechanism. 
Although behaving like deep copies, it is only upon write access, 
that the individual CPU and memory intense copy processes are carried out.
This way, the PXL provides a flexible generalized event container 
with comfortable high-performance duplication functionality, 
meeting the needs of HEP analyses in channels with ambiguous 
event topologies.
*/ 
typedef pxl::Object<pxl::EventViewData> EventView; // NOTICE: EventViews cannot be managed by Copy On Write!
/// This typedef defines a weak pointer for pxl::EventView
typedef pxl::WkPtrSpec<pxl::EventViewData, pxl::EventView> EventViewWkPtr;
/// This typedef defines a reference for pxl::EventView
typedef EventView& EventViewRef;
/// This typedef defines a const reference for pxl::EventView
typedef const EventView& EventViewConstRef;

template<> std::ostream& Object<pxl::EventViewData>::print(int level, std::ostream& os, int pan) const;

} // namespace pxl

iotl__declareObjectTypeProto(pxl::EventView)
iotl__declareDataTypeProto(pxl::EventViewData)

#endif // pxl_pol_EventView_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_pol_AnalysisProcess_hh
#define pxl_pol_AnalysisProcess_hh



namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// pol

/**
This class is designed for assisting the evolution of different 
combinatorial hypotheses of an event according to a certain physics 
process. It is derived from pxl::BasicObjectManagerData and 
pxl::BasicObjectData and the  analyzer can store and manage an 
arbitrary number of event views including physics objects and 
relations.  
*/ 
class AnalysisProcessData : public pxl::BasicObjectManagerData {
  public:
    AnalysisProcessData() :
        BasicObjectManagerData() {}
    AnalysisProcessData(const pxl::AnalysisProcessData& original) :
        BasicObjectManagerData(original) {}

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;  
};

/**
This class is designed as a base class to assist the analyzer in the 
evolution of different combinatorial hypotheses of an event according 
to a certain physics process; data is aggregated in pxl::AnalysisProcessData.
This class provides virtual methods to be called at the beginning and end 
of a job, at the beginning and end of a run, and, of course, at event analysis 
and event finishing time (just as needed in a stand-alone analysis framework, 
for instance). When inheriting from this class, the analyst can
place user code in the according reimplementations of these methods. 
*/ 
class AnalysisProcess : public pxl::Object<pxl::AnalysisProcessData> {
	// NOTICE: AnalysisProcesses cannot be managed by Copy On Write!
  public: 
    AnalysisProcess() :
        pxl::Object<pxl::AnalysisProcessData>() {}
    AnalysisProcess(const pxl::AnalysisProcess& original) :
        pxl::Object<pxl::AnalysisProcessData>(original) {}
    virtual ~AnalysisProcess() {}

    /// This method can be reimplemented to build/destroy 
    /// a static template for user-defined tree creation. \p mode is a freely usable parameter.
    virtual void buildTemplate(int mode = 0) {}

    /// This method can be reimplemented to hold physics analysis code executed at the begin of a computing job 
    /// (as needed for histogram booking etc.).
    /// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::Objects instance (that might carry the reconstructed event data or generator information).  
    virtual void beginJob(const pxl::Objects* input = 0) {}
    /// This method can be reimplemented to hold physics analysis code executed at the begin of a run. 
    /// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::Objects instance (that might carry the reconstructed event data or generator information).  
    virtual void beginRun(const pxl::Objects* input = 0) {}
    /// This method can be reimplemented to hold physics analysis code executed for the actual event analysis. 
    /// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::Objects instance (that might carry the reconstructed event data or generator information).  
    virtual void analyseEvent(const pxl::Objects* input = 0) {}
    /// This method can be reimplemented to hold physics analysis code executed at the end of each event.
    /// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::Objects instance (that might carry the reconstructed event data or generator information).  
    /// By default, this method clears the object owner and deletes all owned objects. 
    virtual void finishEvent(const pxl::Objects* input = 0) { set().clearObjects(); }
    /// This method can be reimplemented to hold physics analysis code executed at the end of a run. 
    /// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::Objects instance (that might carry the reconstructed event data or generator information).  
    virtual void endRun(const pxl::Objects* input = 0) {}
    /// This method can be reimplemented to hold physics analysis code executed at the end of a computing job 
    /// (as needed for histogram storing etc.). 
    /// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::Objects instance (that might carry the reconstructed event data or generator information).  
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

    virtual pxl::WkPtrBase* createSelfWkPtr();

  protected:
    virtual void storeYourSelf(pxl::oStreamer& output) const;
};

/// This typedef defines a weak pointer for pxl::AnalysisProcess
typedef pxl::WkPtrSpec<pxl::AnalysisProcessData, pxl::AnalysisProcess> AnalysisProcessWkPtr;
/// This typedef defines a reference for pxl::AnalysisProcess
typedef AnalysisProcess& AnalysisProcessRef;
/// This typedef defines a const reference for pxl::AnalysisProcess
typedef const AnalysisProcess& AnalysisProcessConstRef;

template<> std::ostream& Object<pxl::AnalysisProcessData>::print(int level, std::ostream& os, int pan) const;

} // namespace pxl

iotl__declareObjectTypeProto(pxl::AnalysisProcess)
iotl__declareDataTypeProto(pxl::AnalysisProcessData)

#endif // pxl_pol_AnalysisProcess_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_pol_AnalysisFork_hh
#define pxl_pol_AnalysisFork_hh




namespace pxl {

// iotl
class iStreamer;
class oStreamer;

// pol

/**
This class is designed for assisting the parallel evolution of 
different physics process hypotheses or the analysis of different  
(instrumental) aspects of an event. To allow easy management of 
different physics process hypotheses of an event, this class 
provides storage and easy access to an arbitrary number of 
processes (pxl::AnalysisProcess derivatives) and further 
analysis forks (pxl::AnalysisFork derivatives).
*/ 
class AnalysisForkData : public pxl::BasicObjectManagerData {
  public: 
    AnalysisForkData() :
        BasicObjectManagerData() {}
    AnalysisForkData(const pxl::AnalysisForkData& original) :
        BasicObjectManagerData(original) {}

  friend class pxl::iStreamer;
  friend class pxl::oStreamer;
};

/**
This class is designed as a base class to assist the analyzer in the 
parallel evolution of different physics process hypotheses or the analysis of different  
(instrumental) aspects of an event; data is aggregated in pxl::AnalysisForkData
*/
class AnalysisFork : public pxl::Object<pxl::AnalysisForkData> {
	// NOTICE: AnalysisForks cannot be managed by Copy On Write!
  public: 
    AnalysisFork() :
        pxl::Object<pxl::AnalysisForkData>() {}
    AnalysisFork(const pxl::AnalysisFork& original) :
        pxl::Object<pxl::AnalysisForkData>(original) {}
    virtual ~AnalysisFork() {}

    /// This method can be reimplemented to build/destroy 
    /// a static template for user-defined tree creation. \p mode is a freely usable parameter.
    /// By default, this method invokes the corresponding method of all managed pxl::AnalysisProcess instances.    
    virtual void buildTemplate(int mode = 0);

    /// This method can be reimplemented to hold physics analysis code executed at the begin of a computing job.  
    /// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::Objects instance (that might carry the reconstructed event data or generator information).  
    /// By default, this method invokes the corresponding method of all managed pxl::AnalysisProcess instances.    
    virtual void beginJob(const pxl::Objects*input = 0);
    /// This method can be reimplemented to hold physics analysis code executed at the begin of a run. 
    /// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::Objects instance (that might carry the reconstructed event data or generator information).  
    /// By default, this method invokes the corresponding method of all managed pxl::AnalysisProcess instances.    
    virtual void beginRun(const pxl::Objects*input = 0);
    /// This method can be reimplemented to hold physics analysis code executed for the actual event analysis. 
    /// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::Objects instance (that might carry the reconstructed event data or generator information).  
    /// By default, this method invokes the corresponding method of all managed pxl::AnalysisProcess instances.    
    virtual void analyseEvent(const pxl::Objects*input = 0);
    /// This method can be reimplemented to hold physics analysis code executed at the end of each event.
    /// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::Objects instance (that might carry the reconstructed event data or generator information).  
    /// By default, this method invokes the corresponding method of all managed pxl::AnalysisProcess instances.    
    virtual void finishEvent(const pxl::Objects*input = 0);
    /// This method can be reimplemented to hold physics analysis code executed at the end of a run. 
    /// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::Objects instance (that might carry the reconstructed event data or generator information).  
    /// By default, this method invokes the corresponding method of all managed pxl::AnalysisProcess instances.    
    virtual void endRun(const pxl::Objects*input = 0);
    /// This method can be reimplemented to hold physics analysis code executed at the end of a computing job 
    /// (as needed for histogram storing etc.). 
    /// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::Objects instance (that might carry the reconstructed event data or generator information).  
    /// By default, this method invokes the corresponding method of all managed pxl::AnalysisProcess instances.    
    virtual void endJob(const pxl::Objects*input = 0);

    virtual pxl::ObjectBase* clone() const;

    virtual pxl::WkPtrBase* createSelfWkPtr();

  protected:
    virtual void storeYourSelf(pxl::oStreamer& output) const;
};

/// This typedef defines a weak pointer for pxl::AnalysisFork
typedef pxl::WkPtrSpec<pxl::AnalysisForkData, pxl::AnalysisFork> AnalysisForkWkPtr;
/// This typedef defines a reference for pxl::AnalysisFork
typedef AnalysisFork& AnalysisForkRef;
/// This typedef defines a const reference for pxl::AnalysisFork
typedef const AnalysisFork& AnalysisForkConstRef;

template<> std::ostream& Object<pxl::AnalysisForkData>::print(int level, std::ostream& os, int pan) const;

} // namespace pxl

iotl__declareObjectTypeProto(pxl::AnalysisFork)
iotl__declareDataTypeProto(pxl::AnalysisForkData)

#endif // pxl_pol_AnalysisFork_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_pol_ParticleFilter_hh
#define pxl_pol_ParticleFilter_hh




namespace pxl {

/** 
This class provides a pT-sorted filter for PXL physics objects (requires the name, pT and |eta| to match).
*/ 
class ParticleFilter : public pxl::Filter<pxl::Particle, double> {
  public:
    /// This constructor defines filter name, minimum pT and maximum |eta| value and runs the filter on the \p objects container. 
    ParticleFilter(const pxl::ObjectOwner& objects, const std::string& name,
                   double ptMin = 0.0, double etaMax = 0.0);

    virtual bool pass(const pxl::Particle& pa) const;
    virtual double sort(const pxl::Particle& pa) const;

  private: 
    std::string _name;
    double _ptMin;
    double _etaMax;
};

} // namespace pxl

#endif // pxl_pol_ParticleFilter_hh
//------------------------------------------------------------------------------
//
//    Physics eXtension Library (PXL)
//    C++ Toolkit for Fourvector Analysis, Relation Management 
//    and Hypothesis Evolution in High Energy Physics
//    Copyright (C) 2006-2007  S. Kappler, G. Mueller, C. Saout
//    E-mail contact: project-PXL@cern.ch
//    
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, 5th Floor, Boston, MA 02110-1301 USA
//
//------------------------------------------------------------------------------
#ifndef pxl_pol_types_hh
#define pxl_pol_types_hh



namespace pxl {
	
/// This typedef provides a PXL-style Iterator for pxl::Objects
typedef Objects::Iterator ObjectIterator;
/// This typedef provides a type selective PXL-style Iterator for pxl::Objects
typedef Objects::TypeIterator<Particle> ParticleIterator;
/// This typedef provides a type selective PXL-style Iterator for pxl::Objects
typedef Objects::TypeIterator<Vertex> VertexIterator;
/// This typedef provides a type selective PXL-style Iterator for pxl::Objects
typedef Objects::TypeIterator<Collision> CollisionIterator;
/// This typedef provides a type selective PXL-style Iterator for pxl::Objects
typedef Objects::TypeIterator<EventView> EventViewIterator;
/// This typedef provides a type selective PXL-style Iterator for pxl::Objects
typedef Objects::TypeIterator<AnalysisProcess> AnalysisProcessIterator;
/// This typedef provides a type selective PXL-style Iterator for pxl::Objects
typedef Objects::TypeIterator<AnalysisFork> AnalysisForkIterator;

/// This typedef provides a PXL-style Iterator for pxl::Relations
typedef Relations::Iterator RelationIterator;
/// This typedef provides a type selective PXL-style Iterator for pxl::Relations
typedef Relations::TypeIterator<ParticleWkPtr> ParticleRelationIterator;
/// This typedef provides a type selective PXL-style Iterator for pxl::Relations
typedef Relations::TypeIterator<VertexWkPtr> VertexRelationIterator;
/// This typedef provides a type selective PXL-style Iterator for pxl::Relations
typedef Relations::TypeIterator<CollisionWkPtr> CollisionRelationIterator;
/// This typedef provides a type selective PXL-style Iterator for pxl::Relations
typedef Relations::TypeIterator<EventViewWkPtr> EventViewRelationIterator;
/// This typedef provides a type selective PXL-style Iterator for pxl::Relations
typedef Relations::TypeIterator<AnalysisProcessWkPtr> AnalysisProcessRelationIterator;
/// This typedef provides a type selective PXL-style Iterator for pxl::Relations
typedef Relations::TypeIterator<AnalysisForkWkPtr> AnalysisForkRelationIterator;

/// This typedef provides a type selective PXL-style Iterator for pxl::Filter
typedef WkPtrOwner<int>::TypeIterator<ParticleWkPtr> ParticleFilterIterator;
/// This typedef provides a type selective PXL-style Iterator for pxl::Filter
typedef WkPtrOwner<int>::TypeIterator<VertexWkPtr> VertexFilterIterator;
/// This typedef provides a type selective PXL-style Iterator for pxl::Filter
typedef WkPtrOwner<int>::TypeIterator<CollisionWkPtr> CollisionFilterIterator;
/// This typedef provides a type selective PXL-style Iterator for pxl::Filter
typedef WkPtrOwner<int>::TypeIterator<EventViewWkPtr> EventViewFilterIterator;
/// This typedef provides a type selective PXL-style Iterator for pxl::Filter
typedef WkPtrOwner<int>::TypeIterator<AnalysisProcessWkPtr> AnalysisProcessFilterIterator;
/// This typedef provides a type selective PXL-style Iterator for pxl::Filter
typedef WkPtrOwner<int>::TypeIterator<AnalysisForkWkPtr> AnalysisForkFilterIterator;

} // namespace pxl

#endif // pxl_pol_types_hh

#endif // pxl_pol_hh

#endif // pxl_hh
