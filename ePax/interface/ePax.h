//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_hh
#define pxl_hh

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_base_hh
#define pxl_base_hh

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

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
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_pcl_VariantBase_hh
#define pxl_pcl_VariantBase_hh

#include <cstdlib>
#include <vector>
#include <string>
#include <stdexcept>

#include <string.h>


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

    inline VariantBase() { memset(&v, 0, sizeof v); }
    inline VariantBase(const pxl::VariantBase& orig) : v(orig.v) {}
    inline VariantBase(const pxl::VariantBase* orig) : v(orig->v) {}

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
        memset(&v, 0, sizeof v);
    }

    /// This method ensures correct duplication of possible internal memory allocations after a copy of the variant, depending on the type given by \p t.
    inline void dup(Type t)
    { 
        if (t == TYPE_STRING) {
            std::size_t len = strlen((char*)v.p) + 1;
            char* copy = new char[len];
            memcpy(copy, v.p, len);
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
    void wrongType(Type tShould, Type tIs) const throw (std::runtime_error);

  private:
    static const TypeInfo& fallbackGetTypeInfo(Type t) throw (std::runtime_error);

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
    memcpy(v.p, arg.c_str(), size);
    ((char*)v.p)[size] = 0;
}

template<>
inline VariantBase::Type VariantBase::findType<std::string>()
{ return TYPE_STRING; }

#undef PCL_ANYBASE_CHECK

} // namespace pxl

#endif // pxl_pcl_VariantBase_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

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
    inline Variant(const pxl::Variant* orig) :
            VariantBase(orig), type(orig->type) { dup(type); }
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
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_pcl_functions_hh
#define pxl_pcl_functions_hh


namespace pxl {

/// @internal This function returns a platform-specific CPU timestamp; internally used for performance tests.
double getCpuTime();

} // namespace pxl

#endif // pxl_pcl_functions_hh

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_ptl_MutableId_hh
#define pxl_ptl_MutableId_hh

// MersenneTwister.h
// Mersenne Twister random number generator -- a C++ class MTRand
// Based on code by Makoto Matsumoto, Takuji Nishimura, and Shawn Cokus
// Richard J. Wagner  v1.0  15 May 2003  rjwagner@writeme.com

// The Mersenne Twister is an algorithm for generating random numbers.  It
// was designed with consideration of the flaws in various other generators.
// The period, 2^19937-1, and the order of equidistribution, 623 dimensions,
// are far greater.  The generator is also fast; it avoids multiplication and
// division, and it benefits from caches and pipelines.  For more information
// see the inventors' web page at http://www.math.keio.ac.jp/~matumoto/emt.html

// Reference
// M. Matsumoto and T. Nishimura, "Mersenne Twister: A 623-Dimensionally
// Equidistributed Uniform Pseudo-Random Number Generator", ACM Transactions on
// Modeling and Computer Simulation, Vol. 8, No. 1, January 1998, pp 3-30.

// Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
// Copyright (C) 2000 - 2003, Richard J. Wagner
// All rights reserved.                          
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
//   1. Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//
//   2. Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//   3. The names of its contributors may not be used to endorse or promote 
//      products derived from this software without specific prior written 
//      permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// The original code included the following notice:
//
//     When you use this, send an email to: matumoto@math.keio.ac.jp
//     with an appropriate reference to your work.
//
// It would be nice to CC: rjwagner@writeme.com and Cokus@math.washington.edu
// when you write.

#ifndef MERSENNETWISTER_H
#define MERSENNETWISTER_H

// Not thread safe (unless auto-initialization is avoided and each thread has
// its own MTRand object)

#include <iostream>
#include <limits.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

class MTRand {
// Data
public:
	typedef unsigned long uint32;  // unsigned integer type, at least 32 bits
	
	enum { N = 624 };       // length of state vector
	enum { SAVE = N + 1 };  // length of array for save()

protected:
	enum { M = 397 };  // period parameter
	
	uint32 state[N];   // internal state
	uint32 *pNext;     // next value to get from state
	int left;          // number of values left before reload needed


//Methods
public:
	MTRand( const uint32& oneSeed );  // initialize with a simple uint32
	MTRand( uint32 *const bigSeed, uint32 const seedLength = N );  // or an array
	MTRand();  // auto-initialize with /dev/urandom or time() and clock()
	
	// Do NOT use for CRYPTOGRAPHY without securely hashing several returned
	// values together, otherwise the generator state can be learned after
	// reading 624 consecutive values.
	
	// Access to 32-bit random numbers
	double rand();                          // real number in [0,1]
	double rand( const double& n );         // real number in [0,n]
	double randExc();                       // real number in [0,1)
	double randExc( const double& n );      // real number in [0,n)
	double randDblExc();                    // real number in (0,1)
	double randDblExc( const double& n );   // real number in (0,n)
	uint32 randInt();                       // integer in [0,2^32-1]
	uint32 randInt( const uint32& n );      // integer in [0,n] for n < 2^32
	double operator()() { return rand(); }  // same as rand()
	
	// Access to 53-bit random numbers (capacity of IEEE double precision)
	double rand53();  // real number in [0,1)
	
	// Access to nonuniform random number distributions
	double randNorm( const double& mean = 0.0, const double& variance = 0.0 );
	
	// Re-seeding functions with same behavior as initializers
	void seed( const uint32 oneSeed );
	void seed( uint32 *const bigSeed, const uint32 seedLength = N );
	void seed();
	
	// Saving and loading generator state
	void save( uint32* saveArray ) const;  // to array of size SAVE
	void load( uint32 *const loadArray );  // from such array
	friend std::ostream& operator<<( std::ostream& os, const MTRand& mtrand );
	friend std::istream& operator>>( std::istream& is, MTRand& mtrand );

protected:
	void initialize( const uint32 oneSeed );
	void reload();
	uint32 hiBit( const uint32& u ) const { return u & 0x80000000UL; }
	uint32 loBit( const uint32& u ) const { return u & 0x00000001UL; }
	uint32 loBits( const uint32& u ) const { return u & 0x7fffffffUL; }
	uint32 mixBits( const uint32& u, const uint32& v ) const
		{ return hiBit(u) | loBits(v); }
	uint32 twist( const uint32& m, const uint32& s0, const uint32& s1 ) const
		{ return m ^ (mixBits(s0,s1)>>1) ^ (-loBit(s1) & 0x9908b0dfUL); }
	static uint32 hash( time_t t, clock_t c );
};


inline MTRand::MTRand( const uint32& oneSeed )
	{ seed(oneSeed); }

inline MTRand::MTRand( uint32 *const bigSeed, const uint32 seedLength )
	{ seed(bigSeed,seedLength); }

inline MTRand::MTRand()
	{ seed(); }

inline double MTRand::rand()
	{ return double(randInt()) * (1.0/4294967295.0); }

inline double MTRand::rand( const double& n )
	{ return rand() * n; }

inline double MTRand::randExc()
	{ return double(randInt()) * (1.0/4294967296.0); }

inline double MTRand::randExc( const double& n )
	{ return randExc() * n; }

inline double MTRand::randDblExc()
	{ return ( double(randInt()) + 0.5 ) * (1.0/4294967296.0); }

inline double MTRand::randDblExc( const double& n )
	{ return randDblExc() * n; }

inline double MTRand::rand53()
{
	uint32 a = randInt() >> 5, b = randInt() >> 6;
	return ( a * 67108864.0 + b ) * (1.0/9007199254740992.0);  // by Isaku Wada
}

inline double MTRand::randNorm( const double& mean, const double& variance )
{
	// Return a real number from a normal (Gaussian) distribution with given
	// mean and variance by Box-Muller method
	double r = sqrt( -2.0 * log( 1.0-randDblExc()) ) * variance;
	double phi = 2.0 * 3.14159265358979323846264338328 * randExc();
	return mean + r * cos(phi);
}

inline MTRand::uint32 MTRand::randInt()
{
	// Pull a 32-bit integer from the generator state
	// Every other access function simply transforms the numbers extracted here
	
	if( left == 0 ) reload();
	--left;
		
	register uint32 s1;
	s1 = *pNext++;
	s1 ^= (s1 >> 11);
	s1 ^= (s1 <<  7) & 0x9d2c5680UL;
	s1 ^= (s1 << 15) & 0xefc60000UL;
	return ( s1 ^ (s1 >> 18) );
}

inline MTRand::uint32 MTRand::randInt( const uint32& n )
{
	// Find which bits are used in n
	// Optimized by Magnus Jonsson (magnus@smartelectronix.com)
	uint32 used = n;
	used |= used >> 1;
	used |= used >> 2;
	used |= used >> 4;
	used |= used >> 8;
	used |= used >> 16;
	
	// Draw numbers until one is found in [0,n]
	uint32 i;
	do
		i = randInt() & used;  // toss unused bits to shorten search
	while( i > n );
	return i;
}


inline void MTRand::seed( const uint32 oneSeed )
{
	// Seed the generator with a simple uint32
	initialize(oneSeed);
	reload();
}


inline void MTRand::seed( uint32 *const bigSeed, const uint32 seedLength )
{
	// Seed the generator with an array of uint32's
	// There are 2^19937-1 possible initial states.  This function allows
	// all of those to be accessed by providing at least 19937 bits (with a
	// default seed length of N = 624 uint32's).  Any bits above the lower 32
	// in each element are discarded.
	// Just call seed() if you want to get array from /dev/urandom
	initialize(19650218UL);
	register int i = 1;
	register uint32 j = 0;
	register int k = ( N > seedLength ? N : seedLength );
	for( ; k; --k )
	{
		state[i] =
			state[i] ^ ( (state[i-1] ^ (state[i-1] >> 30)) * 1664525UL );
		state[i] += ( bigSeed[j] & 0xffffffffUL ) + j;
		state[i] &= 0xffffffffUL;
		++i;  ++j;
		if( i >= N ) { state[0] = state[N-1];  i = 1; }
		if( j >= seedLength ) j = 0;
	}
	for( k = N - 1; k; --k )
	{
		state[i] =
			state[i] ^ ( (state[i-1] ^ (state[i-1] >> 30)) * 1566083941UL );
		state[i] -= i;
		state[i] &= 0xffffffffUL;
		++i;
		if( i >= N ) { state[0] = state[N-1];  i = 1; }
	}
	state[0] = 0x80000000UL;  // MSB is 1, assuring non-zero initial array
	reload();
}


inline void MTRand::seed()
{
	// Seed the generator with an array from /dev/urandom if available
	// Otherwise use a hash of time() and clock() values
	
	// First try getting an array from /dev/urandom
	FILE* urandom = fopen( "/dev/urandom", "rb" );
	if( urandom )
	{
		uint32 bigSeed[N];
		register uint32 *s = bigSeed;
		register int i = N;
		register bool success = true;
		while( success && i-- )
			success = fread( s++, sizeof(uint32), 1, urandom );
		fclose(urandom);
		if( success ) { seed( bigSeed, N );  return; }
	}
	
	// Was not successful, so use time() and clock() instead
	seed( hash( time(NULL), clock() ) );
}


inline void MTRand::initialize( const uint32 seed )
{
	// Initialize generator state with seed
	// See Knuth TAOCP Vol 2, 3rd Ed, p.106 for multiplier.
	// In previous versions, most significant bits (MSBs) of the seed affect
	// only MSBs of the state array.  Modified 9 Jan 2002 by Makoto Matsumoto.
	register uint32 *s = state;
	register uint32 *r = state;
	register int i = 1;
	*s++ = seed & 0xffffffffUL;
	for( ; i < N; ++i )
	{
		*s++ = ( 1812433253UL * ( *r ^ (*r >> 30) ) + i ) & 0xffffffffUL;
		r++;
	}
}


inline void MTRand::reload()
{
	// Generate N new values in state
	// Made clearer and faster by Matthew Bellew (matthew.bellew@home.com)
	register uint32 *p = state;
	register int i;
	for( i = N - M; i--; ++p )
		*p = twist( p[M], p[0], p[1] );
	for( i = M; --i; ++p )
		*p = twist( p[M-N], p[0], p[1] );
	*p = twist( p[M-N], p[0], state[0] );

	left = N, pNext = state;
}


inline MTRand::uint32 MTRand::hash( time_t t, clock_t c )
{
	// Get a uint32 from t and c
	// Better than uint32(x) in case x is floating point in [0,1]
	// Based on code by Lawrence Kirby (fred@genesis.demon.co.uk)

	static uint32 differ = 0;  // guarantee time-based seeds will change

	uint32 h1 = 0;
	unsigned char *p = (unsigned char *) &t;
	for( size_t i = 0; i < sizeof(t); ++i )
	{
		h1 *= UCHAR_MAX + 2U;
		h1 += p[i];
	}
	uint32 h2 = 0;
	p = (unsigned char *) &c;
	for( size_t j = 0; j < sizeof(c); ++j )
	{
		h2 *= UCHAR_MAX + 2U;
		h2 += p[j];
	}
	return ( h1 + differ++ ) ^ h2;
}


inline void MTRand::save( uint32* saveArray ) const
{
	register uint32 *sa = saveArray;
	register const uint32 *s = state;
	register int i = N;
	for( ; i--; *sa++ = *s++ ) {}
	*sa = left;
}


inline void MTRand::load( uint32 *const loadArray )
{
	register uint32 *s = state;
	register uint32 *la = loadArray;
	register int i = N;
	for( ; i--; *s++ = *la++ ) {}
	left = *la;
	pNext = &state[N-left];
}


inline std::ostream& operator<<( std::ostream& os, const MTRand& mtrand )
{
	register const MTRand::uint32 *s = mtrand.state;
	register int i = mtrand.N;
	for( ; i--; os << *s++ << "\t" ) {}
	return os << mtrand.left;
}


inline std::istream& operator>>( std::istream& is, MTRand& mtrand )
{
	register MTRand::uint32 *s = mtrand.state;
	register int i = mtrand.N;
	for( ; i--; is >> *s++ ) {}
	is >> mtrand.left;
	mtrand.pNext = &mtrand.state[mtrand.N-mtrand.left];
	return is;
}

#endif  // MERSENNETWISTER_H

// Change log:
//
// v0.1 - First release on 15 May 2000
//      - Based on code by Makoto Matsumoto, Takuji Nishimura, and Shawn Cokus
//      - Translated from C to C++
//      - Made completely ANSI compliant
//      - Designed convenient interface for initialization, seeding, and
//        obtaining numbers in default or user-defined ranges
//      - Added automatic seeding from /dev/urandom or time() and clock()
//      - Provided functions for saving and loading generator state
//
// v0.2 - Fixed bug which reloaded generator one step too late
//
// v0.3 - Switched to clearer, faster reload() code from Matthew Bellew
//
// v0.4 - Removed trailing newline in saved generator format to be consistent
//        with output format of built-in types
//
// v0.5 - Improved portability by replacing static const int's with enum's and
//        clarifying return values in seed(); suggested by Eric Heimburg
//      - Removed MAXINT constant; use 0xffffffffUL instead
//
// v0.6 - Eliminated seed overflow when uint32 is larger than 32 bits
//      - Changed integer [0,n] generator to give better uniformity
//
// v0.7 - Fixed operator precedence ambiguity in reload()
//      - Added access for real numbers in (0,1) and (0,n)
//
// v0.8 - Included time.h header to properly support time_t and clock_t
//
// v1.0 - Revised seeding to match 26 Jan 2002 update of Nishimura and Matsumoto
//      - Allowed for seeding with arrays of any length
//      - Added access for real numbers in [0,1) with 53-bit resolution
//      - Added access for real numbers from normal (Gaussian) distributions
//      - Increased overall speed by optimizing twist()
//      - Doubled speed of integer [0,n] generation
//      - Fixed out-of-range number generation on 64-bit machines
//      - Improved portability by substituting literal constants for long enum's
//      - Changed license from GNU LGPL to BSD


//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_io_BinaryStream_hh
#define pxl_io_BinaryStream_hh



namespace pxl
{

typedef unsigned char uint8_t;
typedef signed char int8_t;

typedef unsigned short uint16_t;
typedef signed short int16_t;

typedef unsigned int uint32_t;
typedef signed int int32_t;

// io
/**
 This abstract class serves the internal PXL I/O scheme by defining how basic C++ types
 are written to an output buffer.
 */
class OutputStream
{

public:
	virtual void writeChar(char) const = 0;
	virtual void writeUnsignedChar(unsigned char) const = 0;
	virtual void writeFloat(float) const = 0;
	virtual void writeDouble(double) const = 0;
	virtual void writeShort(short) const = 0;
	virtual void writeUnsignedShort(unsigned short) const = 0;
	virtual void writeInt(int) const = 0;
	virtual void writeUnsignedInt(unsigned int) const = 0;
	virtual void writeLong(long) const = 0;
	virtual void writeUnsignedLong(unsigned long) const = 0;
	virtual void writeBool(bool) const = 0;
	virtual void writeString(const std::string&) const = 0;
};

// iotl
/**
 This class serves the internal PXL I/O scheme by implementing how basic C++ types
 are written to an output buffer.
 */

class BufferOutput : public OutputStream
{
	mutable char *_Buffer;
	mutable size_t _BufferSize;
	mutable size_t _BufferPos;

public:

	BufferOutput() :
		_Buffer(0), _BufferSize(0), _BufferPos(0)
	{

	}

	virtual ~BufferOutput()
	{
		destroy();
	}

	void destroy()
	{
		delete [] _Buffer;
		_Buffer = 0;
		_BufferSize = 0;
		_BufferPos = 0;
	}

	void resize(size_t size) const
	{
		char *newBuffer = new char[size];
		if (_BufferPos > 0)
			memcpy(newBuffer, _Buffer, _BufferPos > size ? size : _BufferPos);
		delete[] _Buffer;
		_Buffer = newBuffer;
		_BufferSize = size;
	}

	void reserve(size_t count) const
	{
		size_t needed = _BufferPos + count;
		if (_BufferSize < needed)
		{
			size_t size = _BufferSize + 4096;
			while (size < needed)
				size += 4096;
			resize(size);
		}
	}

	void writeChar(char c) const
	{
		write(&c, 1);
	}

	void writeUnsignedChar(unsigned char c) const
	{
		write(&c, 1);
	}

	void writeString(const std::string& s) const
	{
		unsigned int len = static_cast<unsigned int>(s.size());
		writeUnsignedInt(len);
		write(s.c_str(), len);
	}

	void writeFloat(float f) const
	{
		reserve(4);
		writeLE(&f, 4);
	}

	void writeDouble(double d) const
	{
		reserve(8);
		writeLE(&d, 8);
	}

	void writeInt(int i) const
	{
		reserve(4);
		writeLE(&i, 4);
	}

	void writeUnsignedInt(unsigned int i) const
	{
		reserve(4);
		writeLE(&i, 4);
	}

	void writeLong(long l) const
	{
		reserve(4);
		writeLE(&l, 4);
	}

	void writeUnsignedLong(unsigned long l) const
	{
		reserve(4);
		writeLE(&l, 4);
	}

	void writeShort(short s) const
	{
		reserve(2);
		writeLE(&s, 2);
	}

	void writeUnsignedShort(unsigned short s) const
	{
		reserve(2);
		writeLE(&s, 2);
	}

	void writeBool(bool b) const
	{
		write(&b, 1);
	}

	void writeLE(const void *data, size_t size) const
	{
		reserve(size);
#ifdef PXL_BIG_ENDIAN
		char* src = data + size;
		char* dst = _Buffer + _BufferPos;

		for (size_t i = 0; i < size; i++)
		{
			*dst = *src;
			dst++;
			src--;
		}
#else
		memcpy(_Buffer + _BufferPos, data, size);
#endif
		_BufferPos += size;
	}

	void write(const void *data, size_t size) const
	{
		reserve(size);
		memcpy(_Buffer + _BufferPos, data, size);
		_BufferPos += size;
	}

	const char *getData() const
	{
		return _Buffer;
	}

	size_t getSize() const
	{
		return _BufferPos;
	}
};

// iotl
/**
 This abstract class serves the internal PXL I/O scheme by defining how basic C++ types
 are read from an input buffer.
 */
class InputStream
{

public:
	virtual bool good() const = 0;
	virtual bool readChar(char &) const = 0;
	virtual bool readUnsignedChar(unsigned char &) const = 0;
	virtual bool readFloat(float &) const = 0;
	virtual bool readDouble(double &) const = 0;
	virtual bool readInt(int &) const = 0;
	virtual bool readUnsignedInt(unsigned int &) const = 0;
	virtual bool readLong(long &) const = 0;
	virtual bool readUnsignedLong(unsigned long &) const = 0;
	virtual bool readShort(short &) const = 0;
	virtual bool readUnsignedShort(unsigned short &) const = 0;
	virtual bool readBool(bool &) const = 0;
	virtual bool readString(std::string &) const = 0;
};

// iotl
/**
 This class serves the internal PXL I/O scheme by implementing how basic C++ types
 are read from an input buffer.
 */
class BufferInput : public InputStream
{
	char *_Buffer;
	size_t _BufferSize;
	mutable size_t _BufferPos;

public:

	BufferInput() :
		_Buffer(0), _BufferSize(0), _BufferPos(0)
	{

	}

	virtual ~BufferInput()
	{
		destroy();
	}

	void destroy()
	{
		delete [] _Buffer;
		_Buffer = 0;
		_BufferSize = 0;
		_BufferPos = 0;
	}

	char *data()
	{
		return _Buffer;
	}

	const char *data() const
	{
		return _Buffer;
	}

	void resize(size_t size)
	{
		char *newBuffer = new char[size];
		if (_BufferPos > 0)
			memcpy(newBuffer, _Buffer, _BufferPos > size ? size : _BufferPos);
		delete [] _Buffer;
		_Buffer = newBuffer;
		_BufferSize = size;
	}

	bool good() const
	{
		return ((_BufferSize - _BufferPos) > 0);
	}

	size_t available() const
	{
		return (_BufferSize - _BufferPos);
	}

	bool readChar(char& c) const
	{
		if (available() < 1)
			return false;

		read(&c, 1);
		return true;
	}

	bool readUnsignedChar(unsigned char& c) const
	{
		if (available() < 1)
			return false;

		read(&c, 1);
		return true;
	}

	bool readString(std::string& s) const
	{
		unsigned int size = 0;
		readUnsignedInt(size);
		if (available() < size)
			return false;
		char* temp = 0;
		temp = new char[size+1];
		temp[size]=0;
		read(temp, size);
		s.assign(temp);
		delete [] temp;
		return true;
	}

	bool readFloat(float& f) const
	{
		if (available() < 4)
			return false;

		readLE(&f, 4);
		return true;
	}

	bool readDouble(double& d) const
	{
		if (available() < 8)
			return false;

		readLE(&d, 8);
		return true;
	}

	bool readInt(int& i) const
	{
		if (available() < 4)
			return false;

		readLE(&i, 4);
		return true;
	}

	bool readUnsignedInt(unsigned int& i) const
	{
		if (available() < 4)
			return false;

		readLE(&i, 4);
		return true;
	}

	bool readLong(long& l) const
	{
		if (available() < 4)
			return false;

		readLE(&l, 4);
		return true;
	}

	bool readUnsignedLong(unsigned long& l) const
	{
		if (available() < 4)
			return false;

		readLE(&l, 4);
		return true;
	}

	bool readShort(short& s) const
	{
		if (available() < 2)
			return false;

		readLE(&s, 2);
		return true;
	}

	bool readUnsignedShort(unsigned short& s) const
	{
		if (available() < 2)
			return false;

		readLE(&s, 2);
		return true;
	}

	bool readBool(bool& b) const
	{
		if (available() < 1)
			return false;

		read(&b, 1);
		return true;
	}

	void readLE(void *data, size_t size) const
	{
#ifdef PXL_BIG_ENDIAN
		char* src = _Buffer + _BufferPos;
		char* dst = data + size;

		for (size_t i = 0; i < size; i++)
		{
			*dst = *src;
			dst++;
			src--;
		}
#else
		memcpy(data, _Buffer + _BufferPos, size);
#endif
		_BufferPos += size;
	}

	void read(void *data, size_t size) const
	{
		memcpy(data, _Buffer + _BufferPos, size);
		_BufferPos += size;
	}
};

}

#endif /*pxl_io_BinaryStream_hh*/

namespace pxl
{

/// This typedef is intended to provide a data type for the PXL unique object-id
class MutableId
{
	unsigned char bytes[16];

public:


	MutableId()
	{
		generate();
	}

	MutableId (const InputStream& in)
	{
		deserialize (in);
	}

	MutableId(const char* id);

	explicit MutableId(const std::string& id);

	MutableId(MTRand& rand)
	{
		generate(rand);
	}

	void generate(MTRand& rand)
	{
		for (int i = 0; i < 4; i++)
		{
			MTRand::uint32 value = rand.randInt();
			char *c = (char *)&value;
			bytes[i*4+0] = c[0];
			bytes[i*4+1] = c[1];
			bytes[i*4+2] = c[2];
			bytes[i*4+3] = c[3];
		}

		/* set version 4 (random)*/
		bytes[7] &= ((1 << 4) - 1);
		bytes[7] |= 4 << 4;

		/* set variant (always DCE 1.1 only) */
		bytes[8] &= ((1 << 7) - 1);
		bytes[8] |= 2 << 6;
	}

	void generate()
	{
		static MTRand rand;
		generate(rand);
	}

	static MutableId create()
	{
		MutableId id;
		id.generate();
		return id;
	}
	
	bool operator ==(const MutableId& id) const
	{
		for (int i = 0; i < 16; i++)
		{
			if (bytes[i] != id.bytes[i])
				return false;
		}

		return true;
	}

	bool operator !=(const MutableId& id) const
	{
		for (int i = 0; i < 16; i++)
		{
			if (bytes[i] != id.bytes[i])
				return true;
		}

		return false;
	}

	void reset()
	{
		for (int j = 0; j < 16; j++)
			bytes[j] = 0;
	}

	std::string toString() const;

	bool operator <(const MutableId& op) const
	{
		for (int i = 0; i < 16; i++)
		{
			if (bytes[i] < op.bytes[i])
				return true;
			else if (bytes[i] > op.bytes[i])
				return false;
		}

		return false;
	}

	friend inline std::ostream& operator <<(std::ostream& os,
			const MutableId &uid);

	void serialize(const OutputStream &out) const
	{
		for (size_t i = 0; i < 16; i++)
			out.writeUnsignedChar(bytes[i]);
	}

	void deserialize(const InputStream &in)
	{
		for (size_t i = 0; i < 16; i++)
			in.readUnsignedChar(bytes[i]);
	}

};

inline std::ostream& operator <<(std::ostream& os, const MutableId &id)
{
	static const char *hex[] =
	{ "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "a", "b", "c", "d",
			"e", "f" };

	for (int i = 0; i < 4; i++)
	{
		os << hex[id.bytes[i] >> 4];
		os << hex[id.bytes[i] & 0x0F];
	}
	os << "-";
	for (int i = 4; i < 6; i++)
	{
		os << hex[id.bytes[i] >> 4];
		os << hex[id.bytes[i] & 0x0F];
	}
	os << "-";
	for (int i = 6; i < 8; i++)
	{
		os << hex[id.bytes[i] >> 4];
		os << hex[id.bytes[i] & 0x0F];
	}
	os << "-";
	for (int i = 8; i < 10; i++)
	{
		os << hex[id.bytes[i] >> 4];
		os << hex[id.bytes[i] & 0x0F];
	}
	os << "-";
	for (int i = 10; i < 16; i++)
	{
		os << hex[id.bytes[i] >> 4];
		os << hex[id.bytes[i] & 0x0F];
	}

	return os;
}

} // namespace pxl

#endif // pxl_ptl_MutableId_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_ptl_Id_hh
#define pxl_ptl_Id_hh


namespace pxl {

/// This typedef is intended to provide a read-only data type for the PXL unique object-id
typedef pxl::MutableId Id;

} // namespace pxl

#endif // pxl_ptl_Id_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_ptl_CopyHistory_hh
#define pxl_ptl_CopyHistory_hh

#include <map>
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_ptl_ObjectBase_hh
#define pxl_ptl_ObjectBase_hh


//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_io_Serializable_hh
#define pxl_io_Serializable_hh


namespace pxl
{

// io
/**
 This class is the abstract base class for all objects to be stored in a PXL I/O file.
 It holds the unique ID of each individual object. In addition, a UUID indicating the class
 type must be implemented in derived classes. (De-) Serialization happens via consecutive calls to the
 serialize/deserialize methods of the base classes.
 */
class Serializable
{
public:

	virtual ~Serializable() {}

	/// Id of the object's type.
	virtual const Id& getTypeId() const = 0;
	
	/// Returns the unique ID of the individual object.
	const Id& getId() const
	{
		return _id;
	}

	/// This method serializes the type ID and unique object ID. Derived classes should extend this method by first calling the base class method.
	virtual void serialize(const OutputStream &out) const
	{
		// Serialize ID of the type.
		getTypeId().serialize(out);
		// Serialize UUID.
		_id.serialize(out);
	}

	/// This method deserializes unique object ID. Derived classes should extend this method by first calling the base class method.
	virtual void deserialize(const InputStream &in)
	{
		// Deserialize uuid;
		_id.deserialize(in);
	}
private:
	/// The unique ID.
	Id _id;
};

}

#endif /*pxl_io_Serializable_hh*/
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_ptl_WkPtrBase_hh
#define pxl_ptl_WkPtrBase_hh


namespace pxl
{

// ptl

class Relative;

/** 
 This base class provides common functionalities for all derived PXL weak pointers. 
 */
class WkPtrBase
{
public:
	virtual ~WkPtrBase()
	{
		connect(0);
	}

	/// This virtual method creates a deep copy and returns a C++ pointer to the newly created weak pointer instance.  
	virtual pxl::WkPtrBase* clone() const
	{
		return new pxl::WkPtrBase(*this);
	}

	/// This method returns a C++ pointer of type pxl::Relative to the referenced object  
	inline pxl::Relative* pointer() const
	{
		return _objectRef;
	}
	/// This method returns true, if the referenced object exists.  
	inline bool valid() const
	{
		return _objectRef != 0;
	}

	/// This allows pointer-like tests if the weak pointer is valid.
	inline operator bool()
	{
		return valid();
	}
	/// This arrow operator de-references the weak pointer.   
	inline pxl::Relative* operator->() const
	{
		return access();
	}
	/// compare the referenced object pointers
	inline bool operator==(pxl::WkPtrBase &other) const
	{
		return (_objectRef == other.pointer());
	}
	/// compare the referenced object pointers   
	inline bool operator!=(pxl::WkPtrBase &other) const
	{
		return (_objectRef != other.pointer());
	}

	/// This method attempts a dynamic cast on the referenced object
	static inline WkPtrBase* cast_dynamic(WkPtrBase* orig)
	{
		return orig;
	}

	// safe access to object
	inline Relative* access() const throw (std::runtime_error)
	{
		if (_objectRef)
			return _objectRef;
		throw std::runtime_error("pxl::WkPtrBase::access(): FATAL: The object you intend to access does not exist!");
		return 0;
	}

protected:
	WkPtrBase() :
		_notifyChainIn(0), _notifyChainOut(0), _objectRef(0)
	{
	}

	void notifyDeleted();

	void connect(pxl::Relative* pointer);

	pxl::WkPtrBase* _notifyChainIn;
	pxl::WkPtrBase* _notifyChainOut;

	pxl::Relative* _objectRef;

	friend class Relative;
};

} // namespace pxl

#endif // pxl_ptl_WkPtrBase_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_ptl_Relations_hh
#define pxl_ptl_Relations_hh

#include <set>


namespace pxl
{

class Relative;

/**
 
 The class pxl::Relative owns three instances of pxl::Relations 
 for managing mother, daughter and flat relations to other pxl::Relative derivatives. 
 For decay tree integrity reasons, the class pxl::Relative allows to 
 establish relations only to objects contained in the same object owner 
 (and in the case of both objects not being contained in owners). 
 */

class Relations
{
public:
	typedef std::set<pxl::Relative*>::const_iterator const_iterator;
	typedef std::set<pxl::Relative*>::iterator iterator;

	Relations()
	{
	}

	virtual ~Relations()
	{
	}

	void serialize(const OutputStream &out) const;

	const std::set<pxl::Relative*>& getContainer() const
	{
		return _relatives;
	}

	bool set(pxl::Relative* relative)
	{
		return (_relatives.insert(relative)).second;
	}

	bool erase(pxl::Relative* relative)
	{
		return _relatives.erase(relative);
	}

	bool has(pxl::Relative* relative) const
	{
		return (_relatives.count(relative) > 0 );
	}

	template<class objecttype> int getObjectsOfType(
			std::vector<objecttype*>& relatives) const
	{
		int size = relatives.size();
		for (const_iterator iter = _relatives.begin(); iter!=_relatives.end(); ++iter)
		{
			objecttype* obj = dynamic_cast<objecttype*>(*iter);
			if (obj)
				relatives.push_back(obj);
		}
		return relatives.size()-size;
	}

	int size() const
	{
		return _relatives.size();
	}

	const_iterator begin() const
	{
		return _relatives.begin();
	}

	iterator begin()
	{
		return _relatives.begin();
	}

	const_iterator end() const
	{
		return _relatives.end();
	}

	iterator end()
	{
		return _relatives.end();
	}

	void clearContainer()
	{
		_relatives.clear();
	}

private:
	std::set<pxl::Relative*> _relatives;
};

} // namespace pxl

#endif // pxl_ptl_Relations_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef SOFTRELATIONS_HH_
#define SOFTRELATIONS_HH_



namespace pxl
{
class Relative;
class ObjectOwner;

class SoftRelations
{
public:
	typedef std::multimap<std::string, pxl::Id>::const_iterator const_iterator;
	typedef std::multimap<std::string, pxl::Id>::iterator iterator;

	void serialize(const OutputStream &out) const
	{
		unsigned int size = _relationsMap.size();
		out.writeUnsignedInt(size);
		for (const_iterator iter = _relationsMap.begin(); iter
				!=_relationsMap.end(); ++iter)
		{
			out.writeString(iter->first);
			iter->second.serialize(out);
		}
	}

	void deserialize(const InputStream &in)
	{
		_relationsMap.clear();
		unsigned int size;
		in.readUnsignedInt(size);
		for (unsigned int i=0; i<size; ++i)
		{
			std::string name;
			in.readString(name);
			pxl::Id id(in);
			_relationsMap.insert(std::pair<std::string, pxl::Id>(name, id));
		}
	}

	pxl::Relative* getFirst(const pxl::ObjectOwner& owner, const std::string& type = "") const;

	int getSoftRelatives(std::vector<pxl::Relative*>& vec, const pxl::ObjectOwner& owner, const std::string& type) const;
	int getSoftRelatives(std::vector<pxl::Relative*>& vec, const pxl::ObjectOwner& owner) const;

	template <class objecttype>
	int getSoftRelativesOfType(std::vector<objecttype*>& vec, const pxl::ObjectOwner& owner, const std::string& type) const;

	template <class objecttype>
	int getSoftRelativesOfType(std::vector<objecttype*>& vec, const pxl::ObjectOwner& owner) const;

	int keepSoftRelatives(std::vector<pxl::Relative*>& vec, const std::string& type) const;
	int keepSoftRelatives(std::vector<pxl::Relative*>& vec) const;

	bool has(const pxl::Relative* relative) const;
	
	bool has(const pxl::Relative* relative, const std::string& name) const;

	bool has(const pxl::Id& id) const;
	bool has(const pxl::Id& id, const std::string& name) const;

	bool hasType(const std::string& name) const;

	int count(const std::string& name) const;

	void set(const pxl::Relative* relative, const std::string& type);

	void remove(const pxl::Relative* relative, const std::string& type);

	void remove(const pxl::Relative* relative);

	int size() const
	{
		return _relationsMap.size();
	}

	void clearContainer()
	{
		_relationsMap.clear();
	}

	const std::multimap<std::string, pxl::Id>& getContainer() const
	{
		return _relationsMap;
	}

	const_iterator begin() const
	{
		return _relationsMap.begin();
	}

	iterator begin()
	{
		return _relationsMap.begin();
	}

	const_iterator end() const
	{
		return _relationsMap.end();
	}

	iterator end()
	{
		return _relationsMap.end();
	}

	std::ostream& print(int level = 0, std::ostream& os = std::cout) const;

private:
	//this way round?
	std::multimap<std::string, pxl::Id> _relationsMap;
};

} //namespace pxl

#endif /*SOFTRELATIONS_HH_*/
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_ptl_Layout_hh
#define pxl_ptl_Layout_hh


namespace pxl
{

// ptl

/** 
 This class holds layout information of PXL objects when visualized
 in the Graphical User Interface VisualPxl. For internal use only, 
 methods and data members are self-explanatory.  
 */
class Layout
{
public:
	Layout() :
		_type(0), _style(0), _color(0), _a(0.), _b(0.), _c(0.), _d(0.)
	{
	}
	~Layout()
	{
	}

	inline int getType()
	{
		return _type;
	}
	inline int getStyle()
	{
		return _style;
	}
	inline int getColor()
	{
		return _color;
	}

	inline double getA()
	{
		return _a;
	}
	inline double getB()
	{
		return _b;
	}
	inline double getC()
	{
		return _c;
	}
	inline double getD()
	{
		return _d;
	}

	inline void setType(int v)
	{
		_type = v;
	}
	inline void setStyle(int v)
	{
		_style = v;
	}
	inline void setColor(int v)
	{
		_color = v;
	}

	inline void setA(double v)
	{
		_a = v;
	}
	inline void setB(double v)
	{
		_b = v;
	}
	inline void setC(double v)
	{
		_c = v;
	}
	inline void setD(double v)
	{
		_d = v;
	}

	void serialize(const OutputStream &out) const;
	
	void deserialize(const InputStream &in);

private:
	int _type;
	int _style;
	int _color;

	double _a;
	double _b;
	double _c;
	double _d;
};

} // namespace pxl

#endif // pxl_ptl_Layout_hh

namespace pxl
{

// ptl

class ObjectOwner;

/** 
 This base class provides common functionalities for all derived PXL objects,  
 such as mother-daughter relations, weak pointer concept and related service routines. 
 In addition, it aggregates the PXL unique object-id and layout information for visualization in VisualPxl.
 It has a C++ pointer to the pxl::ObjectOwner it is aggregated in in order to avoid
 mother/daughter relations to outside objects to be established.
 */
class Relative : public Serializable
{
public:

	virtual ~Relative()
	{
		// remove all relations:
		unlinkMothers();
		unlinkDaughters();
		unlinkFlat();
		if (_refWkPtrSpec)
			_refWkPtrSpec->notifyDeleted();
		delete _ptrLayout;
		if (_refObjectOwner)
			std::cerr
					<< "Error in ~Relative(): Relative derivative must be deleted by ObjectOwner to guarantee safe object handling!"
					<< std::endl; //don't throw exception in destructor:
	}

	/// This method returns the PXL unique object-id.
	inline pxl::Id id() const
	{
		return getId();
	}

	static const Id& getStaticTypeId()
	{
		static const Id id("5dee644d-906f-4d8e-aecc-d9a644293260");
		return id;
	}

	virtual const pxl::Id& getTypeId() const
	{
		return getStaticTypeId();
	}

	virtual void serialize(const OutputStream &out) const
	{
		Serializable::serialize(out);

		_softRelations.serialize(out);

		out.writeString(_name);

		if (_ptrLayout)
		{
			out.writeBool(true);
			_ptrLayout->serialize(out);
		}
		else
		{
			out.writeBool(false);
		}
	}

	virtual void deserialize(const InputStream &in)
	{
		Serializable::deserialize(in);

		_softRelations.deserialize(in);

		in.readString(_name);

		bool hasLayout;
		in.readBool(hasLayout);
		if (hasLayout)
		{
			layout()->deserialize(in);
		}
		else
		{
			//necessary?
			delete _ptrLayout;
			_ptrLayout = 0;
		}
	}

	/// This method returns a C++ pointer to the pxl::ObjectOwner it is owned by.  
	inline pxl::ObjectOwner* owner() const
	{
		return _refObjectOwner;
	}

	/// This virtual method creates a deep copy and returns a C++ pointer to the newly created object.  
	virtual pxl::Relative* clone() const
	{
		return new pxl::Relative(*this);
	}

	/// This method grants access to the pxl::Relations instance managing mother relations  
	inline const pxl::Relations& getMotherRelations() const
	{
		return _motherRelations;
	}
	/// This method grants access to the pxl::Relations instance managing daughter relations  
	inline const pxl::Relations& getDaughterRelations() const
	{
		return _daughterRelations;
	}

	/// This method grants access to the pxl::Relations instance managing flat relations  
	inline const pxl::Relations& getFlatRelations() const
	{
		return _flatRelations;
	}

	/// This method establishes a mother relation to the \p target object; please notice, that
	/// only relations between objects owned by the same object owner will be established.   
	void linkMother(pxl::Relative* target) throw(std::runtime_error);
	/// This method establishes a daughter relation to the \p target object; please notice, that
	/// only relations between objects owned by the same object owner will be established.   
	void linkDaughter(pxl::Relative* target) throw(std::runtime_error);
	/// This method establishes a flat relation to the \p target object; please notice, that
	/// only relations between objects owned by the same object owner will be established.
	void linkFlat(pxl::Relative* target) throw(std::runtime_error);

	/// This method removes an existing daughter relation to the \p target object.
	void unlinkMother(pxl::Relative* target);
	/// This method removes an existing daughter relation to the \p target object.
	void unlinkDaughter(pxl::Relative* target);
	/// This method removes an existing daughter relation to the \p target object.
	void unlinkFlat(pxl::Relative* target);

	/// This method removes all existing mother relations.
	void unlinkMothers();
	/// This method removes all existing daughter relations.
	void unlinkDaughters();
	/// This method removes all existing flat relations.
	void unlinkFlat();

	void linkSoft(pxl::Relative* relative, const std::string& type)
	{
		if (relative)
		{
			_softRelations.set(relative, type);
			relative->_softRelations.set(this, type);
		}
	}

	void unlinkSoft(pxl::Relative* relative, const std::string& type)
	{
		if (relative)
		{
			_softRelations.remove(relative, type);
			relative->_softRelations.remove(this, type);
		}
	}

	const pxl::SoftRelations& getSoftRelations() const
	{
		return _softRelations;
	}

	pxl::SoftRelations& setSoftRelations()
	{
		return _softRelations;
	}

	/// This method returns the name.
	inline const std::string& getName() const
	{
		return _name;
	}

	/// This method sets the name to the contents of \p v.
	inline void setName(const std::string& v)
	{
		_name = v;
	}

	/// This method grants access to the layout information provided for visualization in VisualPxl.
	inline pxl::Layout* layout()
	{
		if (!_ptrLayout)
			_ptrLayout = new pxl::Layout;
		return _ptrLayout;
	}

	/// This method recursively invokes its own and the print() methods of all daughter objects.
	/// @param level verbosity level
	/// @param os output _stream, default is std::cout
	/// @param pan print indention
	/// @return output _stream
	std::ostream& printDecayTree(int level = 0, std::ostream& os = std::cout,
			int pan = 1) const;

	/// This virtual method is intended to print out object state information on various verbosity levels.
	/// @param level verbosity level
	/// @param os output _stream, default is std::cout
	/// @param pan print indention
	/// @return output _stream
	virtual std::ostream& print(int level = 1, std::ostream& os = std::cout,
			int pan = 0) const;

	/// This method creates a weak pointer to itself and returns a pxl::WkPtrBase* to the newly created weak pointer instance. 
	virtual pxl::WkPtrBase* createSelfWkPtr()
	{
		throw std::runtime_error("pxl::ObjectBase::createSelfWkPtr(): ATTENTION! Inheriting class must reimplement this virtual method.");
		return 0;
	}
	
	inline const pxl::Relative& operator=(const pxl::Relative& pa)
	{
		if (pa._ptrLayout)
			_ptrLayout = new pxl::Layout;
		_name = pa._name;
		
		return *this;
	}
	
protected:
	Relative() :
		pxl::Serializable(), _refWkPtrSpec(0), _refObjectOwner(0),
				_motherRelations(), _daughterRelations(), _name("default"),
				_ptrLayout(0)
	{
	}

	explicit Relative(const Relative& original) :
		pxl::Serializable(), _refWkPtrSpec(0), _refObjectOwner(0),
				_motherRelations(), _daughterRelations(),
				_name(original._name), _ptrLayout(0)
	{
		this->init(original);
	}

	explicit Relative(const Relative* original) :
		pxl::Serializable(), _refWkPtrSpec(0), _refObjectOwner(0),
				_motherRelations(), _daughterRelations(),
				_name(original->_name), _ptrLayout(0)
	{
		this->init(*original);
	}

	inline void init(const Relative& original)
	{
		if (original._ptrLayout)
			_ptrLayout = new pxl::Layout;
	}

	std::ostream& printPan1st(std::ostream& os, int pan) const;
	std::ostream& printPan(std::ostream& os, int pan) const;

private:
	pxl::WkPtrBase* _refWkPtrSpec;
	pxl::ObjectOwner* _refObjectOwner;

	pxl::Relations _motherRelations;
	pxl::Relations _daughterRelations;
	pxl::Relations _flatRelations;

	pxl::SoftRelations _softRelations;

	std::string _name;

	pxl::Layout* _ptrLayout;

	friend class WkPtrBase;
	friend class ObjectOwner;
};

} // namespace pxl

// operators
std::ostream& operator<<(std::ostream& cxxx, const pxl::Relative& obj);

#endif // pxl_ptl_ObjectBase_hh

namespace pxl {

// ptl

/// This typedef is intended to hold the origin information of copied objects in pxl::ObjectOwner instances. 
typedef std::map<pxl::Id, pxl::Relative*> CopyHistory;

} // namespace pxl

#endif // pxl_ptl_CopyHistory_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_ptl_Index_hh
#define pxl_ptl_Index_hh


namespace pxl {

// ptl

/// This typedef is intended to hold the index information (string-object associations) in pxl::ObjectOwner instances. 
typedef std::map<std::string, pxl::Relative*> Index;


} // namespace pxl

#endif // pxl_ptl_Index_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_ptl_WkPtrSpec_hh
#define pxl_ptl_WkPtrSpec_hh



namespace pxl
{

// ptl

/** 
 This class template represents a weak pointer to PXL objects of \p objecttype that aggregate data of \p datatype.   
 */
template<class objecttype> class weak_ptr : public pxl::WkPtrBase
{
public:
	weak_ptr() :
		pxl::WkPtrBase()
	{
	}
	weak_ptr(objecttype* ptr) :
		pxl::WkPtrBase()
	{
		pxl::WkPtrBase::connect(ptr);
	}
	weak_ptr(objecttype& object) :
		pxl::WkPtrBase()
	{
		pxl::WkPtrBase::connect(&object);
	}
	weak_ptr(const pxl::weak_ptr<objecttype>& original) :
		pxl::WkPtrBase()
	{
		pxl::WkPtrBase::connect((objecttype*) original._objectRef);
	}
	weak_ptr(const pxl::weak_ptr<objecttype>* original) :
		pxl::WkPtrBase()
	{
		pxl::WkPtrBase::connect((objecttype*) original->_objectRef);
	}

	virtual ~weak_ptr()
	{
		pxl::WkPtrBase::connect(0);
	}

	/// This virtual method creates a deep copy and returns a C++ pointer to the newly created weak pointer instance.  
	virtual pxl::WkPtrBase* clone() const
	{
		return new pxl::weak_ptr<objecttype>(*this);
	}

	/// This assignment operator causes the weak pointer to reference the object referenced by \p pptr.
	inline void operator=(const pxl::weak_ptr<objecttype>& pptr)
	{
		connect(pptr._objectRef);
	}
	/// This assignment operator causes the weak pointer to reference the object.
	inline void operator=(objecttype& object)
	{
		connect(&object);
	}
	/// This assignment operator causes the weak pointer to reference the object pointed to by \p objectptr.
	inline void operator=(objecttype* objectptr)
	{
		connect(objectptr);
	}

	// methods to grant object & data access
	/// This method provides direct access to the referenced object.
	inline objecttype& object() const
	{
		return *access();
	}

	/// This arrow operator de-references the weak pointer.   
	inline objecttype* operator->() const
	{
		return access();
	}

	/// This arrow operator de-references the weak pointer.   
	inline objecttype* ptr() const
	{
		return dynamic_cast<objecttype*>(_objectRef);
	}

	inline operator objecttype*() const
	{
		return dynamic_cast<objecttype*>(_objectRef);
	}

	/// This method attempts a dynamic cast on the referenced object
	static weak_ptr<objecttype>* cast_dynamic(WkPtrBase* orig) throw (std::runtime_error)
	{
		objecttype *object = dynamic_cast<objecttype*>(orig->pointer());
		if (!object)
			return 0;
		// FIXME: This is crude but required:
		if (PXL_UNLIKELY(reinterpret_cast<void*>(object)
				!= reinterpret_cast<void*>(orig->pointer())))
			throw std::runtime_error("pxl::WkPtrSpec::cast_dynamic(): Unsupported multiple inheritance configuration.");
		return reinterpret_cast<pxl::weak_ptr<objecttype>*>(orig);
	}

	// safe access to object
	inline objecttype* access() const throw (std::runtime_error)
	{
		if (_objectRef)
			return (objecttype*)_objectRef;
		throw std::runtime_error("pxl::WkPtrSpec::access(): FATAL: The object you intend to access does not exist!");
		return 0;
	}

};

template<class objecttype> objecttype& operator*(pxl::weak_ptr<objecttype>& wkPtr)
{
	return wkPtr.object();
}

template<class objecttype> const objecttype& operator*(
		const pxl::weak_ptr<objecttype>& wkPtr)
{
	return wkPtr.object();
}

} // namespace pxl


#endif // pxl_ptl_WkPtrSpec_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_ptl_ObjectOwner_hh
#define pxl_ptl_ObjectOwner_hh




namespace pxl
{

// ptl

/** 
 This class is an active container for pxl::Relative derivatives, such as 
 pxl::Vertex, pxl::Particle or pxl::EventView; 
 it has the ownership and deletion responsability for the contained objects. 
 The method template create() can be used to create derivatives 
 of pxl::Relative within object owners. The method set() can be 
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
 can be used to directly access objects by their index-ids or object-ids.
 The pxl::ObjectOwner extends the functionality of the contained STL vector. It provides a selective iterator, the class template 
 pxl::ObjectOwner::TypeIterator, that ignores all objects other than the 
 specialized data type. 
 */
class ObjectOwner
{
public:
	ObjectOwner() :
		_container(), _copyHistory(), _index(), _uuidSearchMap()
	{
	}
	/// This copy constructor performs a deep copy of object 
	/// owner \p original and all contained objects with (redirected) relations.
	/// A copy history keeps track of originals and copies
	/// and the findCopyOf() method allows quick access to the copies. 
	ObjectOwner(const pxl::ObjectOwner& original) :
		_container(), _copyHistory(), _index(), _uuidSearchMap()
	{
		this->init(original);
	}

	/// This copy constructor performs a deep copy of object 
	/// owner \p original and all contained objects with (redirected) relations.
	/// A copy history keeps track of originals and copies
	/// and the findCopyOf() method allows quick access to the copies. 
	ObjectOwner(const pxl::ObjectOwner* original) :
		_container(), _copyHistory(), _index(), _uuidSearchMap()
	{
		this->init(*original);
	}
	/// This destructor deletes all contained objects.
	virtual ~ObjectOwner()
	{
		pxl::ObjectOwner::clearContainer();
	}

	virtual void serialize(const OutputStream &out) const;

	virtual void deserialize(const InputStream &in);

	/// This method template creates a new instance of \p objecttype;
	/// objecttype must be a class inheriting from pxl::Relative;
	/// the newly created instance is owned and will be deleted by this object owner. 
	template<class objecttype> objecttype* create()
	{
		objecttype* pitem = new objecttype;
		pitem->_refObjectOwner = this;
		_container.push_back(static_cast<pxl::Relative*>(pitem));
		_uuidSearchMap.insert(std::pair<pxl::Id, pxl::Relative*>(pitem->getId(), pitem));
		return pitem;
	}

	/// This method template creates a copy of \p original by invoking the copy constructor of \p objecttype; 
	/// \p objecttype must be a class inheriting from pxl::Relative;
	/// the newly created instance is owned and will be deleted by this object owner. 
	template<class objecttype> objecttype* create(
	const objecttype* original)
	{
		objecttype* pitem = new objecttype(*original);
		pitem->_refObjectOwner = this;
		_container.push_back(static_cast<pxl::Relative*>(pitem));
		_uuidSearchMap.insert(std::pair<pxl::Id, pxl::Relative*>(pitem->getId(), pitem));
		return pitem;
	}

	/// This method template creates a new \p objecttype instance by invoking a \p ctrtype overloaded constructor; 
	/// \p objecttype must be a class inheriting from pxl::Relative;
	/// the newly created instance is owned and will be deleted by this object owner. 
	template<class objecttype, class ctrtype> objecttype* create(
	const ctrtype& original)
	{
		objecttype* pitem = new objecttype(*original);
		pitem->_refObjectOwner = this;
		_container.push_back(static_cast<pxl::Relative*>(pitem));
		_uuidSearchMap.insert(std::pair<pxl::Id, pxl::Relative*>(pitem->getId(), pitem));
		return pitem;
	}

	/// This method inserts \p item in the container of this object owner and takes deletion responsability.
	void set(pxl::Relative* item);
	/// This method deletes \p item.
	void remove(pxl::Relative* item);
	/// This method returns true if \p item is owned by this object owner.
	bool has(const pxl::Relative* item) const;

	/// This method clears the object owner and deletes all owned objects. 
	void clearContainer();

	/// This method searches the index for the index-id \p idx and returns a dynamically casted 
	/// C++ pointer of type \p objecttype* to the corresponding object; 
	/// in case idx is not found a null pointer is returned.
	template<class objecttype> inline objecttype* findObject(
	const std::string& idx) const // goes via Index & casts

	{
		pxl::Index::const_iterator it = _index.find(idx);
		if (it!=_index.end())
		return dynamic_cast<objecttype*>(it->second);
		return 0;
	}

	inline pxl::Relative* findObject(const std::string& idx) const
	{
		pxl::Index::const_iterator it = _index.find(idx);
		if (it!=_index.end())
		return it->second;
		return 0;
	}

	pxl::Relative* getById(pxl::Id id) const
	{
		std::map<pxl::Id, pxl::Relative*>::const_iterator found = _uuidSearchMap.find(id);
		if ( found != _uuidSearchMap.end() )
		return found->second;
		return 0;
	}

	/// This method searches the copy history to locate the copy of \p original and 
	/// returns a dynamically casted C++ pointer of type \p objecttype* to the corresponding copy; 
	/// in case no copy can be traced a null pointer is returned.
	template<class objecttype> objecttype* findCopyOf(
	const pxl::Relative* original) const // goes via CopyHistory & casts

	{
		pxl::CopyHistory::const_iterator it = _copyHistory.find(original->id());
		if (it!=_copyHistory.end())
		return dynamic_cast<objecttype*>(it->second);
		return 0;
	}

	/// This method provides direct access to the copy history (created by the copy constructor). 
	inline const pxl::CopyHistory& getCopyHistory() const
	{
		return _copyHistory;
	}
	/// This method clears the copy history  (created by the copy constructor). 
	inline void clearCopyHistory()
	{
		_copyHistory.clear();
	}

	/// This method registers the object \p obj with the index-id \p idx in the index and returns true in case of success; 
	/// in case the index string is present, by default an error message is given. A bool can be set which allows overwriting.
	/// Please notice that \p obj must be owned by this object owner and \p idx must not be a zero length string.  
	bool setIndex(const std::string& idx, pxl::Relative* obj, bool overwrite = false);

	/// This method provides direct read access to the index. 
	inline const pxl::Index& getIndex() const
	{
		return _index;
	}
	/// This method removes the index entry with index-id \p idx; please notice: it does not remove the object itself. 
	inline void removeIndex(const std::string& idx)
	{
		_index.erase(idx);
	}
	/// This method clears the index; please notice: it does not remove the objects themselves.
	inline void clearIndex()
	{
		_index.clear();
	}

	/// This method allows read access to the contained STL vector of Relative pointers to, e.g., use STL algorithms.
	const std::vector<pxl::Relative*>& getObjects() const
	{
		return _container;
	}

	/// Typedef for standard const_iterator.
	typedef std::vector<pxl::Relative*>::const_iterator const_iterator;
	typedef std::vector<pxl::Relative*>::iterator iterator;

	/// This returns the const iterator to the first element of the contained vector.  
	inline const_iterator begin() const
	{
		return _container.begin();
	}

	/// This returns the iterator to the first element of the contained vector.  
	inline iterator begin()
	{
		return _container.begin();
	}

	/// This returns the const iterator to the end of the contained vector.
	inline const_iterator end() const
	{
		return _container.end();
	}

	/// This returns the iterator to the end of the contained vector.
	inline iterator end()
	{
		return _container.end();
	}

	/// Returns the number of elements the ObjectOwner holds.
	inline unsigned int size() const
	{
		return _container.size();
	}

	/// Fills into the passed vector weak pointers to the objects of the type specified by the template argument.
	template<class objecttype> int getObjectsOfType(std::vector<objecttype*>& vec) const
	{
		int size = vec.size();
		for (pxl::ObjectOwner::const_iterator iter = begin(); iter!=end(); ++iter)
		{
			objecttype* obj = dynamic_cast<objecttype*>(*iter);
			if (obj!=0)
			vec.push_back(obj);
		}
		return vec.size()-size;
	}

	// - - - - - - - - - - - - - - - - - - - - - - - - - - -
	/// For STL-style iteration on selective class: iterator class template; 
	/// this iterator behaves like a normal STL iterator but ignores all objects 
	/// that cannot be interpreted as type objecttype (tested by dynamic casts).
	/// Use in STL-style, except that either begin gets the objecttype class as template argument,
	/// or a constructor of the TypeIterator is used.
	template<class objecttype> class TypeIterator
	{
	public:
		/// Copy constructor.	
		TypeIterator(const TypeIterator& other) :
		_containerRef(other._containerRef), _iter(other._iter)
		{
		}

		/// Constructor from ObjectOwner instance.
		TypeIterator(const pxl::ObjectOwner* container) :
		_containerRef(container), _iter(container->begin())
		{
			if ( dynamic_cast<objecttype*>(*_iter)==0)
			(*this)++;
		}

		const TypeIterator operator++(int)
		{
			TypeIterator orig = *this;
			if (_iter!=_containerRef->end())
			do
			_iter++;
			while (_iter!=_containerRef->end()
			&& dynamic_cast<objecttype*>(*_iter)==0);
			return orig;
		}

		const TypeIterator& operator++()
		{
			if (_iter!=_containerRef->end())
			do
			_iter++;
			while (_iter!=_containerRef->end()
			&& dynamic_cast<objecttype*>(*_iter)==0);
			return *this;
		}

		inline objecttype* operator*()
		{
			return _iter==_containerRef->end() ? 0
			: dynamic_cast<objecttype*>(*_iter);
		}

		inline bool operator==(pxl::ObjectOwner::const_iterator iter)
		{
			return (_iter==iter);
		}

		inline bool operator!=(pxl::ObjectOwner::const_iterator iter)
		{
			return (_iter!=iter);
		}

	private:
		const pxl::ObjectOwner* _containerRef;
		pxl::ObjectOwner::const_iterator _iter;
	};

	/// This templated method provides an STL-style begin()-method to
	/// initialise the TypeIterator.
	template<class objecttype> const TypeIterator<objecttype> begin() const
	{
		TypeIterator<objecttype> it(this);
		return it;
	}

private:
	void init(const pxl::ObjectOwner& original);

	std::vector<pxl::Relative*> _container;
	pxl::CopyHistory _copyHistory;
	pxl::Index _index;
	std::map<pxl::Id, pxl::Relative*> _uuidSearchMap;

};

} // namespace pxl

#endif // pxl_ptl_ObjectOwner_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_ptl_Filter_hh
#define pxl_ptl_Filter_hh

#include <algorithm>


namespace pxl
{

template<class comparetype> class Comparator
{
public:
	virtual bool operator()(const comparetype&, const comparetype&) = 0;
	virtual bool operator()(const comparetype*, const comparetype*) = 0;	
	virtual ~Comparator()
	{
	}
};

template<class objecttype> class FilterCriterion
{
public:
	virtual bool operator()(const objecttype&) const = 0;
	virtual ~FilterCriterion()
	{
	}
};

/** 
 This class template provides a sorted filter for PXL physics objects;
 it handles objects of type \p objecttype sorted by \p sorttype. The user can write
 his own filters by public inheritance from this class and reimplementation
 of the methods pass() and sort().  
 */
template<class objecttype, class compare> class Filter
{
public:
	virtual ~Filter()
	{
	}

	/// This method applies the filter by running over the \p objects container and fills
	/// the passed vector with pointers to the objects passing the filter criteria.
	virtual int apply(const pxl::ObjectOwner* objects,
			std::vector<objecttype*>& fillVector,
			const FilterCriterion<objecttype>& criterion)
	{
		int size = fillVector.size();

		// fill map:
		for (pxl::ObjectOwner::TypeIterator<objecttype> iter(objects); iter
				!=objects->end(); ++iter)
		{
			if (!criterion(**iter))
				continue;
			fillVector.push_back(*iter);
		}
		//<typename std::vector<objecttype*>::iterator, compare>
		compare comp;
		std::sort(fillVector.begin(), fillVector.end(), comp);
		return fillVector.size()-size;
	}

};

} // namespace pxl

#endif // pxl_ptl_Filter_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_pol_UserRecord_hh
#define pxl_pol_UserRecord_hh




namespace pxl
{
/**
 This class is intented to aggregate information complementary to data members in form
 of string-variant pairs. 
 All PXL physics objects own user records and provide methods for quick 
 access to individual user record entries.
 */
class UserRecord
{

private:
	class DataSocket
	{
public:
		DataSocket() :
			_references(1)
		{
		}
		DataSocket(const DataSocket& original) :
			_references(1), _data(original._data)
		{
		}
		DataSocket(const DataSocket* original) :
			_references(1), _data(original->_data)
		{
		}
		virtual ~DataSocket()
		{
		}

		// for deep copies
		virtual DataSocket* clone() const
		{
			return new DataSocket(this);
		}

		// methods to grant data access
		inline std::map<std::string, pxl::Variant>* getData()
		{
			return &_data;
		}
		inline void setData(const std::map<std::string, pxl::Variant>* object)
		{
			_data = *object;
		}

		unsigned int _references;
		std::map<std::string, pxl::Variant> _data;

	}; //class Datasocket

public:
	typedef std::map<std::string, pxl::Variant> PtlMap;
	typedef std::pair<std::string, pxl::Variant> StlPair;
	typedef std::map<std::string, pxl::Variant>::const_iterator const_iterator;
	typedef std::map<std::string, pxl::Variant>::iterator iterator;

	UserRecord()
	{
		_dataSocket = new DataSocket;
	}
	UserRecord(const UserRecord& original)
	{
		_dataSocket = original._dataSocket;
		_dataSocket->_references++;
	}
	UserRecord(const UserRecord* original)
	{
		_dataSocket = original->_dataSocket;
		_dataSocket->_references++;
	}
	~UserRecord()
	{
		dropDataSocket();
	}

	void serialize(const OutputStream &out) const;
	void deserialize(const InputStream &in);

	/// This assignment operator acts directly on the aggregated data.
	inline UserRecord& operator=(const UserRecord& original)
	{
		dropDataSocket();
		_dataSocket = original._dataSocket;
		_dataSocket->_references++;
		return *this;
	}

	/// This method template inserts (or replaces) the user record indetified by \p key.
	template<typename datatype> void set(const std::string& key,
			const datatype& item)
	{
		pxl::Variant& value = findOrAlloc(key);
		if (PXL_UNLIKELY(value.getType() != Variant::TYPE_NULL))
			value.clear();

		value.template init<datatype>();
		value.template set<datatype>(item);
	}

	/// This method template searches and returns the user record item identified by \p key; \p defaultitem is returned in case the key is not found. 
	template<typename datatype> datatype find(const std::string& key,
			const datatype& defaultitem) const
	{
		const pxl::Variant* value = findOrReturn(key);
		if (!value)
			return defaultitem;
		return value->template get<datatype>();
	}

	/// This method template searches and returns the user record item indetified by \p key; a pxl::Exception is thrown in case the key is not found. 
	template<typename datatype> datatype find(const std::string& key) const
			throw (std::runtime_error)
	{
		const pxl::Variant* value = findOrReturn(key);
		if (!value)
			throw std::runtime_error("pxl::UserRecord::find(...): key not found and no default item provided");
		return value->template get<datatype>();
	}

	/// This method templates checks if the user record entry identified by key is present.
	bool check(const std::string& key) const
	{
		const pxl::Variant* value = findOrReturn(key);
		if (!value)
			return false;
		return true;
	}

	/// This method template checks if user record entry identified by \p key is present.
	/// If yes, its value is put into the passed \p item. 
	template<typename datatype> bool check(const std::string& key,
			datatype& item) const
	{
		const pxl::Variant* value = findOrReturn(key);
		if (!value)
			return false;
		item = value->template get<datatype>();
		return true;
	}

	inline void clear()
	{
		setContainer()->clear();
	}

	inline void erase(const std::string& key)
	{
		setContainer()->erase(key);
	}

	/// This method grants read access to the aggregated data.
	inline const std::map<std::string, pxl::Variant>* getContainer() const
	{
		return _dataSocket->getData();
	}

	inline const_iterator begin() const
	{
		return getContainer()->begin();
	}

	inline const_iterator end() const
	{
		return getContainer()->end();
	}

	inline int size() const
	{
		return getContainer()->size();
	}

	std::ostream
			& print(int level = 0, std::ostream& os = std::cout, int pan = 0) const;

private:
	DataSocket* _dataSocket;

	/// This method grants write access to the aggregated data; 
	/// if necessary, the copy-on-write mechanism performs a deep copy of the aggregated data first. 
	inline std::map<std::string, pxl::Variant>* setContainer()
	{
		if (_dataSocket->_references > 1)
		{
			_dataSocket->_references--;
			_dataSocket = new DataSocket(*_dataSocket);
		}
		return _dataSocket->getData();
	}

	inline void dropDataSocket()
	{
		if (_dataSocket->_references-- == 1)
			delete _dataSocket;
	}

	pxl::Variant& findOrAlloc(const std::string &key)
	{
		iterator insertPos = setContainer()->lower_bound(key);
		if (insertPos == getContainer()->end() || insertPos->first != key)
			return setContainer()->insert(insertPos, StlPair(key, pxl::Variant()))->second;
		else
			return insertPos->second;
	}

	const pxl::Variant* findOrReturn(const std::string &key) const
	{
		const_iterator found = getContainer()->find(key);
		if (found == getContainer()->end())
			return 0;
		return &found->second;
	}

};

} // namespace pxl

#endif // pxl_pol_UserRecord_hh

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_io_hh
#define pxl_io_hh

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_io_ObjectFactory_hh
#define pxl_io_ObjectFactory_hh



// iotl
/**
 This class serves the PXL I/O scheme by managing
 the relation of classes to UUIDs.
 */

namespace pxl
{

class ObjectFactory
{
public:
	
	class ProducerInterface
	{
	public:
		virtual Serializable *create () const = 0;
	};
	
	template <class T>
	class ProducerTemplate : public ProducerInterface
	{
	public:

		ProducerTemplate (const Id& id)
		{
			ObjectFactory::registerClass (id, this);
		}
		
		virtual Serializable *create () const
		{
			return new T();
		}
	};

private:
	
	ObjectFactory()
	{
	}
	
	std::map<Id,const ProducerInterface *> _Producers;

public:
	static ObjectFactory& instance()
	{
		static ObjectFactory f;
		return f;
	}
		
	static Serializable *create (const Id& id)
	{
		std::map<Id,const ProducerInterface *>::iterator result;
		result = instance()._Producers.find (id);
		if (result == instance()._Producers.end())
			return 0;
		else
			return (*result).second->create();
	}
	
	static void registerClass (const Id& id, const ProducerInterface* producer)
	{
		instance()._Producers[id] = producer;
	}
};

}

#endif //pxl_io_ObjectFactory_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_io_ChunkWriter_hh
#define pxl_io_ChunkWriter_hh

#include <fstream>
#include <zlib.h>


namespace pxl
{

// io
/**
 This class implemenents methods for writing to PXL I/O files.
 PXL I/O allows the storage of complete physics events and information chunks.
 Each event or information chunk makes up a section in the output file. 
 Each section consists of a header, and a number of blocks which can be compressed individually.
 The compression is incorporated via zlib.
 The entry point for the standard user is the class OutputFile. 
 */
class ChunkWriter
{
public:
	ChunkWriter(std::ostream& stream) :
		_stream(stream), _nBytes(0) 
	{
	}

	~ChunkWriter()
	{
	}

	/// This method writes the current block to the output file stream.
	bool write()
	{
		return write("");
	}

	/// Begin a new event, optionally pass information string.
	bool newEvent(const std::string& info = "")
	{
		return newFileSection(info, 'E');
	}

	/// Begin a new information chunk, optionally pass information string.
	bool newInformationChunk(const std::string& info = "")
	{
		return newFileSection(info, 'I');
	}

	/// Writes a new file section, indicating the content by the section marker char, and the passed information string.
	bool newFileSection(const std::string& info, char cSectionMarker);

	/// Writes a new block marker.
	inline bool newBlock()
	{
		// begin block marker:
		return writeFlag('B');
	}

	/// Writes an end-of-event marker and the number of bytes stored in the event.
	bool endEvent();

	/// End information chunk.
	inline bool endInformationChunk()
	{
		return endEvent();
	}

	/// This method writes the current stream. An information string, and a compression mode char (allowed values between '0' and '9') are passed.
	bool write(std::string info, char compressionMode = '6') throw(std::runtime_error);

	const OutputStream& getOutputStream()
	{
		return _buffer;
	}
	
protected:
	/// Write char flag.
	inline bool writeFlag(char cEvtMarker)
	{
		_stream.write(&cEvtMarker, 1);
		_nBytes+=1;
		return true;
	}
	
private:
	std::ostream& _stream;
	BufferOutput _buffer;
	pxl::int32_t _nBytes;
};
}
#endif /*pxl_io_ChunkWriter_hh*/
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_io_ChunkReader_hh
#define pxl_io_ChunkReader_hh



#define iotl__iStreamer__lengthUnzipBuffer 65536

namespace pxl
{

namespace skipMode
{

enum skipMode
{
	off = 0,
	on
};

}
// iotl
/**
 This class implemenents various methods for reading from PXL I/O and
 bases on the event structure defined in the ChunkWriter class. 
 The entry point for the user is the class InputFile.
 */

class ChunkReader
{
public:
	typedef skipMode::skipMode skipMode;

	enum statusFlag
	{
		preHeader = 0,
		evPreBlock,
		infoPreBlock,
		preBlock
	};

	enum readMode
	{
		all = 0,
		event,
		infoChunk
	};

	enum infoMode
	{
		ignore = 0,
		evaluate
	};

	ChunkReader(std::istream& stream) :
		_stream(stream), _status(preHeader)
	{
	}

	~ChunkReader()
	{
	}

	/// Reads in the next event header. 
	bool readHeader(readMode mode, skipMode skip, infoMode checkInfo,
			const std::string& infoCondition);

	/// Reads in the next block.
	bool readBlock(skipMode skip, infoMode checkInfo,
			const std::string& infoCondition) throw(std::runtime_error);

	/// This methods skips an event/information chunk.
	bool skip();

	/// This method goes back one event or information chunk.
	bool previous();

	/// Reads in the header of the next event. False is returned if not successful.
	bool next(skipMode skip = skipMode::on, infoMode checkInfo = ignore,
			const std::string& infoCondition = "")
	{
		return readHeader(all, skip, checkInfo, infoCondition);
	}

	/// Reads in the header of the next event. False is returned if not successful.
	bool nextEvent(skipMode skip = skipMode::on, infoMode checkInfo = ignore,
			const std::string& infoCondition = "")
	{
		return readHeader(event, skip, checkInfo, infoCondition);
	}

	/// Reads in the header of the next event. False is returned if not successful.
	bool nextInformationChunk(skipMode skip = skipMode::on,
			infoMode checkInfo = ignore, const std::string& infoCondition = "")
	{
		return readHeader(infoChunk, skip, checkInfo, infoCondition);
	}

	/// Reads the next block and puts data into the input stream. False is returned if not successful.
	bool nextBlock(skipMode skip = skipMode::on, infoMode checkInfo = ignore,
			const std::string& infoCondition = "")
	{
		return readBlock(skip, checkInfo, infoCondition);
	}

	/// Access to the data read in the individual blocks.
	inline const InputStream& getInputStream()
	{
		return _buffer;
	}

	bool isInformationChunk()
	{
		return ( (_status == preHeader && _stream.peek()=='I') || _status
				== infoPreBlock );
	}

	bool isEvent()
	{
		return ( (_status == preHeader && _stream.peek()=='E') || _status
				== evPreBlock );
	}

	bool isBlock()
	{
		return (_stream.peek()=='B');
	}

	bool isEnd()
	{
		return (_stream.peek()=='e');
	}

	/// Method used internally to get the status, indicating the position in the I/O file.
	inline const statusFlag getStatus() const
	{
		return _status;
	}

	void endEvent()
	{
		if (_status!=preHeader)
			while (nextBlock())
				;
	}

protected:
	/// Helper method to perform the unzipping.
	int unzipEventData(int nBytes) throw(std::runtime_error);

	/// Reads in a char from file and returns this.
	inline char nextBlockId()
	{
		char identifier;
		_stream.read(&identifier, 1);
		return identifier;
	}

private:
	std::istream& _stream;
	BufferInput _buffer;
	/// Status flag. 0 at end of event, 1 at end of block.
	statusFlag _status;

};

} //namespace pxl

#endif /*pxl_iotl_ChunkReader_hh*/
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_io_InputHandler_hh
#define pxl_io_InputHandler_hh


//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_pol_Event_hh
#define pxl_pol_Event_hh



namespace pxl
{

class Event : public pxl::Serializable
{
public:
	Event() :
		pxl::Serializable()
	{
	}

	Event(const Event& event) :
		_objects(event._objects), _userRecords(event._userRecords)
	{
	}

	virtual ~Event()
	{
	}

	virtual const pxl::Id& getTypeId() const
	{
		return getStaticTypeId();
	}

	static const Id& getStaticTypeId();

	virtual void serialize(const OutputStream &out) const
	{
		pxl::Serializable::serialize(out);
		_objects.serialize(out);
		_userRecords.serialize(out);
	}

	virtual void deserialize(const InputStream &in)
	{
		pxl::Serializable::deserialize(in);
		_objects.deserialize(in);
		_userRecords.deserialize(in);
	}

	/// This method template creates a new instance of \p objecttype;
	/// objecttype must be a class inheriting from pxl::Relative;
	/// the newly created instance is owned and will be deleted by the object owner. 
	template<class datatype> datatype* create()
	{
		return _objects.create<datatype>();
	}

	/// This method template creates a new \p objecttype instance by invoking a \p ctrtype overloaded constructor; 
	/// \p objecttype must be a class inheriting from pxl::Relative;
	/// the newly created instance is owned and will be deleted by the object owner. 
	template<class datatype, class ctrdatatype> datatype* create(
			const ctrdatatype& ori)
	{
		return _objects.create<datatype,ctrdatatype>(ori);
	}

	// crateIndexed
	/// This method template acts like create() and registers the newly created instance under \p idx in the index.
	template<class datatype> datatype* createIndexed(const std::string& idx)
	{
		datatype* obj = _objects.create<datatype>();
		setIndex(idx, obj);
		return obj;
	}

	/// This method template acts like create() and registers the newly created instance under \p idx in the index.
	template<class datatype, class ctrdatatype> datatype* createIndexed(
			const ctrdatatype& ori, const std::string& idx)
	{
		datatype* obj = _objects.create<datatype,ctrdatatype>(ori);
		setIndex(idx, obj);
		return obj;
	}

	/// This method inserts \p obj in the container of the object owner and takes deletion responsability.
	inline void setObject(pxl::Relative* obj)
	{
		_objects.set(obj);
	}

	/// This method inserts \p obj with the index-id \p idx in the container of the object owner and takes deletion responsability.
	inline void setObject(pxl::Relative* obj, const std::string& idx)
	{
		_objects.set(obj);
		setIndex(idx, obj);
	}

	/// This method registers the object \p obj with the index-id \p idx in the index and returns true in case of success;
	/// please notice, that obj must be owned by this object owner and \p idx must not be a zero length string.  
	inline bool setIndex(const std::string& idx, pxl::Relative* obj)
	{
		return _objects.setIndex(idx, obj);
	}

	
	inline pxl::ObjectOwner& getObjectOwner()
	{
		return _objects;
	}

	/// This method provides const access to the object owner.
	inline const pxl::ObjectOwner& getObjectOwner() const
	{
		return _objects;
	}

	template<class objecttype> inline void getObjectsOfType(
			std::vector<objecttype*>& vec) const
	{
		_objects.getObjectsOfType<objecttype>(vec);
	}
	
	template<class objecttype, class comparator> inline void getObjectsOfType(
			std::vector<objecttype*>& vec, const FilterCriterion<objecttype>& criterion) const
	{
	        pxl::Filter<objecttype, comparator> filter; 
		filter.apply(&_objects, vec, criterion);
		//_objects.getObjectsOfType<objecttype>(vec);
	}        

	inline const std::vector<pxl::Relative*>& getObjects() const
	{
		return _objects.getObjects();
	}
	
	/// This method deletes the object \p obj.
	inline void removeObject(pxl::Relative* obj)
	{
		_objects.remove(obj);
	}

	/// This method clears the object owner and deletes all owned objects. 
	inline void clearObjects()
	{
		_objects.clearContainer();
	}

	/// This method searches the index for the index-id \p idx and returns a dynamically casted 
	/// C++ pointer of type \p objecttype* to the corresponding object; 
	/// in case idx is not found a null pointer is returned.
	template<class objecttype> inline objecttype* findObject(
			const std::string idx) const
	{
		return _objects.findObject<objecttype>(idx);
	}

	/// This method searches the copy history to locate the copy of \p original and 
	/// returns a dynamically casted C++ pointer of type \p objecttype* to the corresponding copy; 
	/// in case no copy can be traced a null pointer is returned.
	template<class objecttype> inline objecttype* findCopyOf(
			const pxl::Relative* original) const
	{
		return _objects.findCopyOf<objecttype>(original);
	}

	/// This method provides direct access to the index. 
	inline const pxl::Index& getIndex() const
	{
		return _objects.getIndex();
	}

	/// This method removes the index entry with index-id \p idx; please notice: it does not remove the object itself. 
	inline void removeIndex(const std::string& idx)
	{
		_objects.removeIndex(idx);
	}

	/// This method clears the index; please notice: it does not remove the objects themself.
	inline void clearIndex()
	{
		_objects.clearIndex();
	}
	
	/// This method sets the user record entry identified by \p key to \p item.
	template<typename datatype> inline void setUserRecord(
			const std::string& key, const datatype& item)
	{
		_userRecords.template set<datatype>(key, item);
	}

	/// This method removes the user record entry identified by \p key.
	inline void removeUserRecord(const std::string& key)
	{
		_userRecords.erase(key);
	}

	inline void clearUserRecords()
	{
		_userRecords.clear();
	}

	/// This method searches the user record entry identified by \p key; \p defaultitem is returned in case key is not found.
	template<typename datatype> inline datatype findUserRecord(
			const std::string& key, const datatype& defaultitem) const
	{
		return _userRecords.template find<datatype>(key, defaultitem);
	}

	/// This method searches the user record entry identified by \p key; an exception is thrown in case key is not found.
	template<typename datatype> inline datatype findUserRecord(
			const std::string& key) const throw(std::runtime_error)
	{
		return _userRecords.template find<datatype>(key);
	}

	/// This method checks if the user record entry identified by \p key is present.
	template<typename datatype> inline bool checkUserRecord(
			const std::string& key) const
	{
		return _userRecords.template check<datatype>(key);
	}

	/// This method checks if the user record entry identified by \p key is present. If yes, \p item is set to the according value.
	template<typename datatype> inline bool checkUserRecord(
			const std::string& key, datatype& item) const
	{
		return _userRecords.template check<datatype>(key, item);
	}
	
	std::ostream& print(int level=1, std::ostream& os=std::cout, int pan=1) const;

private:
	pxl::ObjectOwner _objects;
	pxl::UserRecord _userRecords;
};

} // namespace pxl

#endif // pxl_pol_Event_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_io_InformationChunk_hh
#define pxl_io_InformationChunk_hh



namespace pxl
{

class InformationChunk : public pxl::Serializable
{
public:

	virtual const pxl::Id& getTypeId() const
	{
		return getStaticTypeId();
	}

	static const Id& getStaticTypeId()
	{
		static const Id id("7be73236-5038-4988-ba8e-9f65a26c4e72");
		return id;
	}

	virtual void serialize(const OutputStream &out) const
	{
		pxl::Serializable::serialize(out);
		_userRecords.serialize(out);
	}

	virtual void deserialize(const InputStream &in)
	{
		pxl::Serializable::deserialize(in);
		_userRecords.deserialize(in);
	}

	void setName(const std::string& name)
	{
		_name=name;
	}

	inline const std::string& getName() const
	{
		return _name;
	}

	/// This method provides access to the user records.
	inline const pxl::UserRecord& getUserRecord() const
	{
		return _userRecords;
	}

	/// This method sets the user record entry identified by \p key to \p item.
	template<typename datatype> inline void setUserRecord(
			const std::string& key, const datatype& item)
	{
		_userRecords.template set<datatype>(key, item);
	}

	/// This method removes the user record entry identified by \p key.
	inline void removeUserRecord(const std::string& key)
	{
		_userRecords.erase(key);
	}

	inline void clearUserRecords()
	{
		_userRecords.clear();
	}

	/// This method searches the user record entry identified by \p key; \p defaultitem is returned in case key is not found.
	template<typename datatype> inline datatype findUserRecord(
			const std::string& key, const datatype& defaultitem) const
	{
		return _userRecords.template find<datatype>(key, defaultitem);
	}

	/// This method searches the user record entry identified by \p key; an exception is thrown in case key is not found.
	template<typename datatype> inline datatype findUserRecord(
			const std::string& key) const throw(std::runtime_error)
	{
		return _userRecords.template find<datatype>(key);
	}

	/// This method checks if the user record entry identified by \p key is present.
	template<typename datatype> inline bool checkUserRecord(
			const std::string& key) const
	{
		return _userRecords.template check<datatype>(key);
	}

	/// This method checks if the user record entry identified by \p key is present. If yes, \p item is set to the according value.
	template<typename datatype> inline bool checkUserRecord(
			const std::string& key, datatype& item) const
	{
		return _userRecords.template check<datatype>(key, item);
	}

private:
	std::string _name;
	pxl::UserRecord _userRecords;

};

} //namespace pxl

#endif /*pxl_io_InformationChunk_hh*/

namespace pxl
{
// io
/**
 This abstract class offers the basic functionality for reading the PXL physics event structure.
 Derived classes can handle concrete I/O operations.
 */
class InputHandler
{
public:
	typedef skipMode::skipMode skipMode;

	InputHandler()
	{
	}

	virtual ~InputHandler()
	{
	}

	virtual ChunkReader& getChunkReader() = 0;

	/// This method reads in the header of the next event or information chunk.
	bool next()
	{
		if (getChunkReader().next())
			return true;
		return false;
	}

	bool nextEvent(skipMode doSkip = skipMode::on)
	{
		if (getChunkReader().nextEvent(doSkip))
			return true;
		return false;
	}

	bool nextInformationChunk(skipMode doSkip = skipMode::on)
	{
		if (getChunkReader().nextInformationChunk(doSkip))
			return true;
		return false;
	}

	/// This method reads in the next event if the information condition is fulfilled. Else, false is returned.
	bool nextIf(const std::string& info, skipMode doSkip = skipMode::on)
	{
		if (getChunkReader().next(doSkip, ChunkReader::evaluate, info))
			return true;
		return false;
	}

	/// This method reads in the next event if the information condition is fulfilled. Else, false is returned.
	bool nextEventIf(const std::string& info, skipMode doSkip = skipMode::on)
	{
		if (getChunkReader().nextEvent(doSkip, ChunkReader::evaluate, info))
			return true;
		return false;
	}

	/// This method reads in the next event if the information condition is fulfilled. Else, false is returned.
	bool nextInformationChunkIf(const std::string& info,
			skipMode doSkip = skipMode::on)
	{
		if (getChunkReader().nextInformationChunk(doSkip, ChunkReader::evaluate, info))
			return true;
		return false;
	}

	/// A pxl::Event is passed to this method and filled with the current event.
	bool readEvent(pxl::Event* event)
	{
		if (getChunkReader().nextBlock() && getChunkReader().getInputStream().good())
		{
			pxl::Id id(getChunkReader().getInputStream());
			if (id==event->getStaticTypeId())
			{
				event->deserialize(getChunkReader().getInputStream());
				return true;
			}
		}
		return false;
	}

	/// A pxl::Event is passed to this method and filled with the current event.
	bool readEventIf(pxl::Event* event, const std::string& blockInfo,
			skipMode doSkip = skipMode::on)
	{
		if (readBlockIf(blockInfo, doSkip))
		{
			pxl::Id id(getChunkReader().getInputStream());
			if (id==event->getStaticTypeId())
			{
				event->deserialize(getChunkReader().getInputStream());
				return true;
			}
		}
		return false;
	}

	/// A pxl::InformationChunk is passed to this method and filled if the next block contains an information chunk.
	bool readInformationChunk(pxl::InformationChunk* chunk)
	{
		if (getChunkReader().nextBlock() && getChunkReader().getInputStream().good())
		{
			pxl::Id id(getChunkReader().getInputStream());
			if (id==chunk->getStaticTypeId())
			{
				chunk->deserialize(getChunkReader().getInputStream());
				return true;
			}
		}
		return false;
	}

	/// A pxl::InformationChunk is passed to this method and filled if the info condition is fulfilled and there is an information chunk.
	bool readInformationChunkIf(pxl::InformationChunk* event,
			const std::string& blockInfo, skipMode doSkip = skipMode::on)
	{
		if (readBlockIf(blockInfo, doSkip))
		{
			pxl::Id id(getChunkReader().getInputStream());
			if (id==event->getStaticTypeId())
			{
				event->deserialize(getChunkReader().getInputStream());
				return true;
			}
		}
		return false;
	}

	inline bool isInformationChunk()
	{
		return getChunkReader().isInformationChunk();
	}

	inline bool isEvent()
	{
		return getChunkReader().isEvent();
	}

	/// This methods skips one event or information chunk.
	bool skip()
	{
		return getChunkReader().skip();
	}

	/// This method goes to the previous event.
	bool previous()
	{
		return getChunkReader().previous();
	}

	/// With this method, n events can be skipped in forward or backward direction.
	int skipEvents(int n)
	{
		int skipped = 0;
		for (; n<0 && getChunkReader().previous(); ++n)
			--skipped;
		for (; n>0 && getChunkReader().skip(); --n)
			++skipped;
		return skipped;
	}

	/// This method reads in the next block.
	bool readBlock()
	{
		return getChunkReader().nextBlock();
	}

	/// This method reads in the next block if the info condition is fulfilled.
	bool readBlockIf(const std::string& blockInfo,
			skipMode doSkip = skipMode::on)
	{
		bool success = false;
		while (!success && getChunkReader().getStatus()!=0)
		{
			success = getChunkReader().nextBlock(doSkip, ChunkReader::evaluate, blockInfo);
		}
		return success;
	}

	/// This method explicitly reads an object of type objecttype. Caution: This method should only be used if the type of the following object is known by hard.
	template<class objecttype> bool readObject(objecttype* obj) throw(std::runtime_error);

	/// This method fills the objects from the read-in block into the passed vector. The number of added objects is returned.
	int readObjects(std::vector<pxl::Serializable*>& objects);

	/// This method fills the objects from the read-in block into the passed pxl::Event. The number of added objects is returned.
	int readObjects(pxl::Event* event);

};

} //namespace pxl

#endif //pxl_io_InputHandler_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_iotl_OutputHandler_hh
#define pxl_iotl_OutputHandler_hh



namespace pxl
{

// io
/**
 This abstract class allows the user an easy handling of the PXL general output. Methods to write event headers,
 the PXL event class, and information chunks are offered.
 */
class OutputHandler
{
public:

	OutputHandler() :
		_newEvent(true)
	{
	}

	virtual ~OutputHandler()
	{
	}

	virtual pxl::ChunkWriter& getChunkWriter() = 0;

	/// This method writes the passed pxl::Event to the output file and finishes the current event.
	bool writeEvent(pxl::Event* event, const std::string& eventInfo = "",
			const std::string& blockInfo = "")
	{
		if (_newEvent)
		{
			getChunkWriter().newEvent(eventInfo);
			_newEvent = false;
		}
		getChunkWriter().newBlock();
		event->serialize(getChunkWriter().getOutputStream());
		getChunkWriter().write(blockInfo);
		_newEvent = true;
		return getChunkWriter().endEvent();
	}

	/// This method writes the passed pxl::InformationChunk to the output file and finishes the current event.
	bool writeInformationChunk(pxl::InformationChunk* infoChunk,
			const std::string& chunkInfo = "")
	{
		if (_newEvent)
		{
			getChunkWriter().newInformationChunk(chunkInfo);
			_newEvent = false;
		}
		getChunkWriter().newBlock();
		infoChunk->serialize(getChunkWriter().getOutputStream());
		getChunkWriter().write();
		_newEvent = true;
		return getChunkWriter().endInformationChunk();
	}

	/// This method queues the passed object for later writing to the output file.
	template<class objecttype> void streamObject(objecttype* obj)
	{
		obj->serialize(getChunkWriter().getOutputStream());
	}

	/// Use this method to write an information string describing the new event. Otherwise, this method need not necessarily be used.
	bool newEvent(const std::string& info)
	{
		if (!_newEvent)
		{
			std::cerr
					<< "Error in OutputFile::newEvent: Finish the current event first."
					<< std::endl;
			return false;
		}
		getChunkWriter().newEvent(info);
		_newEvent = false;
		return true;
	}

	/// Use this method to write out a block to file. This method is not needed if you use the writeEvent-method.
	bool writeStream(const std::string& info = "")
	{
		if (_newEvent)
		{
			getChunkWriter().newEvent("");
			_newEvent = false;
		}
		getChunkWriter().newBlock();
		return getChunkWriter().write(info);
	}

	/// Use this method to write out a block to disk and finish the current event.
	bool writeEvent(const std::string& info = "")
	{
		if (_newEvent)
		{
			getChunkWriter().newEvent("");
			_newEvent = false;
		}
		getChunkWriter().newBlock();
		getChunkWriter().write(info);
		_newEvent = true;
		return getChunkWriter().endEvent();
	}

private:
	bool _newEvent;
};

}//namespace pxl

#endif /*pxl_iotl_OutputHandler_hh*/
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_io_InputFile_hh
#define pxl_io_InputFile_hh



namespace pxl
{
// io
/**
 This class offers an easy handling of the PXL I/O. Various methods to access
 the content of I/O files are offered.
 */
class InputFile : public InputHandler
{
public:

	InputFile() :
		InputHandler(), _stream(), _reader(_stream)
	{
	}

	InputFile(const std::string& filename) :
		InputHandler(), _stream(filename.c_str(), std::ios::binary),
				_reader(_stream)
	{
	}

	virtual void open(const std::string& filename)
	{
		_stream.open(filename.c_str(), std::ios::binary);
	}

	virtual void close()
	{
		_stream.close();
	}

	virtual ~InputFile()
	{
		_stream.close();
	}

	virtual pxl::ChunkReader& getChunkReader()
	{
		return _reader;
	}

private:
	std::ifstream _stream;
	pxl::ChunkReader _reader;
};

} //namespace pxl

#endif /*pxl_io_InputFile_hh*/
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_iotl_OutputFile_hh
#define pxl_iotl_OutputFile_hh



namespace pxl
{

// io
/**
 This class allows the user an easy handling of the PXL output to file. Methods to write event headers,
 the PXL event class, and information chunks are offered by inheritance from the OutputHandler.
 */

class OutputFile : public OutputHandler
{
public:

	OutputFile(const std::string& filename) :
		OutputHandler(), _stream(filename.c_str(), std::ios::binary),
				_writer(_stream)
	{
	}

	virtual ~OutputFile()
	{
		_stream.close();
	}
	
	virtual void open(const std::string& filename)
	{
		_stream.open(filename.c_str(), std::ios::binary);
	}
	
	virtual void close()
	{
		_stream.close();
	}

	virtual pxl::ChunkWriter& getChunkWriter()
	{
		return _writer;
	}

private:
	std::ofstream _stream;
	pxl::ChunkWriter _writer;
};

}//namespace pxl

#endif /*pxl_io_OutputFile_hh*/
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_io_GenericInputHandler_hh
#define pxl_io_GenericInputHandler_hh



namespace pxl
{
// io
/**
 This class offers a generic handling of the PXL I/O. Various methods to access
 the PXL I/O content are offered.
 */
class GenericInputHandler : public InputHandler
{
public:

	GenericInputHandler(pxl::ChunkReader& reader) :
		InputHandler(),	_reader(&reader)
	{
	}

	virtual ~GenericInputHandler()
	{
	}

	virtual pxl::ChunkReader& getChunkReader() throw (std::runtime_error)
	{
		if (!_reader)
			throw std::runtime_error("GenericInputHandler::getChunkReader(): ChunkReader pointer invalid.");			
		return *_reader;
	}
	
	virtual void setChunkReader(pxl::ChunkReader* reader)
	{
		_reader=reader;
	}

private:
	pxl::ChunkReader* _reader;
};

} //namespace pxl

#endif /*pxl_io_GenericInputHandler_hh*/
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_iotl_GenericOutputHandler_hh
#define pxl_iotl_GenericOutputHandler_hh



namespace pxl
{

// io
/**
 This class allows the user an easy handling of the PXL output using any ChunkWriter. Methods to write event headers,
 the PXL event class, and information chunks are offered by inheritance from the OutputHandler.
 */

class GenericOutputHandler : public OutputHandler
{
public:
	GenericOutputHandler(pxl::ChunkWriter& writer) :
		OutputHandler(), _writer(&writer)
	{
	}

	virtual ~GenericOutputHandler()
	{
	}

	virtual pxl::ChunkWriter& getChunkWriter() throw(std::runtime_error)
	{
		if (!_writer)
			throw std::runtime_error("GenericOutputHandler::getChunkWriter(): ChunkWriter pointer invalid.");
		return *_writer;
	}

	virtual void setChunkWriter(pxl::ChunkWriter* writer)
	{
		_writer=writer;
	}

private:
	pxl::ChunkWriter* _writer;
};

}//namespace pxl

#endif /*pxl_iotl_GenericOutputHandler_hh*/

#endif // pxl_io_hh

#endif // pxl_base_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_pol_hh
#define pxl_pol_hh

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_pol_Basic3Vector_hh
#define pxl_pol_Basic3Vector_hh

#include <cmath>


namespace pxl
{
// pol

#define EPSILON 1.0e-9
#ifndef M_PI
#define M_PI        3.14159265358979323846
#endif

/** 
 This class provides a simple threevector with basic algebra. The methods provided are self-explanatory.
 */
class Basic3Vector
{
public:
    Basic3Vector() :
        _x(0), _y(0), _z(0)
    {
    }
    Basic3Vector(const Basic3Vector& orig) :
        _x(orig._x), _y(orig._y), _z(orig._z)
    {
    }
    Basic3Vector(const Basic3Vector* orig) :
        _x(orig->_x), _y(orig->_y), _z(orig->_z)
    {
    }
    Basic3Vector(double x, double y, double z) :
        _x(x), _y(y), _z(z)
    {
    }

	virtual void serialize(const OutputStream &out) const
	{
		out.writeDouble(_x);
		out.writeDouble(_y);
		out.writeDouble(_z);
	}

	virtual void deserialize(const InputStream &in)
	{
	    in.readDouble(_x);
	    in.readDouble(_y);
	    in.readDouble(_z);
	}
    
    inline void setX(double x)
    {
        _x = x;
    }
    inline void setY(double y)
    {
        _y = y;
    }
    inline void setZ(double z)
    {
        _z = z;
    }

    inline void setRhoPhi(double perp, double phi)
    {
        _x = std::cos(phi) * perp;
        _y = std::sin(phi) * perp;
    }
    inline void setRhoPhiZ(double perp, double phi, double z)
    {
        setRhoPhi(perp, phi);
        _z = z;
    }
    inline void setRThetaPhi(double r, double theta, double phi)
    {
        setRhoPhiZ(std::cos(theta) * r, phi, std::sin(theta) * r);
    }

    inline double getX() const
    {
        return _x;
    }
    inline double getY() const
    {
        return _y;
    }
    inline double getZ() const
    {
        return _z;
    }

    inline bool isNullPerp() const
    {
        return _x > -EPSILON && _x< EPSILON && _y> -EPSILON && _y < EPSILON;
    }
    inline bool isNull() const
    {
        return isNullPerp() && _z > -EPSILON && _z < EPSILON;
    }

    inline double getPerp2() const
    {
        return _x*_x + _y*_y;
    }
    inline double getPerp() const
    {
        return std::sqrt(getPerp2());
    }
    inline double getPhi() const
    {
        return isNullPerp() ? 0.0 : std::atan2(_y, _x);
    }

    inline double getMag2() const
    {
        return _x*_x + _y*_y + _z*_z;
    }
    inline double getMag() const
    {
        return std::sqrt(getMag2());
    }

    inline double getCosTheta() const
    {
        double mag = getMag();
        return mag < EPSILON ? 1.0 : _z / mag;
    }
    inline double getCos2Theta() const
    {
        double mag2 = getMag2();
        return mag2 < EPSILON ? 1.0 : _z*_z / mag2;
    }

    inline double getTheta() const
    {
        return isNull() ? 0.0 : std::atan2(getPerp(), _z);
    }

    inline double deltaRho(const Basic3Vector* fv) const
    {
        double dDtheta = deltaTheta(fv);
        double dDphi = deltaPhi(fv);
        return std::sqrt(dDtheta*dDtheta + dDphi*dDphi);
    }

    inline double deltaPhi(const Basic3Vector* fv) const
    {
        double dDphi = getPhi() - fv->getPhi();
        while (dDphi > M_PI)
            dDphi -= 2 * M_PI;
        while (dDphi < -M_PI)
            dDphi += 2 * M_PI;
        return dDphi;
    }

    inline double deltaTheta(const Basic3Vector* fv) const
    {
        double dDtheta = getTheta() - fv->getTheta();
        while (dDtheta > M_PI)
            dDtheta -= 2 * M_PI;
        while (dDtheta < -M_PI)
            dDtheta += 2 * M_PI;
        return dDtheta;
    }

    inline const pxl::Basic3Vector& operator=(const pxl::Basic3Vector& vec)
    {
        _x = vec._x;
        _y = vec._y;
        _z = vec._z;
        return *this;
    }

    //FIXME: add functions for pointers, eg pxl::Basic3Vector* add(pxl::Basic3Vector* vec);
    inline const pxl::Basic3Vector& operator+=(const pxl::Basic3Vector& vec)
    {
        _x += vec._x;
        _y += vec._y;
        _z += vec._z;
        return *this;
    }
    inline const pxl::Basic3Vector& operator-=(const pxl::Basic3Vector& vec)
    {
        _x -= vec._x;
        _y -= vec._y;
        _z -= vec._z;
        return *this;
    }

private:
    double _x;
    double _y;
    double _z;

};

#undef EPSILON

// non-member operators
bool const operator==(const pxl::Basic3Vector& obj1, const pxl::Basic3Vector& obj2);
bool const operator!=(const pxl::Basic3Vector& obj1, const pxl::Basic3Vector& obj2);

} // namespace pxl

#endif // pxl_pol_Basic3Vector_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_pol_Basic4Vector_hh
#define pxl_pol_Basic4Vector_hh



namespace pxl
{
// pol

/** 
 This class provides a simple Lorentz-fourvector with basic algebra. The methods provided are self-explanatory.
 */
class Basic4Vector : public Basic3Vector
{
public:
    Basic4Vector() :
        Basic3Vector(0, 0, 0), _t(0)
    {
    }
    Basic4Vector(const Basic3Vector &orig, double t = 0) :
        Basic3Vector(orig), _t(t)
    {
    }
    Basic4Vector(const Basic4Vector &orig) :
        Basic3Vector(orig), _t(orig._t)
    {
    }
    Basic4Vector(const Basic4Vector* orig) :
        Basic3Vector(orig), _t(orig->_t)
    {
    }
    Basic4Vector(double x, double y, double z, double t = 0) :
        Basic3Vector(x, y, z), _t(t)
    {
    }

	virtual void serialize(const OutputStream &out) const
	{
		pxl::Basic3Vector::serialize(out);
		out.writeDouble(_t);
	}

	virtual void deserialize(const InputStream &in)
	{
		pxl::Basic3Vector::deserialize(in);
		in.readDouble(_t);
	}
    
    // setX inherited
    // setY inherited
    // setZ inherited
    inline void setT(double t)
    {
        _t = t;
    }

    inline void setPx(double px)
    {
        setX(px);
    }
    inline void setPy(double py)
    {
        setY(py);
    }
    inline void setPz(double pz)
    {
        setZ(pz);
    }

    inline void setE(double e)
    {
        _t = e;
    }
    inline void setMass(double m)
    {
        _t = std::sqrt(m*m + getMag2());
    }

    // getX inherited
    // getY inherited
    // getZ inherited
    inline double getT() const
    {
        return _t;
    }

    inline double getPx() const
    {
        return getX();
    }
    inline double getPy() const
    {
        return getY();
    }
    inline double getPz() const
    {
        return getZ();
    }

    inline double getE() const
    {
        return _t;
    }

    inline double getMass2() const
    {
        return _t*_t - getMag2();
    }

    inline double getMass() const
    {
        double m2 = getMass2();
        return m2 < 0.0 ? 0.0 : std::sqrt(m2);
    }
    // getPerp inherited
    inline double getPt() const
    {
        return getPerp();
    }

    // getPhi inherited
    // getTheta inherited
    // deltaRho inherited
    // deltaPhi inherited
    // deltaTheta inherited

    inline double getEta() const
    {
        return -std::log(std::tan(getTheta()*0.5));
    }

    inline double getEt2() const
    {
        double pt2 = getPerp2();
        return pt2 == 0.0 ? 0.0 : _t*_t * pt2 / getMag2();
    }
    inline double getEt() const
    {
        return std::sqrt(getEt2());
    }

    inline double deltaR(const Basic4Vector* fv) const
    {
        double dDeta = deltaEta(fv);
        double dDphi = deltaPhi(fv);
        return std::sqrt(dDeta*dDeta + dDphi*dDphi);
    }

    inline double deltaEta(const Basic4Vector* fv) const
    {
        return getEta() - fv->getEta();
    }

    inline const pxl::Basic4Vector& operator=(const pxl::Basic3Vector& vec)
    {
        Basic3Vector::operator=(vec); return *this;}
    inline const pxl::Basic4Vector& operator=(const pxl::Basic4Vector& vec)
    {   Basic3Vector::operator=(vec); _t = vec._t; return *this;}

    inline const pxl::Basic4Vector& operator+=(const pxl::Basic4Vector& vec)
    {   Basic3Vector::operator+=(vec); _t += vec._t; return *this;}
    inline const pxl::Basic4Vector& operator-=(const pxl::Basic4Vector& vec)
    {   Basic3Vector::operator-=(vec); _t -= vec._t; return *this;}

private:
    double _t;

};

// non-member operators
bool const operator==(const pxl::Basic4Vector& obj1, const pxl::Basic4Vector& obj2);
bool const operator!=(const pxl::Basic4Vector& obj1, const pxl::Basic4Vector& obj2);

} // namespace pxl


#endif // pxl_pol_Basic4Vector_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_pol_BasicObject_hh
#define pxl_pol_BasicObject_hh




namespace pxl
{

// pol

/** 
 This class provides common functionalities of PXL physics objects like 
 data members for storing an object name and flags for status, Monte-Carlo mode and 
 object locking; more specific information, such as b-tags, jet cone sizes or energy 
 corrections, for instance, can be stored in the so-called user records (see pxl::UserRecord). 
 An integer workflag facilitates tagging of individual objects. 
 */
class Object : public pxl::Relative
{
public:
	Object() :
		pxl::Relative(), _locked(0), _workflag(0), _userRecords()
	{
	}

	explicit Object(const Object* original) :
		pxl::Relative(original), _locked(original->_locked),
				_workflag(original->_workflag),
				_userRecords(original->_userRecords)
	{
	}

	virtual const pxl::Id& getTypeId() const
	{
		return getStaticTypeId();
	}

	static const Id& getStaticTypeId()
	{
		static const Id id("3b3a2442-04f6-400e-8e30-1de2dbc8d628");
		return id;
	}

	virtual void serialize(const OutputStream &out) const
	{
		pxl::Relative::serialize(out);
		out.writeBool(_locked);
		out.writeInt(_workflag);
		_userRecords.serialize(out);
	}

	virtual void deserialize(const InputStream &in)
	{
		pxl::Relative::deserialize(in);
		in.readBool(_locked);
		in.readInt(_workflag);
		_userRecords.deserialize(in);
	}

	/// This method returns the value of the lock flag.
	inline bool getLocked() const
	{
		return _locked;
	}

	/// This method returns the value of the workflag.
	inline int getWorkflag() const
	{
		return _workflag;
	}

	/// This method sets the value of the lock flag to \p v.
	inline void setLocked(bool v)
	{
		_locked = v;
	}

	/// This method sets the value of the workflag to \p v.
	inline void setWorkflag(int v)
	{
		_workflag = v;
	}

	/// This method provides access to the user records.
	inline const pxl::UserRecord& getUserRecord() const
	{
		return _userRecords;
	}

	/// This method sets the user record entry identified by \p key to \p item.
	template<typename datatype> inline void setUserRecord(
			const std::string& key, const datatype& item)
	{
		_userRecords.template set<datatype>(key, item);
	}

	/// This method removes the user record entry identified by \p key.
	inline void removeUserRecord(const std::string& key)
	{
		_userRecords.erase(key);
	}

	inline void clearUserRecords()
	{
		_userRecords.clear();
	}

	/// This method searches the user record entry identified by \p key; \p defaultitem is returned in case key is not found.
	template<typename datatype> inline datatype findUserRecord(
			const std::string& key, const datatype& defaultitem) const
	{
		return _userRecords.template find<datatype>(key, defaultitem);
	}

	/// This method searches the user record entry identified by \p key; an exception is thrown in case key is not found.
	template<typename datatype> inline datatype findUserRecord(
			const std::string& key) const throw(std::runtime_error)
	{
		return _userRecords.template find<datatype>(key);
	}

	/// This method checks if the user record entry identified by \p key is present.
	inline bool checkUserRecord(const std::string& key) const
	{
		return _userRecords.check(key);
	}

	/// This method checks if the user record entry identified by \p key is present. If yes, \p item is set to the according value.
	template<typename datatype> inline bool checkUserRecord(
			const std::string& key, datatype& item) const
	{
		return _userRecords.template check<datatype>(key, item);
	}

	inline const pxl::Object& operator=(const pxl::Object& pa)
	{
		pxl::Relative::operator=(pa);
		_locked = pa._locked;
		_workflag = pa._workflag;
		_userRecords = pa._userRecords;
		
		return *this;
	}

	/// This virtual method is intended to print out object state information on various verbosity levels.
	/// @param level verbosity level
	/// @param os output _stream, default is std::cout
	/// @param pan print indention
	/// @return output _stream
	virtual std::ostream& print(int level = 1, std::ostream& os = std::cout,
			int pan = 0) const;

	virtual std::ostream& printContent(int level = 1,
			std::ostream& os = std::cout, int pan = 0) const;

	virtual pxl::WkPtrBase* createSelfWkPtr()
	{
		return new pxl::weak_ptr<Object>(this);
	}

private:
	bool _locked;
	int _workflag;

	pxl::UserRecord _userRecords;
};

///// This typedef defines a weak pointer for pxl::Object
typedef pxl::weak_ptr<Object> ObjectWkPtr;

} // namespace pxl

#endif // pxl_pol_BasicObject_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_pol_BasicManager_hh
#define pxl_pol_BasicManager_hh




namespace pxl
{

// pol

/** 
 This class the functionality of the pxl::BasicObjectData class by providing an object owner (see pxl::ObjectOwner) and 
 corresponding service methods. This way, physics objects (like instances of the classes 
 pxl::Particle, pxl::Vertex or pxl::Collision as well as other arbitrary pxl::Relative derivatives can be 
 aggregated and managed. 
 */
class ObjectManager : public pxl::Object
{
public:
	ObjectManager() :
		Object(), _objects()
	{
	}
	/// This copy constructor performs a deep copy of \p original 
	/// with all contained objects and their (redirected) relations.
	ObjectManager(const pxl::ObjectManager& original) :
		Object(original), _objects(original._objects)
	{
	}
	/// This copy constructor performs a deep copy of \p original 
	/// with all contained objects and their (redirected) relations.
	ObjectManager(const pxl::ObjectManager* original) :
		Object(original), _objects(original->_objects)
	{
	}

	virtual const pxl::Id& getTypeId() const
	{
		return getStaticTypeId();
	}

	static const Id& getStaticTypeId()
	{
		static const Id id("86cab8f4-6c08-477d-a4e9-bd718d6f899f");
		return id;
	}

	virtual void serialize(const OutputStream &out) const
	{
		pxl::Object::serialize(out);
		_objects.serialize(out);
	}

	virtual void deserialize(const InputStream &in)
	{
		pxl::Object::deserialize(in);
		_objects.deserialize(in);
	}

	// create
	/// This method template creates a new instance of \p objecttype;
	/// objecttype must be a class inheriting from pxl::Relative;
	/// the newly created instance is owned and will be deleted by the object owner. 
	template<class datatype> datatype* create()
	{
		return _objects.create<datatype>();
	}

	/// This method template creates a new \p objecttype instance by invoking a \p ctrtype overloaded constructor; 
	/// \p objecttype must be a class inheriting from pxl::Relative;
	/// the newly created instance is owned and will be deleted by the object owner. 
#ifndef SWIG
	template<class datatype, class ctrdatatype> datatype* create(
			const ctrdatatype& ori)
	{
		return _objects.create<datatype,ctrdatatype>(ori);
	}
#endif
	// crateIndexed
	/// This method template acts like create() and registers the newly created instance under \p idx in the index.
	template<class datatype> datatype* createIndexed(const std::string& idx)
	{
		datatype* obj = _objects.create<datatype>();
		setIndex(idx, obj);
		return obj;
	}

	/// This method template acts like create() and registers the newly created instance under \p idx in the index.
	template<class datatype, class ctrdatatype> datatype* createIndexed(
			const ctrdatatype& ori, const std::string& idx)
	{
		datatype* obj = _objects.create<datatype,ctrdatatype>(ori);
		setIndex(idx, obj);
		return obj;
	}

	/// This method inserts \p obj in the container of the object owner and takes deletion responsability.
	inline void setObject(pxl::Relative* obj)
	{
		_objects.set(obj);
	}

	/// This method inserts \p obj with the index-id \p idx in the container of the object owner and takes deletion responsability.
	inline void setObject(pxl::Relative* obj, const std::string& idx)
	{
		_objects.set(obj);
		setIndex(idx, obj);
	}

	/// This method registers the object \p obj with the index-id \p idx in the index and returns true in case of success;
	/// please notice, that obj must be owned by this object owner and \p idx must not be a zero length string.  
	inline bool setIndex(const std::string& idx, pxl::Relative* obj)
	{
		return _objects.setIndex(idx, obj);
	}

	/// This method provides access to the object owner.
	inline pxl::ObjectOwner& getObjectOwner()
	{
		return _objects;
	}

	inline const pxl::ObjectOwner& getObjectOwner() const
	{
		return _objects;
	}

	inline const std::vector<pxl::Relative*>& getObjects() const
	{
		return _objects.getObjects();
	}

	template<class objecttype> inline void getObjectsOfType(
			std::vector<objecttype*>& vec) const
	{
		_objects.getObjectsOfType<objecttype>(vec);
	}

	template<class objecttype, class comparator> inline void getObjectsOfType(
			std::vector<objecttype*>& vec, const FilterCriterion<objecttype>& criterion) const
	{
	        pxl::Filter<objecttype, comparator> filter; 
		filter.apply(&_objects, vec, criterion);
		//_objects.getObjectsOfType<objecttype>(vec);
	}   
	/*
	template<> inline void getObjectsOfType(std::vector<Particle*>& vec, const FilterCriterion<Particle>& criterion)
	const
	{
	   pxl::Filter<Particle, ParticlePtComparator> filter;
	   filter.apply(&_objects, vec, criterion);
	}     */

	/// This method deletes the object \p obj.
	inline void removeObject(pxl::Relative* obj)
	{
		_objects.remove(obj);
	}

	/// This method clears the object owner and deletes all owned objects. 
	inline void clearObjects()
	{
		_objects.clearContainer();
	}

	/// This method searches the index for the index-id \p idx and returns a dynamically casted 
	/// C++ pointer of type \p objecttype* to the corresponding object; 
	/// in case idx is not found a null pointer is returned.
	template<class objecttype> inline objecttype* findObject(
			const std::string idx) const
	{
		return _objects.findObject<objecttype>(idx);
	}

	/// This method searches the copy history to locate the copy of \p original and 
	/// returns a dynamically casted C++ pointer of type \p objecttype* to the corresponding copy; 
	/// in case no copy can be traced a null pointer is returned.
	template<class objecttype> inline objecttype* findCopyOf(
			const pxl::Relative* original) const
	{
		return _objects.findCopyOf<objecttype>(original);
	}

	/// This method provides direct access to the copy history (created by the copy constructor). 
	inline const pxl::CopyHistory& getCopyHistory() const
	{
		return _objects.getCopyHistory();
	}

	/// This method clears the copy history  (created by the copy constructor). 
	inline void clearCopyHistory()
	{
		_objects.clearCopyHistory();
	}

	/// This method provides direct access to the index. 
	inline const pxl::Index& getIndex() const
	{
		return _objects.getIndex();
	}

	/// This method removes the index entry with index-id \p idx; please notice: it does not remove the object itself. 
	inline void removeIndex(const std::string& idx)
	{
		_objects.removeIndex(idx);
	}

	/// This method clears the index; please notice: it does not remove the objects themself.
	inline void clearIndex()
	{
		_objects.clearIndex();
	}

	virtual pxl::WkPtrBase* createSelfWkPtr()
	{
		return new pxl::weak_ptr<ObjectManager>(this);
	}

private:
	pxl::ObjectOwner _objects;
};

///// This typedef defines a weak pointer for pxl::ObjectManager
typedef pxl::weak_ptr<pxl::ObjectManager> ObjectManagerWkPtr;

} // namespace pxl


#endif // pxl_pol_BasicManager_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_pol_Vertex_hh
#define pxl_pol_Vertex_hh



//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_pol_CommonVertex_hh
#define pxl_pol_CommonVertex_hh

namespace pxl
{
/**
 * This is the common, pure virtual interface class for vertices.
 */

class CommonVertex
{
public:
	//getters for basic vector quantities.
	virtual double getX() const = 0;
	virtual double getY() const = 0;
	virtual double getZ() const = 0;

	//setters for basic fourvector quantities
	virtual void setX(double x) = 0;
	virtual void setY(double y) = 0;
	virtual void setZ(double z) = 0;

	//setters for basic fourvector quantities
	virtual void setXYZ(double x, double y, double z) = 0;
	virtual void addXYZ(double x, double y, double z) = 0;

};

}

#endif /*pxl_pol_CommonVertex_hh*/

namespace pxl
{

// pol
/**
 This class allows to store threevector and further properties of the decay vertex; see also pxl::BasicObjectData.
 */
class Vertex : public pxl::Object, public pxl::CommonVertex
{
public:
	Vertex() :
		Object(), _vector()
	{
	}

	Vertex(const Vertex& original) :
		Object(original), _vector(original._vector)
	{
	}
	Vertex(const Vertex* original) :
		Object(original), _vector(original->_vector)
	{
	}

	virtual const pxl::Id& getTypeId() const
	{
		return getStaticTypeId();
	}

	static const Id& getStaticTypeId()
	{
		static const Id id("80da63d3-c838-466d-9d5b-cddb6110f0e3");
		return id;
	}

	virtual void serialize(const OutputStream &out) const
	{
		pxl::Object::serialize(out);
		_vector.serialize(out);
	}

	virtual void deserialize(const InputStream &in)
	{
		pxl::Object::deserialize(in);
		_vector.deserialize(in);
	}

	/// This method grants read access to the vector.
    inline const pxl::Basic3Vector& getVector() const
	{
		return _vector;
	}

	inline double getX() const
	{
		return _vector.getX();
	}
	inline double getY() const
	{
		return _vector.getY();
	}
	inline double getZ() const
	{
		return _vector.getZ();
	}

	inline void setX(double x)
	{
		_vector.setX(x);
	}
	inline void setY(double y)
	{
		_vector.setY(y);
	}
	inline void setZ(double z)
	{
		_vector.setZ(z);
	}

	inline void setXYZ(double x, double y, double z)
	{
		_vector.setX(x);
		_vector.setY(y);
		_vector.setZ(z);
	}
	inline void setVector(const pxl::Basic3Vector& vector)
	{
		_vector = vector;
	}

	inline void addXYZ(double x, double y, double z)
	{
		_vector.setX(x + _vector.getX());
		_vector.setY(y + _vector.getY());
		_vector.setZ(z + _vector.getZ());
	}

	inline void addVector(const pxl::Basic3Vector& vector)
	{
		_vector+=vector;
	}

	inline void addVertex(const pxl::Vertex* vx)
	{
		_vector += vx->getVector();
	}

	/// This method adds the vector of \p vxd.
	inline const pxl::Vertex& operator+=(const pxl::Vertex& vx)
	{
		_vector += vx._vector;
		return *this;
	}
	/// This method subtracts the vector of \p vxd.
	inline const pxl::Vertex& operator-=(const pxl::Vertex& vx)
	{
		_vector -= vx._vector;
		return *this;
	}

	virtual pxl::Relative* clone() const
	{
		return new pxl::Vertex(*this);
	}

	virtual std::ostream& print(int level = 1, std::ostream& os = std::cout,
			int pan = 0) const;

private:
	pxl::Basic3Vector _vector;

};

// non-member operators
bool const operator==(const pxl::Vertex& obj1, const pxl::Vertex& obj2);
bool const operator!=(const pxl::Vertex& obj1, const pxl::Vertex& obj2);

// typedefs

/// This typedef defines a weak pointer for pxl::Vertex
typedef pxl::weak_ptr<pxl::Vertex> VertexWkPtr;

} // namespace pxl

#endif // pxl_pol_Vertex_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_pol_Particle_hh
#define pxl_pol_Particle_hh


//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_pol_CommonParticle_hh
#define pxl_pol_CommonParticle_hh

namespace pxl
{
/**
 * This is the common, pure virtual interface class for particles.
 * The physical representation, (px, py, pz, E) or (pt, eta, phi, m/E), can
 * be chosen by the concrete implementation.
 */

class CommonParticle
{
public:
	//getters for basic fourvector quantities in (px, py, pz, E)-representation
	virtual double getPx() const = 0;
	virtual double getPy() const = 0;
	virtual double getPz() const = 0;
	virtual double getE() const = 0;

	//getters for basic fourvector quantities in (pt, eta, phi, mass)-representation
	virtual double getPt() const = 0;
	virtual double getEta() const = 0;
	virtual double getPhi() const = 0;
	virtual double getMass() const = 0;
	
	virtual double getCharge() const = 0;

	//setters for basic fourvector quantities
	virtual void setP4(double px, double py, double pz, double e) = 0;
	virtual void addP4(double px, double py, double pz, double e) = 0;
	virtual void setCharge(double q) = 0;

};

}

#endif /*pxl_pol_CommonParticle_hh*/

namespace pxl
{
// pol
/**
 This class allows to store Lorentz-fourvector and further properties of particles or reconstructed objects
 such as charge, particle-id plus the inherited properties of pxl::BasicObjectData.
 */
class Particle : public pxl::Object, public pxl::CommonParticle
{
public:
	Particle() :
		Object(), _charge(0), _particleId(0)
	{
	}

	explicit Particle(const Particle& original) :
		Object(original), _vector(original._vector), _charge(original._charge),
				_particleId(original._particleId)
	{
	}

	explicit Particle(const Particle* original) :
		Object(original), _vector(original->_vector), _charge(original->_charge),
				_particleId(original->_particleId)
	{
	}

	virtual const pxl::Id& getTypeId() const
	{
		return getStaticTypeId();
	}

	static const Id& getStaticTypeId()
	{
		static const Id id("c5515a0d-36bf-4076-bf33-e14343cf5a88");
		return id;
	}

	virtual void serialize(const OutputStream &out) const
	{
		pxl::Object::serialize(out);
		_vector.serialize(out);
		out.writeDouble(_charge);
		out.writeInt(_particleId);
	}

	virtual void deserialize(const InputStream &in)
	{
		pxl::Object::deserialize(in);
		_vector.deserialize(in);
		in.readDouble(_charge);
		in.readInt(_particleId);
	}

	/// This method grants read access to the vector.
        inline const pxl::Basic4Vector& getVector() const
	{
		return _vector;
	}

	/// This method returns the particle charge.
	inline double getCharge() const
	{
		return _charge;
	}
	/// This method sets the particle charge to v.
	inline void setCharge(double v)
	{
		_charge = v;
	}

	/// This method returns the particle-id.
	inline int getParticleId() const
	{
		return _particleId;
	}
	/// This method sets the particle-id to \p v.
	inline void setParticleId(int v)
	{
		_particleId = v;
	}

	inline const pxl::Particle& operator=(const pxl::Particle& pa)
	{
		pxl::Object::operator=(pa);
		_vector = pa._vector;
		_charge = pa._charge;
		_particleId = pa._particleId;
		return *this;
	}

	/// This method adds vector and charge of \p pad.
	inline const pxl::Particle& operator+=(const pxl::Particle& pa)
	{
		_vector += pa._vector;
		_charge += pa._charge;
		return *this;
	}

	/// This method subtracts vector and charge of of \p pad.
	inline const pxl::Particle& operator-=(const pxl::Particle& pa)
	{
		_vector -= pa._vector;
		_charge += pa._charge;
		return *this;
	}

	virtual pxl::Relative* clone() const
	{
		return new pxl::Particle(*this);
	}

	inline double getPx() const
	{
		return _vector.getPx();
	}

	inline double getPy() const
	{
		return _vector.getPy();
	}

	inline double getPz() const
	{
		return _vector.getPz();
	}

	inline double getE() const
	{
		return _vector.getE();
	}

	inline double getMass() const
	{
		return _vector.getMass();
	}

	inline double getPt() const
	{
		return _vector.getPt();
	}

	inline double getEta() const
	{
		return _vector.getEta();
	}

	inline double getEt() const
	{
		return _vector.getEt();
	}

	inline double getPhi() const
	{
		return _vector.getPhi();
	}

	inline double getTheta() const
	{
		return _vector.getTheta();
	}

	inline void setP4(double px, double py, double pz, double e)
	{
		_vector.setPx(px);
		_vector.setPy(py);
		_vector.setPz(pz);
		_vector.setE(e);
	}

	inline void setP4(const pxl::Basic4Vector& vector)
	{
		_vector = vector;
	}

	inline void addP4(double px, double py, double pz, double e)
	{
		_vector.setPx(px + _vector.getPx());
		_vector.setPy(py + _vector.getPy());
		_vector.setPz(pz + _vector.getPz());
		_vector.setE(e + _vector.getE());
	}

	void setP4FromDaughters();

	void setP4FromDaughtersRecursive();

	inline void addP4(const pxl::Basic4Vector& vector)
	{
		_vector+=vector;
	}

	inline void addP4(const pxl::Particle* particle)
	{
		_vector+=particle->getVector();
	}

	inline void addParticle(const pxl::Particle* pa)
	{
		_vector += pa->getVector();
		_charge += pa->getCharge();
	}

	virtual std::ostream& print(int level = 1, std::ostream& os = std::cout,
			int pan = 0) const;

	virtual pxl::WkPtrBase* createSelfWkPtr()
	{
		return new pxl::weak_ptr<Particle>(this);
	}

private:
	pxl::Basic4Vector _vector;
	double _charge;
	int _particleId;

};

// non-member operators
bool const operator==(const pxl::Particle& obj1, const pxl::Particle& obj2);
bool const operator!=(const pxl::Particle& obj1, const pxl::Particle& obj2);

// typedefs
/**
 This typedef represents particles and reconstructed
 objects such as muons, electrons, photons, jets; data is aggregated in pxl::ParticleData.
 */
typedef pxl::weak_ptr<pxl::Particle> ParticleWkPtr;

} // namespace pxl


#endif // pxl_pol_Particle_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_pol_Collision_hh
#define pxl_pol_Collision_hh


namespace pxl
{
/**
 This class represents individual interactions in  multicollision events.
 It allows the separation of different collisions as they occur
 at high-rate hadron colliders by providing the relation management
 necessary to associate pxl::Vertex or pxl::Particle objects, for instance.
 */

class Collision : public pxl::Object
{
public:
	Collision() :
		Object()
	{
	}
	/// This copy constructor provides a deep copy of the event container \p original with all data members,
	/// physics objects, and their (redirected) relations.
	Collision(const pxl::Collision& original) :
		Object(original)
	{
	}
	/// This copy constructor provides a deep copy of the event container \p original with all data members,
	/// physics objects, and their (redirected) relations.
	Collision(const pxl::Collision* original) :
		Object(original)
	{
	}

	virtual pxl::WkPtrBase* createSelfWkPtr()
	{
		return new pxl::weak_ptr<Collision>(this);
	}

	virtual const pxl::Id& getTypeId() const
	{
		return getStaticTypeId();
	}

	static const Id& getStaticTypeId()
	{
		static const Id id("59b2f95c-5142-4970-844f-226ebbc57a99");
		return id;
	}

	virtual void serialize(const OutputStream &out) const
	{
		pxl::Object::serialize(out);
	}

	virtual void deserialize(const InputStream &in)
	{
		pxl::Object::deserialize(in);
	}

	std::ostream& print(int level=1, std::ostream& os=std::cout, int pan=0) const;

};

/// This typedef defines a weak pointer for pxl::Collision
typedef pxl::weak_ptr<pxl::Collision> CollisionWkPtr;

} // namespace pxl

#endif // pxl_pol_Collision_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_pol_EventView_hh
#define pxl_pol_EventView_hh



namespace pxl
{

// pol
/**
 By inheritance from pxl::ObjectManager, 
 this class is capable of holding the complete information of one 
 multicollision event with decay trees, spatial vertex information, 
 four-momenta as well as additional event-related reconstruction data 
 in the user records. Physics objects (i.e. instances of the classes pxl::Particle, 
 pxl::Vertex or pxl::Collision) as well as other arbitrary pxl::Relative 
 derivatives can be aggregated and managed.
 The name 'event view'  arises from the fact that it is 
 intended to represent a distinct view of an event (e.g. connecting 
 particles to the decay tree according to one out of a
 number of hypotheses, applying different jet energy corrections, etc.). 
 To facilitate the development of numerous 
 parallel or subsequent event views, as needed for hypothesis evolution, 
 for instance, this class features a copy constructor, 
 which provides a deep copy of the event container with all data members, 
 physics objects, and their (redirected) relations. 
 This way, the PXL provides a flexible generalized event container  
 meeting the needs of HEP analyses in channels with ambiguous 
 event topologies.
 */
class EventView : public pxl::ObjectManager
{
public:
	EventView() :
		ObjectManager()
	{
	}
	/// This copy constructor provides a deep copy of the event container \p original with all data members, 
	/// physics objects, and their (redirected) relations. 
	EventView(const pxl::EventView& original) :
		ObjectManager(original)
	{
	}
	/// This copy constructor provides a deep copy of the event container \p original with all data members, 
	/// physics objects, and their (redirected) relations.
	EventView(const pxl::EventView* original) :
		ObjectManager(original)
	{
	}

	virtual pxl::WkPtrBase* createSelfWkPtr()
	{
		return new pxl::weak_ptr<EventView>(this);
	}

	virtual const pxl::Id& getTypeId() const
	{
		return getStaticTypeId();
	}

	static const Id& getStaticTypeId()
	{
		static const Id id("c8db3cce-dc4b-421e-882a-83e213c9451f");
		return id;
	}
	
	virtual void serialize(const OutputStream &out) const
	{
		pxl::ObjectManager::serialize(out);
	}

	virtual void deserialize(const InputStream &in)
	{
		pxl::ObjectManager::deserialize(in);
	}

	virtual pxl::Relative* clone() const
	{
		return new pxl::EventView(*this);
	}
	
	std::ostream& print(int level=0, std::ostream& os=std::cout, int pan=1) const;

};

/// This typedef defines a weak pointer for pxl::EventView
typedef pxl::weak_ptr<pxl::EventView> EventViewWkPtr;
} // namespace pxl


#endif // pxl_pol_EventView_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_pol_AnalysisProcess_hh
#define pxl_pol_AnalysisProcess_hh


namespace pxl
{
// pol

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
class AnalysisProcess : public ObjectManager
{
public:
	AnalysisProcess() :
		pxl::ObjectManager()
	{
	}
	AnalysisProcess(const pxl::AnalysisProcess& original) :
		pxl::ObjectManager(original)
	{
	}
	AnalysisProcess(const pxl::AnalysisProcess* original) :
		pxl::ObjectManager(original)
	{
	}
	virtual ~AnalysisProcess()
	{
	}

	inline virtual const pxl::Id& getTypeId() const
	{
		return getStaticTypeId();
	}

	static const Id& getStaticTypeId()
	{
		static const Id id("36c128e0-b14a-4f35-a317-d972d28f1802");
		return id;
	}

	virtual void serialize(const OutputStream &out) const
	{
		pxl::ObjectManager::serialize(out);
	}

	virtual void deserialize(const InputStream &in)
	{
		pxl::ObjectManager::deserialize(in);
	}

	/// This method can be reimplemented to build/destroy 
	/// a static template for user-defined tree creation. \p mode is a freely usable parameter.
	virtual void buildTemplate(int mode = 0)
	{
	}

	/// This method can be reimplemented to hold physics analysis code executed at the begin of a computing job 
	/// (as needed for histogram booking etc.).
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).  
	virtual void beginJob(const pxl::ObjectOwner* input = 0)
	{
	}
	/// This method can be reimplemented to hold physics analysis code executed at the begin of a run. 
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).  
	virtual void beginRun(const pxl::ObjectOwner* input = 0)
	{
	}
	/// This method can be reimplemented to hold physics analysis code executed for the actual event analysis. 
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).  
	virtual void analyseEvent(const pxl::ObjectOwner* input = 0)
	{
	}
	/// This method can be reimplemented to hold physics analysis code executed at the end of each event.
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).  
	/// By default, this method clears the object owner and deletes all owned objects. 
	virtual void finishEvent(const pxl::ObjectOwner* input = 0)
	{
		clearObjects();
	}
	/// This method can be reimplemented to hold physics analysis code executed at the end of a run. 
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).  
	virtual void endRun(const pxl::ObjectOwner* input = 0)
	{
	}
	/// This method can be reimplemented to hold physics analysis code executed at the end of a computing job 
	/// (as needed for histogram storing etc.). 
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).  
	virtual void endJob(const pxl::ObjectOwner* input = 0)
	{
	}

	virtual pxl::Relative* clone() const;

	virtual pxl::WkPtrBase* createSelfWkPtr()
	{
		return new pxl::weak_ptr<AnalysisProcess>(this);
	}

	std::ostream& print(int level=1, std::ostream& os=std::cout, int pan=0) const;
};

static ObjectFactory::ProducerTemplate<pxl::AnalysisProcess>
		_AnalysisProcessProducer(pxl::AnalysisProcess::getStaticTypeId());

/// This typedef defines a weak pointer for pxl::AnalysisProcess
typedef pxl::weak_ptr<pxl::AnalysisProcess> AnalysisProcessWkPtr;

} // namespace pxl


#endif // pxl_pol_AnalysisProcess_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_pol_AnalysisFork_hh
#define pxl_pol_AnalysisFork_hh


namespace pxl
{

// pol

/**
 This class is designed as a base class to assist the analyzer in the 
 parallel evolution of different physics process hypotheses or the analysis of different  
 (instrumental) aspects of an event; data is aggregated in pxl::AnalysisForkData
 */
class AnalysisFork : public pxl::ObjectManager
{
public:
	AnalysisFork() :
		pxl::ObjectManager()
	{
	}
	AnalysisFork(const pxl::AnalysisFork& original) :
		pxl::ObjectManager(original)
	{
	}
	AnalysisFork(const pxl::AnalysisFork* original) :
		pxl::ObjectManager(original)
	{
	}
	virtual ~AnalysisFork()
	{
	}

	inline virtual const pxl::Id& getTypeId() const
	{
		return getStaticTypeId();
	}

	static const Id& getStaticTypeId()
	{
		static const Id id("91b6a6ec-4ecf-490f-ba92-47d20e42bc16");
		return id;
	}

	virtual void serialize(const OutputStream &out) const
	{
		pxl::ObjectManager::serialize(out);
	}

	virtual void deserialize(const InputStream &in)
	{
		pxl::ObjectManager::deserialize(in);
	}

	/// This method can be reimplemented to hold physics analysis code executed at the begin of a computing job.  
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).  
	/// By default, this method invokes the corresponding method of all managed pxl::AnalysisProcess instances.    
	virtual void beginJob(const pxl::ObjectOwner* input = 0);
	/// This method can be reimplemented to hold physics analysis code executed at the begin of a run. 
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).  
	/// By default, this method invokes the corresponding method of all managed pxl::AnalysisProcess instances.    
	virtual void beginRun(const pxl::ObjectOwner* input = 0);
	/// This method can be reimplemented to hold physics analysis code executed for the actual event analysis. 
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).  
	/// By default, this method invokes the corresponding method of all managed pxl::AnalysisProcess instances.    
	virtual void analyseEvent(const pxl::ObjectOwner* input = 0);
	/// This method can be reimplemented to hold physics analysis code executed at the end of each event.
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).  
	/// By default, this method invokes the corresponding method of all managed pxl::AnalysisProcess instances.    
	virtual void finishEvent(const pxl::ObjectOwner* input = 0);
	/// This method can be reimplemented to hold physics analysis code executed at the end of a run. 
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).  
	/// By default, this method invokes the corresponding method of all managed pxl::AnalysisProcess instances.    
	virtual void endRun(const pxl::ObjectOwner* input = 0);
	/// This method can be reimplemented to hold physics analysis code executed at the end of a computing job 
	/// (as needed for histogram storing etc.). 
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).  
	/// By default, this method invokes the corresponding method of all managed pxl::AnalysisProcess instances.    
	virtual void endJob(const pxl::ObjectOwner* input = 0);

	virtual pxl::Relative* clone() const;

	virtual std::ostream& print(int level=1, std::ostream& os=std::cout, int pan=0) const;

	virtual pxl::WkPtrBase* createSelfWkPtr()
	{
		return new pxl::weak_ptr<AnalysisFork>(this);
	}

};

static ObjectFactory::ProducerTemplate<pxl::AnalysisFork>
		_AnalysisForkProducer(pxl::AnalysisFork::getStaticTypeId());
//
/// This typedef defines a weak pointer for pxl::AnalysisFork
typedef pxl::weak_ptr<pxl::AnalysisFork> AnalysisForkWkPtr;

} // namespace pxl

#endif // pxl_pol_AnalysisFork_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_pol_ParticleFilter_hh
#define pxl_pol_ParticleFilter_hh




namespace pxl
{

/** 
 This class provides a pT-sorted filter for PXL physics objects (requires the name, pT and |eta| to match).
 */
class ParticlePtComparator : public Comparator<pxl::Particle>
{
public:
	virtual bool operator()(const pxl::Particle& p1, const pxl::Particle& p2)
	{
		return (p1.getPt()>p2.getPt());
	}
	virtual bool operator()(const pxl::Particle* p1, const pxl::Particle* p2)
	{
		return (p1->getPt()>p2->getPt());
	}
	
};
typedef ParticlePtComparator PtComparator;

class ParticlePtEtaNameCriterion : public FilterCriterion<pxl::Particle>
{
public:
	ParticlePtEtaNameCriterion(const std::string& name, double ptMin = 0.0,
			double etaMax = 0.0) :
		_name(name), _ptMin(ptMin), _etaMax(etaMax)
	{
	}

	virtual bool operator()(const pxl::Particle& pa) const
	{
		if ((_name != "" && pa.getName() != _name) || (_ptMin > 0.0
				&& pa.getPt() < _ptMin) || (_etaMax > 0.0
				&& std::fabs(pa.getEta()) > _etaMax))
			return false;
		return true;
	}

private:
	std::string _name;
	double _ptMin;
	double _etaMax;
};

class ParticlePtCriterion : public FilterCriterion<pxl::Particle>
{
public:
	ParticlePtCriterion(double ptMin = 0.0) :
		_ptMin(ptMin)
	{
	}

	virtual bool operator()(const pxl::Particle& pa) const
	{
		if ((_ptMin > 0.0 && pa.getPt() < _ptMin))
			return false;
		return true;
	}

private:
	double _ptMin;
};

class ParticleNameCriterion : public FilterCriterion<pxl::Particle>
{
public:
	ParticleNameCriterion(const std::string& name) :
		_name(name)
	{
	}

	virtual bool operator()(const pxl::Particle& pa) const
	{
		if (_name != "" && pa.getName() != _name) return false;
		return true;
	}

private:
	std::string _name;
};


typedef pxl::Filter<pxl::Particle, ParticlePtComparator> ParticleFilter;

} // namespace pxl

#endif // pxl_pol_ParticleFilter_hh

#endif // pxl_pol_hh

#endif // pxl_hh
