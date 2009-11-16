//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef pxl_hh
#define pxl_hh

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_HH
#define PXL_BASE_HH

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_BASIC3VECTOR_HH
#define PXL_BASE_BASIC3VECTOR_HH

#include <cmath>
#include <cfloat>
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_IO_STREAM_HH
#define PXL_IO_STREAM_HH

#include <string>
#include <iostream>

#include <string.h>

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_MACROS_HH
#define PXL_BASE_MACROS_HH

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

#ifdef _MSC_VER
	#ifdef PXL_EXPORT
		#define PXL_DLL_EXPORT __declspec( dllexport )
	#else
		#define PXL_DLL_EXPORT __declspec( dllimport )
	#endif // PXL_EXPORT
#else
	#define PXL_DLL_EXPORT 
#endif // _MSC_VER

#ifndef PXL_VERSION
	#define PXL_VERSION "trunk"
#endif

#endif // PXL_BASE_MACROS_HH

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
class PXL_DLL_EXPORT OutputStream
{

public:
	virtual ~OutputStream()
	{
	}
	
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

class PXL_DLL_EXPORT BufferOutput : public OutputStream
{
	mutable char *_Buffer;
	mutable size_t _BufferSize;
	mutable size_t _BufferPos;

public:

	BufferOutput() :
		_Buffer(0), _BufferSize(0), _BufferPos(0)
	{
		// optimal size to be determined
		resize(1048576);
	}
	
	BufferOutput(size_t bufferSize) :
		_Buffer(0), _BufferSize(0), _BufferPos(0)
	{
		resize(bufferSize);
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
	
	void clear()
	{
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
			size_t size = _BufferSize*2;
			while (size < needed)
				size *= 2;
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
	
private:
	BufferOutput(const BufferOutput& original)
	{
	}
	
	BufferOutput& operator= (const BufferOutput& other)
	{
		return *this;
	}
};

// iotl
/**
 This abstract class serves the internal PXL I/O scheme by defining how basic C++ types
 are read from an input buffer.
 */
class PXL_DLL_EXPORT InputStream
{

public:
	virtual ~InputStream()
	{
	}
	
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
class PXL_DLL_EXPORT BufferInput : public InputStream
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
		if (_BufferSize > 0)
			memcpy(newBuffer, _Buffer, _BufferSize > size ? size : _BufferSize);
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
	
private:
	BufferInput(const BufferInput& original)
	{
	}
	
	BufferInput& operator= (const BufferInput& other)
	{
		return *this;
	}
};

}

#endif /*PXL_IO_STREAM_HH*/

namespace pxl
{

/**
 This class provides a simple threevector with basic algebra. The methods provided are self-explanatory.
 Note that theta is defined as 0 in z-direction, as usual in collider physics.
 */
class PXL_DLL_EXPORT Basic3Vector
{
public:
	Basic3Vector()
	{
		_v[0] = 0;
		_v[1] = 0;
		_v[2] = 0;
	}
	Basic3Vector(const Basic3Vector& orig)
	{
		_v[0] = orig._v[0];
		_v[1] = orig._v[1];
		_v[2] = orig._v[2];
	}
	explicit Basic3Vector(const Basic3Vector* orig)
	{
		_v[0] = orig->_v[0];
		_v[1] = orig->_v[1];
		_v[2] = orig->_v[2];
	}
	Basic3Vector(double x, double y, double z)
	{
		_v[0] = x;
		_v[1] = y;
		_v[2] = z;
	}

	Basic3Vector(const double x[3])
	{
		_v[0] = x[0];
		_v[1] = x[1];
		_v[2] = x[2];
	}

	virtual ~Basic3Vector()
	{
	}

	virtual void serialize(const OutputStream &out) const
	{
		out.writeDouble(_v[0]);
		out.writeDouble(_v[1]);
		out.writeDouble(_v[2]);
	}

	virtual void deserialize(const InputStream &in)
	{
		in.readDouble(_v[0]);
		in.readDouble(_v[1]);
		in.readDouble(_v[2]);
	}

	inline void setX(double x)
	{
		_v[0] = x;
	}
	inline void setY(double y)
	{
		_v[1] = y;
	}
	inline void setZ(double z)
	{
		_v[2] = z;
	}

	inline void setComponent(int i, double val)
	{
		_v[i] = val;
	}

	inline double getComponent(int i) const
	{
		return _v[i];
	}
	
	inline double getElement(int i) const
	{
		return _v[i];
	}
	
	inline double* getArray()
	/// returns a pointer to the C Array double[3] for the three
	/// components of the vector
	{
		return _v;
	}
	
	const double* getConstArray() const
	/// returns a pointer to the C Array double[3] for the three
	/// components of the vector
	{
		return _v;
	}
	
	void setCArray(double val[3])
	{
		_v[0] = val[0];
		_v[1] = val[1];
		_v[2] = val[2];
	}

	void setRhoPhi(double perp, double phi);

	void setRhoPhiZ(double perp, double phi, double z);

	void setRThetaPhi(double r, double theta, double phi);

	inline double getX() const
	{
		return _v[0];
	}
	inline double getY() const
	{
		return _v[1];
	}
	inline double getZ() const
	{
		return _v[2];
	}

	bool isNullPerp() const;

	bool isNull() const;

	double getPerp2() const;

	double getPerp() const;

	double getPhi() const;

	double getMag2() const;

	double getMag() const;

	double getCosTheta() const;

	double getCos2Theta() const;

	double getTheta() const;

	///Returns unit vector in spherical coordinates
	Basic3Vector getETheta() const;

	///Returns unit vector in spherical coordinates
	Basic3Vector getEPhi() const;

	double deltaRho(const Basic3Vector* fv) const;

	double deltaPhi(const Basic3Vector* fv) const;

	double deltaTheta(const Basic3Vector* fv) const;

	const Basic3Vector& operator=(const Basic3Vector& vec);

	const Basic3Vector& operator+=(const Basic3Vector& vec);

	const Basic3Vector& operator-=(const Basic3Vector& vec);

	const Basic3Vector& operator*=(double skalar);

	const Basic3Vector& operator/=(double skalar);

	Basic3Vector operator+(const Basic3Vector& vec);
	Basic3Vector operator-(const Basic3Vector& vec);

	Basic3Vector operator/(double skalar) const;

	// Scalar product
	double operator*( const Basic3Vector& vec) const;

	bool isUnityVector() const;
	
	void normalize();

private:
	double _v[3];

};

// non-member operators
PXL_DLL_EXPORT bool const operator==(const Basic3Vector& obj1,
		const Basic3Vector& obj2);
PXL_DLL_EXPORT bool const operator!=(const Basic3Vector& obj1,
		const Basic3Vector& obj2);

PXL_DLL_EXPORT Basic3Vector operator*(const double skalar,
		const Basic3Vector& vec);
PXL_DLL_EXPORT Basic3Vector operator*(const Basic3Vector& vec,
		const double skalar);

} // namespace pxl

#endif // PXL_BASE_BASIC3VECTOR_HH
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_BASICNVECTOR_HH
#define PXL_BASE_BASICNVECTOR_HH

#include <stdexcept>
#include <vector>

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_IO_SERIALIZABLE_HH
#define PXL_IO_SERIALIZABLE_HH

#include <sstream>
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_ID_HH
#define PXL_BASE_ID_HH

// Random.h
// Mersenne Twister random number generator -- a C++ class Random
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

// Parts of this file are modified beginning in 29.10.09 for adaption in
// PXL.


#ifndef PXL_BASE_RANDOM_HH
#define PXL_BASE_RANDOM_HH

// Not thread safe (unless auto-initialization is avoided and each thread has
// its own Random object)
#include <limits.h>
#include <stdio.h>
#include <time.h>

//necessary for win32
#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif


namespace pxl{

class PXL_DLL_EXPORT Random{
public:
	typedef unsigned long uint32;  // unsigned integer type, at least 32 bits
	enum { N = 624 };              // length of state vector
	enum { SAVE = N + 1 };         // length of array for save()

protected:
	enum { M = 397 };  // period parameter
	uint32 state[N];   // internal state
	uint32 *pNext;     // next value to get from state
	int left;          // number of values left before reload needed


//Methods
public:
	/// initialize with a simple uint32
	Random( const uint32& oneSeed );
	// initialize with an array
	Random( uint32 *const bigSeed, uint32 const seedLength = N );
	/// auto-initialize with /dev/urandom or time() and clock()
	/// Do NOT use for CRYPTOGRAPHY without securely hashing several returned
	/// values together, otherwise the generator state can be learned after
	/// reading 624 consecutive values.
	Random();
	// Access to 32-bit random numbers
	double rand();                          ///< real number in [0,1]
	double rand( const double& n );         ///< real number in [0,n]
	double randExc();                       ///< real number in [0,1)
	double randExc( const double& n );      ///< real number in [0,n)
	double randDblExc();                    ///< real number in (0,1)
	double randDblExc( const double& n );   ///< real number in (0,n)
	uint32 randInt();                       ///< integer in [0,2^32-1]
	uint32 randInt( const uint32& n );      ///< integer in [0,n] for n < 2^32
	double operator()() { return rand(); }  ///< same as rand()

	/// Access to 53-bit random numbers (capacity of IEEE double precision)
	double rand53();  // real number in [0,1)
	///Exponential distribution in (0,inf)
	double randExponential();
	/// Normal distruibuted random number
	double randNorm( const double& mean = 0.0, const double& variance = 0.0 );
	/// Uniform distribution in [min, max]
	double randUniform(double min, double max);
	/// Rayleigh distributed random number	
	double randRayleigh(double sigma);
	/// Power-Law distribution, not possible for index == -1
	double randPowerLaw(double index, double min, double max);
	/// Broken power-law distribution 
	double randBrokenPowerLaw(double index1, double index2, double breakpoint, double min, double max );
	/// Random point on a unit-sphere
	Basic3Vector randUnitVectorOnSphere();
	// Re-seeding functions with same behavior as initializers
	void seed( const uint32 oneSeed );
	void seed( uint32 *const bigSeed, const uint32 seedLength = N );
	void seed();

	// Saving and loading generator state
	void save( uint32* saveArray ) const;  // to array of size SAVE
	void load( uint32 *const loadArray );  // from such array
	friend std::ostream& operator<<( std::ostream& os, const Random& mtrand );
	friend std::istream& operator>>( std::istream& is, Random& mtrand );

protected:
	void initialize( const uint32 oneSeed );
	void reload();
	uint32 hiBit( const uint32& u ) const { return u & 0x80000000UL; }
	uint32 loBit( const uint32& u ) const { return u & 0x00000001UL; }
	uint32 loBits( const uint32& u ) const { return u & 0x7fffffffUL; }
	uint32 mixBits( const uint32& u, const uint32& v ) const
		{ return hiBit(u) | loBits(v); }

#ifdef _MSC_VER
	#pragma warning( push )
	#pragma warning( disable : 4146 )
#endif
	uint32 twist( const uint32& m, const uint32& s0, const uint32& s1 ) const
		{ return m ^ (mixBits(s0,s1)>>1) ^ (-loBit(s1) & 0x9908b0dfUL); }

#ifdef _MSC_VER
	#pragma warning( pop )
#endif

	static uint32 hash( time_t t, clock_t c );
};

inline Random::Random( const uint32& oneSeed )
	{ seed(oneSeed); }

inline Random::Random( uint32 *const bigSeed, const uint32 seedLength )
	{ seed(bigSeed,seedLength); }

inline Random::Random()
	{ seed(); }

inline double Random::rand()
	{ return double(randInt()) * (1.0/4294967295.0); }

inline double Random::rand( const double& n )
	{ return rand() * n; }

inline double Random::randExc()
	{ return double(randInt()) * (1.0/4294967296.0); }

inline double Random::randExc( const double& n )
	{ return randExc() * n; }

inline double Random::randDblExc()
	{ return ( double(randInt()) + 0.5 ) * (1.0/4294967296.0); }

inline double Random::randDblExc( const double& n )
	{ return randDblExc() * n; }

inline double Random::rand53()
{
	uint32 a = randInt() >> 5, b = randInt() >> 6;
	return ( a * 67108864.0 + b ) * (1.0/9007199254740992.0);  // by Isaku Wada
}

inline double Random::randNorm( const double& mean, const double& variance )
{
	// Return a real number from a normal (Gaussian) distribution with given
	// mean and variance by Box-Muller method
	double r = sqrt( -2.0 * log( 1.0-randDblExc()) ) * variance;
	double phi = 2.0 * 3.14159265358979323846264338328 * randExc();
	return mean + r * cos(phi);
}

inline double Random::randUniform(double min, double max)
{
	return min + (max-min) * this->rand();
}

	///returns a rayleigh distributed value
inline double Random::randRayleigh(double sigma) 
{
	return sigma * sqrt(-2.0 * log(1-this->rand()));
}

inline Basic3Vector Random::randUnitVectorOnSphere()
{
	double z = this->randUniform(-1.0,1.0);
	double t = this->randUniform(-1.0*M_PI,M_PI);
	double r  = sqrt(1-z*z);
	return Basic3Vector(r*cos(t),r*sin(t),z);
}

inline double Random::randPowerLaw(double index, double min, double max)
{
	if ((min < 0) || (max < min))
	{
		throw std::runtime_error("Power law distribution only possible for 0 <= min <= max");
	}
	//check for index -1!
	double myrand = this->rand();
	if ((std::abs(index +1.0) )< DBL_EPSILON)
	{
		throw std::runtime_error("Power law distribution only for index !=-1");
	}
	double part1 = pow(max, index + 1);
	double part2 = pow(min, index + 1);
	double ex = 1 / (index + 1);
	return pow((part1 - part2) * myrand + part2, ex);
}

inline double Random::randExponential()
{
	double dum;
	do
	{
		dum = this->rand();
	} while (dum < DBL_EPSILON);
	return -1.0 * log(dum);
}

inline double Random::randBrokenPowerLaw(double index1, double index2, double breakpoint, double min, double max )
{
	if (min >= breakpoint)
	{
		return this->randPowerLaw(index2, min, max);
	}
	else if (max <= breakpoint)
	{
		return this->randPowerLaw(index2, min, max);
	}
	else
	{
		double decide = this->rand();
		double intPL1 = (pow(breakpoint, index1 + 1) - pow(min, index1 + 1))
				/ (index1 + 1);
		double intPL2 = pow(breakpoint, index1 - index2)
				* (pow(max, index2 + 1) - pow(breakpoint, index2 + 1))
				/ (index2 + 1);
		if (decide > intPL1 / (intPL1 + intPL2))
			return this->randPowerLaw(index2, breakpoint, max);
		else
			return this->randPowerLaw(index1, min, breakpoint);
	}

}

inline Random::uint32 Random::randInt()
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

inline Random::uint32 Random::randInt( const uint32& n )
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


inline void Random::seed( const uint32 oneSeed )
{
	// Seed the generator with a simple uint32
	initialize(oneSeed);
	reload();
}


inline void Random::seed( uint32 *const bigSeed, const uint32 seedLength )
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


inline void Random::seed()
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
			success = fread( s++, sizeof(uint32), 1, urandom ) != 0;
		fclose(urandom);
		if( success ) { seed( bigSeed, N );  return; }
	}

	// Was not successful, so use time() and clock() instead
	seed( hash( time(NULL), clock() ) );
}


inline void Random::initialize( const uint32 seed )
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


inline void Random::reload()
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


inline Random::uint32 Random::hash( time_t t, clock_t c )
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


inline void Random::save( uint32* saveArray ) const
{
	register uint32 *sa = saveArray;
	register const uint32 *s = state;
	register int i = N;
	for( ; i--; *sa++ = *s++ ) {}
	*sa = left;
}


inline void Random::load( uint32 *const loadArray )
{
	register uint32 *s = state;
	register uint32 *la = loadArray;
	register int i = N;
	for( ; i--; *s++ = *la++ ) {}
	left = *la;
	pNext = &state[N-left];
}


inline std::ostream& operator<<( std::ostream& os, const Random& mtrand )
{
	register const Random::uint32 *s = mtrand.state;
	register int i = mtrand.N;
	for( ; i--; os << *s++ << "\t" ) {}
	return os << mtrand.left;
}


inline std::istream& operator>>( std::istream& is, Random& mtrand )
{
	register Random::uint32 *s = mtrand.state;
	register int i = mtrand.N;
	for( ; i--; is >> *s++ ) {}
	is >> mtrand.left;
	mtrand.pNext = &mtrand.state[mtrand.N-mtrand.left];
	return is;
}

}
#endif  // PXL_BASE_RANDOM_HH 




namespace pxl {

/// This typedef is intended to provide a data type for the PXL unique object-id
class PXL_DLL_EXPORT Id
{
	unsigned char bytes[16];

public:


	Id()
	{
		generate();
	}

	Id (const InputStream& in)
	{
		deserialize (in);
	}

	Id(const char* id);

	explicit Id(const std::string& id);

	Id(Random& rand)
	{
		generate(rand);
	}

	void generate(Random& rand)
	{
		for (int i = 0; i < 4; i++)
		{
			Random::uint32 value = rand.randInt();
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
		static Random rand;
		generate(rand);
	}

	static Id create()
	{
		Id id;
		id.generate();
		return id;
	}

	bool operator ==(const Id& id) const
	{
		for (int i = 0; i < 16; i++)
		{
			if (bytes[i] != id.bytes[i])
				return false;
		}

		return true;
	}

	bool operator !=(const Id& id) const
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

	bool operator <(const Id& op) const
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
			const Id &uid);

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

inline std::ostream& operator <<(std::ostream& os, const Id &id)
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

#endif // PXL_BASE_ID_HH


namespace pxl
{

// io
/**
 This class is the abstract base class for all objects to be stored in a PXL I/O file.
 It holds the unique ID of each individual object. In addition, a UUID indicating the class
 type must be implemented in derived classes. (De-) Serialization happens via consecutive calls to the
 serialize/deserialize methods of the base classes.
 */
class PXL_DLL_EXPORT Serializable
{
public:

	Serializable()
	{
	}

	///Copy constructor. A copied object gets a new unique ID.
	Serializable(const Serializable& original) : _id()
	{
	}

	//Nothing needs to be done, a UUID is assigned already
	Serializable& operator=(const Serializable& original)
	{
		if (this != &original)
			_id = Id::create();
		return *this;
	}

	virtual ~Serializable()
	{
	}

	/// Id of the object's type.
	virtual const Id& getTypeId() const = 0;

	/// Returns the unique ID of the individual object.
	const Id& getId() const
	{
		return _id;
	}

	virtual Serializable* clone() const = 0;

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

	virtual std::ostream& print(int level=1, std::ostream& os=std::cout, int pan=1) const
	{
		 os << "Serializable [" << getId() <<"]"<< std::endl;
		 return os;
	}
	virtual const std::string toString() const
	{
		std::ostringstream ss;
		this->print(0,ss);
		return ss.str();
	}
private:
	/// The unique ID.
	Id _id;
};

}

#endif /*PXL_IO_SERIALIZABLE_HH*/
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_IO_OBJECT_FACTORY_HH
#define PXL_IO_OBJECT_FACTORY_HH

#include <map>


// io
/**
 This class serves the PXL I/O scheme by managing
 the relation of classes to UUIDs.
 */

namespace pxl {

class ObjectProducerInterface;

class PXL_DLL_EXPORT ObjectFactory
{
private:
	
	ObjectFactory();
	
	std::map<Id,const ObjectProducerInterface *> _Producers;

public:

	static ObjectFactory& instance();
		
	static Serializable *create (const Id& id);
	
	static void registerClass (const Id& id, const ObjectProducerInterface* producer);
};

class ObjectProducerInterface
{
public:
	virtual ~ObjectProducerInterface()
	{
	}
	
	virtual Serializable *create () const = 0;
};

template <class T>
class ObjectProducerTemplate : public ObjectProducerInterface
{
public:

	ObjectProducerTemplate (const Id& id)
	{
		ObjectFactory::registerClass (id, this);
	}
	
	virtual Serializable *create () const
	{
		return new T();
	}
};

} // namespace pxl

#endif // PXL_IO_OBJECT_FACTORY_HH

namespace pxl
{

/**
 This class provides a vector of arbitrary length for data storage or algebra in R**N
 */
class PXL_DLL_EXPORT BasicNVector : public Serializable
{
public:

	BasicNVector() : Serializable(), _data(NULL), _size1(0), _alienarray(false)
	{
	}

	BasicNVector(size_t size) : Serializable(), _data(NULL), _size1(size), _alienarray(false)
	{
		_data = new double[size];
		std::fill_n(_data, size, 0 );
	}

	BasicNVector(const BasicNVector& orig): Serializable(orig)
	{
		_size1 = orig.getSize();
		_data = new double[_size1];
		_alienarray = orig._alienarray;
		for (size_t i=0;i<_size1;i++)
		{
			_data[i] = orig.getElement(i);
		}
	}

	BasicNVector(const Basic3Vector& vec) : Serializable()
	{
		_size1 = 3;
		_data = new double[_size1];
		_alienarray = false;
		_data[0] = vec.getX();
		_data[1] = vec.getY();
		_data[2] = vec.getZ();
	}

	explicit BasicNVector(const BasicNVector* orig): Serializable(*orig)
	{
		_size1 = orig->getSize();
		_data = new double[_size1];
		_alienarray = orig->_alienarray;
		for (size_t i=0;i<_size1;i++)
		{
			_data[i] = orig->getElement(i);
		}
	}

	BasicNVector(size_t size, double *data): Serializable()

	{
		_size1 = size;
		_data = data;
		_alienarray = true;
	}

	virtual ~BasicNVector()
	{
		if (_data && (!_alienarray))
		delete[] _data;
	}

	/// Uses the given array as data array without copying it
	/// If the Vector is deleted, the array is also deleted!
	void use(size_t size, double *data);

	/// resets the size (and allocates memory if not yet done)
	void setSize(size_t size)
	{
		if (!_data)
		{
			_data = new double[size];
			std::fill_n(_data, size, 0 );
		}
		_size1 = size;

	}

	size_t getSize() const
	{
		return _size1;
	}

	double getElement(size_t i) const
	{
		if (i<_size1)
		return _data[i];
		else
		throw std::runtime_error("Index out of range");
	}

	void setElement(size_t i,double value)
	{
		if (i<_size1)
		_data[i] = value;
		else
		throw std::runtime_error("Index out of range");
	}

	virtual std::ostream& print(int level=1, std::ostream& os=std::cout, int pan=1) const;

	double* getArray()
	{
		return _data;
	}

	const double* getConstArray() const
	{
		return _data;
	}

	virtual Serializable* clone() const
	{
		return new BasicNVector(*this);
	}

	static const Id& getStaticTypeId()
	{
		static const Id id("b95fb401-e0a7-87b6-1737-9ce54aae6f29");
		return id;
	}

	virtual const Id& getTypeId() const
	{
		return getStaticTypeId();
	}

	virtual void serialize(const OutputStream &out) const;

	virtual void deserialize(const InputStream &in);

	const BasicNVector& operator=(const BasicNVector& vec);
	BasicNVector operator+(const BasicNVector& vec);
	BasicNVector operator-(const BasicNVector& vec);
	const BasicNVector& operator+=(const BasicNVector& vec);
	const BasicNVector& operator-=(const BasicNVector& vec);
	const BasicNVector& operator*=(double skalar);
	const BasicNVector& operator/=(double skalar);
	const BasicNVector operator/(double skalar) const;

	double& operator()(size_t i);
	//// Scalar product
	double operator*( const BasicNVector& vec) const;

private:
	double* _data;
	size_t _size1;
	bool _alienarray;

};
// non-member operators
PXL_DLL_EXPORT bool const operator==(const BasicNVector& obj1,
		const BasicNVector& obj2);
PXL_DLL_EXPORT bool const operator!=(const BasicNVector& obj1,
		const BasicNVector& obj2);

PXL_DLL_EXPORT BasicNVector operator*(double skalar, const BasicNVector& vec);
PXL_DLL_EXPORT BasicNVector operator*(const BasicNVector& vec, double skalar);

} // namespace pxl

#endif // PXL_BASE_BASICNVECTOR_HH
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_BASICMATRIX_HH
#define PXL_BASE_BASICMATRIX_HH


namespace pxl
{

/**
 * Storage Order
 */
typedef enum
{
	ROWMAJOR, COLUMNMAJOR
} StorageOrder;


/**
 This class provides a Matrix of arbitrary dimensions for data storage
 or algebra in R**NxM

 The Matrix is stored as BLAS compatible array of doubles. The default storage order is the C default of Row-Major
 order, meaning, that a Matrix is stored row by row
 double M[] = { 1,2,3,
 4,5,6
 }
 The Same matrix in Column Order would be stored as
 double M[] = {1,4,2,5,3,6};

 size1, usually defined by index i, denotes the row in row-based order,
 respectivelty the column in column-based order. In row-based order
 element i,j is hence element i*size2+j of the array, and in column-based
 order element i,j is also element i*size2+j of the array.
 **/
class PXL_DLL_EXPORT BasicMatrix : public Serializable
{
public:
	BasicMatrix() : Serializable() ,_size1(0), _size2(0), _data(NULL), _storageType(ROWMAJOR), _alienarray(false)
{
	}
	BasicMatrix(size_t size1, size_t size2) : Serializable(), _size1(size1), _size2(size2), _storageType(ROWMAJOR), _alienarray(false)
	{
		_data = new double[size1*size2];
		std::fill_n(_data, size1*size2, 0 );
	}

	BasicMatrix(size_t size1, size_t size2, StorageOrder storage) : Serializable(), _size1(size1), _size2(size2), _storageType(storage), _alienarray(false)
	{
		_data = new double[size1*size2];
		std::fill_n(_data, size1*size2, 0 );
	}

	BasicMatrix(const BasicMatrix& orig) : Serializable(orig)
	{
		_size1 = orig._size1;
		_size2 = orig._size2;
		_data = new double[_size1*_size2];
		_alienarray = orig._alienarray;
		_storageType = orig._storageType;
		for (size_t i=0;i<_size1*_size2;i++)
		{
				_data[i] = orig.getConstArray()[i];
		}
	}

	explicit BasicMatrix(const BasicMatrix* orig) : Serializable(*orig)
	{
		_size1 = orig->_size1;
		_size2 = orig->_size2;
		_data = new double[_size1*_size2];
		_alienarray = orig->_alienarray;
		_storageType = orig->_storageType;
		for (size_t i=0;i<_size1*_size2;i++)
		{
				_data[i] = orig->getConstArray()[i];
		}
	}

	BasicMatrix(size_t size1, size_t size2, double *data, StorageOrder storage) : Serializable(), _size1(size1), _size2(size2), _data(data), _storageType(storage), _alienarray(true)
	{
	}

	virtual ~BasicMatrix()
	{
		if ((_data) && (!_alienarray))
		delete[] _data;
	}

	void use(size_t size1, size_t size2, double *data);

	void setRowBasedStorage()
	{
		_storageType = ROWMAJOR;
	}

	void setColumnBasedStorage()
	{
		_storageType = COLUMNMAJOR;
	}

	bool isRowBasedStorage() const;
	bool isColumnBasedStorage() const;

	size_t getSize1() const
	{
		return _size1;
	}

	size_t getSize2() const
	{
		return _size2;
	}

	size_t getNumberOfRows() const
	{
		if (_storageType==ROWMAJOR)
			return _size1;
		else
			return _size2;
	}

	size_t getNumberOfColumns() const
	{
		if (_storageType==ROWMAJOR)
			return _size2;
		else
			return _size1;
	}

	// resizes the array, data loss!
	void resize(size_t i, size_t j);

	// only reshapes the array and preserves the data
	void reshape(size_t i, size_t j);

	double getElement(size_t i,size_t j) const
	{
		if (i >= _size1)
		{
			if (this->isRowBasedStorage())
			throw std::runtime_error("Index exceeds number of rows");
			else
			throw std::runtime_error("Index exceeds number of columns");
		}
		if (j >= _size2)
		{
			if (this->isColumnBasedStorage())
			throw std::runtime_error("Index exceeds number of rows");
			else
			throw std::runtime_error("Index exceeds number of columns");
		}

		return _data[i*_size2+j];
	}

	void setElement(size_t i,size_t j, double val)
	{
		if (i >= _size1)
		{
			if (this->isRowBasedStorage())
			throw std::runtime_error("Index exceeds number of rows");
			else
			throw std::runtime_error("Index exceeds number of columns");
		}
		if (j >= _size2)
		{
			if (this->isColumnBasedStorage())
			throw std::runtime_error("Index exceeds number of rows");
			else
			throw std::runtime_error("Index exceeds number of columns");
		}

		_data[i*_size2+j] = val;
	}

	double* getArray()
	{
		return _data;
	}

	const double* getConstArray() const
	{
		return _data;
	}

	const BasicMatrix& operator=(const BasicMatrix& M);
	BasicMatrix operator+(const BasicMatrix& M);
	BasicMatrix operator-(const BasicMatrix& M);
	const BasicMatrix& operator+=(const BasicMatrix& M);
	const BasicMatrix& operator-=(const BasicMatrix& M);
	const BasicMatrix& operator*=(double skalar);
	const BasicMatrix& operator/=(double skalar);
	const BasicMatrix operator/(double skalar) const;

	double& operator()(size_t i, size_t j);

	virtual std::ostream& print(int level=1, std::ostream& os=std::cout, int pan=1) const;

	virtual Serializable* clone() const
	{
		return new BasicMatrix(*this);
	}

	virtual void serialize(const OutputStream &out) const;

	virtual void deserialize(const InputStream &in);

	static const Id& getStaticTypeId()
	{
		static const Id id("85ebd9cb-08ae-6c7b-2257-45604aae6f55");
		return id;
	}

	virtual const Id& getTypeId() const
	{
		return getStaticTypeId();
	}
private:
	size_t _size1;
	size_t _size2;
	double* _data;
	StorageOrder _storageType;
	bool _alienarray; // true if working with an not self allocated array - if so, memory is not freed!


};
PXL_DLL_EXPORT bool const operator==(const BasicMatrix& obj1, const BasicMatrix& obj2);
PXL_DLL_EXPORT bool const operator!=(const BasicMatrix& obj1, const BasicMatrix& obj2);

PXL_DLL_EXPORT BasicMatrix operator*(double skalar, const BasicMatrix& vec);
PXL_DLL_EXPORT BasicMatrix operator*(const BasicMatrix& vec, double skalar);

// Fundamental Matrix Vector Products

PXL_DLL_EXPORT BasicNVector operator*(const BasicMatrix& M, const BasicNVector& vec);
PXL_DLL_EXPORT Basic3Vector operator*(const BasicMatrix& M, const Basic3Vector& vec);



} // namespace pxl

#endif // PXL_BASE_BASICMATRIX_HH
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_LORENTZVECTOR_HH
#define PXL_BASE_LORENTZVECTOR_HH



namespace pxl
{
// pol

/**
 This class provides a simple Lorentz-fourvector with basic algebra. The methods provided are self-explanatory.
 */
class PXL_DLL_EXPORT LorentzVector : public Basic3Vector
{
public:
	LorentzVector() : Basic3Vector(0, 0, 0), _t(0)
	{
	}

	LorentzVector(const Basic3Vector &orig, double t = 0) :
	Basic3Vector(orig), _t(t)
	{
	}
	LorentzVector(const LorentzVector &orig) :
	Basic3Vector(orig), _t(orig._t)
	{

	}
	explicit LorentzVector(const LorentzVector* orig) :
	Basic3Vector(orig), _t(orig->_t)
	{
	}

	LorentzVector(double x, double y, double z, double t = 0) :
	Basic3Vector(x, y, z), _t(t)
	{
	}

	virtual ~LorentzVector()
	{
	}

	virtual void serialize(const OutputStream &out) const
	{
		Basic3Vector::serialize(out);
		out.writeDouble(_t);
	}

	virtual void deserialize(const InputStream &in)
	{
		Basic3Vector::deserialize(in);
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
	
	double getP() const
	{
		return getMag();
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

	inline double deltaR(const LorentzVector* fv) const
	{
		double dDeta = deltaEta(fv);
		double dDphi = deltaPhi(fv);
		return std::sqrt(dDeta*dDeta + dDphi*dDphi);
	}

	inline double deltaEta(const LorentzVector* fv) const
	{
		return getEta() - fv->getEta();
	}

	inline const LorentzVector& operator=(const Basic3Vector& vec)
	{
		Basic3Vector::operator=(vec); return *this;
	}

	inline const LorentzVector& operator=(const LorentzVector& vec)
	{
		Basic3Vector::operator=(vec); _t = vec._t; return *this;
	}

	inline const LorentzVector& operator+=(const LorentzVector& vec)
	{
		Basic3Vector::operator+=(vec); _t += vec._t; return *this;
	}

	inline const LorentzVector& operator-=(const LorentzVector& vec)
	{
		Basic3Vector::operator-=(vec); _t -= vec._t; return *this;
	}

private:
	double _t;

};

// non-member operators
PXL_DLL_EXPORT bool const operator==(const LorentzVector& obj1,
		const LorentzVector& obj2);
PXL_DLL_EXPORT bool const operator!=(const LorentzVector& obj1,
		const LorentzVector& obj2);

} // namespace pxl


#endif // PXL_BASE_LORENTZVECTOR_HH
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_BASICCONTAINER_HH
#define PXL_BASE_BASICCONTAINER_HH


//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_USERRECORD_HH
#define PXL_BASE_USERRECORD_HH


//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_VARIANT_HH
#define PXL_BASE_VARIANT_HH

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_VARIANTBASE_HH
#define PXL_BASE_VARIANTBASE_HH

#include <cstdlib>
#include <memory>



//#include <boost/shared_ptr.hpp>
//#include "pxl/io/Serializable.hh"

namespace pxl
{

/**
 This base class serves the PXL variant data type.
 It does not contain the information about the type store, this task
 is delegated to the caller. This way memory can be saved, if, for instance,
 table-like structures are to be stored, where the type information is stored
 in a separate header structure.
 */
// Important notice: pxl::VariantBase is not a fully opaque type, you have to
// take care of calling clear() and dup() after copying or freeing yourself!
class PXL_DLL_EXPORT VariantBase
{
public:
	/// This enum represents the possible value types that can be stored by the PXL variant data type.
	enum Type
	{
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
//		TYPE_SERIALIZABLEPTR, // The Serializable Pointer
		TYPE_PTR,
		TYPE_USER
	};

	inline VariantBase()
	{	memset(&v, 0, sizeof v);}
	inline VariantBase(const VariantBase& orig) : v(orig.v)
	{}
	inline VariantBase(const VariantBase* orig) : v(orig->v)
	{}

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
		{
			delete[] (char*)v.p;
		}
		//else if (t == TYPE_SERIALIZABLEPTR)
		//{
		//	delete (boost::shared_ptr<Serializable>*)v.p;
		//}
		else if (PXL_UNLIKELY(t >= TYPE_USER))
		{
			const TypeInfo& tInfo = getTypeInfo(t);
			if (tInfo.clear)
			tInfo.clear(&v);
		}
		memset(&v, 0, sizeof v);
	}

	/// This method ensures correct duplication of possible internal memory allocations after a copy of the variant, depending on the type given by \p t.
	inline void dup(Type t)
	{
		if (t == TYPE_STRING)
		{
			std::size_t len = strlen((char*)v.p) + 1;
			char* copy = new char[len];
			memcpy(copy, v.p, len);
			v.p = copy;
		}
		//else if (t == TYPE_SERIALIZABLEPTR)
		//{
		//	boost::shared_ptr<Serializable>* copy(new boost::shared_ptr<Serializable>(* static_cast<boost::shared_ptr<Serializable>*>(v.p)));
		//	v.p = copy;

		//}
		else if (PXL_UNLIKELY(t >= TYPE_USER))
		{
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
	union Value
	{
		bool b;
		char c;
		unsigned char uc;
		short s;
		unsigned short us;
		int i;
		unsigned int ui;
		long l;
		unsigned long ul;
		float f;
		double d;
		void *p;
	}v;

	/// This subclass serves as internal type information the PXL variant data type.
	class TypeInfo
	{
	public:
		inline TypeInfo(const char* name) :
		clear(0), dup(0), name(name)
		{}

		void (*clear)(Value* v);
		void (*dup)(Value* v);
		const char* name;
	};

	/// This method returns internal type information for the type given by \p t .
	static inline const TypeInfo& getTypeInfo(Type t)
	{
		if (PXL_UNLIKELY((std::size_t)t >= getTypes().size()))
		return fallbackGetTypeInfo(t);
		return getTypes()[t];
	}

	/// This method throws a detailed exception about a type mismatch between \p tShould and \p tIs .
	void wrongType(Type tShould, Type tIs) const throw (std::runtime_error);

private:
	static const TypeInfo& fallbackGetTypeInfo(Type t) throw (std::runtime_error);

	static std::vector<TypeInfo> & getTypes();
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
PCL_ANYBASE_SIMPLE(DOUBLE, double, d)
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

//template<>
//inline pxl::Serializable* VariantBase::get(Type t) const
//{
//	PCL_ANYBASE_CHECK(SERIALIZABLEPTR, t)
//	return &**((boost::shared_ptr<pxl::Serializable>*)(v.p));
//}

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

//template<>
//inline void VariantBase::set(Type t, pxl::Serializable* ser)
//{
//	PCL_ANYBASE_CHECK(SERIALIZABLEPTR, t)
//	if (v.p)
//	delete[] (boost::shared_ptr<pxl::Serializable>*)v.p;
//
//	v.p = new boost::shared_ptr<pxl::Serializable>(ser);
//}

template<>
inline VariantBase::Type VariantBase::findType<std::string>()
{	return TYPE_STRING;}

//template<>
//inline VariantBase::Type VariantBase::findType<pxl::Serializable*>()
//{	return TYPE_SERIALIZABLEPTR;}

#undef PCL_ANYBASE_CHECK

} // namespace pxl

#endif // PXL_BASE_VARIANTBASE_HH

namespace pxl
{

/**
 This class represents the PXL variant data type.
 It carries a generic data holder along with type information about
 the type of the currently stored value.
 */
class PXL_DLL_EXPORT Variant : public VariantBase
{
public:
	inline Variant() : type(TYPE_NULL)
	{}
	inline Variant(Type t) : type(t)
	{}
	inline Variant(const Variant& orig) :
	VariantBase(orig), type(orig.type)
	{
		dup(type);
	}

	inline Variant(const Variant* orig) :
	VariantBase(orig), type(orig->type)
	{
		dup(type);
	}

	inline ~Variant()
	{
		if (type != TYPE_NULL) VariantBase::clear(type);
	}

	inline Variant& operator=(const Variant& orig)
	{
		if (type != TYPE_NULL)
		{
			VariantBase::clear(type);
			type = orig.type;
		}
		else if (PXL_UNLIKELY(type != orig.type))
		wrongType(orig.type, type);
		v = orig.v;
		dup(type);
		return *this;
	}

	/// This method returns the type id of the currently stored value.
	inline Type getType() const
	{	return type;}
	/// This method returns the type name of the currently stored value.
	inline const char *getTypeName() const
	{	return getTypeInfo(type).name;}

	/// This method sets the type of the value to \p t and throws a pxl::Exception if a type is already assigned.
	inline void setType(Type t)
	{
		if (PXL_UNLIKELY(type != TYPE_NULL))
		wrongType(TYPE_NULL, type);
		type = t;
	}

	/// This method clears the contents and type of the currently stored value.
	inline void clear()
	{	VariantBase::clear(type); type = TYPE_NULL;}

	/// This method initialises the variant with an empty value of the type that matches the template instantiation.
	template<typename datatype>
	inline void init()
	{	setType(findType<datatype>());}

	/// This method returns the value of the variant or throws a Exception if the currently stored type doesn't match the template instantiation.
	template<typename datatype>
	inline datatype get() const
	{	return VariantBase::get<datatype>(type);}

	/// This method set the value of the variant or throws a Exception if the currently assigned type doesn't match the template instantiation.
	template<typename datatype>
	inline void set(datatype arg)
	{	VariantBase::set<datatype>(type, arg);}

protected:
	Type type;
};

} // namespace pxl

#endif // PXL_BASE_VARIANT_HH


namespace pxl
{
/**
 This class is intented to aggregate information complementary to data members in form
 of string-variant pairs.
 All PXL physics objects own user records and provide methods for quick
 access to individual user record entries.
 */
class PXL_DLL_EXPORT UserRecord
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
		inline std::map<std::string, Variant>* getData()
		{
			return &_data;
		}
		inline void setData(const std::map<std::string, Variant>* object)
		{
			_data = *object;
		}

		unsigned int _references;
		std::map<std::string, Variant> _data;

	}; //class Datasocket

public:
	typedef std::map<std::string, Variant> PtlMap;
	typedef std::pair<std::string, Variant> StlPair;
	typedef std::map<std::string, Variant>::const_iterator const_iterator;
	typedef std::map<std::string, Variant>::iterator iterator;

	UserRecord()
	{
		_dataSocket = new DataSocket;
	}
	UserRecord(const UserRecord& original)
	{
		_dataSocket = original._dataSocket;
		_dataSocket->_references++;
	}
	explicit UserRecord(const UserRecord* original)
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
		Variant& value = findOrAlloc(key);
		if (PXL_UNLIKELY(value.getType() != Variant::TYPE_NULL))
			value.clear();

		value.template init<datatype>();
		value.template set<datatype>(item);
	}

	/// This method template searches and returns the user record item identified by \p key; \p defaultitem is returned in case the key is not found.
	template<typename datatype> datatype find(const std::string& key,
			const datatype& defaultitem) const
	{
		const Variant* value = findOrReturn(key);
		if (!value)
			return defaultitem;
		return value->template get<datatype>();
	}

	/// This method template searches and returns the user record item indetified by \p key; a pxl::Exception is thrown in case the key is not found.
	template<typename datatype> datatype find(const std::string& key) const
			throw (std::runtime_error)
	{
		const Variant* value = findOrReturn(key);
		if (!value)
			throw std::runtime_error("pxl::UserRecord::find(...): key not found and no default item provided");
		return value->template get<datatype>();
	}

	/// This method templates checks if the user record entry identified by key is present.
	bool check(const std::string& key) const
	{
		const Variant* value = findOrReturn(key);
		if (!value)
			return false;
		return true;
	}

	/// This method template checks if user record entry identified by \p key is present.
	/// If yes, its value is put into the passed \p item.
	template<typename datatype> bool check(const std::string& key,
			datatype& item) const
	{
		const Variant* value = findOrReturn(key);
		if (!value)
			return false;
		item = value->template get<datatype>();
		return true;
	}
	
	/// This method template checks if a user record entry identified by \p key is present,
	/// and changes it to the passed value in case it is present. If not, an exception is thrown.
	template<typename datatype> void change(const std::string& key,
			datatype item) throw (std::runtime_error)
	{
		iterator found = setContainer()->find(key);
		if (found == getContainer()->end())
			throw std::runtime_error("pxl::UserRecord::change(...): UserRecord entry not found");
		
		Variant& value = found->second;
		if (value.getType() != VariantBase::findType<datatype>())
			throw std::runtime_error("pxl::UserRecord::change(...): UserRecord entry of wrong type");
				
		value.template set<datatype>(item);
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
	inline const std::map<std::string, Variant>* getContainer() const
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

	inline size_t size() const
	{
		return getContainer()->size();
	}

	std::ostream
			& print(int level = 0, std::ostream& os = std::cout, int pan = 0) const;
	const std::string toString() const
	{
		std::ostringstream ss;
		this->print(0,ss);
		return ss.str();
	}

private:
	DataSocket* _dataSocket;

	/// This method grants write access to the aggregated data;
	/// if necessary, the copy-on-write mechanism performs a deep copy of the aggregated data first.
	inline std::map<std::string, Variant>* setContainer()
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

	Variant& findOrAlloc(const std::string &key)
	{
		iterator insertPos = setContainer()->lower_bound(key);
		if (insertPos == getContainer()->end() || insertPos->first != key)
			return setContainer()->insert(insertPos, StlPair(key, Variant()))->second;
		else
			return insertPos->second;
	}

	const Variant* findOrReturn(const std::string &key) const
	{
		const_iterator found = getContainer()->find(key);
		if (found == getContainer()->end())
			return 0;
		return &found->second;
	}

};

} // namespace pxl

#endif // PXL_BASE_USERRECORD_HH


namespace pxl
{

/// The Index for the basicContainer - should be merged with the
/// Index.hh someday??
typedef std::map<std::string, Serializable*> ContainerIndex;

/**
 A Container for all objects that inherit from pxl::Serializable. This should be a rudiment version of the pxl::ObjectOwner.
 */
class PXL_DLL_EXPORT BasicContainer : public Serializable
{
public:

	BasicContainer() :
		Serializable()
	{
	}

	BasicContainer(const BasicContainer& basicContainer) :
		_container(), _index(), _uuidSearchMap(),
				_userRecords(basicContainer._userRecords)
	{
		this->init(basicContainer);
	}

	explicit BasicContainer(const BasicContainer* basicContainer) :
		_container(), _index(), _uuidSearchMap(),
				_userRecords(basicContainer->_userRecords)
	{
		this->init(*basicContainer);
	}

	BasicContainer& operator=(const BasicContainer& original)
	{
		if (this != &original)
		{
			clearContainer();
			this->init(original);
		}
		return *this;
	}

	void init(const BasicContainer& basicContainer);

	/// This destructor deletes all contained objects.
	virtual ~BasicContainer()
	{
		this->clearContainer();
	}

	/// This virtual method creates a deep copy and returns a C++ pointer to the newly created object.
	virtual Serializable* clone() const
	{
		return new BasicContainer(*this);
	}

	/// This method allows read access to the contained STL vector of
	/// Relative pointers to, e.g., use STL algorithms.
	const std::vector<Serializable*>& getObjects() const
	{
		return _container;
	}

	/// This method template creates a new instance of \p objecttype;
	/// objecttype must be a class inheriting from pxl::Serializable;
	/// the newly created instance is owned and will be deleted by this object owner.
	template<class objecttype> objecttype* create()
	{
		objecttype* pitem = new objecttype;
		_container.push_back(static_cast<Serializable*>(pitem));
		_uuidSearchMap.insert(std::pair<Id, Serializable*>(pitem->getId(), pitem));
		return pitem;
	}

	/// This method template creates a copy of \p original by invoking the
	/// copy constructor of \p objecttype; / \p objecttype must be a class
	/// inheriting from pxl::Serializable; / the newly created instance
	/// is owned and will be deleted by this object owner.
	template<class objecttype> objecttype* create(
	const objecttype* original)
	{
		objecttype* pitem = new objecttype(*original);
		_container.push_back(static_cast<Serializable*>(pitem));
		_uuidSearchMap.insert(std::pair<Id, Serializable*>(pitem->getId(), pitem));
		return pitem;
	}

	/// This method template creates a new \p objecttype instance by
	/// invoking a \p ctrtype overloaded constructor; / \p objecttype must
	/// be a class inheriting from pxl::Serializable / the newly created
	/// instance is owned and will be deleted by this object owner.
	template<class objecttype, class ctrtype> objecttype* create(
	const ctrtype& original)
	{
		objecttype* pitem = new objecttype(*original);
		_container.push_back(static_cast<Serializable*>(pitem));
		_uuidSearchMap.insert(std::pair<Id, Serializable*>(pitem->getId(), pitem));
		return pitem;
	}

	/// This method inserts \p item in the container of this object owner
	///and takes deletion responsability.
	void setObject(Serializable* item);

	/// This method deletes \p item.
	void remove(Serializable* item);

	/// This method takes \p item from the container.
	void take(Serializable* item);

	/// This method returns true if \p item is owned by this object owner.
	bool has(const Serializable* item) const;

	/// This method clears the object owner and deletes all owned objects.
	void clearContainer();

	/// Typedef for standard const_iterator.
	typedef std::vector<Serializable*>::const_iterator const_iterator;
	typedef std::vector<Serializable*>::iterator iterator;

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

	/// Returns the number of elements the BasicContainer holds.
	inline size_t size() const
	{
		return _container.size();
	}

	/// This method searches the index for the index-id \p idx and returns a
	///dynamically casted / C++ pointer of type \p objecttype* to the
	///corresponding object; / in case idx is not found a null pointer is
	///returned.
	template<class objecttype> inline objecttype* findObject(
	const std::string& idx) const // goes via Index & casts

	{
		ContainerIndex::const_iterator it = _index.find(idx);
		if (it!=_index.end())
		return dynamic_cast<objecttype*>(it->second);
		return 0;
	}

	inline Serializable* findObject(const std::string& idx) const
	{
		ContainerIndex::const_iterator it = _index.find(idx);
		if (it!=_index.end())
		return it->second;
		return 0;
	}

	Serializable* getById(Id id) const
	{
		std::map<Id, Serializable*>::const_iterator found = _uuidSearchMap.find(id);
		if ( found != _uuidSearchMap.end() )
		return found->second;
		return 0;
	}

	/// Fills into the passed vector weak pointers to the objects of the
	/// type specified by the template argument.
	template<class objecttype> size_t
	getObjectsOfType(std::vector<objecttype*>& vec) const
	{
		size_t size = vec.size();
		for (BasicContainer::const_iterator iter = begin(); iter!=end(); ++iter)
		{
			objecttype* obj = dynamic_cast<objecttype*>(*iter);
			if (obj!=0)
			vec.push_back(obj);
		}
		return vec.size()-size;
	}

	virtual const Id& getTypeId() const
	{
		return getStaticTypeId();
	}

	static const Id& getStaticTypeId()
	{
		static const Id id("38ebda57-df6f-a577-b811-0bba49745b09");
		return id;
	}

	virtual void serialize(const OutputStream &out) const;

	virtual void deserialize(const InputStream &in);

	//////////////////////////////////////////////
	//User Record
	//////////////////////////////////////////////

	/// This method provides access to the user records.
	inline const UserRecord& getUserRecord() const
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

	/// This method provides direct read access to the index.
	inline const ContainerIndex& getIndex() const
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

private:
	std::vector<Serializable*> _container;
	ContainerIndex _index;
	std::map<Id, Serializable*> _uuidSearchMap;
	UserRecord _userRecords;
};

} // namespace pxl

#endif // PXL_BASE_BASICCONTAINER_HH

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_EVENT_HH
#define PXL_BASE_EVENT_HH


//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_OBJECT_OWNER
#define PXL_BASE_OBJECT_OWNER

#include <algorithm>

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_OBJECTBASE_HH
#define PXL_BASE_OBJECTBASE_HH


//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_WKPTRBASE_HH
#define PXL_BASE_WKPTRBASE_HH


namespace pxl {

// ptl

class Relative;

/** 
 This base class provides common functionalities for all derived PXL weak pointers. 
 */
class PXL_DLL_EXPORT WkPtrBase
{
public:
	virtual ~WkPtrBase()
	{
		connect(0);
	}

	/// This virtual method creates a deep copy and returns a C++ pointer to the newly created weak pointer instance.  
	virtual WkPtrBase* clone() const
	{
		return new WkPtrBase(*this);
	}

	/// This method returns a C++ pointer of type Relative to the referenced object  
	inline Relative* pointer() const
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
	inline Relative* operator->() const
	{
		return access();
	}
	/// compare the referenced object pointers
	inline bool operator==(WkPtrBase &other) const
	{
		return (_objectRef == other.pointer());
	}
	/// compare the referenced object pointers   
	inline bool operator!=(WkPtrBase &other) const
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
		throw std::runtime_error("WkPtrBase::access(): FATAL: The object you intend to access does not exist!");
		return 0;
	}

protected:
	WkPtrBase() :
		_notifyChainIn(0), _notifyChainOut(0), _objectRef(0)
	{
	}

	void notifyDeleted();

	void connect(Relative* pointer);

	WkPtrBase* _notifyChainIn;
	WkPtrBase* _notifyChainOut;

	Relative* _objectRef;

	friend class Relative;
};

} // namespace pxl

#endif // PXL_BASE_WK_PTR_BASE_HH
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_RELATIONS_HH
#define PXL_BASE_RELATIONS_HH

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

class PXL_DLL_EXPORT Relations
{
public:
	typedef std::set<Relative*>::const_iterator const_iterator;
	typedef std::set<Relative*>::iterator iterator;

	Relations()
	{
	}

	void serialize(const OutputStream &out) const;

	const std::set<Relative*>& getContainer() const
	{
		return _relatives;
	}

	bool set(Relative* relative)
	{
		return (_relatives.insert(relative)).second;
	}

	bool erase(Relative* relative)
	{
		return (_relatives.erase(relative) > 0);
	}

	bool has(Relative* relative) const
	{
		return (_relatives.count(relative) > 0 );
	}
	
	Relative* getFirst() const
	{
		if (_relatives.begin() != _relatives.end())
			return (*_relatives.begin());
		return 0;
	}

	template<class objecttype> size_t getObjectsOfType(
			std::vector<objecttype*>& relatives) const
	{
		size_t size = relatives.size();
		for (const_iterator iter = _relatives.begin(); iter!=_relatives.end(); ++iter)
		{
			objecttype* obj = dynamic_cast<objecttype*>(*iter);
			if (obj)
				relatives.push_back(obj);
		}
		return relatives.size()-size;
	}

	size_t size() const
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
	std::set<Relative*> _relatives;
};

} // namespace pxl

#endif // PXL_BASE_RELATIONS_HH
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_SOFTRELATIONS_HH
#define PXL_BASE_SOFTRELATIONS_HH



namespace pxl
{
class ObjectOwner;
class BasicContainer;

class PXL_DLL_EXPORT SoftRelations
{
public:
	typedef std::multimap<std::string, Id>::const_iterator const_iterator;
	typedef std::multimap<std::string, Id>::iterator iterator;

	void serialize(const OutputStream &out) const
	{
		size_t size = _relationsMap.size();
		out.writeUnsignedInt((unsigned int)size);
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
		for (size_t i=0; i<size; ++i)
		{
			std::string name;
			in.readString(name);
			Id id(in);
			_relationsMap.insert(std::pair<std::string, Id>(name, id));
		}
	}

	Serializable* getFirst(const ObjectOwner& owner, const std::string& type = "") const;
	Serializable* getFirst(const BasicContainer& owner, const std::string& type = "") const;

	int getSoftRelatives(std::vector<Serializable*>& vec, const ObjectOwner& owner, const std::string& type) const;
	int getSoftRelatives(std::vector<Serializable*>& vec, const ObjectOwner& owner) const;

	int getSoftRelatives(std::vector<Serializable*>& vec, const BasicContainer& owner, const std::string& type) const;
	int getSoftRelatives(std::vector<Serializable*>& vec, const BasicContainer& owner) const;

	template <class objecttype>
	int getSoftRelativesOfType(std::vector<objecttype*>& vec, const ObjectOwner& owner, const std::string& type) const;

	template <class objecttype>
	int getSoftRelativesOfType(std::vector<objecttype*>& vec, const ObjectOwner& owner) const;

	int keepSoftRelatives(std::vector<Serializable*>& vec, const std::string& type) const;
	int keepSoftRelatives(std::vector<Serializable*>& vec) const;

	bool has(const Serializable* relative) const;

	bool has(const Serializable* relative, const std::string& name) const;

	bool has(const Id& id) const;
	bool has(const Id& id, const std::string& name) const;

	bool hasType(const std::string& name) const;

	int count(const std::string& name) const;

	void set(const Serializable* relative, const std::string& type);

	void remove(const Serializable* relative, const std::string& type);

	void remove(const Serializable* relative);

	size_t size() const
	{
		return _relationsMap.size();
	}

	void clearContainer()
	{
		_relationsMap.clear();
	}

	const std::multimap<std::string, Id>& getContainer() const
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

	std::ostream& print(int level = 0, std::ostream& os = std::cout, int pan=1) const;
	const std::string toString() const
	{
		std::ostringstream ss;
		this->print(0,ss);
		return ss.str();
	}

private:
	//this way round?
	std::multimap<std::string, Id> _relationsMap;
};

} //namespace pxl

#endif /*PXL_BASE_SOFTRELATIONS_HH*/
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_LAYOUT_HH
#define PXL_BASE_LAYOUT_HH


namespace pxl
{

// ptl

/** 
 This class holds layout information of PXL objects when visualized
 in the Graphical User Interface VisualPxl. For internal use only, 
 methods and data members are self-explanatory.  
 */
class PXL_DLL_EXPORT Layout
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

#endif // PXL_BASE_LAYOUT_HH

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
class PXL_DLL_EXPORT Relative : public Serializable
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
	inline Id id() const
	{
		return getId();
	}

	static const Id& getStaticTypeId()
	{
		static const Id id("5dee644d-906f-4d8e-aecc-d9a644293260");
		return id;
	}

	virtual const Id& getTypeId() const
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
	inline ObjectOwner* owner() const
	{
		return _refObjectOwner;
	}

	/// This virtual method creates a deep copy and returns a C++ pointer to the newly created object.  
	virtual Serializable* clone() const
	{
		return new Relative(*this);
	}

	/// This method grants access to the Relations instance managing mother relations  
	inline const Relations& getMotherRelations() const
	{
		return _motherRelations;
	}
	/// This method grants access to the pxl::Relations instance managing daughter relations  
	inline const Relations& getDaughterRelations() const
	{
		return _daughterRelations;
	}

	/// This method grants access to the pxl::Relations instance managing flat relations  
	inline const Relations& getFlatRelations() const
	{
		return _flatRelations;
	}
	
	/// This method returns the first entry of the mother relations.
	/// In case the collection is empty, 0 is returned.
	Relative* getMother() const
	{
		return _motherRelations.getFirst();
	}
	
	size_t numberOfDaughters() const
	{
		return _daughterRelations.size();
	}
	
	size_t numberOfMothers() const
	{
		return _motherRelations.size();
	}

	/// This method establishes a mother relation to the \p target object; please notice, that
	/// only relations between objects owned by the same object owner will be established.   
	void linkMother(Relative* target) throw(std::runtime_error);
	/// This method establishes a daughter relation to the \p target object; please notice, that
	/// only relations between objects owned by the same object owner will be established.   
	void linkDaughter(Relative* target) throw(std::runtime_error);
	/// This method establishes a flat relation to the \p target object; please notice, that
	/// only relations between objects owned by the same object owner will be established.
	void linkFlat(Relative* target) throw(std::runtime_error);

	/// This method removes an existing daughter relation to the \p target object.
	void unlinkMother(Relative* target);
	/// This method removes an existing daughter relation to the \p target object.
	void unlinkDaughter(Relative* target);
	/// This method removes an existing daughter relation to the \p target object.
	void unlinkFlat(Relative* target);

	/// This method removes all existing mother relations.
	void unlinkMothers();
	/// This method removes all existing daughter relations.
	void unlinkDaughters();
	/// This method removes all existing flat relations.
	void unlinkFlat();

	void linkSoft(Relative* relative, const std::string& type)
	{
		if (relative)
		{
			_softRelations.set(relative, type);
			relative->_softRelations.set(this, type);
		}
	}

	void unlinkSoft(Relative* relative, const std::string& type)
	{
		if (relative)
		{
			_softRelations.remove(relative, type);
			relative->_softRelations.remove(this, type);
		}
	}

	const SoftRelations& getSoftRelations() const
	{
		return _softRelations;
	}

	SoftRelations& setSoftRelations()
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
	inline Layout* layout()
	{
		if (!_ptrLayout)
		_ptrLayout = new Layout;
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
	virtual WkPtrBase* createSelfWkPtr()
	{
		throw std::runtime_error("pxl::ObjectBase::createSelfWkPtr(): ATTENTION! Inheriting class must reimplement this virtual method.");
		return 0;
	}

protected:
	/// Default constructor.
	Relative() :
	Serializable(), _refWkPtrSpec(0), _refObjectOwner(0),
	_name("default"), _ptrLayout(0)
	{
	}

	/// Copy constructor. Relations are not copied.
	Relative(const Relative& original) :
	Serializable(), _refWkPtrSpec(0), _refObjectOwner(0),
	_name(original._name), _ptrLayout(0)
	{
	}

	/// Copy constructor. Relations are not copied.
	explicit Relative(const Relative* original) :
	Serializable(), _refWkPtrSpec(0), _refObjectOwner(0),
	_name(original->_name), _ptrLayout(0)
	{
	}

	std::ostream& printPan1st(std::ostream& os, int pan) const;
	std::ostream& printPan(std::ostream& os, int pan) const;

private:
	WkPtrBase* _refWkPtrSpec;
	ObjectOwner* _refObjectOwner;

	Relations _motherRelations;
	Relations _daughterRelations;
	Relations _flatRelations;

	SoftRelations _softRelations;

	std::string _name;

	Layout* _ptrLayout;

	friend class WkPtrBase;
	friend class ObjectOwner;

	Relative& operator=(const Relative& original)
	{
		return *this;
	}
};

}
// namespace pxl

// operators
PXL_DLL_EXPORT std::ostream& operator<<(std::ostream& cxxx, const pxl::Relative& obj);

#endif // PXL_BASE_OBJECTBASE_HH
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_WEAK_PTR_HH
#define PXL_BASE_WEAK_PTR_HH



namespace pxl
{

// ptl

/** 
 This class template represents a weak pointer to PXL objects of \p objecttype that aggregate data of \p datatype.   
 */
template<class objecttype> class weak_ptr : public WkPtrBase
{
public:
	weak_ptr() :
		WkPtrBase()
	{
	}
	weak_ptr(objecttype* ptr) :
		WkPtrBase()
	{
		WkPtrBase::connect(ptr);
	}
	weak_ptr(objecttype& object) :
		WkPtrBase()
	{
		WkPtrBase::connect(&object);
	}
	weak_ptr(const weak_ptr<objecttype>& original) :
		WkPtrBase()
	{
		WkPtrBase::connect((objecttype*) original._objectRef);
	}
	explicit weak_ptr(const weak_ptr<objecttype>* original) :
		WkPtrBase()
	{
		WkPtrBase::connect((objecttype*) original->_objectRef);
	}

	virtual ~weak_ptr()
	{
		WkPtrBase::connect(0);
	}

	/// This virtual method creates a deep copy and returns a C++ pointer to the newly created weak pointer instance.  
	virtual WkPtrBase* clone() const
	{
		return new weak_ptr<objecttype>(*this);
	}

	/// This assignment operator causes the weak pointer to reference the object referenced by \p pptr.
	inline void operator=(const weak_ptr<objecttype>& pptr)
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
			throw std::runtime_error("WkPtrSpec::cast_dynamic(): Unsupported multiple inheritance configuration.");
		return reinterpret_cast<weak_ptr<objecttype>*>(orig);
	}

	// safe access to object
	inline objecttype* access() const throw (std::runtime_error)
	{
		if (_objectRef)
			return (objecttype*)_objectRef;
		throw std::runtime_error("WkPtrSpec::access(): FATAL: The object you intend to access does not exist!");
		return 0;
	}

};

template<class objecttype> objecttype& operator*(weak_ptr<objecttype>& wkPtr)
{
	return wkPtr.object();
}

template<class objecttype> const objecttype& operator*(
		const weak_ptr<objecttype>& wkPtr)
{
	return wkPtr.object();
}

} // namespace pxl


#endif // PXL_BASE_WEAK_PTR_HH

namespace pxl
{

class ObjectOwner;

// - - - - - - - - - - - - - - - - - - - - - - - - - - -
/// For STL-style iteration on selective class: iterator class template;
/// this iterator behaves like a normal STL iterator but ignores all objects
/// that cannot be interpreted as type objecttype (tested by dynamic casts).
/// Use in STL-style, except that either begin gets the objecttype class as template argument,
/// or a constructor of the TypeIterator is used.
template<class objecttype> class ObjectOwnerTypeIterator
{
	typedef std::vector<Relative*>::const_iterator const_iterator;
	typedef std::vector<Relative*>::iterator iterator;

public:
	/// Copy constructor.
	ObjectOwnerTypeIterator(const ObjectOwnerTypeIterator& other);

	/// Constructor from ObjectOwner instance.
	ObjectOwnerTypeIterator(const ObjectOwner& container);

	/// Constructor from ObjectOwner instance.
	ObjectOwnerTypeIterator(const ObjectOwner* container);

	const ObjectOwnerTypeIterator<objecttype> operator++(int);

	const ObjectOwnerTypeIterator<objecttype>& operator++();

	inline objecttype* operator*();

	inline bool operator==(const_iterator iter);

	inline bool operator!=(const_iterator iter);

private:
	const ObjectOwner* _containerRef;
	const_iterator _iter;
};

/**
 This class is an active container for Relative derivatives, such as
 Vertex, Particle or EventView;
 it has the ownership and deletion responsability for the contained objects.
 The method template create() can be used to create derivatives
 of Relative within object owners. The method set() can be
 used to manually add objects, has() tests object ownerships, and
 remove() explicitely removes objects from the owner and deletes them.
 The copy constructor of the class ObjectOwner also produces copies of
 the contained objects, and re-establishes corresponding mother-daughter relations
 amongst the copied objects. For the convenience of a quick and targeted object access, the
 newly created owner carries a so-called copy history for mapping original and
 copied objects. This information is used by the findCopyOf() method:
 provided a reference to the original object, this method returns a pointer to the copied object.
 A further, powerful tool for targeted object access is the so-called index, which
 allows to map objects to unique string identifiers, the index-id. The method findObject()
 can be used to directly access objects by their index-ids or object-ids.
 The ObjectOwner extends the functionality of the contained STL vector. It provides a selective iterator, the class template
 ObjectOwner::TypeIterator, that ignores all objects other than the
 specialized data type.
 */
class PXL_DLL_EXPORT ObjectOwner
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
	ObjectOwner(const ObjectOwner& original) :
		_container(), _copyHistory(), _index(), _uuidSearchMap()
	{
		this->init(original);
	}

	/// This copy constructor performs a deep copy of object
	/// owner \p original and all contained objects with (redirected) relations.
	/// A copy history keeps track of originals and copies
	/// and the findCopyOf() method allows quick access to the copies.
	explicit ObjectOwner(const ObjectOwner* original) :
		_container(), _copyHistory(), _index(), _uuidSearchMap()
	{
		this->init(*original);
	}
	/// This destructor deletes all contained objects.
	virtual ~ObjectOwner()
	{
		this->clearContainer();
	}

	ObjectOwner& operator=(const ObjectOwner& original)
	{
		if (this != &original)
		{
			this->clearContainer();
			this->init(original);
		}
		return *this;
	}

	virtual void serialize(const OutputStream &out) const;

	virtual void deserialize(const InputStream &in);

	/// This method template creates a new instance of \p objecttype;
	/// objecttype must be a class inheriting from Relative;
	/// the newly created instance is owned and will be deleted by this object owner.
	template<class objecttype> objecttype* create()
	{
		objecttype* pitem = new objecttype;
		pitem->_refObjectOwner = this;
		_container.push_back(static_cast<Relative*>(pitem));
		_uuidSearchMap.insert(std::pair<Id, Relative*>(pitem->getId(), pitem));
		return pitem;
	}

	/// This method template creates a copy of \p original by invoking the copy constructor of \p objecttype;
	/// \p objecttype must be a class inheriting from Relative;
	/// the newly created instance is owned and will be deleted by this object owner.
	template<class objecttype> objecttype* create(
	const objecttype* original)
	{
		objecttype* pitem = new objecttype(*original);
		pitem->_refObjectOwner = this;
		_container.push_back(static_cast<Relative*>(pitem));
		_uuidSearchMap.insert(std::pair<Id, Relative*>(pitem->getId(), pitem));
		return pitem;
	}

	/// This method template creates a new \p objecttype instance by invoking a \p ctrtype overloaded constructor;
	/// \p objecttype must be a class inheriting from Relative;
	/// the newly created instance is owned and will be deleted by this object owner.
	template<class objecttype, class ctrtype> objecttype* create(
	const ctrtype& original)
	{
		objecttype* pitem = new objecttype(*original);
		pitem->_refObjectOwner = this;
		_container.push_back(static_cast<Relative*>(pitem));
		_uuidSearchMap.insert(std::pair<Id, Relative*>(pitem->getId(), pitem));
		return pitem;
	}

	/// This method inserts \p item in the container of this object owner and takes deletion responsability.
	void set(Relative* item);
	/// This method deletes \p item.
	void remove(Relative* item);
	/// This method takes \p item from the object owner
	void take(Relative* item);
	/// This method returns true if \p item is owned by this object owner.
	bool has(const Relative* item) const;

	/// This method clears the object owner and deletes all owned objects.
	void clearContainer();

	/// This method searches the index for the index-id \p idx and returns a dynamically casted
	/// C++ pointer of type \p objecttype* to the corresponding object;
	/// in case idx is not found a null pointer is returned.
	template<class objecttype> inline objecttype* findObject(
	const std::string& idx) const // goes via Index & casts

	{
		std::map<std::string, Relative*>::const_iterator it = _index.find(idx);
		if (it!=_index.end())
		return dynamic_cast<objecttype*>(it->second);
		return 0;
	}

	inline Relative* findObject(const std::string& idx) const
	{
		std::map<std::string, Relative*>::const_iterator it = _index.find(idx);
		if (it!=_index.end())
		return it->second;
		return 0;
	}

	Relative* getById(Id id) const
	{
		std::map<Id, Relative*>::const_iterator found = _uuidSearchMap.find(id);
		if ( found != _uuidSearchMap.end() )
		return found->second;
		return 0;
	}

	/// This method searches the copy history to locate the copy of \p original and
	/// returns a dynamically casted C++ pointer of type \p objecttype* to the corresponding copy;
	/// in case no copy can be traced a null pointer is returned.
	template<class objecttype> objecttype* findCopyOf(
	const Relative* original) const // goes via CopyHistory & casts

	{
		std::map<Id, Relative*>::const_iterator it = _copyHistory.find(original->id());
		if (it!=_copyHistory.end())
		return dynamic_cast<objecttype*>(it->second);
		return 0;
	}

	/// This method provides direct access to the copy history (created by the copy constructor).
	inline const std::map<Id, Relative*>& getCopyHistory() const
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
	bool setIndex(const std::string& idx, Relative* obj, bool overwrite = false);

	/// This method provides direct read access to the index.
	inline const std::map<std::string, Relative*>& getIndex() const
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
	const std::vector<Relative*>& getObjects() const
	{
		return _container;
	}

	/// Typedef for standard const_iterator.
	typedef std::vector<Relative*>::const_iterator const_iterator;
	typedef std::vector<Relative*>::iterator iterator;

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
	inline size_t size() const
	{
		return _container.size();
	}

	/// Fills into the passed vector weak pointers to the objects of the type specified by the template argument.
	template<class objecttype> size_t getObjectsOfType(std::vector<objecttype*>& vec) const
	{
		size_t size = vec.size();
		for (ObjectOwner::const_iterator iter = begin(); iter!=end(); ++iter)
		{
			objecttype* obj = dynamic_cast<objecttype*>(*iter);
			if (obj!=0)
			vec.push_back(obj);
		}
		return vec.size()-size;
	}



	/// This templated method provides an STL-style begin()-method to
	/// initialise the TypeIterator.
	template<class objecttype> const ObjectOwnerTypeIterator<objecttype> begin() const
	{
		ObjectOwnerTypeIterator<objecttype> it(this);
		return it;
	}
	
	void sort(int (*comp)(Relative*, Relative*))
	{
		std::sort(_container.begin(), _container.end(), comp);
	}
		

private:
	void init(const ObjectOwner& original);

	std::vector<Relative*> _container;
	std::map<Id, Relative*> _copyHistory;
	std::map<std::string, Relative*> _index;
	std::map<Id, Relative*> _uuidSearchMap;

};

/// Copy constructor.
template <class objecttype>
ObjectOwnerTypeIterator<objecttype>::ObjectOwnerTypeIterator(const ObjectOwnerTypeIterator& other) :
	_containerRef(other._containerRef), _iter(other._iter)
{
}

/// Constructor from ObjectOwner instance.
template <class objecttype>
ObjectOwnerTypeIterator<objecttype>::ObjectOwnerTypeIterator(const ObjectOwner& container) :
	_containerRef(&container), _iter(container.begin())
{
	if ( _iter!=_containerRef->end() && dynamic_cast<objecttype*>(*_iter)==0)
	(*this)++;
}

/// Constructor from ObjectOwner instance.
template <class objecttype>
ObjectOwnerTypeIterator<objecttype>::ObjectOwnerTypeIterator(const ObjectOwner* container) :
_containerRef(container), _iter(container->begin())
{
	if ( _iter!=_containerRef->end() && dynamic_cast<objecttype*>(*_iter)==0)
	(*this)++;
}

template <class objecttype>
const ObjectOwnerTypeIterator<objecttype> ObjectOwnerTypeIterator<objecttype>::operator++(int)
{
	ObjectOwnerTypeIterator orig = *this;
	if (_iter!=_containerRef->end())
	do
	_iter++;
	while (_iter!=_containerRef->end()
	&& dynamic_cast<objecttype*>(*_iter)==0);
	return orig;
}

template <class objecttype>
const ObjectOwnerTypeIterator<objecttype>& ObjectOwnerTypeIterator<objecttype>::operator++()
{
	if (_iter!=_containerRef->end())
	do
	_iter++;
	while (_iter!=_containerRef->end()
	&& dynamic_cast<objecttype*>(*_iter)==0);
	return *this;
}

template <class objecttype>
inline objecttype* ObjectOwnerTypeIterator<objecttype>::operator*()
{
	return _iter==_containerRef->end() ? 0
	: dynamic_cast<objecttype*>(*_iter);
}

template <class objecttype>
inline bool ObjectOwnerTypeIterator<objecttype>::operator==(const_iterator iter)
{
	return (_iter==iter);
}

template <class objecttype>
inline bool ObjectOwnerTypeIterator<objecttype>::operator!=(const_iterator iter)
{
	return (_iter!=iter);
}


} // namespace pxl

#endif // PXL_BASE_OBJECT_OWNER


namespace pxl
{

class PXL_DLL_EXPORT Event : public Serializable
{
public:
	Event() :
		Serializable()
	{
	}

	Event(const Event& event) :
		_objects(event._objects), _userRecords(event._userRecords)
	{
	}

	explicit Event(const Event* event) :
		_objects(event->_objects), _userRecords(event->_userRecords)
	{
	}

	virtual ~Event()
	{
	}

	virtual const Id& getTypeId() const
	{
		return getStaticTypeId();
	}

	static const Id& getStaticTypeId();

	virtual void serialize(const OutputStream &out) const
	{
		Serializable::serialize(out);
		_objects.serialize(out);
		_userRecords.serialize(out);
	}

	virtual void deserialize(const InputStream &in)
	{
		Serializable::deserialize(in);
		_objects.deserialize(in);
		_userRecords.deserialize(in);
	}


	/// This virtual method creates a deep copy and returns a C++ pointer to the newly created object.
	virtual Serializable* clone() const
	{
		return new Event(*this);
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
	inline void setObject(Relative* obj)
	{
		_objects.set(obj);
	}

	/// This method inserts \p obj with the index-id \p idx in the container of the object owner and takes deletion responsability.
	inline void setObject(Relative* obj, const std::string& idx)
	{
		_objects.set(obj);
		setIndex(idx, obj);
	}

	/// This method registers the object \p obj with the index-id \p idx in the index and returns true in case of success;
	/// please notice, that obj must be owned by this object owner and \p idx must not be a zero length string.
	inline bool setIndex(const std::string& idx, Relative* obj)
	{
		return _objects.setIndex(idx, obj);
	}

	inline ObjectOwner& getObjectOwner()
	{
		return _objects;
	}

	/// This method provides const access to the object owner.
	inline const ObjectOwner& getObjectOwner() const
	{
		return _objects;
	}

	template<class objecttype> inline void getObjectsOfType(
			std::vector<objecttype*>& vec) const
	{
		_objects.getObjectsOfType<objecttype>(vec);
	}

	inline const std::vector<Relative*>& getObjects() const
	{
		return _objects.getObjects();
	}

	/// This method deletes the object \p obj.
	inline void removeObject(Relative* obj)
	{
		_objects.remove(obj);
	}

	/// This method takes the object \p obj from the object owner.
	inline void takeObject(Relative* obj)
	{
		_objects.take(obj);
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
			const Relative* original) const
	{
		return _objects.findCopyOf<objecttype>(original);
	}

	/// This method provides direct access to the index.
	inline const std::map<std::string, Relative*>& getIndex() const
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

	//////////////////////////////////////////////
	//User Record
	//////////////////////////////////////////////

	/// This method provides access to the user records.
	inline const UserRecord& getUserRecord() const
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
	inline bool checkUserRecord(
			const std::string& key) const
	{
		return _userRecords.check(key);
	}

	/// This method checks if the user record entry identified by \p key is present. If yes, \p item is set to the according value.
	template<typename datatype> inline bool checkUserRecord(
			const std::string& key, datatype& item) const
	{
		return _userRecords.template check<datatype>(key, item);
	}
	
		
	/// This method checks if the user record entry identified by \p key is present. If yes, \p item is set to the according value.
	template<typename datatype> 
	void changeUserRecord(const std::string& key, datatype item)
	{
		_userRecords.template change<datatype>(key, item);
	}

	virtual std::ostream
			& print(int level=1, std::ostream& os=std::cout, int pan=1) const;

private:
	ObjectOwner _objects;
	UserRecord _userRecords;
};

} // namespace pxl

#endif // PXL_BASE_EVENT_HH
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_FILTER_HH
#define PXL_BASE_FILTER_HH



namespace pxl
{

template<class comparetype> class Comparator
{
public:
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
 it handles objects of type \p objecttype sorted by the criterion \p compare. The user can write
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
	static size_t apply(const ObjectOwner& objects,
			std::vector<objecttype*>& fillVector,
			const FilterCriterion<objecttype>& criterion)
	{
		size_t size = fillVector.size();

		// fill map:
		for (ObjectOwnerTypeIterator<objecttype> iter(objects); iter
				!=objects.end(); ++iter)
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


/** 
 This class template provides a sorted filter for PXL physics objects;
 it handles objects of type \p objecttype sorted by the criterion \p compare. The user can write
 his own filters by public inheritance from this class and reimplementation
 of the methods pass() and sort().  
 */
template<class objecttype, class criteriontype, class compare> class Filter2
{
public:
	
	Filter2(const ObjectOwner& objects,
			std::vector<objecttype*>& fillVector,
			const criteriontype& criterion)
	{
		this->apply(objects, fillVector, criterion);
	}
	
	virtual ~Filter2()
	{
	}

	/// This method applies the filter by running over the \p objects container and fills
	/// the passed vector with pointers to the objects passing the filter criteria.
	virtual int apply(const ObjectOwner& objects,
			std::vector<objecttype*>& fillVector,
			const FilterCriterion<objecttype>& criterion)
	{
		int size = fillVector.size();

		// fill map:
		for (ObjectOwnerTypeIterator<objecttype> iter(objects); iter
				!=objects.end(); ++iter)
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

#endif // PXL_BASE_FILTER_HH
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_FUNCTIONS_HH
#define PXL_BASE_FUNCTIONS_HH


namespace pxl {

/// @internal This function returns a platform-specific CPU timestamp; internally used for performance tests.
double getCpuTime();

} // namespace pxl

#endif // PXL_BASE_FUNCTIONS_HH
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_INFORMATIONCHUNK_HH
#define PXL_BASE_INFORMATIONCHUNK_HH



namespace pxl
{

class PXL_DLL_EXPORT InformationChunk : public Serializable
{
public:

	virtual const Id& getTypeId() const
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
		Serializable::serialize(out);
		_userRecords.serialize(out);
	}

	virtual void deserialize(const InputStream &in)
	{
		Serializable::deserialize(in);
		_userRecords.deserialize(in);
	}

	/// This virtual method creates a deep copy and returns a C++ pointer to the newly created object.
	virtual Serializable* clone() const
	{
		return new InformationChunk(*this);
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
	inline const UserRecord& getUserRecord() const
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
	inline bool checkUserRecord(
			const std::string& key) const
	{
		return _userRecords.check(key);
	}

	/// This method checks if the user record entry identified by \p key is present. If yes, \p item is set to the according value.
	template<typename datatype> inline bool checkUserRecord(
			const std::string& key, datatype& item) const
	{
		return _userRecords.template check<datatype>(key, item);
	}

private:
	std::string _name;
	UserRecord _userRecords;

};

} //namespace pxl

#endif /*PXL_BASE_INFORMATIONCHUNK_HH*/
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_OBJECT_HH
#define PXL_BASE_OBJECT_HH



namespace pxl {

/** 
 This class provides common functionalities of PXL physics objects like 
 data members for storing an object name and flags for status, Monte-Carlo mode and 
 object locking; more specific information, such as b-tags, jet cone sizes or energy 
 corrections, for instance, can be stored in the so-called user records (see UserRecord). 
 An integer workflag facilitates tagging of individual objects. 
 */
class PXL_DLL_EXPORT Object : public Relative
{
public:
	Object() :
	Relative(), _locked(0), _workflag(0), _userRecords()
	{
	}

	Object(const Object& original) :
	Relative(original), _locked(original._locked),
	_workflag(original._workflag),
	_userRecords(original._userRecords)
	{
	}

	explicit Object(const Object* original) :
	Relative(*original), _locked(original->_locked),
	_workflag(original->_workflag),
	_userRecords(original->_userRecords)
	{
	}

	virtual const Id& getTypeId() const
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
		Relative::serialize(out);
		out.writeBool(_locked);
		out.writeInt(_workflag);
		_userRecords.serialize(out);
	}

	virtual void deserialize(const InputStream &in)
	{
		Relative::deserialize(in);
		in.readBool(_locked);
		in.readInt(_workflag);
		_userRecords.deserialize(in);
	}

	/// This virtual method creates a deep copy and returns a C++ pointer to the newly created object.  
	virtual Serializable* clone() const
	{
		return new Object(*this);
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
	inline const UserRecord& getUserRecord() const
	{
		return _userRecords;
	}

	inline void setUserRecord(const UserRecord& userRecord)
	{
		_userRecords = userRecord;
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
	
	/// This method checks if the user record entry identified by \p key is present. If yes, \p item is set to the according value.
	template<typename datatype> 
	void changeUserRecord(const std::string& key, datatype item)
	{
		_userRecords.template change<datatype>(key, item);
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

	virtual WkPtrBase* createSelfWkPtr()
	{
		return new weak_ptr<Object>(this);
	}

private:
	bool _locked;
	int _workflag;

	UserRecord _userRecords;
	
	Object& operator=(const Object& original)
	{
		return *this;
	}
};

///// This typedef defines a weak pointer for Object
typedef weak_ptr<Object> ObjectWkPtr;

}
// namespace pxl

#endif // PXL_BASE_OBJECT_HH

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_BASE_BASICMANAGER_HH
#define PXL_BASE_BASICMANAGER_HH




namespace pxl
{

// pol

/**
 This class the functionality of the pxl::BasicObjectData class by providing an object owner (see pxl::ObjectOwner) and
 corresponding service methods. This way, physics objects (like instances of the classes
 pxl::Particle, pxl::Vertex or pxl::Collision as well as other arbitrary pxl::Relative derivatives can be
 aggregated and managed.
 */
class PXL_DLL_EXPORT ObjectManager : public Object
{
public:
	ObjectManager() :
		Object(), _objects()
	{
	}
	/// This copy constructor performs a deep copy of \p original
	/// with all contained objects and their (redirected) relations.
	ObjectManager(const ObjectManager& original) :
		Object(original), _objects(original._objects)
	{
	}
	/// This copy constructor performs a deep copy of \p original
	/// with all contained objects and their (redirected) relations.
	explicit ObjectManager(const ObjectManager* original) :
		Object(original), _objects(original->_objects)
	{
	}

	virtual const Id& getTypeId() const
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
		Object::serialize(out);
		_objects.serialize(out);
	}

	virtual void deserialize(const InputStream &in)
	{
		Object::deserialize(in);
		_objects.deserialize(in);
	}

	/// This virtual method creates a deep copy and returns a C++ pointer to the newly created object.
	virtual Serializable* clone() const
	{
		return new ObjectManager(*this);
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
	inline void setObject(Relative* obj)
	{
		_objects.set(obj);
	}

	/// This method inserts \p obj with the index-id \p idx in the container of the object owner and takes deletion responsability.
	inline void setObject(Relative* obj, const std::string& idx)
	{
		_objects.set(obj);
		setIndex(idx, obj);
	}

	/// This method registers the object \p obj with the index-id \p idx in the index and returns true in case of success;
	/// please notice, that obj must be owned by this object owner and \p idx must not be a zero length string.
	inline bool setIndex(const std::string& idx, Relative* obj)
	{
		return _objects.setIndex(idx, obj);
	}

	/// This method provides access to the object owner.
	inline ObjectOwner& getObjectOwner()
	{
		return _objects;
	}

	inline const ObjectOwner& getObjectOwner() const
	{
		return _objects;
	}

	inline const std::vector<Relative*>& getObjects() const
	{
		return _objects.getObjects();
	}

	template<class objecttype> inline void getObjectsOfType(
			std::vector<objecttype*>& vec) const
	{
		_objects.getObjectsOfType<objecttype>(vec);
	}

	/// This method deletes the object \p obj.
	inline void removeObject(Relative* obj)
	{
		_objects.remove(obj);
	}

	/// This method takes the object \p obj from the object owner.
	inline void takeObject(Relative* obj)
	{
		_objects.take(obj);
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
			const Relative* original) const
	{
		return _objects.findCopyOf<objecttype>(original);
	}

	/// This method provides direct access to the copy history (created by the copy constructor).
	inline const std::map<Id, Relative*>& getCopyHistory() const
	{
		return _objects.getCopyHistory();
	}

	/// This method clears the copy history  (created by the copy constructor).
	inline void clearCopyHistory()
	{
		_objects.clearCopyHistory();
	}

	/// This method provides direct access to the index.
	inline const std::map<std::string, Relative*>& getIndex() const
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

	virtual WkPtrBase* createSelfWkPtr()
	{
		return new weak_ptr<ObjectManager>(this);
	}

private:
	ObjectOwner _objects;

	ObjectManager& operator=(const ObjectManager& original)
	{
		return *this;
	}
};

///// This typedef defines a weak pointer for ObjectManager
typedef weak_ptr<ObjectManager> ObjectManagerWkPtr;

} // namespace pxl


#endif // PXL_BASE_BASICMANAGER_HH

#endif // PXL_BASE_HH
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_IO_HH
#define PXL_IO_HH

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_IO_CHUNK_READER_HH
#define PXL_IO_CHUNK_READER_HH

#include <fstream>


#define iotl__iStreamer__lengthUnzipBuffer 65536

namespace pxl
{

namespace skipSpace
{

enum skipMode
{
	off = 0,
	on
};

}

/**
 This class implemenents various methods for reading from PXL I/O and
 bases on the event structure defined in the ChunkWriter class. 
 The entry point for the user is the class InputFile.
 */

using namespace skipSpace;

class PXL_DLL_EXPORT ChunkReader
{
public:

	/// The status flag can take one out four values, either before an
	/// event (pre-header), or before a block, then either within an information chunk,
	/// or an event, or something unknown
	enum statusFlag
	{
		preHeader = 0,
		evPreBlock,
		infoPreBlock,
		preBlock
	};

	/// The readMode indicates if anything is read, just events, or only information chunks
	enum readMode
	{
		all = 0,
		event,
		infoChunk
	};

	/// The infoMode can be passed as a flag to indicate that a string condition
	/// must be fulfilled to read a block, an event, or an information chunk
	enum infoMode
	{
		ignore = 0,
		evaluate
	};

	/// The fileMode flag says if the read-in istream is seekable, or not. It needs to be modified by
	/// implementing classes. 
	enum fileMode
	{
		nonSeekable = 0,
		seekable
	};

	ChunkReader(std::istream& stream, fileMode seekMode = seekable) :
	_stream(stream), _status(preHeader), _sectionCount(0), _seekMode(seekMode)
	{
		_inputBuffer =	new unsigned char[iotl__iStreamer__lengthUnzipBuffer];
		_outputBuffer = new unsigned char[iotl__iStreamer__lengthUnzipBuffer];
	}

	~ChunkReader()
	{
		delete[] _inputBuffer;
		delete[] _outputBuffer;
	}
	
	void reset()
	{
		_sectionCount = 0;
		_status = preHeader;
		_buffer.destroy();
	}

	unsigned long getSectionCount()
	{
		return _sectionCount;
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
	bool next(skipMode skip = on, infoMode checkInfo = ignore,
			const std::string& infoCondition = "")
	{
		return readHeader(all, skip, checkInfo, infoCondition);
	}

	/// Reads in the header of the next event. False is returned if not successful.
	bool nextEvent(skipMode skip = on, infoMode checkInfo = ignore,
			const std::string& infoCondition = "")
	{
		return readHeader(event, skip, checkInfo, infoCondition);
	}

	/// Reads in the header of the next event. False is returned if not successful.
	bool nextInformationChunk(skipMode skip = on,
			infoMode checkInfo = ignore, const std::string& infoCondition = "")
	{
		return readHeader(infoChunk, skip, checkInfo, infoCondition);
	}

	/// Reads the next block and puts data into the input stream. False is returned if not successful.
	bool nextBlock(skipMode skip = on, infoMode checkInfo = ignore,
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
	
	inline void setStatus(statusFlag flag)
	{
		_status=flag;
	}

	void endEvent()
	{
		if (_status!=preHeader)
		while (nextBlock())
		;
	}
	
	/// Returns the size of the associated file.
	size_t getSize() const;
	
	/// Returns the current position in the associated file.
	size_t getPosition() const
	{
		return _stream.tellg();
	}
	
	bool eof() const
	{
		return _stream.peek()==EOF;
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
	ChunkReader(const ChunkReader& original) : _stream(original._stream)
	{
	}
	
	ChunkReader& operator= (const ChunkReader& other)
	{
		return *this;
	}
		
		
	std::istream& _stream;
	BufferInput _buffer;
	/// Status flag. 0 at end of event, 1 at end of block.
	statusFlag _status;
	unsigned long _sectionCount;
	fileMode _seekMode;
	
	unsigned char* _inputBuffer;
	unsigned char* _outputBuffer;
};

} //namespace pxl

#endif // PXL_IO_CHUNK_READER_HH
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_IO_CHUNKWRITER_HH
#define PXL_IO_CHUNKWRITER_HH



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
class PXL_DLL_EXPORT ChunkWriter
{
public:
	ChunkWriter(std::ostream& stream, char compressionMode = '6') :
		_stream(stream), _nBytes(0), _compressionMode(compressionMode)
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
	bool write(std::string info) throw(std::runtime_error);

	const BufferOutput& getOutputStream()
	{
		return _buffer;
	}
	
	void setCompressionMode(char compressionMode)
	{
		if ('0'<= compressionMode && compressionMode <= '9')
			_compressionMode = compressionMode;
	}
	
	void setCompressionMode(int compressionMode)
	{
		if (0<= compressionMode && compressionMode <= 9)
		_compressionMode = '0' + compressionMode;
	}
	
protected:
	/// Write char flag.
	bool writeFlag(char cEvtMarker);
	
private:
	ChunkWriter(const ChunkWriter& original) : _stream(original._stream)
	{
	}
	
	ChunkWriter& operator= (const ChunkWriter& other)
	{
		return *this;
	}
		
	std::ostream& _stream;
	BufferOutput _buffer;
	pxl::int32_t _nBytes;
	char _compressionMode;
};
}
#endif /*PXL_IO_CHUNKWRITER_HH*/
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_IO_GENERICINPUTHANDLER_HH
#define PXL_IO_GENERICINPUTHANDLER_HH


//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_IO_INPUTHANDLER_HH
#define PXL_IO_INPUTHANDLER_HH



namespace pxl
{
// io
/**
 This abstract class offers the basic functionality for reading the PXL physics event structure.
 Derived classes can handle concrete I/O operations.
 */
using namespace skipSpace;

class PXL_DLL_EXPORT InputHandler
{
public:

	InputHandler() : _objectCount(0)
	{
	}

	virtual ~InputHandler()
	{
	}

	virtual ChunkReader& getChunkReader() = 0;

	/// This method returns the number of file sections from the current ChunkReader.
	/// File sections are counted regardless of their actual content (InformationChunk, Event, etc). 
	unsigned long getSectionCount()
	{
		return getChunkReader().getSectionCount();
	}

	/// This method reads in the header of the next file section (e.g. event or information chunk)
	bool nextFileSection()
	{
		if (getChunkReader().next())
		return true;
		return false;
	}

	/// This method reads in the header of the next file section. If an information chunk occurs,
	/// either false is returned, or the file section is skipped, depending on the skip mode.
	/// In case of an event, success is reported.
	bool nextEvent(skipMode doSkip = on)
	{
		if (getChunkReader().nextEvent(doSkip))
		return true;
		return false;
	}

	/// This method reads in the header of the next file section. If an event occurs,
	/// either false is returned, or the file section is skipped, depending on the skip mode.
	/// In case of an information chunk, success is reported.
	bool nextInformationChunk(skipMode doSkip = on)
	{
		if (getChunkReader().nextInformationChunk(doSkip))
		return true;
		return false;
	}

	/// This method reads in the next file section if the information condition is fulfilled. Else, false is returned.
	bool nextFileSectionIf(const std::string& info, skipMode doSkip = on)
	{
		if (getChunkReader().next(doSkip, ChunkReader::evaluate, info))
		return true;
		return false;
	}

	/// This method reads in the next file section if it is an event and
	/// if the information condition is fulfilled. Else, false is returned.
	bool nextEventIf(const std::string& info, skipMode doSkip = on)
	{
		if (getChunkReader().nextEvent(doSkip, ChunkReader::evaluate, info))
		return true;
		return false;
	}

	/// This method reads in the next file section if it is an information chunk and
	/// if the information condition is fulfilled. Else, false is returned.
	bool nextInformationChunkIf(const std::string& info,
			skipMode doSkip = on)
	{
		if (getChunkReader().nextInformationChunk(doSkip, ChunkReader::evaluate, info))
		return true;
		return false;
	}

	/// Use this method in case the file contains pxl::Events after a next-statement. 
	/// A pxl::Event is passed to this method and filled with the current event.
	/// If the file section header has not been read or other objects occur,
	/// false is returned.
	bool readEvent(Event* event);

	/// Use this method in case the file contains pxl::Events after a next-statement. 
	/// A pxl::Event is passed to this method and filled with the current event.
	/// If the file section header has not been read or other objects occur,
	/// false is returned.
	bool readEventIf(Event* event, const std::string& blockInfo,
			skipMode doSkip = on);

	/// Use this method in case the file contains pxl::BasicContainers after a next-statement. 
	/// A pxl::BasicContainer is passed to this method and filled with the current container.
	/// If the file section header has not been read or other objects occur,
	/// false is returned.
	bool readBasicContainer(BasicContainer* basicContainer);

	/// Use this method in case the file contains pxl::BasicContainers after a next-statement. 
	/// A pxl::BasicContainer is passed to this method and filled with the current container.
	/// If the file section header has not been read, the info condition is not fulfilled,
	/// or other objects occur, false is returned.
	bool readBasicContainerIf(BasicContainer* basicContainer, const std::string& blockInfo,
			skipMode doSkip = on);

	/// A pxl::InformationChunk is passed to this method and filled if the next block contains an information chunk.
	/// Else, false is returned.
	bool readInformationChunk(InformationChunk* chunk);

	/// A pxl::InformationChunk is passed to this method and filled if the info condition is fulfilled and there is an information chunk.
	/// Else, false is returned.
	bool readInformationChunkIf(InformationChunk* event,
			const std::string& blockInfo, skipMode doSkip = on);

	/// Returns if the current file section is an information chunk. 
	inline bool isInformationChunk()
	{
		return getChunkReader().isInformationChunk();
	}

	/// Returns if the current file section is an event.
	inline bool isEvent()
	{
		return getChunkReader().isEvent();
	}

	/// This methods skips one file section.
	bool skip()
	{
		return getChunkReader().skip();
	}

	/// This method goes to the previous file section.
	bool previous()
	{
		return getChunkReader().previous();
	}

	/// With this method, n file sections can be skipped in forward or backward direction.
	/// The number of actually skipped file sections is returned (positive or negative).
	int skipFileSections(int n);

	/// With this method, n file sections can be skipped in forward or backward direction.
	/// The number of actually skipped file sections is returned (positive or negative).
	/// WARNING: This method is depreciated, use skipFileSections instead
	int skipEvents(int n)
	{
		std::cerr << "WARNING: This method InputHandler::skipEvents is depreciated, use skipFileSections instead" << std::endl;
		return skipFileSections(n);
	}
	
	/// seek to the desired file section (with 0 being the first file section)
	bool seekToFileSection(int index);
	
	/// seek to the desired file section (with 0 being the first file section
	/// WARNING: This method is depreciated, use seekToFileSection instead
	bool seekToEvent(int index)
	{
		std::cerr << "WARNING: This method InputHandler::seekToEvent is depreciated, use seekToFileSection instead" << std::endl;
		return seekToFileSection(index);
	}
	
	
	/// This method reads in the next block.
	bool readBlock()
	{
		return getChunkReader().nextBlock();
	}

	/// This method reads in the next block if the info condition is fulfilled.
	bool readBlockIf(const std::string& blockInfo,
			skipMode doSkip = on);

	/// This method explicitly reads an object of type objecttype. 
	/// Caution: This method should only be used if the type of the following object is known by hard.
	/// Else use the readNextObject method. This method may be deprecated in the future.
	template<class objecttype> bool readObject(objecttype* obj) throw(std::runtime_error);

	/// Seek to the desired object (with 0 being the first object).
	/// Current implementation inefficient due to missing file index.
	Serializable* seekToObject(size_t index) throw(std::runtime_error);
	
	/// This method reads in the next object from the file, regardless of file section boundaries.
	/// In case there are no more objects to be read, a zero pointer is returned.
	/// Attention: This method returns an object which was created with new. The user takes
	/// deletion responsibility.
	Serializable* readNextObject() throw(std::runtime_error);

	/// This method reads in the previous object from the file, regardless of file section boundaries.
	/// Attention: This method returns an object which was created with new. The user takes
	/// deletion responsibility.	
	Serializable* readPreviousObject() throw(std::runtime_error);

	/// This method fills the objects from the read-in block into the passed vector. The number of added objects is returned.
	/// Attention: The objects in the vector are created with new, the user takes deletion responsibility.
	int readObjects(std::vector<Serializable*>& objects) throw(std::runtime_error);

	/// This method fills the objects from the read-in block into the passed pxl::BasicContainer. The number of added objects is returned.
	/// Deletion responsibility is taken by the BasicContainer.
	int readObjects(BasicContainer* container) throw(std::runtime_error);

	/// This method fills the objects from the read-in block into the passed pxl::Event. The number of added objects is returned.
	/// Only derivatives of pxl::Relative are filled, other items are skipped. 
	/// Deletion responsibility is taken by the Event.
	int readObjects(Event* event) throw(std::runtime_error);
	
	/// Returns the number of read objects
	size_t objectCount() const
	{
		return _objectCount;
	}
	
	/// Resets InputHandler in case a new file/stream is opened.
	void reset()
	{
		_objectCount = 0;
		getChunkReader().reset();
	}
	
	/// Returns the size of the associated file.
	size_t getSize()
	{
		return getChunkReader().getSize();
	}
	
	/// Returns the current position in the associated file.
	size_t getPosition()
	{
		return getChunkReader().getPosition();
	}
	

private:
	InputHandler(const InputHandler& original)
	{
	}
	
	InputHandler& operator= (const InputHandler& other)
	{
		return *this;
	}
		
	size_t _objectCount;
	
};

}
//namespace pxl

#endif //PXL_IO_INPUTHANDLER_HH

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

	GenericInputHandler(ChunkReader& reader) :
		InputHandler(),	_reader(&reader)
	{
	}

	virtual ~GenericInputHandler()
	{
	}

	virtual ChunkReader& getChunkReader() throw (std::runtime_error)
	{
		if (!_reader)
			throw std::runtime_error("GenericInputHandler::getChunkReader(): ChunkReader pointer invalid.");			
		return *_reader;
	}
	
	virtual void setChunkReader(ChunkReader* reader)
	{
		_reader=reader;
	}

private:
	GenericInputHandler(const GenericInputHandler& original)
	{
	}
	
	GenericInputHandler& operator= (const GenericInputHandler& other)
	{
		return *this;
	}
	
	ChunkReader* _reader;
};

} //namespace pxl

#endif /*PXL_IO_GENERICINPUTHANDLER_HH*/
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_IO_GENERICOUTPUTHANDLER_HH
#define PXL_IO_GENERICOUTPUTHANDLER_HH


//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_IO_OUTPUTHANDLER_HH
#define PXL_IO_OUTPUTHANDLER_HH



namespace pxl
{

// io
/**
 This abstract class allows the user an easy handling of the PXL general output. Methods to write event headers,
 the PXL event class, and information chunks are offered.
 */
class PXL_DLL_EXPORT OutputHandler
{
public:
        /// Creates an OutputHandler with a maximum size after which a block is written.
	/// Set the maximum size to 0 (or smaller) to not let this happen automatically.
	OutputHandler(size_t maxSize = 1048576, size_t maxNObjects = 1000);

	virtual ~OutputHandler();

	virtual ChunkWriter& getChunkWriter() = 0;

		
	/// This method queues the passed object for later writing to the output file.
	void streamObject(const Serializable* obj)
	{
		obj->serialize(getChunkWriter().getOutputStream());
		_nObjects++;
		if ((_maxSize>0 && getOutputStream().getSize() > _maxSize) || (_maxNObjects>0 &&_nObjects > _maxNObjects))
			writeFileSection();
	}
	
	/// This method streams the passed pxl::Event to the output file.
	/// A file section is finished if the given maximum section size is reached.
	void writeEvent(const Event* event)
	{
		event->serialize(getChunkWriter().getOutputStream());
		_nObjects++;
//		if ((_maxSize>0 && getOutputStream().getSize() > _maxSize) || (_maxNObjects>0 &&_nObjects > _maxNObjects))
//			writeStream();
		// Restore backwards compatible behaviour for VISPA. Will be changed if necessary 
		// functionality in PXL
		writeFileSection();
		
	}
	
	/// This method writes the passed pxl::InformationChunk to the output file.
	void writeInformationChunk(const InformationChunk* infoChunk)
	{
		infoChunk->serialize(getChunkWriter().getOutputStream());
		_nObjects++;
//		if ((_maxSize>0 && getOutputStream().getSize() > _maxSize) || (_maxNObjects>0 &&_nObjects > _maxNObjects))
//			writeStream();
		// Restore backwards compatible behaviour for VISPA. Will be changed if necessary 
		// functionality in PXL
		writeFileSection();
	}
	
	/// This method writes the passed pxl::BasicContainer to the output file.
	void writeBasicContainer(const BasicContainer* basicContainer)
	{
		basicContainer->serialize(getChunkWriter().getOutputStream());
		_nObjects++;
// 		if ((_maxSize>0 && getOutputStream().getSize() > _maxSize) || (_maxNObjects>0 &&_nObjects > _maxNObjects))
// 			writeStream();
		// Restore backwards compatible behaviour for VISPA. Will be changed if necessary 
		// functionality in PXL
		writeFileSection();
	}

	/// return the associated OutputStream
	const BufferOutput& getOutputStream()
	{
		return getChunkWriter().getOutputStream();
	}
	
	/// Use this method to write an information string describing the new file section (same as event). Otherwise, this method need not necessarily be used.
	bool newFileSection(const std::string& info);
	
	/// Use this method to write an information string describing the new event. Otherwise, this method need not necessarily be used.
	/// WARNING: The use of this method is depreciated.
	bool newEvent(const std::string& info)
	{
		return newFileSection(info);
	}

	/// Use this method to write out a block to file. This method is not needed if you use the writeFileSection-method.
	bool writeStream(const std::string& info = "");

	/// Use this method to write out a block to disk and finish the current file section.
	bool writeFileSection(const std::string& info = "");
	
	/// Use this method to write out a block to disk and finish the current event.
	/// WARNING: The use of this method is depreciated.
	bool writeEvent(const std::string& info = "")
	{
		return writeFileSection(info);
	}
	
	/// Finish section in case stream not written to file
	void finish()
	{
		if (getOutputStream().getSize() > 0)
			writeFileSection();
	}
	
	void setMaxNObjects(size_t maxNObjects)
	{
		_maxNObjects = maxNObjects;
	}
	
	size_t getMaxNObjects() const
	{
		return _maxNObjects;
	}
	
	void setMaxSize(size_t maxSize)
	{
		_maxSize = maxSize;
	}

	size_t getMaxSize() const
	{
		return _maxSize;
	}
	

private:
	OutputHandler(const OutputHandler& original)
	{
	}
	
	OutputHandler& operator= (const OutputHandler& other)
	{
		return *this;
	}
			
	size_t _maxSize;
	bool _newEvent;
	size_t _maxNObjects;
	size_t _nObjects;
};

}
//namespace pxl

#endif /*PXL_IO_OUTPUTHANDLER_HH*/

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
	GenericOutputHandler(ChunkWriter& writer) :
		OutputHandler(), _writer(&writer)
	{
	}

	virtual ~GenericOutputHandler()
	{
	}

	virtual ChunkWriter& getChunkWriter() throw(std::runtime_error)
	{
		if (!_writer)
			throw std::runtime_error("GenericOutputHandler::getChunkWriter(): ChunkWriter pointer invalid.");
		return *_writer;
	}

	virtual void setChunkWriter(ChunkWriter* writer)
	{
		_writer=writer;
	}

private:
	GenericOutputHandler(const GenericOutputHandler& original)
	{
	}
	
	GenericOutputHandler& operator= (const GenericOutputHandler& other)
	{
		return *this;
	}
		
		
	ChunkWriter* _writer;
};

}//namespace pxl

#endif /*PXL_IO_GENERICOUTPUTHANDLER_HH*/
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_IO_INPUTFILE_HH
#define PXL_IO_INPUTFILE_HH



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
				_reader(_stream, ChunkReader::seekable)
	{
	}

	virtual void open(const std::string& filename)
	{
		if (_stream.is_open())
			_stream.close();
		_stream.clear();
		reset();
		_stream.open(filename.c_str(), std::ios::binary);
		_stream.seekg (0, std::ios::beg);
		_reader.setStatus(ChunkReader::preHeader);
	}

	virtual void close()
	{
		if (_stream.is_open())
			_stream.close();
		reset();
		_reader.setStatus(ChunkReader::preHeader);
	}

	virtual ~InputFile()
	{
		if (_stream.is_open())
			_stream.close();
	}

	virtual ChunkReader& getChunkReader()
	{
		return _reader;
	}

	virtual bool good()
	{
		return _stream.good();
	}
	
	virtual bool eof()
	{
		return _stream.eof();
	}
	
	virtual bool bad()
	{
		return _stream.bad();
	}
	
private:
	InputFile(const InputFile& original) : _stream(),
				_reader(_stream, ChunkReader::seekable)
	{
	}
	
	InputFile& operator= (const InputFile& other)
	{
		return *this;
	}
	
	std::ifstream _stream;
	ChunkReader _reader;
};

} //namespace pxl

#endif /*PXL_IO_INPUTFILE_HH*/
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_IO_OUTPUTFILE_HH
#define PXL_IO_OUTPUTFILE_HH



namespace pxl
{

// io
/**
 This class allows the user an easy handling of the PXL output to file. Methods to write event headers,
 the PXL event class, and information chunks are offered by inheritance from the OutputHandler.
 */

class PXL_DLL_EXPORT OutputFile : public OutputHandler
{
public:

	OutputFile(const std::string& filename, size_t maxBlockSize = 1048576, size_t maxNObjects = 1000);
	virtual ~OutputFile();
	virtual void open(const std::string& filename);
	virtual void close();
	virtual ChunkWriter& getChunkWriter();
		
	void setCompressionMode(char compressionMode);
	void setCompressionMode(int compressionMode);
	
private:
	
	OutputFile(const OutputFile& original);

	OutputFile& operator= (const OutputFile& other);
	
	std::ofstream _stream;
	ChunkWriter _writer;
};

}//namespace pxl

#endif /*PXL_IO_OUTPUTFILE_HH*/

#endif // PXL_IO_HH
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_HEP_HH
#define PXL_HEP_HH

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_HEP_VERTEX_HH
#define PXL_HEP_VERTEX_HH



//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_HEP_COMMON_VERTEX_HH
#define PXL_HEP_COMMON_VERTEX_HH

namespace pxl
{
/**
 * This is the common, pure virtual interface class for vertices.
 */

class PXL_DLL_EXPORT CommonVertex
{
public:
	virtual ~CommonVertex()
	{
	}
		
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

#endif /*PXL_HEP_COMMON_VERTEX_HH*/

namespace pxl
{

// pol
/**
 This class allows to store threevector and further properties of the decay vertex; see also pxl::BasicObjectData.
 */
class PXL_DLL_EXPORT Vertex : public Object, public CommonVertex
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
	explicit Vertex(const Vertex* original) :
		Object(original), _vector(original->_vector)
	{
	}

	virtual const Id& getTypeId() const
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
		Object::serialize(out);
		_vector.serialize(out);
	}

	virtual void deserialize(const InputStream &in)
	{
		Object::deserialize(in);
		_vector.deserialize(in);
	}

	/// This method grants read access to the vector.
	inline const Basic3Vector& getVector() const
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
	inline void setVector(const Basic3Vector& vector)
	{
		_vector = vector;
	}

	inline void addXYZ(double x, double y, double z)
	{
		_vector.setX(x + _vector.getX());
		_vector.setY(y + _vector.getY());
		_vector.setZ(z + _vector.getZ());
	}

	inline void addVector(const Basic3Vector& vector)
	{
		_vector+=vector;
	}

	inline void addVertex(const Vertex* vx)
	{
		_vector += vx->getVector();
	}

	/// This method adds the vector of \p vxd.
	inline const Vertex& operator+=(const Vertex& vx)
	{
		_vector += vx._vector;
		return *this;
	}
	/// This method subtracts the vector of \p vxd.
	inline const Vertex& operator-=(const Vertex& vx)
	{
		_vector -= vx._vector;
		return *this;
	}

	virtual Serializable* clone() const
	{
		return new Vertex(*this);
	}

	virtual std::ostream& print(int level = 1, std::ostream& os = std::cout,
			int pan = 0) const;

private:
	Basic3Vector _vector;

	Vertex& operator=(const Vertex& original)
	{
		return *this;
	}
};

// non-member operators
PXL_DLL_EXPORT bool const operator==(const Vertex& obj1, const Vertex& obj2);
PXL_DLL_EXPORT bool const operator!=(const Vertex& obj1, const Vertex& obj2);

// typedefs

/// This typedef defines a weak pointer for pxl::Vertex
typedef weak_ptr<Vertex> VertexWkPtr;

} // namespace pxl

#endif // PXL_HEP_VERTEX_HH
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_HEP_PARTICLE_HH
#define PXL_HEP_PARTICLE_HH


//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_HEP_COMMON_PARTICLE_HH
#define PXL_HEP_COMMON_PARTICLE_HH


namespace pxl
{
/**
 * This is the common, pure virtual interface class for particles.
 * The physical representation, (px, py, pz, E) or (pt, eta, phi, m/E), can
 * be chosen by the concrete implementation.
 */
class CommonParticle;
	
template<class objecttype>
class PXL_DLL_EXPORT ParticleHelper
{
public:
	ParticleHelper() : object(0), particleRef(0)
	{
	}

	virtual ~ParticleHelper()
	{
		this->transferBack();
		delete object;
	}

	void connect(CommonParticle& part)
	{
		particleRef = &part;
	}


	objecttype* set()
	{
		if (!object)
			object = new objecttype;
		this->setObjectData();
		return object;
	}

private:

	void transferBack()
	{
	}

	void setObjectData()
	{
	}

	objecttype* object;
	CommonParticle* particleRef;

};
	
class PXL_DLL_EXPORT CommonParticle
{
public:
	virtual ~CommonParticle()
	{
	}
		
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
	
	/// get auto_ptr to helper class which enables using setters
	template<class objecttype>
	std::auto_ptr<ParticleHelper<objecttype> > getAs()
	{
		ParticleHelper<objecttype>* temp = new ParticleHelper<objecttype>;
		temp->connect(*this);
		return std::auto_ptr<ParticleHelper<objecttype> >(temp);
	}

};


}

#endif /*PXL_HEP_COMMON_PARTICLE_HH*/

namespace pxl
{
// pol
/**
 This class allows to store Lorentz-fourvector and further properties of particles or reconstructed objects
 such as charge, particle-id plus the inherited properties of BasicObjectData.
 */
class PXL_DLL_EXPORT Particle : public Object, public CommonParticle
{
public:
	Particle() :
		Object(), _charge(0), _particleId(0)
	{
	}

	Particle(const Particle& original) :
		Object(original), _vector(original._vector), _charge(original._charge),
				_particleId(original._particleId)
	{
	}

	explicit Particle(const Particle* original) :
		Object(original), _vector(original->_vector),
				_charge(original->_charge), _particleId(original->_particleId)
	{
	}

	virtual const Id& getTypeId() const
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
		Object::serialize(out);
		_vector.serialize(out);
		out.writeDouble(_charge);
		out.writeInt(_particleId);
	}

	virtual void deserialize(const InputStream &in)
	{
		Object::deserialize(in);
		_vector.deserialize(in);
		in.readDouble(_charge);
		in.readInt(_particleId);
	}

	/// This method grants read access to the vector.
	inline const LorentzVector& getVector() const
	{
		return _vector;
	}

	/// This method grants write access to the vector.
	inline LorentzVector& setVector()
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

	/// This method adds vector and charge of \p pad.
	inline const Particle& operator+=(const Particle& pa)
	{
		_vector += pa._vector;
		_charge += pa._charge;
		return *this;
	}

	/// This method subtracts vector and charge of of \p pad.
	inline const Particle& operator-=(const Particle& pa)
	{
		_vector -= pa._vector;
		_charge += pa._charge;
		return *this;
	}

	virtual Serializable* clone() const
	{
		return new Particle(*this);
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

	inline double getP() const
	{
		return _vector.getP();
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

	inline void setP4(const LorentzVector& vector)
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

	inline void addP4(const LorentzVector& vector)
	{
		_vector+=vector;
	}

	inline void addP4(const Particle* particle)
	{
		_vector+=particle->getVector();
	}

	inline void addParticle(const Particle* pa)
	{
		_vector += pa->getVector();
		_charge += pa->getCharge();
	}

	virtual std::ostream& print(int level = 1, std::ostream& os = std::cout,
			int pan = 0) const;

	virtual WkPtrBase* createSelfWkPtr()
	{
		return new weak_ptr<Particle>(this);
	}

private:
	LorentzVector _vector;
	double _charge;
	int _particleId;
	
	Particle& operator=(const Particle& original)
	{
		return *this;
	}

};

// non-member operators
PXL_DLL_EXPORT bool const operator==(const Particle& obj1, const Particle& obj2);
PXL_DLL_EXPORT bool const operator!=(const Particle& obj1, const Particle& obj2);

// typedefs
/**
 This typedef represents particles and reconstructed
 objects such as muons, electrons, photons, jets; data is aggregated in ParticleData.
 */
typedef weak_ptr<Particle> ParticleWkPtr;

} // namespace pxl

#endif // PXL_HEP_PARTICLE_HH

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_HEP_COLLISION_HH
#define PXL_HEP_COLLISION_HH


namespace pxl
{
/**
 This class represents individual interactions in  multicollision events.
 It allows the separation of different collisions as they occur
 at high-rate hadron colliders by providing the relation management
 necessary to associate pxl::Vertex or pxl::Particle objects, for instance.
 */

class PXL_DLL_EXPORT Collision : public Object
{
public:
	Collision() :
		Object()
	{
	}
	/// This copy constructor provides a deep copy of the event container \p original with all data members,
	/// hep objects, and their (redirected) relations.
	Collision(const Collision& original) :
		Object(original)
	{
	}
	/// This copy constructor provides a deep copy of the event container \p original with all data members,
	/// hep objects, and their (redirected) relations.
	explicit Collision(const Collision* original) :
		Object(original)
	{
	}

	virtual WkPtrBase* createSelfWkPtr()
	{
		return new weak_ptr<Collision>(this);
	}

	virtual const Id& getTypeId() const
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
		Object::serialize(out);
	}

	virtual void deserialize(const InputStream &in)
	{
		Object::deserialize(in);
	}

	/// This virtual method creates a deep copy and returns a C++ pointer to the newly created object.  
	virtual Serializable* clone() const
	{
		return new Collision(*this);
	}
	
	virtual std::ostream& print(int level=1, std::ostream& os=std::cout, int pan=0) const;
	
private:
	Collision& operator=(const Collision& original)
	{
		return *this;
	}
};

/// This typedef defines a weak pointer for pxl::Collision
typedef weak_ptr<Collision> CollisionWkPtr;

} // namespace pxl

#endif // PXL_HEP_COLLISION_HH
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_HEP_EVENT_VIEW_hh
#define PXL_HEP_EVENT_VIEW_hh



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
 hep objects, and their (redirected) relations. 
 This way, the PXL provides a flexible generalized event container  
 meeting the needs of HEP analyses in channels with ambiguous 
 event topologies.
 */
class PXL_DLL_EXPORT EventView : public ObjectManager
{
public:
	EventView() :
		ObjectManager()
	{
	}
	/// This copy constructor provides a deep copy of the event container \p original with all data members, 
	/// hep objects, and their (redirected) relations. 
	EventView(const EventView& original) :
		ObjectManager(original)
	{
	}
	/// This copy constructor provides a deep copy of the event container \p original with all data members, 
	/// hep objects, and their (redirected) relations.
	explicit EventView(const EventView* original) :
		ObjectManager(original)
	{
	}

	virtual WkPtrBase* createSelfWkPtr()
	{
		return new weak_ptr<EventView>(this);
	}

	virtual const Id& getTypeId() const
	{
		return getStaticTypeId();
	}

	static const Id& getStaticTypeId();
	
	virtual void serialize(const OutputStream &out) const
	{
		ObjectManager::serialize(out);
	}

	virtual void deserialize(const InputStream &in)
	{
		ObjectManager::deserialize(in);
	}

	virtual Serializable* clone() const
	{
		return new EventView(*this);
	}
	
	virtual std::ostream& print(int level=0, std::ostream& os=std::cout, int pan=1) const;

private:
		
	EventView& operator=(const EventView& original)
	{
		return *this;
	}
};

/// This typedef defines a weak pointer for pxl::EventView
typedef weak_ptr<EventView> EventViewWkPtr;

} // namespace pxl


#endif // PXL_HEP_EVENT_VIEW_hh
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_HEP_ANALYSIS_PROCESS_HH
#define PXL_HEP_ANALYSIS_PROCESS_HH


namespace pxl
{
// pol

/**
 This class is designed as a base class to assist the analyzer in the
 evolution of different combinatorial hypotheses of an event according
 to a certain hep process; data is aggregated in pxl::AnalysisProcessData.
 This class provides virtual methods to be called at the beginning and end
 of a job, at the beginning and end of a run, and, of course, at event analysis
 and event finishing time (just as needed in a stand-alone analysis framework,
 for instance). When inheriting from this class, the analyst can
 place user code in the according reimplementations of these methods.
 */
class PXL_DLL_EXPORT AnalysisProcess : public ObjectManager
{
public:
	AnalysisProcess() :
		ObjectManager()
	{
	}
	AnalysisProcess(const AnalysisProcess& original) :
		ObjectManager(original)
	{
	}
	explicit AnalysisProcess(const AnalysisProcess* original) :
		ObjectManager(original)
	{
	}
	virtual ~AnalysisProcess()
	{
	}

	inline virtual const Id& getTypeId() const
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
		ObjectManager::serialize(out);
	}

	virtual void deserialize(const InputStream &in)
	{
		ObjectManager::deserialize(in);
	}

	/// This method can be reimplemented to build/destroy
	/// a static template for user-defined tree creation. \p mode is a freely usable parameter.
	virtual void buildTemplate(int mode = 0)
	{
	}

	/// This method can be reimplemented to hold hep analysis code executed at the begin of a computing job
	/// (as needed for histogram booking etc.).
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).
	virtual void beginJob(const ObjectOwner* input = 0)
	{
	}
	/// This method can be reimplemented to hold hep analysis code executed at the begin of a run.
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).
	virtual void beginRun(const ObjectOwner* input = 0)
	{
	}
	/// This method can be reimplemented to hold hep analysis code executed for the actual event analysis.
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).
	virtual void analyseEvent(const ObjectOwner* input = 0)
	{
	}
	/// This method can be reimplemented to hold hep analysis code executed at the end of each event.
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).
	/// By default, this method clears the object owner and deletes all owned objects.
	virtual void finishEvent(const ObjectOwner* input = 0)
	{
		clearObjects();
	}
	/// This method can be reimplemented to hold hep analysis code executed at the end of a run.
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).
	virtual void endRun(const ObjectOwner* input = 0)
	{
	}
	/// This method can be reimplemented to hold hep analysis code executed at the end of a computing job
	/// (as needed for histogram storing etc.).
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).
	virtual void endJob(const ObjectOwner* input = 0)
	{
	}

	virtual Serializable* clone() const;

	virtual WkPtrBase* createSelfWkPtr()
	{
		return new weak_ptr<AnalysisProcess>(this);
	}

	virtual std::ostream& print(int level=1, std::ostream& os=std::cout, int pan=0) const;

private:
	AnalysisProcess& operator=(const AnalysisProcess& original)
	{
		return *this;
	}
};

/// This typedef defines a weak pointer for pxl::AnalysisProcess
typedef weak_ptr<AnalysisProcess> AnalysisProcessWkPtr;

} // namespace pxl


#endif // PXL_HEP_ANALYSIS_PROCESS_HH
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_HEP_ANALYSIS_FORK_HH
#define PXL_HEP_ANALYSIS_FORK_HH


namespace pxl
{

// pol

/**
 This class is designed as a base class to assist the analyzer in the
 parallel evolution of different hep process hypotheses or the analysis of different
 (instrumental) aspects of an event; data is aggregated in pxl::AnalysisForkData
 */
class PXL_DLL_EXPORT AnalysisFork : public ObjectManager
{
public:
	AnalysisFork() :
		ObjectManager()
	{
	}
	AnalysisFork(const AnalysisFork& original) :
		ObjectManager(original)
	{
	}
	explicit AnalysisFork(const AnalysisFork* original) :
		ObjectManager(original)
	{
	}
	virtual ~AnalysisFork()
	{
	}

	inline virtual const Id& getTypeId() const
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
		ObjectManager::serialize(out);
	}

	virtual void deserialize(const InputStream &in)
	{
		ObjectManager::deserialize(in);
	}

	/// This method can be reimplemented to hold hep analysis code executed at the begin of a computing job.
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).
	/// By default, this method invokes the corresponding method of all managed pxl::AnalysisProcess instances.
	virtual void beginJob(const ObjectOwner* input = 0);
	/// This method can be reimplemented to hold hep analysis code executed at the begin of a run.
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).
	/// By default, this method invokes the corresponding method of all managed pxl::AnalysisProcess instances.
	virtual void beginRun(const ObjectOwner* input = 0);
	/// This method can be reimplemented to hold hep analysis code executed for the actual event analysis.
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).
	/// By default, this method invokes the corresponding method of all managed pxl::AnalysisProcess instances.
	virtual void analyseEvent(const ObjectOwner* input = 0);
	/// This method can be reimplemented to hold hep analysis code executed at the end of each event.
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).
	/// By default, this method invokes the corresponding method of all managed pxl::AnalysisProcess instances.
	virtual void finishEvent(const ObjectOwner* input = 0);
	/// This method can be reimplemented to hold hep analysis code executed at the end of a run.
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).
	/// By default, this method invokes the corresponding method of all managed pxl::AnalysisProcess instances.
	virtual void endRun(const ObjectOwner* input = 0);
	/// This method can be reimplemented to hold hep analysis code executed at the end of a computing job
	/// (as needed for histogram storing etc.).
	/// The optional parameter \p input is a const pointer to a pxl::ObjectOwner or pxl::ObjectOwner instance (that might carry the reconstructed event data or generator information).
	/// By default, this method invokes the corresponding method of all managed pxl::AnalysisProcess instances.
	virtual void endJob(const ObjectOwner* input = 0);

	virtual Serializable* clone() const;

	virtual std::ostream& print(int level=1, std::ostream& os=std::cout, int pan=0) const;

	virtual WkPtrBase* createSelfWkPtr()
	{
		return new weak_ptr<AnalysisFork>(this);
	}
private:
	AnalysisFork& operator=(const AnalysisFork& original)
	{
		return *this;
	}
};

//
/// This typedef defines a weak pointer for pxl::AnalysisFork
typedef weak_ptr<AnalysisFork> AnalysisForkWkPtr;

} // namespace pxl

#endif // PXL_HEP_ANALYSIS_FORK_HH
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_HEP_PARTICLE_FILTER_HH
#define PXL_HEP_PARTICLE_FILTER_HH



namespace pxl
{

/** 
 This class provides a pT-sorted filter for PXL hep objects (requires the name, pT and |eta| to match).
 */
class ParticlePtComparator : public Comparator<Particle>
{
public:
	virtual bool operator()(const Particle* p1, const Particle* p2)
	{
		return (p1->getPt()>p2->getPt());
	}
};

class ParticlePtEtaNameCriterion : public FilterCriterion<Particle>
{
public:
	ParticlePtEtaNameCriterion(const std::string& name, double ptMin = 0.0,
			double etaMax = 0.0) :
		_name(name), _ptMin(ptMin), _etaMax(etaMax)
	{
	}

	virtual bool operator()(const Particle& pa) const
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

class ParticlePtCriterion : public FilterCriterion<Particle>
{
public:
	ParticlePtCriterion(double ptMin = 0.0) :
		_ptMin(ptMin)
	{
	}

	virtual bool operator()(const Particle& pa) const
	{
		if ((_ptMin > 0.0 && pa.getPt() < _ptMin))
			return false;
		return true;
	}

private:
	double _ptMin;
};

typedef Filter<Particle, ParticlePtComparator> ParticleFilter;

} // namespace pxl

#endif // PXL_HEP_PARTICLE_FILTER_HH

#endif // PXL_HEP_HH
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef  astro_hh
#define  astro_hh

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_ASTRO_ASTRO_BASIC_OBJECT_HH
#define PXL_ASTRO_ASTRO_BASIC_OBJECT_HH

//#include <time.h>
//#include <string>
//#include <cmath>
//#include <vector>
//#include <stdexcept>
//#include <cfloat>


namespace pxl
{
class BasicContainer;
/**
 This Class contains basic astrophysical objects, observed at a
 direction at a given time. The direction is stored as a 3Vector. For conveniece, access to the 3Vector is provided in two different flavours of spherical coordinates: Latitude (0 at the equator) and Longitude, and Zenith (pi/2 at the 'equator') and Azimuth.
 Errors on the observation are stored as errors on the angles.
 */

class PXL_DLL_EXPORT AstroBasicObject : public Serializable
{
public:
	AstroBasicObject() : Serializable(), _time(0),_vector(0,0,0),_dphi(0),_dtheta(0)
	{
	}

	AstroBasicObject(const AstroBasicObject& original) :
	Serializable(), _time(original._time),
	_vector(original._vector),_dphi(original._dphi), _dtheta(original._dtheta)
	{
	}

	explicit AstroBasicObject(const AstroBasicObject* original) :
	Serializable(), _time(original->_time),
	_vector(original->_vector),_dphi(original->_dphi), _dtheta(original->_dtheta)
	{
	}

	virtual ~AstroBasicObject()
	{
	}

	void setAzimuthZenith(double azimuth,double zenith);

	void setAzimuth(double azimuth);

	void setZenith(double zenith);

	inline void setAzimuthError(double deltaphi)
	{
		_dphi = deltaphi;
	}

	inline const double getAzimuthError()
	{
		return _dphi;
	}

	inline const double getZenithError() const
	{
		return _dtheta;
	}

	inline void setZenithError(double deltatheta)
	{
		_dtheta = deltatheta;
	}
	////////////////////////////////////////////////////////////////////////

	void setLongitudeLatitude(double longitude,double latitude);

	void setLongitude(double longitude);

	void setLatitude(double latitude);

	inline void setLongitudeError(double deltaphi)
	{
		_dphi = deltaphi;
	}

	inline const double getLatitudeError() const
	{
		return _dtheta;
	}

	inline void setLatitudeError(double deltatheta)
	{
		_dtheta = deltatheta;
	}

	inline const double getLongitudeError()
	{
		return _dphi;
	}

	inline void setTime(const time_t &time)
	{
		_time = time;
	}

	double getAzimuth() const;

	double getZenith() const;

	double getLongitude() const;

	double getLatitude() const;

	inline const Basic3Vector& getVector() const
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
	inline void setVector(const Basic3Vector& vector)
	{
		_vector = vector;
	}

	inline time_t getTime() const
	{
		return _time;
	}

	/// Returns the angular distance to given object
	double angularDistanceTo(const AstroBasicObject &obj) const;
	double angularDistanceTo(const AstroBasicObject *obj) const;

	virtual const Id& getTypeId() const
	{
		return getStaticTypeId();
	}

	static const Id& getStaticTypeId()
	{
		static const Id id("033dfe86-1a17-767a-8dd9-2877496dab5f");
		return id;
	}

	virtual void serialize(const OutputStream &out) const;

	virtual void deserialize(const InputStream &in);

	/// This virtual method creates a deep copy and returns a C++ pointer to the newly created object.
	virtual Serializable* clone() const
	{
		return new AstroBasicObject(*this);
	}

private:
	time_t _time;
	Basic3Vector _vector;
	double _dphi;
	double _dtheta;
};

}

#endif	 // ----- #ifndef ASTROBASICOBJECT_INC	-----
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_ASTRO_ASTRO_OBJECT_HH
#define PXL_ASTRO_ASTRO_OBJECT_HH


namespace pxl
{
/**
 AstroBasicObject with UserRecord and SoftRelations.
 */

class PXL_DLL_EXPORT AstroObject : public AstroBasicObject
{
public:
	AstroObject() : AstroBasicObject()
	{
	}

	/// Copy constructor. SoftRelations are not copied since they only exist
	/// between individual objects
	AstroObject(const AstroObject& original) :
	AstroBasicObject(original)
	{
		_userRecords = original._userRecords;
	}

	/// Copy constructor. SoftRelations are not copied since they only exist
	/// between individual objects
	explicit AstroObject(const AstroObject* original) :
	AstroBasicObject(original)
	{
		_userRecords = original->_userRecords;
	}

	AstroObject(const AstroBasicObject& original) :
	AstroBasicObject(original)
	{
	}

	explicit AstroObject(const AstroBasicObject* original) :
	AstroBasicObject(original)
	{
	}

	virtual ~AstroObject()
	{
	}

	virtual void serialize(const OutputStream &out) const;

	virtual void deserialize(const InputStream &in);

	virtual const Id& getTypeId() const
	{
		return getStaticTypeId();
	}

	static const Id& getStaticTypeId()
	{
		static const Id id("115a29d6-ecfe-37bc-b647-937d49705165");
		return id;
	}

	/// This virtual method creates a deep copy and returns a C++ pointer to the newly created object.
	virtual Serializable* clone() const
	{
		return new AstroObject(*this);
	}

	void linkSoft(AstroObject* astroobject, const std::string& type);
	void linkSoft(AstroObject& astroobject, const std::string& type);

	void unlinkSoft(AstroObject* astroobject, const std::string& type);
	void unlinkSoft(AstroObject& astroobject, const std::string& type);

	const SoftRelations& getSoftRelations() const
	{
		return _softRelations;
	}

	SoftRelations& setSoftRelations()
	{
		return _softRelations;
	}

	//////////////////////////////////////////////
	//User Record
	//////////////////////////////////////////////

	/// This method provides access to the user records.
	inline const UserRecord& getUserRecord() const
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

private:
	UserRecord _userRecords;
	SoftRelations _softRelations;

	// Make AstroObject private since it's not clear how copied
	// SoftRelations should be managed
	AstroObject& operator=(const AstroObject& original)
	{
		return *this;
	}
};

}
#endif
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_ASTRO_UHECR_HH
#define PXL_ASTRO_UHECR_HH


namespace pxl
{
/**
 A class for Ultra-High-Energy Cosmic Rays, with energy, mass and charge
 */

class PXL_DLL_EXPORT UHECR : public AstroObject
{
public:

	UHECR() : AstroObject(), _energy(0.), _denergy(0.), _mass(0.), _charge(0.)
	{
	}

	UHECR(const UHECR& original) :
	AstroObject(original), _energy(original._energy), _denergy(original._denergy), _mass(original._mass), _charge(original._charge)
	{
	}

	explicit UHECR(const UHECR* original) :
	AstroObject(original), _energy(original->_energy), _denergy(original->_denergy), _mass(original->_mass), _charge(original->_charge)
	{
	}

	virtual ~UHECR()
	{
	}

	virtual void serialize(const OutputStream &out) const;

	virtual void deserialize(const InputStream &in);

	virtual const Id& getTypeId() const
	{
		return getStaticTypeId();
	}

	static const Id& getStaticTypeId()
	{
		static const Id id("a4d8433f-fa90-06e9-9c28-8d784974c1f9");
		return id;
	}

	/// This virtual method creates a deep copy and returns a C++ pointer to the newly created object.
	virtual Serializable* clone() const
	{
		return new UHECR(*this);
	}

	/// Prints out some Informations about the CR
	virtual std::ostream& print(int level=1, std::ostream& os=std::cout, int pan=1) const;

	inline double getEnergyError() const
	{
		return _denergy;
	}

	inline double getEnergy() const
	{
		return _energy;
	}

	inline double getCharge() const
	{
		return _charge;
	}

	inline double getMass() const
	{
		return _mass;
	}

	inline void setEnergy(double energy)
	{
		_energy= energy;
	}

	inline void setEnergyError(double denergy)
	{
		_denergy= denergy;
	}

	inline void setMass(double mass)
	{
		_mass = mass;
	}

	inline void setCharge(double charge)
	{
		_charge = charge;
	}

private:
	double _energy;
	double _denergy;
	double _mass;
	double _charge;

	UHECR& operator=(const UHECR& original)
	{
		return *this;
	}
};
}
//namespace pxl

#endif
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_ASTRO_REGION_OF_INTEREST_HH
#define PXL_ASTRO_REGION_OF_INTEREST_HH


namespace pxl
{
/**
 RegionOfInterest are interesting regions in the sky, which can be
 named and related to other AstroObjects with SoftRelations
 */

class PXL_DLL_EXPORT RegionOfInterest: public AstroObject
{
public:

	RegionOfInterest() : AstroObject(), _coneRadius(0.0), _name()
	{
	}




	RegionOfInterest(const RegionOfInterest& original) :
	AstroObject(original), _coneRadius(original._coneRadius), _name(original._name)
	{
	}

	explicit RegionOfInterest(const RegionOfInterest* original) :
	AstroObject(original), _coneRadius(original->_coneRadius), _name(original->_name)
	{
	}

	RegionOfInterest(const AstroBasicObject& original,double coneRadius=0.0, std::string name="") :
	AstroObject(original), _coneRadius(coneRadius), _name(name)
	{
	}

	explicit RegionOfInterest(const AstroBasicObject* original, double coneRadius=0.0,std::string name = "") :
	AstroObject(original), _coneRadius(coneRadius), _name(name)
	{
	}

	virtual ~RegionOfInterest()
	{
	}

	virtual void serialize(const OutputStream &out) const;

	virtual void deserialize(const InputStream &in);

	virtual const Id& getTypeId() const
	{
		return getStaticTypeId();
	}

	static const Id& getStaticTypeId()
	{
		static const Id id("201e9d8b-561a-c349-bb65997c49e757a6");
		return id;
	}

	/// This virtual method creates a deep copy and returns a C++ pointer to the newly created object.
	virtual Serializable* clone() const
	{
		return new RegionOfInterest(*this);
	}

	/// Prints out some Informations about the CR
	virtual std::ostream& print(int level=1, std::ostream& os=std::cout, int pan=1) const;

	inline double getConeRadius() const
	{
		return _coneRadius;
	}

	inline std::string getName() const
	{
		return _name;
	}

	inline void setConeRadius(double coneRadius)
	{
		_coneRadius= coneRadius;
	}

	inline void setName(std::string name)
	{
		_name = name;
	}

	bool objectInCone(const AstroBasicObject& obj) const;
	bool objectInCone(const AstroBasicObject* obj) const;

	bool linkIfObjectInCone(AstroObject& obj,std::string name="") ;
	bool linkIfObjectInCone(AstroObject* obj,std::string name="") ;

private:
	double _coneRadius;
	std::string _name;

	RegionOfInterest& operator=(const RegionOfInterest& original)
	{
		return *this;
	}
};
}
//namespace pxl

#endif
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifndef PXL_ASTRO_SOURCE_HH
#define PXL_ASTRO_SOURCE_HH


namespace pxl
{
/**
 UHECRSource are objects in the sky, which can be
 named and related to other AstroObjects with SoftRelations,but have a distance and a  luminosity
 */

class PXL_DLL_EXPORT UHECRSource: public AstroObject
{
public:

	UHECRSource() : AstroObject(), _luminosity(0.0), _distance(0.0), _name()
	{
	}

	UHECRSource(const UHECRSource& original) :
	AstroObject(original), _luminosity(original._luminosity), _distance(original._distance),
		_name(original._name)
	{
	}

	explicit UHECRSource(const UHECRSource* original) :
	AstroObject(original),_luminosity(original->_luminosity), _distance(original->_distance),
		_name(original->_name)
	{
	}

	virtual ~UHECRSource()
	{
	}

	virtual void serialize(const OutputStream &out) const;

	virtual void deserialize(const InputStream &in);

	virtual const Id& getTypeId() const
	{
		return getStaticTypeId();
	}

	static const Id& getStaticTypeId()
	{
		static const Id id("0855b033-8b1a-c7e0-c816e0044b041bff");
		return id;
	}

	/// This virtual method creates a deep copy and returns a C++ pointer to the newly created object.
	virtual Serializable* clone() const
	{
		return new UHECRSource(*this);
	}

	/// Prints out some Informations about the CR
	virtual std::ostream& print(int level=1, std::ostream& os=std::cout, int pan=1) const;

	inline double getLuminosity() const
	{
		return _luminosity;
	}

	inline double getDistance() const
	{
		return _distance;
	}

	inline std::string getName() const
	{
		return _name;
	}

	inline void setDistance(double distance)
	{
		_distance = distance;
	}

	inline void setLuminosity(double luminosity)
	{
		_luminosity = luminosity;
	}

	inline void setName(std::string name)
	{
		_name = name;
	}


private:
	double _luminosity;
	double _distance;
	std::string _name;

	UHECRSource& operator=(const UHECRSource& original)
	{
		return *this;
	}
};
}
//namespace pxl

#endif

#endif   // ----- #ifndef ASTRO_INC  ----- 

//#include "pxl/modules.hh"
//#include "pxl/scripting.hh"
	 
#endif // pxl_hh
