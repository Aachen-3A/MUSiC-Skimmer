//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#include <string>

#include "MUSiCProject/Pxl/interface/Pxl.h"

#ifndef EPSILON
#define EPSILON 1.0e-9
#endif

#ifndef M_PI
#define M_PI        3.14159265358979323846
#endif

namespace pxl
{

bool const operator==(const Basic3Vector& obj1, const Basic3Vector& obj2)
{
	return obj1.getX() == obj2.getX() && obj1.getY() == obj2.getY()
			&& obj1.getZ() == obj2.getZ();
}

bool const operator!=(const Basic3Vector& obj1, const Basic3Vector& obj2)
{
	return obj1.getX() != obj2.getX() || obj1.getY() != obj2.getY()
			|| obj1.getZ() != obj2.getZ();
}

bool Basic3Vector::isUnityVector() const
{
	if ((std::abs(_v[0] * _v[0] + _v[1] * _v[1] + _v[2] * _v[2] - 1.0))
			<= DBL_EPSILON)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void Basic3Vector::setRhoPhi(double perp, double phi)
{
	_v[0] = std::cos(phi) * perp;
	_v[1] = std::sin(phi) * perp;
}

void Basic3Vector::setRhoPhiZ(double perp, double phi, double z)
{
	setRhoPhi(perp, phi);
	_v[2] = z;
}

void Basic3Vector::setRThetaPhi(double r, double theta, double phi)
{
	setRhoPhiZ(std::sin(theta) * r, phi, std::cos(theta) * r);
}

bool Basic3Vector::isNullPerp() const
{
	return _v[0] > -EPSILON && _v[0] < EPSILON && _v[1] > -EPSILON && _v[1]
			< EPSILON;
}

bool Basic3Vector::isNull() const
{
	return isNullPerp() && _v[2] > -EPSILON && _v[2] < EPSILON;
}

double Basic3Vector::getPerp2() const
{
	return _v[0] * _v[0] + _v[1] * _v[1];
}
double Basic3Vector::getPerp() const
{
	return std::sqrt(getPerp2());
}

double Basic3Vector::getPhi() const
{
	return isNullPerp() ? 0.0 : std::atan2(_v[1], _v[0]);
}
double Basic3Vector::getMag2() const
{
	return _v[0] * _v[0] + _v[1] * _v[1] + _v[2] * _v[2];
}
double Basic3Vector::getMag() const
{
	return std::sqrt(getMag2());
}

double Basic3Vector::getCosTheta() const
{
	double mag = getMag();
	return mag < EPSILON ? 1.0 : _v[2] / mag;
}

double Basic3Vector::getCos2Theta() const
{
	double mag2 = getMag2();
	return mag2 < EPSILON ? 1.0 : _v[2] * _v[2] / mag2;
}

double Basic3Vector::getTheta() const
{
	return isNull() ? 0.0 : std::atan2(getPerp(), _v[2]);
}

double Basic3Vector::deltaRho(const Basic3Vector* fv) const
{
	double dDtheta = deltaTheta(fv);
	double dDphi = deltaPhi(fv);
	return std::sqrt(dDtheta * dDtheta + dDphi * dDphi);
}

double Basic3Vector::deltaPhi(const Basic3Vector* fv) const
{
	double dDphi = getPhi() - fv->getPhi();
	while (dDphi > M_PI)
		dDphi -= 2 * M_PI;
	while (dDphi < -M_PI)
		dDphi += 2 * M_PI;
	return dDphi;
}

double Basic3Vector::deltaTheta(const Basic3Vector* fv) const
{
	double dDtheta = getTheta() - fv->getTheta();
	while (dDtheta > M_PI)
		dDtheta -= 2 * M_PI;
	while (dDtheta < -M_PI)
		dDtheta += 2 * M_PI;
	return dDtheta;
}

Basic3Vector operator*(double skalar, const Basic3Vector& vec)
{
	Basic3Vector out;
	out.setX(vec.getX() * skalar);
	out.setY(vec.getY() * skalar);
	out.setZ(vec.getZ() * skalar);
	return out;
}

Basic3Vector operator*(const Basic3Vector& vec, double skalar)
{
	Basic3Vector out;
	out.setX(vec.getX() * skalar);
	out.setY(vec.getY() * skalar);
	out.setZ(vec.getZ() * skalar);
	return out;
}

Basic3Vector Basic3Vector::getETheta() const
{
	Basic3Vector out;
	out.setX(cos(this->getTheta()) * cos(this->getPhi()));
	out.setY(cos(this->getTheta()) * sin(this->getPhi()));
	out.setZ(-1 * sin(this->getTheta()));
	return out;
}

Basic3Vector Basic3Vector::getEPhi() const
{
	Basic3Vector out;
	out.setX(-1 * sin(this->getPhi()));
	out.setY(cos(this->getPhi()));
	out.setZ(0);
	return out;
}

const Basic3Vector& Basic3Vector::operator=(const Basic3Vector& vec)
{
	_v[0] = vec._v[0];
	_v[1] = vec._v[1];
	_v[2] = vec._v[2];
	return *this;
}

const Basic3Vector& Basic3Vector::operator+=(const Basic3Vector& vec)
{
	_v[0] += vec._v[0];
	_v[1] += vec._v[1];
	_v[2] += vec._v[2];
	return *this;
}

const Basic3Vector& Basic3Vector::operator-=(const Basic3Vector& vec)
{
	_v[0] -= vec._v[0];
	_v[1] -= vec._v[1];
	_v[2] -= vec._v[2];
	return *this;
}

const Basic3Vector& Basic3Vector::operator*=(double skalar)
{
	_v[0] *= skalar;
	_v[1] *= skalar;
	_v[2] *= skalar;
	return *this;
}

const Basic3Vector& Basic3Vector::operator/=(double skalar)
{
	_v[0] /= skalar;
	_v[1] /= skalar;
	_v[2] /= skalar;
	return *this;
}

Basic3Vector Basic3Vector::operator/(double skalar) const
{
	Basic3Vector out;
	out.setX(_v[0] / skalar);
	out.setY(_v[1] / skalar);
	out.setZ(_v[2] / skalar);
	return out;
}

// Scalar product
double Basic3Vector::operator*(const Basic3Vector& vec) const
{
	return vec._v[0] * _v[0] + vec._v[1] * _v[1] + vec._v[2] * _v[2];
}

Basic3Vector Basic3Vector::operator+(const Basic3Vector& vec)
{
	Basic3Vector out = *this;
	out += vec;
	return out;
}

Basic3Vector Basic3Vector::operator-(const Basic3Vector& vec)
{

	Basic3Vector out = *this;
	out -= vec;
	return out;
}

void Basic3Vector::normalize()
{
	*this/=this->getMag();
}

} // namespace pxl
#undef EPSILON
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------


namespace pxl
{

void BasicContainer::init(const BasicContainer& original)
{
	std::map<Id, Serializable*> copyHistory;

	for (const_iterator iter = original._container.begin(); iter
			!= original._container.end(); ++iter)
	{

		Serializable* pOld = *iter;
		Serializable* pNew = pOld->clone();

		setObject(pNew);
		copyHistory.insert(std::pair<Id, Serializable*>(pOld->getId(), pNew));
	}

	// redirect index:
	for (ContainerIndex::const_iterator iter = original._index.begin(); iter
	!=original._index.end(); ++iter)
	{

		Serializable* pOld = iter->second;

		Serializable* pNew = 0;
		std::map<Id, Serializable*>::const_iterator found = copyHistory.find(pOld->getId());
		if (found != copyHistory.end())
		pNew = found->second;

		if (pNew)
		_index.insert(ContainerIndex::const_iterator::value_type(iter->first,
				pNew));
	}
}

void BasicContainer::clearContainer()
{
	for (iterator iter = _container.begin(); iter != _container.end(); iter++)
	{
		delete (*iter);
	}
	_container.clear();
	_index.clear();
	_uuidSearchMap.clear();
}

void BasicContainer::setObject(Serializable* item)
{
	_container.push_back(item);
	_uuidSearchMap.insert(std::pair<Id, Serializable*>(item->getId(), item));
}

void BasicContainer::remove(Serializable* item)
{
	for (ContainerIndex::const_iterator iter = _index.begin(); iter != _index.end(); iter++)
	{
		if (item == iter->second)
		_index.erase(iter->first);
	}

	_uuidSearchMap.erase(item->getId());

	for (iterator iter = _container.begin(); iter != _container.end(); iter++)
	{
		if (item == (*iter))
		{
			delete *iter;
			_container.erase(iter);
			break;
		}
	}
}

void BasicContainer::take(Serializable* item)
{
	for (ContainerIndex::const_iterator iter = _index.begin(); iter != _index.end(); iter++)
	{
		if (item == iter->second)
		_index.erase(iter->first);
	}

	_uuidSearchMap.erase(item->getId());

	for (iterator iter = _container.begin(); iter != _container.end(); iter++)
	{
		if (item == (*iter))
		{
			_container.erase(iter);
			break;
		}
	}
}

bool BasicContainer::has(const Serializable* item) const
{
	return _uuidSearchMap.find(item->getId()) != _uuidSearchMap.end();
}

void BasicContainer::serialize(const OutputStream &out) const
{
	Serializable::serialize(out);
	//write length of vector
	out.writeUnsignedInt(size());
	//serialize container objects
	for (const_iterator iter = begin(); iter!=end(); ++iter)
	{
		(*iter)->serialize(out);
	}
	//serialize Container Index
	out.writeInt(_index.size());
	for (ContainerIndex::const_iterator iter = _index.begin(); iter != _index.end(); ++iter)
	{
		out.writeString(iter->first);
		(iter->second->getId()).serialize(out);

	}
	//serialize user record
	_userRecords.serialize(out);
}

void BasicContainer::deserialize(const InputStream &in)
{
	Serializable::deserialize(in);
	std::map<Id, Serializable*> objIdMap;
	unsigned int size = 0;
	in.readUnsignedInt(size);
	//deserialize content
	for (size_t i=0; i<size; ++i)
	{
		Id typeId (in);
		Serializable* object = dynamic_cast<Serializable*>(ObjectFactory::instance().create(typeId));
		if (object)
		{
			object->deserialize(in);
			setObject(object);
			//fill temporary map to store id's
			objIdMap.insert(std::pair<Id, Serializable*>(object->getId(), object));
		}
		else
			throw std::runtime_error("BasicContainer::deserialize(const InputStream &in) : unknown object in container!");
	}

	//read index
	in.readUnsignedInt(size);
	for (size_t i=0; i<size; ++i)
	{
		std::string name;
		in.readString(name);
		Id id (in);
		_index.insert(std::pair<std::string, Serializable*>(name, objIdMap.find(id)->second));
	}
	//deserialize user record
	_userRecords.deserialize(in);

}

static ObjectProducerTemplate<BasicContainer> _BasicContainerProducer(BasicContainer::getStaticTypeId());

} // namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#include <iostream>
namespace pxl
{

static ObjectProducerTemplate<BasicMatrix> _BasicMatrixProducer(
		BasicMatrix::getStaticTypeId());

void BasicMatrix::serialize(const OutputStream &out) const
{
	Serializable::serialize(out);
	out.writeUnsignedInt(_size1);
	out.writeUnsignedInt(_size1);
	out.writeUnsignedInt(_storageType);
	for (size_t i=0; i<_size1*_size2;i++)
	{
	out.writeDouble(_data[i]);
	}
}

void BasicMatrix::deserialize(const InputStream &in)
{
	Serializable::deserialize(in);
	unsigned int dummy;
	in.readUnsignedInt(dummy);
	_size1 = dummy;
	in.readUnsignedInt(dummy);
	_size2 = dummy;
	in.readUnsignedInt(dummy);
	_storageType = (StorageOrder)dummy;
	if (_data && (!_alienarray))
		delete[] _data;
	_data = new double[_size1*_size2];
	_alienarray = false;

	for (size_t i=0; i<_size1*_size2;i++)
	{
		in.readDouble(_data[i]);
	}
}

void BasicMatrix::use(size_t size1, size_t size2, double *data)
{
	if ((_data) && (!_alienarray))
		delete[] _data;
	_alienarray = true;
	_data = data;
	_size1 = size1;
	_size2 = size2;
}


bool BasicMatrix::isRowBasedStorage() const
{
	switch(_storageType)
	{
		case ROWMAJOR:
			return true;
			break;
		case COLUMNMAJOR:
			return false;
			break;
		default:
			throw std::runtime_error("Something went very wrong! Matrix neither in Row nor Column Storage!");
	}
}


bool BasicMatrix::isColumnBasedStorage() const
{
	switch(_storageType)
	{
		case ROWMAJOR:
			return false;
			break;
		case COLUMNMAJOR:
			return true;
			break;
		default:
			throw std::runtime_error("Something went very wrong! Matrix neither in Row nor Column Storage!");
	}
}


void BasicMatrix::resize(size_t i, size_t j)
{
	if (_data)
		delete _data;
	_data = new double[i*j];
	_size1 = i;
	_size2 = j;
}


void BasicMatrix::reshape(size_t i, size_t j)
{
	if ((i*j)!=(_size1*_size2))
	{
			throw std::runtime_error("Number of elements doens't match!");
	}
	_size1 = i;
	_size2 = j;
}


const BasicMatrix& BasicMatrix::operator=(const BasicMatrix& M)
{
		_size1 = M.getSize1();
		_size2 = M.getSize2();
		_data = new double[_size1 * _size2];
		if (M.isColumnBasedStorage())
		{
			this->setColumnBasedStorage();
		}
		else
		{
			this->setRowBasedStorage();
		}
		//const double* vecdata= vec.getConstArray();
		//memcpy(_data,vecdata,_size1);
		for (size_t i=0; i<_size1*_size2;i++)
		{
			_data[i] = M.getConstArray()[i];
		}
		return *this;
}


const BasicMatrix& BasicMatrix::operator+=(const BasicMatrix& M)
{
	if ((M.getSize1() != _size1) || (M.getSize2() != _size2))
	{
		throw std::runtime_error("Size Mismatch! Only Matrices of same size can be added.");
	}
	if (this->isColumnBasedStorage() == M.isColumnBasedStorage())
	{
		for (size_t i=0; i<_size1;i++)
		{
			for (size_t j=0; j<_size2;j++)
			{
				(_data)[i*_size2+j] += M.getElement(i,j);
			}
		}
	}
	else
	{ // different storage scheme!
	for (size_t i=0; i<_size1;i++)
		{
			for (size_t j=0; j<_size2;j++)
			{
				(_data)[i*_size2 + j] += M.getElement(j,i);
			}
		}
	}
	return *this;
}


const BasicMatrix& BasicMatrix::operator-=(const BasicMatrix& M)
{
	if ((M.getSize1() != _size1) || (M.getSize2() != _size2))
	{
		throw std::runtime_error("Size Mismatch! Only Matrices of same size can be added.");
	}
	if (this->isColumnBasedStorage() == M.isColumnBasedStorage())
	{
		for (size_t i=0; i<_size1;i++)
		{
			for (size_t j=0; j<_size2;j++)
			{
				(_data)[i*_size2+j] -= M.getElement(i,j);
			}
		}
	}
	else
	{ // different storage scheme!
	for (size_t i=0; i<_size1;i++)
		{
			for (size_t j=0; j<_size2;j++)
			{
				(_data)[i*_size2+j] -= M.getElement(j,i);
			}
		}
	}
	return *this;
}


BasicMatrix BasicMatrix::operator+(const BasicMatrix& M)
{
	BasicMatrix out = *this;
	out += M;
	return out;
}


BasicMatrix BasicMatrix::operator-(const BasicMatrix& M )
{
	BasicMatrix out = *this;
	out -= M;
	return out;
}


const BasicMatrix& BasicMatrix::operator*=(double skalar)
{
	for (size_t i=0; i<_size1*_size2;i++)
	{
		(_data)[i] *= skalar;
	}
	return *this;
}


const BasicMatrix& BasicMatrix::operator/=(double skalar)
{
	for (size_t i=0; i<_size1*_size2;i++)
	{
		(_data)[i] /= skalar;
	}
	return *this;
}


const BasicMatrix BasicMatrix::operator/(double skalar) const
{
	BasicMatrix out = *this;
	out /= skalar;
	return out;
}

double& BasicMatrix::operator()(size_t i, size_t j)
{
	if ((i >= _size1) || (j >= _size2))
		throw std::runtime_error("Index out of range");
	else
	{
		return _data[i*_size2 + j];
	}
}

std::ostream& BasicMatrix::print(int level, std::ostream& os, int pan) const
{
	os.precision(2);
	os << std::scientific;
	if (this->isRowBasedStorage())
	{
		for (size_t i=0; i<_size1; i++)
		{
			for (size_t j=0; j<_size2; j++)
			{
				os << this->getElement(i,j) << "  ";
			}
			os << std::endl;
		}
	}
	else
	{ // Column Based
		for (size_t j=0; j<_size2; j++)
		{
			for (size_t i=0; i<_size1; i++)
			{
				os << this->getElement(i,j) << "  ";
			}
			os << std::endl;
		}
	}
	return os;
}

BasicMatrix operator*(const BasicMatrix &M, double skalar)
{
	BasicMatrix R(M);
	R*=skalar;
	return R;
}

BasicMatrix operator*(double skalar, const BasicMatrix &M)
{
	BasicMatrix R(M);
	R*=skalar;
	return R;
}

bool const operator==(const BasicMatrix& obj1, const BasicMatrix& obj2)
{
	if ((obj1.getSize1() != obj2.getSize1()) || (obj1.getSize2() != obj2.getSize2()))
	{
		return false;
	}
	else
	{
		bool result = true;
		size_t iter = 0;
		while (result && (iter < obj1.getSize1()*obj1.getSize2()))
		{
			result = (obj1.getConstArray()[iter] == obj2.getConstArray()[iter]);
			iter++;
		}
		return result;
	}
}

bool const operator!=(const BasicMatrix& obj1, const BasicMatrix& obj2)
{
	return !(obj1==obj2);
}

BasicNVector operator*(const BasicMatrix& M, const BasicNVector& vec)
{

	if (vec.getSize() != M.getNumberOfColumns())
	{
		throw std::runtime_error("Size Mismatch in matrix-vector multiplication");
	}

	BasicNVector out(M.getNumberOfRows());
	for (size_t i=0;i<M.getNumberOfRows();i++)
	{
		for (size_t j=0;j<M.getNumberOfColumns();j++)
		{
			if (M.isRowBasedStorage())
			{
				out(i)+=vec.getElement(j)*M.getElement(i,j);
			}
			else
			{
				out(i)+=vec.getElement(j)*M.getElement(j,i);
			}
		}
	}
	return out;
}

Basic3Vector operator*(const BasicMatrix& M, const Basic3Vector& vec)
{

	if (M.getNumberOfColumns()!=3)
	{
		throw std::runtime_error("Size Mismatch in matrix-vector multiplication");
	}

	Basic3Vector out;
	for (size_t i=0;i<M.getNumberOfRows();i++)
	{
		for (size_t j=0;j<M.getNumberOfColumns();j++)
		{
			if (M.isRowBasedStorage())
			{
				out.getArray()[i]+=vec.getElement(j)*M.getElement(i,j);
			}
			else
			{
				out.getArray()[i]+=vec.getElement(j)*M.getElement(j,i);
			}
		}
	}
	return out;
}



//
//BasicNVector operator*(const BasicMatrix& M, const BasicNVector& vec);



} // namespace pxl



//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------



namespace pxl
{

static ObjectProducerTemplate<BasicNVector> _BasicNVectorProducer(
		BasicNVector::getStaticTypeId());

void BasicNVector::serialize(const OutputStream &out) const
{
	Serializable::serialize(out);
	out.writeUnsignedInt(_size1);
	for (size_t i = 0; i < _size1; i++)
	{
		out.writeDouble(_data[i]);
	}
}

void BasicNVector::deserialize(const InputStream &in)
{
	Serializable::deserialize(in);

	unsigned int size;
	in.readUnsignedInt(size);
	_size1 = size;
	if ((_data) && (!_alienarray))
		delete[] _data;
	_data = new double[_size1];
	_alienarray = false;
	for (size_t i = 0; i < _size1; i++)
	{
		in.readDouble(_data[i]);
	}
}

void BasicNVector::use(size_t size, double *data)
{
	if ((_data) && (!_alienarray))
		delete[] _data;
	_alienarray = true;
	_data = data;
	_size1 = size;
}

double BasicNVector::operator*(const BasicNVector& vec) const
{
	if (vec.getSize() != _size1)
	{
		throw std::runtime_error(
				"Size Mismatch! There is only a skalar product defined for vectors of same size.");
	}
	else
	{
		double result = 0;
		for (size_t i = 0; i < _size1; i++)
		{
			result += (_data)[i] * vec.getElement(i);
		}
		return result;
	}
}

double& BasicNVector::operator()(size_t i)
{
	if (i < _size1)
	{
		return _data[i];
	}
	else
	{
		throw std::runtime_error("Index out of range");
	}
}

const BasicNVector BasicNVector::operator/(double skalar) const
{
	BasicNVector out = *this;
	out /= skalar;
	return out;
}

const BasicNVector& BasicNVector::operator/=(double skalar)
{
	for (size_t i = 0; i < _size1; i++)
	{
		(_data)[i] /= skalar;
	}
	return *this;
}

const BasicNVector& BasicNVector::operator-=(const BasicNVector& vec)
{
	if (vec.getSize() != _size1)
	{
		throw std::runtime_error(
				"Size Mismatch! Only vectors of same size can be substracted.");

	}
	else
	{
		for (size_t i = 0; i < _size1; i++)
		{
			(_data)[i] -= vec.getElement(i);
		}
		return *this;
	}
}

const BasicNVector& BasicNVector::operator*=(double skalar)
{
	for (size_t i = 0; i < _size1; i++)
	{
		(_data)[i] *= skalar;
	}
	return *this;
}

const BasicNVector& BasicNVector::operator+=(const BasicNVector& vec)
{
	if (vec.getSize() != _size1)
	{
		throw std::runtime_error(
				"Size Mismatch! Only vectors of same size can be added.");

	}
	else
	{
		for (size_t i = 0; i < _size1; i++)
		{
			(_data)[i] += vec.getElement(i);
		}
		return *this;
	}
}

BasicNVector BasicNVector::operator+(const BasicNVector& vec)
{
	BasicNVector out = *this;
	out += vec;
	return out;
}

BasicNVector BasicNVector::operator-(const BasicNVector& vec)
{

	BasicNVector out = *this;
	out -= vec;
	return out;
}

const BasicNVector& BasicNVector::operator=(const BasicNVector& vec)
{
	_size1 = vec.getSize();
	_data = new double[_size1];

	//const double* vecdata= vec.getConstArray();
	//memcpy(_data,vecdata,_size1);
	for (size_t i = 0; i < vec.getSize(); i++)
	{
		_data[i] = vec.getElement(i);
	}
	return *this;
}

BasicNVector operator*(double skalar, const BasicNVector& vec)
{
	BasicNVector out = vec;
	out *= skalar;
	return out;
}

BasicNVector operator*(const BasicNVector& vec, double skalar)
{
	BasicNVector out = vec;
	out *= skalar;
	return out;
}

bool const operator==(const BasicNVector& obj1, const BasicNVector& obj2)
{
	if (obj1.getSize() != obj2.getSize())
	{
		return false;
	}
	else
	{
		bool result = true;
		size_t iter = 0;
		while (result && (iter < obj1.getSize()))
		{
			result = (obj1.getElement(iter) == obj2.getElement(iter));
			iter++;
		}
		return result;
	}
}

bool const operator!=(const BasicNVector& obj1, const BasicNVector& obj2)
{
	return !(obj1==obj2);
}

std::ostream& BasicNVector::print(int level, std::ostream& os, int pan) const
{
	os.precision(2);
	os << std::scientific;
	for (size_t i=0; i<_size1; i++)
	{
			os << this->getElement(i) << std::endl;
	}
	return os;
}



} // namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------


namespace pxl
{

std::ostream& Event::print(int level, std::ostream& os, int pan) const
{
	os << "Event." << std::endl;

	if (level>0)
		_userRecords.print(level, os, pan);

	for (ObjectOwner::const_iterator iter = _objects.begin(); iter
			!=_objects.end(); ++iter)
	{
		if ((*iter)->getMotherRelations().size() == 0)
			(*iter)->printDecayTree(0, os, pan);
	}
	return os;
}

const Id& Event::getStaticTypeId()
{
	static const Id id("c95f7434-2481-45e2-91dc-9baff0669bb3");
	return id;
}

static ObjectProducerTemplate<Event>
		_EventProducer(Event::getStaticTypeId());

} // namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#ifdef WIN32
#else
	#include <sys/time.h>
#endif


namespace pxl {

double getCpuTime()
{
#ifdef WIN32
	return timeGetTime() / 1000.;
#else
	struct timeval  tv;
	gettimeofday( &tv, NULL );
	return ( (double) tv.tv_sec + (double) tv.tv_usec / 1000000.0 );
#endif
}

} // namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------


namespace pxl {

Id::Id(const char* id)
{
	// reset
	for (int j = 0; j < 16; j++)
		bytes[j] = 0;

	// read 32 char = 16 bytes
	unsigned char first = 0, second = 0;
	int byte = 0;
	const char* source = id;
	unsigned char* target = &first;

	while (*source != 0 && byte < 16)
	{
		// find next a valid character
		if (*source >= '0' && *source <= '9')
		{
			*target = *source -'0';
		}
		else if (*source >= 'a' && *source <= 'f')
		{
			*target = *source -'a' + 10;
		}
		else if (*source >= 'A' && *source <= 'F')
		{
			*target = *source -'A' + 10;
		}
		else
		{
			source++;
			continue;
		}

		// 
		if (target == &first)
			target = &second;
		else
		{
			bytes[byte] = ((first << 4) | second);
			byte++;
			target = &first;
		}

		source++;
	}
}

Id::Id(const std::string& id)
{
	unsigned char first, second;
	int byte = 0;
	unsigned char* target = &first;
	std::string::const_iterator source = id.begin();
	std::string::const_iterator end = id.end();

	while (source != end && byte < 16)
	{
		if (*source >= '0' && *source <= '9')
		{
			*target = *source -'0';
		}
		else if (*source >= 'a' && *source <= 'f')
		{
			*target = *source -'a' + 10;
		}
		else if (*source >= 'A' && *source <= 'F')
		{
			*target = *source -'A' + 10;
		}
		else
		{
			source++;
			continue;
		}

		if (target == &first)
			target = &second;
		else
		{
			bytes[byte] = ((first << 4) | second);
			byte++;
			target = &first;
		}

		source++;
	}
}


std::string Id::toString() const
{
	static const char hex[] =
	{ '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd',
			'e', 'f' };

	std::string id;
	id.reserve(36);

	for (int i = 0; i < 4; i++)
	{
		id += hex[bytes[i] >> 4];
		id += hex[bytes[i] & 0x0F];
	}
	id += '-';
	for (int i = 4; i < 6; i++)
	{
		id += hex[bytes[i] >> 4];
		id += hex[bytes[i] & 0x0F];
	}
	id += '-';
	for (int i = 6; i < 8; i++)
	{
		id += hex[bytes[i] >> 4];
		id += hex[bytes[i] & 0x0F];
	}
	id += '-';
	for (int i = 8; i < 10; i++)
	{
		id += hex[bytes[i] >> 4];
		id += hex[bytes[i] & 0x0F];
	}
	id += '-';
	for (int i = 10; i < 16; i++)
	{
		id += hex[bytes[i] >> 4];
		id += hex[bytes[i] & 0x0F];
	}

	return id;
}

} // namespace pxl

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------


static pxl::ObjectProducerTemplate<pxl::InformationChunk> _InformationChunkProducer(pxl::InformationChunk::getStaticTypeId());
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------


namespace pxl
{

void Layout::serialize(const OutputStream &out) const
{
	out.writeInt(_type);
	out.writeInt(_style);
	out.writeInt(_color);
	out.writeDouble(_a);
	out.writeDouble(_b);
	out.writeDouble(_c);
	out.writeDouble(_d);
}
void Layout::deserialize(const InputStream &in)
{
	in.readInt(_type);
	in.readInt(_style);
	in.readInt(_color);
	in.readDouble(_a);
	in.readDouble(_b);
	in.readDouble(_c);
	in.readDouble(_d);
}

} //namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------



namespace pxl
{

bool const operator==(const LorentzVector& obj1, const LorentzVector& obj2)
{
    return obj1.getX() == obj2.getX() && obj1.getY() == obj2.getY() && obj1.getZ() == obj2.getZ() && obj1.getE()
            == obj2.getE();
}

bool const operator!=(const LorentzVector& obj1, const LorentzVector& obj2)
{
    return obj1.getX() != obj2.getX() || obj1.getY() != obj2.getY() || obj1.getZ() != obj2.getZ() || obj1.getE()
            != obj2.getE();
}

} // namespace pxl

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------


std::ostream& pxl::Object::print(int level, std::ostream& os, int pan) const
{
	printPan1st(os, pan);
	os << "pxl::Object with name: " << getName() << "\n";
	if (level>0) printContent(level, os, pan);
	return os;
}

std::ostream& pxl::Object::printContent(int level, std::ostream& os, int pan) const
{
	os << "Workflag: " << _workflag << ", locked: " << _locked;
	os << std::endl;
	_userRecords.print(level, os, pan);
	return os;
}

static pxl::ObjectProducerTemplate<pxl::Object> _ObjectProducer(pxl::Object::getStaticTypeId());
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------


static pxl::ObjectProducerTemplate<pxl::ObjectManager>
		_ObjectManagerProducer(pxl::ObjectManager::getStaticTypeId());
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------



namespace pxl
{

void ObjectOwner::init(const ObjectOwner& original)
{
	// copy objects: loop in STL style
	for (const_iterator iter = original._container.begin(); iter
			!= original._container.end(); iter++)
	{

		Relative* pOld = *iter;
		Relative* pNew = dynamic_cast<Relative*>(pOld->clone());

		set(pNew);

		std::map<Id, Relative*>::iterator insertPos =
				_copyHistory.lower_bound(pOld->id());
		if (insertPos == _copyHistory.end() || insertPos->first != pOld->id())
			_copyHistory.insert(insertPos,
					std::map<Id, Relative*>::iterator::value_type(pOld->id(), pNew));
		else
			insertPos->second = pNew;
	}

	// FIXME: possibly inefficient, might be done all in one loop
	// redirect relations: loop in PTL style
	for (const_iterator iter = original._container.begin(); iter
			!=original._container.end(); iter++)
	{
		Relative* pOld = *iter;
		Relative* pNew = 0;
		std::map<Id, Relative*>::const_iterator found = _copyHistory.find(pOld->id());
		if (found != _copyHistory.end())
			pNew = found->second;

		// mother relations
		for (Relations::const_iterator iter = pOld->getMotherRelations().begin(); iter!=pOld->getMotherRelations().end(); ++iter)
		{
			Relative* pOldRel = *iter;
			Relative* pNewRel = 0;
			std::map<Id, Relative*>::const_iterator foundRel =
					_copyHistory.find(pOldRel->id());
			if (foundRel != _copyHistory.end())
				pNewRel = foundRel->second;

			if (pOldRel)
			{
				if (pNewRel)
					pNew->linkMother(pNewRel);
				else
					// FIXME: cerr again?
					std::cerr
							<< "ObjectOwner::ObjectOwner(...): WARNING: some original objects had relations to objects of other owners."
							<< std::endl;
			}
			else
				std::cerr
						<< "ObjectOwner::ObjectOwner(...): WARNING: some originally related objects no longer exist."
						<< std::endl;
		}

		// daughter relations
		// have been set automatically above
	}

	// redirect index:
	for (std::map<std::string, Relative*>::const_iterator iter = original._index.begin(); iter
			!=original._index.end(); ++iter)
	{

		Relative* pOld = iter->second;

		Relative* pNew = 0;
		std::map<Id, Relative*>::const_iterator found = _copyHistory.find(pOld->id());
		if (found != _copyHistory.end())
			pNew = found->second;

		if (pNew)
			_index.insert(std::map<std::string, Relative*>::const_iterator::value_type(iter->first,
					pNew));
		else
			std::cerr
					<< "pxl::ObjectOwner::ObjectOwner(...): WARNING: some original indices pointed to objects of other owners."
					<< std::endl;
	}
}

void ObjectOwner::clearContainer()
{
	for (iterator iter = _container.begin(); iter != _container.end(); iter++)
	{
		(*iter)->_refObjectOwner=0;
		delete (*iter);
	}
	_container.clear();
	_copyHistory.clear();
	_index.clear();
	_uuidSearchMap.clear();
}

void ObjectOwner::set(Relative* item)
{
	item->_refObjectOwner = this;
	_container.push_back(item);
	_uuidSearchMap.insert(std::pair<Id, Relative*>(item->getId(), item));
}

void ObjectOwner::remove(Relative* item)
{
	// search & remove possible indices (multiple occurrences possible!)
	for (std::map<std::string, Relative*>::const_iterator iter = _index.begin(); iter != _index.end(); iter++)
	{
		if (item == iter->second)
			_index.erase(iter->first);
	}

	// search & remove possible copy history
	for (std::map<Id, Relative*>::const_iterator iter = _copyHistory.begin(); iter
	!= _copyHistory.end(); iter++)
	{

		// FIXME: inefficient
		if (item == iter->second)
		{
			_copyHistory.erase(iter->first);
			break; // multiple occurrences *not* possible!
		}
	}

	_uuidSearchMap.erase(item->getId());

	item->_refObjectOwner=0;
	for (iterator iter = _container.begin(); iter != _container.end(); iter++)
	{

		if (item == (*iter))
		{
			delete *iter;
			_container.erase(iter);
			break;
		}
	}

}

void ObjectOwner::take(Relative* item)
{
	// search & remove possible indices (multiple occurrences possible!)
	for (std::map<std::string, Relative*>::const_iterator iter = _index.begin(); iter != _index.end(); iter++)
	{
		if (item == iter->second)
			_index.erase(iter->first);
	}

	// search & remove possible copy history
	for (std::map<Id, Relative*>::const_iterator iter = _copyHistory.begin(); iter
	!= _copyHistory.end(); iter++)
	{

		// FIXME: inefficient
		if (item == iter->second)
		{
			_copyHistory.erase(iter->first);
			break; // multiple occurrences *not* possible!
		}
	}

	_uuidSearchMap.erase(item->getId());

	item->_refObjectOwner=0;
	for (iterator iter = _container.begin(); iter != _container.end(); iter++)
	{

		if (item == (*iter))
		{
			_container.erase(iter);
			break;
		}
	}

}

bool ObjectOwner::has(const Relative* item) const
{
	return item->_refObjectOwner == this;
}

bool ObjectOwner::setIndex(const std::string& idx, Relative* obj,
bool overwrite)
{
	if (!idx.length() || !has(obj))
	{
		if (!idx.length())
		std::cerr << "Error in setting index: id has zero length!"
		<< std::endl;
		else
		std::cerr
		<< "Error in setting index: Object does not belong to ObjectOwner."
		<< std::endl;
		return false;
	}

	std::map<std::string, Relative*>::iterator insertPos = _index.lower_bound(idx);
	if (insertPos == _index.end() || insertPos->first != idx)
	_index.insert(insertPos, std::map<std::string, Relative*>::iterator::value_type(idx, obj));
	else
	{
		if (overwrite)
		insertPos->second = obj;
		else
		{
			std::cerr << "Warning in setting Index: Key "<< idx
			<< " already present and bool 'overwrite' set to false."
			<< std::endl;
			return false;
		}
	}
	return true;
}


void ObjectOwner::serialize(const OutputStream &out) const
{
	// contents of vector
	out.writeUnsignedInt(size());
	for (const_iterator iter = begin(); iter!=end(); ++iter)
	{
		(*iter)->serialize(out);

		//serialize relations explicitly
		(*iter)->getMotherRelations().serialize(out);
		(*iter)->getDaughterRelations().serialize(out);
		(*iter)->getFlatRelations().serialize(out);
	}

	// index
	out.writeInt(_index.size());
	for (std::map<std::string, Relative*>::const_iterator iter = _index.begin(); iter != _index.end(); ++iter)
	{
		out.writeString(iter->first);
		(iter->second->id()).serialize(out);
	}

}

void ObjectOwner::deserialize(const InputStream &in)
{
	/* Algorithm:
	 * a) deserialize all objects. those deserialize themselves, and their relations
	 * as the related objects' uuids.
	 * b) create temprorary id-object map within the same loop.
	 * (no need (!) for orphan relations in this new algorithm)
	 * c) recreate relations (fetching objects from map)
	 * d) fill index (fetching objects from map)
	 * [e) shall the CopyHistory be serialised/deserialised? this is hard since
	 * one would need to have all ObjectOwners in memory, ie the whole event, and do some
	 * explicit stuff afterwards.]
	 * no error handling at the moment
	 */

	std::map<Id, Relative*> objIdMap;
	std::multimap<Relative*, Id> daughterRelationsMap;
	std::multimap<Relative*, Id> motherRelationsMap;
	std::multimap<Relative*, Id> flatRelationsMap;

	unsigned int size = 0;
	in.readUnsignedInt(size);
	for (unsigned int i=0; i<size; ++i)
	{
		Id typeId (in);
		// Contained object must be a Relative derivative
		Relative* object = dynamic_cast<Relative*>(ObjectFactory::instance().create(typeId));
		object->deserialize(in);
		set(object);
		objIdMap.insert(std::pair<Id, Relative*>(object->id(), object));

		int msize = 0;
		in.readInt(msize);
		for (int j=0; j<msize;++j)
		{
			Id id (in);
			motherRelationsMap.insert(std::pair<Relative*, Id>(object, id));
		}

		int dsize = 0;
		in.readInt(dsize);
		for (int j=0; j<dsize;++j)
		{
			Id id (in);
			daughterRelationsMap.insert(std::pair<Relative*, Id>(object, id));
		}
		int fsize = 0;
		in.readInt(fsize);
		for (int j=0; j<fsize;++j)
		{
			Id id (in);
			flatRelationsMap.insert(std::pair<Relative*, Id>(object, id));
		}

	}

	for (std::multimap<Relative*, Id>::const_iterator iter = daughterRelationsMap.begin();
	iter!=daughterRelationsMap.end(); ++iter)
	{
		Relative* target = objIdMap.find(iter->second)->second;
		iter->first->linkDaughter(target);
	}

	for (std::multimap<Relative*, Id>::const_iterator iter = flatRelationsMap.begin();
	iter!=flatRelationsMap.end(); ++iter)
	{
		Relative* target = objIdMap.find(iter->second)->second;
		iter->first->linkFlat(target);
	}

	in.readUnsignedInt(size);

	for (unsigned int i=0; i<size; ++i)
	{
		std::string name;
		in.readString(name);
		Id id (in);
		_index.insert(std::pair<std::string, Relative*>(name, objIdMap.find(id)->second));
	}
}

} // namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------


namespace pxl
{

void Relations::serialize(const OutputStream &out) const
{
	out.writeInt(size());
	for (const_iterator iter = begin(); iter != end(); ++iter)
	{
		(*iter)->getId().serialize(out);
	}
}

}
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#include <stdexcept>


namespace pxl
{

void Relative::linkMother(Relative* target) throw(std::runtime_error)
{
	if (target->_refObjectOwner != this->_refObjectOwner)
		throw std::runtime_error("pxl::ObjectBase::linkDaughter(...): WARNING: mother and daughter have not the same object holder!");

	this->_motherRelations.set(target);
	target->_daughterRelations.set(this);
}

void Relative::linkDaughter(Relative* target) throw(std::runtime_error)
{
	if (target->_refObjectOwner != this->_refObjectOwner)
		throw std::runtime_error("pxl::ObjectBase::linkMother(...): WARNING: mother and daughter have not the same object holder!");

	this->_daughterRelations.set(target);
	target->_motherRelations.set(this);
}

void Relative::linkFlat(Relative* target) throw(std::runtime_error)
{
	if (target->_refObjectOwner != this->_refObjectOwner)
		throw std::runtime_error("pxl::ObjectBase::linkFlat(...): WARNING: potential relatives have not the same object holder!");

	this->_flatRelations.set(target);
	target->_flatRelations.set(this);
}

void Relative::unlinkMother(Relative* target)
{
	this->_motherRelations.erase(target);
	target->_daughterRelations.erase(this);
}

void Relative::unlinkDaughter(Relative* target)
{
	this->_daughterRelations.erase(target);
	target->_motherRelations.erase(this);
}

void Relative::unlinkFlat(Relative* target)
{
	this->_flatRelations.erase(target);
	target->_flatRelations.erase(this);
}

void Relative::unlinkMothers()
{
	for (Relations::const_iterator iter = _motherRelations.begin(); iter
			!=_motherRelations.end(); ++iter)
	{
		(*iter)->_daughterRelations.erase(this);
	}

	_motherRelations.clearContainer();
}

void Relative::unlinkDaughters()
{
	for (Relations::const_iterator iter = _daughterRelations.begin(); iter
			!=_daughterRelations.end(); ++iter)
	{

		(*iter)->_motherRelations.erase(this);
	}

	_daughterRelations.clearContainer();
}

void Relative::unlinkFlat()
{
	for (Relations::const_iterator iter = _flatRelations.begin(); iter
			!=_flatRelations.end(); ++iter)
	{

		(*iter)->_flatRelations.erase(this);
	}

	_flatRelations.clearContainer();
}

std::ostream& Relative::printDecayTree(int level, std::ostream& os, int pan) const
{
	int daughters = 0;

	print(level, os, pan);

	for (Relations::const_iterator iter = _daughterRelations.begin(); iter
			!=_daughterRelations.end(); ++iter)
	{

		(*iter)->printDecayTree(level, os, pan + 1);
		daughters++;

	}

	if (daughters && pan > 1)
	{
		for (int p = 0; p < pan; p++)
			os << "|  ";
		os << "*" << std::endl;
	}

	return os;
}

std::ostream& Relative::print(int level, std::ostream& os, int pan) const
{
	return printPan1st(os, pan) << "pxl::ObjectBase with id: " << id()
			<< std::endl;
}

std::ostream& Relative::printPan1st(std::ostream& os, int pan) const
{
	for (int p = 0; p < pan - 2; p++)
		os << "|  ";
	if (pan - 1 > 0)
		os << "+--";
	if (pan)
		os << "{ ";
		
	return os;
}

std::ostream& Relative::printPan(std::ostream& os, int pan) const
{
	for (int p = 0; p < pan - 1; p++)
		os << "|  ";
	if (pan)
		os << "| ";
	return os;
}

} // namespace pxl

std::ostream& operator<<(std::ostream& cxxx, const pxl::Relative& obj)
{
	return obj.print(0, cxxx, 0);
}
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------



namespace pxl
{
Serializable* SoftRelations::getFirst(const ObjectOwner& owner,
		const std::string& type) const
{
	if (type == "")
	{
		const_iterator found = _relationsMap.begin();
		if (found != _relationsMap.end())
		{
			return owner.getById(found->second);
		}
		else
			return 0;
	}
	else
	{
		const_iterator found = _relationsMap.find(type);
		if (found != _relationsMap.end())
		{
			return owner.getById(found->second);
		}
		else
			return 0;
	}

}

// Gets the relatives, if living in the same basic container
int SoftRelations::getSoftRelatives(std::vector<Serializable*>& vec,
		const BasicContainer& owner, const std::string& type) const
{
	int size = vec.size();

	std::pair<const_iterator, const_iterator> iterators =
			_relationsMap.equal_range(type);
	for (const_iterator iter = iterators.first; iter != iterators.second; ++iter)
	{
		Serializable* relative = owner.getById(iter->second);
		if (relative != 0)
			vec.push_back(relative);
	}

	return vec.size() - size;
}

int SoftRelations::getSoftRelatives(std::vector<Serializable*>& vec,
		const BasicContainer& owner) const
{
	int size = vec.size();
	for (const_iterator iter = _relationsMap.begin(); iter
			!= _relationsMap.end(); ++iter)
	{
		Serializable* relative = owner.getById(iter->second);
		if (relative != 0)
			vec.push_back(relative);
	}
	return vec.size() - size;
}

Serializable* SoftRelations::getFirst(const BasicContainer& owner,
		const std::string& type) const
{
	if (type == "")
	{
		const_iterator found = _relationsMap.begin();
		if (found != _relationsMap.end())
		{
			return owner.getById(found->second);
		}
		else
			return 0;
	}
	else
	{
		const_iterator found = _relationsMap.find(type);
		if (found != _relationsMap.end())
		{
			return owner.getById(found->second);
		}
		else
			return 0;
	}

}

int SoftRelations::getSoftRelatives(std::vector<Serializable*>& vec,
		const ObjectOwner& owner, const std::string& type) const
{
	int size = vec.size();

	std::pair<const_iterator, const_iterator> iterators =
			_relationsMap.equal_range(type);
	for (const_iterator iter = iterators.first; iter != iterators.second; ++iter)
	{
		Serializable* relative = owner.getById(iter->second);
		if (relative != 0)
			vec.push_back(relative);
	}

	return vec.size() - size;
}

int SoftRelations::getSoftRelatives(std::vector<Serializable*>& vec,
		const ObjectOwner& owner) const
{
	int size = vec.size();
	for (const_iterator iter = _relationsMap.begin(); iter
			!= _relationsMap.end(); ++iter)
	{
		Serializable* relative = owner.getById(iter->second);
		if (relative != 0)
			vec.push_back(relative);
	}
	return vec.size() - size;
}

template<class objecttype>
int SoftRelations::getSoftRelativesOfType(std::vector<objecttype*>& vec,
		const ObjectOwner& owner, const std::string& type) const
{
	int size = vec.size();
	std::pair<const_iterator, const_iterator> iterators =
			_relationsMap.equal_range(type);
	for (const_iterator iter = iterators.first; iter != iterators.second; ++iter)
	{
		Serializable* relative = owner.getById(iter->second);
		if (relative != 0)
		{
			objecttype* object = dynamic_cast<objecttype*> (relative);
			if (object != 0)
				vec.push_back(object);
		}
	}
	return vec.size() - size;
}

template<class objecttype>
int SoftRelations::getSoftRelativesOfType(std::vector<objecttype*>& vec,
		const ObjectOwner& owner) const
{
	int size = vec.size();
	for (const_iterator iter = _relationsMap.begin(); iter
			!= _relationsMap.end(); ++iter)
	{
		Serializable* relative = owner.getById(iter->second);
		if (relative != 0)
		{
			objecttype* object = dynamic_cast<objecttype*> (relative);
			if (object != 0)
				vec.push_back(object);
		}
	}
	return vec.size() - size;
}

int SoftRelations::keepSoftRelatives(std::vector<Serializable*>& vec) const
{
	std::vector<Serializable*> keepItems;
	for (const_iterator iter = _relationsMap.begin(); iter
			!= _relationsMap.end(); ++iter)
	{
		for (std::vector<Serializable*>::const_iterator itVec = vec.begin(); itVec
				!= vec.end(); ++itVec)
		{
			if ((*itVec)->getId() == iter->second)
				keepItems.push_back(*itVec);
		}
	}
	vec.swap(keepItems);
	return vec.size();
}

int SoftRelations::keepSoftRelatives(std::vector<Serializable*>& vec,
		const std::string& type) const
{
	std::vector<Serializable*> keepItems;
	for (const_iterator iter = _relationsMap.begin(); iter
			!= _relationsMap.end(); ++iter)
	{
		if (iter->first == type)
		{
			for (std::vector<Serializable*>::const_iterator itVec = vec.begin(); itVec
					!= vec.end(); ++itVec)
			{
				if ((*itVec)->getId() == iter->second)
					keepItems.push_back(*itVec);
			}
		}
	}
	vec.swap(keepItems);
	return vec.size();
}

bool SoftRelations::has(const Serializable* relative) const
{
	for (const_iterator iter = _relationsMap.begin(); iter
			!= _relationsMap.end(); ++iter)
	{
		if (relative->getId() == iter->second)
			return true;
	}
	return false;
}

bool SoftRelations::has(const Serializable* relative, const std::string& type) const
{
	std::pair<const_iterator, const_iterator> iterators =
			_relationsMap.equal_range(type);
	for (const_iterator iter = iterators.first; iter != iterators.second; ++iter)
	{
		if (relative->getId() == iter->second)
			return true;
	}
	return false;
}

bool SoftRelations::has(const Id& id) const
{
	for (const_iterator iter = _relationsMap.begin(); iter
			!= _relationsMap.end(); ++iter)
	{
		if (id == iter->second)
			return true;
	}
	return false;
}

bool SoftRelations::has(const Id& id, const std::string& type) const
{
	std::pair<const_iterator, const_iterator> iterators =
			_relationsMap.equal_range(type);
	for (const_iterator iter = iterators.first; iter != iterators.second; ++iter)
	{
		if (id == iter->second)
			return true;
	}
	return false;
}

bool SoftRelations::hasType(const std::string& name) const
{
	if (_relationsMap.count(name) > 0)
		return true;
	return false;
}

void SoftRelations::set(const Serializable* relative, const std::string& type)
{
	_relationsMap.insert(std::pair<std::string, Id>(type, relative->getId()));
}

int SoftRelations::count(const std::string& name) const
{
	return _relationsMap.count(name);
}

void SoftRelations::remove(const Serializable* relative)
{
	for (iterator iter = _relationsMap.begin(); iter != _relationsMap.end(); ++iter)
	{
		if (iter->second == relative->getId())
		{
			_relationsMap.erase(iter);
			break;
		}
	}
}

void SoftRelations::remove(const Serializable* relative,
		const std::string& type)
{
	std::pair<iterator, iterator> iterators = _relationsMap.equal_range(type);
	for (iterator iter = iterators.first; iter != iterators.second; ++iter)
	{
		if (iter->second == relative->getId())
		{
			_relationsMap.erase(iter);
			break;
		}
	}
}

std::ostream& SoftRelations::print(int level, std::ostream& os, int pan) const
{
	os << "SoftRelations of size " << _relationsMap.size() << "\n";
	for (const_iterator iter = _relationsMap.begin(); iter
			!= _relationsMap.end(); ++iter)
	{
		os << "--> ('" << iter->first << "', " << iter->second << ") \n";
	}
	return os;
}

} //namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

namespace pxl
{

void UserRecord::serialize(const OutputStream &out) const
{
	out.writeUnsignedInt(getContainer()->size());
	for (const_iterator iter = getContainer()->begin(); iter
			!= getContainer()->end(); ++iter)
	{
		out.writeString(iter->first);
		Variant::Type type = iter->second.getType();

		char cType = ' ';
		switch (type)
		{
		case Variant::TYPE_BOOL:
			cType = 'b';
			out.writeChar(cType);
			out.writeBool(iter->second.get<bool> ());
			break;
		case Variant::TYPE_CHAR:
			cType = 'c';
			out.writeChar(cType);
			out.writeChar(iter->second.get<char> ());
			break;
		case Variant::TYPE_UCHAR:
			cType = 'C';
			out.writeChar(cType);
			out.writeUnsignedChar(iter->second.get<unsigned char> ());
			break;
		case Variant::TYPE_INT:
			cType = 'i';
			out.writeChar(cType);
			out.writeInt(iter->second.get<int> ());
			break;
		case Variant::TYPE_UINT:
			cType = 'I';
			out.writeChar(cType);
			out.writeUnsignedInt(iter->second.get<unsigned int> ());
			break;
		case Variant::TYPE_SHORT:
			cType = 'o';
			out.writeChar(cType);
			out.writeShort(iter->second.get<short> ());
			break;
		case Variant::TYPE_USHORT:
			cType = 'O';
			out.writeChar(cType);
			out.writeUnsignedShort(iter->second.get<unsigned short> ());
			break;
		case Variant::TYPE_LONG:
			cType = 'l';
			out.writeChar(cType);
			out.writeLong(iter->second.get<long> ());
			break;
		case Variant::TYPE_ULONG:
			cType = 'L';
			out.writeChar(cType);
			out.writeUnsignedLong(iter->second.get<unsigned long> ());
			break;
		case Variant::TYPE_DOUBLE:
			cType = 'd';
			out.writeChar(cType);
			out.writeDouble(iter->second.get<double> ());
			break;
		case Variant::TYPE_FLOAT:
			cType = 'f';
			out.writeChar(cType);
			out.writeFloat(iter->second.get<float> ());
			break;
		case Variant::TYPE_STRING:
			cType = 's';
			out.writeChar(cType);
			out.writeString(iter->second.get<std::string> ());
			break;
		case Variant::TYPE_USER:
			cType = 'u';
			out.writeChar(cType);
			break;
		//case Variant::TYPE_SERIALIZABLEPTR:
		//	cType = 'S';
		//	out.writeChar(cType);
		//	Serializable *f;
		//	f = iter->second.get<Serializable*> ();
		//	f->serialize(out);
		//	break;

		case Variant::TYPE_PTR:
			cType = 'p';
			out.writeChar(cType);
			break;
		default:
			out.writeChar(cType);
			std::cerr << "Type not handled in pxl::Variant I/O." << std::endl;
		}

	}
}

void UserRecord::deserialize(const InputStream &in)
{
	unsigned int size = 0;
	in.readUnsignedInt(size);
	for (unsigned int j = 0; j < size; ++j)
	{
		std::string name;
		in.readString(name);
		char cType;
		in.readChar(cType);
		//FIXME: temporary solution here - could also use static lookup-map,
		//but leave this unchanged until decided if to switch to new UR implementation.
		switch (cType)
		{
		case 'b':
		{
			bool b;
			in.readBool(b);
			set<bool> (name, b);
		}
			break;
		case 'c':
		{
			char c;
			in.readChar(c);
			set<char> (name, c);
		}
			break;
		case 'C':
		{
			unsigned char c;
			in.readUnsignedChar(c);
			set<unsigned char> (name, c);
		}
			break;

		case 'i':
		{
			int ii;
			in.readInt(ii);
			set<int> (name, ii);
		}
			break;
		case 'I':
		{
			unsigned int ui;
			in.readUnsignedInt(ui);
			set<unsigned int> (name, ui);
		}
			break;
		case 'o':
		{
			short s;
			in.readShort(s);
			set<short> (name, s);
		}
			break;
		case 'O':
		{
			unsigned short us;
			in.readUnsignedShort(us);
			set<unsigned short> (name, us);
		}
			break;
		case 'l':
		{
			long l;
			in.readLong(l);
			set<long> (name, l);
		}
			break;
		case 'L':
		{
			unsigned long ul;
			in.readUnsignedLong(ul);
			set<unsigned long> (name, ul);
		}
			break;
		case 'd':
		{
			double d;
			in.readDouble(d);
			set<double> (name, d);
		}
			break;
		case 'f':
		{
			float f;
			in.readFloat(f);
			set<float> (name, f);
		}
			break;
		case 's':
		{
			std::string ss;
			in.readString(ss);
			set<std::string> (name, ss);
		}
			break;
		//case 'S':
		//{
		//	Id id(in);
		//	Serializable* obj = ObjectFactory::instance().create(id);
		//	obj->deserialize(in);
		//	set<Serializable*> (name, obj);
		//	break;
		//}
		case 'u':
		case 'p':
		default:
			std::cerr << "Type " << cType
					<< " not handled in pxl::Variant I/O." << std::endl;
		}
	}
}

std::ostream& UserRecord::print(int level, std::ostream& os, int pan) const
{
	os << "UserRecord size " << size() << "\n";
	for (const_iterator iter = getContainer()->begin(); iter
			!= getContainer()->end(); ++iter)
	{
		os << "-->";
		Variant::Type type = iter->second.getType();

		os << " ('" << iter->first << "', ";
		switch (type)
		{
		case Variant::TYPE_BOOL:
			os << iter->second.get<bool> ();
			break;
		case Variant::TYPE_CHAR:
			os << iter->second.get<char> ();
			break;
		case Variant::TYPE_UCHAR:
			os << iter->second.get<unsigned char> ();
			break;
		case Variant::TYPE_INT:
			os << iter->second.get<int> ();
			break;
		case Variant::TYPE_UINT:
			os << iter->second.get<unsigned int> ();
			break;
		case Variant::TYPE_SHORT:
			os << iter->second.get<short> ();
			break;
		case Variant::TYPE_USHORT:
			os << iter->second.get<unsigned short> ();
			break;
		case Variant::TYPE_LONG:
			os << iter->second.get<long> ();
			break;
		case Variant::TYPE_ULONG:
			os << iter->second.get<unsigned long> ();
			break;
		case Variant::TYPE_DOUBLE:
			os << iter->second.get<double> ();
			break;
		case Variant::TYPE_FLOAT:
			os << iter->second.get<float> ();
			break;
		case Variant::TYPE_STRING:
			os << iter->second.get<std::string> ();
			break;
		case Variant::TYPE_USER:
			break;
		case Variant::TYPE_PTR:
			os << iter->second.get<void*> ();
			break;
		default:
			;
		}
		os << "), type: " << iter->second.getTypeName() << ".\n";
	}

	return os;
}

} //namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#include <cstddef>
#include <sstream>
#include <vector>


namespace pxl
{

std::vector<VariantBase::TypeInfo> &VariantBase::getTypes()
{
	static std::vector<VariantBase::TypeInfo> types;
	return types;
}

/// This constant array serves the PXL variant data type.
static const char *typeNames[] =
{ // Note: must match VariantBase::Type
		"null", "bool", "char", "uchar", "short", "ushort", "int", "uint",
				"long", "ulong", "float", "double", "string", 
				//"Serializable",
				"ptr" };

const VariantBase::TypeInfo& VariantBase::fallbackGetTypeInfo(Type t)
		throw (std::runtime_error)
{
	if (PXL_UNLIKELY(getTypes().begin() == getTypes().end()))
	{
		for (unsigned int i = 0; i < sizeof typeNames / sizeof *typeNames; i++)
			getTypes().push_back(typeNames[i]);
	}

	if ((std::size_t) t >= getTypes().size())
	{
		std::ostringstream ss;
		ss << "VariantBase::fallbackGetTypeInfo: " << "Variant type "
				<< (std::size_t) t << " undefined.";
		throw std::runtime_error(ss.str());
	}

	return getTypes()[t];
}

void VariantBase::wrongType(Type tShould, Type tIs) const
		throw (std::runtime_error)
{
	const TypeInfo& is = getTypeInfo(tIs);
	const TypeInfo& should = getTypeInfo(tShould);

	std::ostringstream ss;
	ss << "VariantBase::wrongType: " << "Trying to access Variant of type \""
			<< is.name << "\" as \"" << should.name << "\".";

	throw std::runtime_error(ss.str());

}

} // namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------


namespace pxl
{

void WkPtrBase::notifyDeleted()
{
	_objectRef = 0;
	if (_notifyChainOut)
		_notifyChainOut->notifyDeleted();
	_notifyChainIn = 0;
	_notifyChainOut = 0;
}

void WkPtrBase::connect(Relative* pointer)
{
	// disconnect:
	if (_objectRef)
	{
		if (_objectRef->_refWkPtrSpec == this)
			_objectRef->_refWkPtrSpec = _notifyChainOut;
		if (_notifyChainIn && _notifyChainOut)
		{
			_notifyChainIn->_notifyChainOut = _notifyChainOut;
			_notifyChainOut->_notifyChainIn = _notifyChainIn;
		}
		else
		{
			if (_notifyChainIn)
				_notifyChainIn->_notifyChainOut = 0;
			if (_notifyChainOut)
				_notifyChainOut->_notifyChainIn = 0;
		}
	}
	_notifyChainOut = 0;
	_notifyChainIn = 0;

	// connect
	if (pointer)
	{
		_notifyChainIn = 0;
		_notifyChainOut = pointer->_refWkPtrSpec;
		if (_notifyChainOut)
			_notifyChainOut->_notifyChainIn = this;
		pointer->_refWkPtrSpec = this;
	}

	_objectRef = pointer;
}

} // namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#include <zlib.h>


namespace pxl
{

bool ChunkReader::skip()
{
	if (_stream.peek()==EOF)
		return false;
	
	//skip event header
	if (_status == preHeader)
	{
		++_sectionCount;
		_stream.ignore(1);
		// read info size
		pxl::int32_t infoSize = 0;
		_stream.read((char *)&infoSize, 4);
		_stream.ignore(infoSize);
	}

	//skip all blocks
	while (nextBlockId()=='B' && !_stream.eof() )
	{
		// read info size
		pxl::int32_t infoSize = 0;
		_stream.read((char *)&infoSize, 4);
		_stream.ignore(infoSize);
		_stream.ignore(1);

		// read chunk size
		pxl::int32_t chunkSize = 0;
		_stream.read((char *)&chunkSize, 4);
		_stream.ignore(chunkSize);
	}

	_stream.ignore(4);
	_status = preHeader;

	return true;
}

bool ChunkReader::previous()
{
	if (_seekMode == nonSeekable)
		return false;

	if (_status != preHeader)
	{
		endEvent();
		previous();
	}
	std::streampos pos = _stream.tellg();

	pxl::int32_t eventSize;
	pos -= 4;
	if (pos<0)
		return false;
	_stream.seekg(pos);
	_stream.read((char*)&eventSize, 4);

	pos -= eventSize;
	if (pos<0)
		return false;
	_stream.seekg(pos);
	_status = preHeader;

	--_sectionCount;

	return true;
}

/// Reads next block from file to stream. If mode is -1, the information condition string is evaluated,
/// i.e. the block is read only if the string equals the one in the input.
bool ChunkReader::readBlock(skipMode skip, infoMode checkInfo,
		const std::string& infoCondition) throw(std::runtime_error)
{
	//if event header not read, return
	if (_status == preHeader)
		return false;
	
	//return false if end of file
	if (_stream.peek()==EOF)
		return false;

	//check if beginning of block
	char id = nextBlockId();

	if (id!='B')
	{
		if (id=='e') //end of event
		{
			_status = preHeader;
			_stream.ignore(4);
			return false;
		}
		else
		{
			std::stringstream ss;
			ss << "pxl::ChunkReader::readBlock(): Unknown char identifier: " << id;
			throw std::runtime_error(ss.str());
		}
	}

	pxl::int32_t infoSize = 0;
	_stream.read((char *)&infoSize, 4);

	bool readStream = true;

	if (checkInfo == evaluate)
	{
		char* infoBuffer = new char[infoSize+1];
		infoBuffer[infoSize]=0;
		_stream.read(infoBuffer, infoSize);
		std::string info;
		info.assign(infoBuffer);
		//the mode is set to -2 if the info condition is not fulfilled.
		//rest of block must be skipped and false be returned.
		if (infoCondition!=info)
			readStream = false;
		delete infoBuffer;
	}
	else
		_stream.ignore(infoSize);

	char compressionMode;
	_stream.read(&compressionMode, 1);

	// read chunk size
	pxl::int32_t chunkSize = 0;
	_stream.read((char *)&chunkSize, 4);
	
	if (_stream.bad() || _stream.eof())
		return false;

	if (readStream == false)
	{
		_stream.ignore(chunkSize);
		if (skip == on)
			return readBlock(skip, checkInfo, infoCondition);
		else
			return false;
	}
	else
	{
		// read chunk into buffer
		if (compressionMode==' ')
		{
			_buffer.destroy();
			_buffer.resize(chunkSize);
			_stream.read(_buffer.data(), chunkSize);
			if (_stream.bad() || _stream.eof() )
				return false;
		}
		else if (compressionMode=='Z')
		{
			_buffer.destroy();
			unzipEventData(chunkSize);
		}
		else
		{
			throw std::runtime_error("pxl::ChunkReader::readBlock(): Invalid compression mode.");
		}
	}
	return true;
}

bool ChunkReader::readHeader(readMode mode, skipMode doSkip,
		infoMode checkInfo, const std::string& infoCondition)
{
	// if position is not before the header, end the previous event
	endEvent();
	
	++_sectionCount;
	
	if (_stream.peek()==EOF || _stream.bad() )
		return false;

	switch (mode)
	{

	// we want to read an info chunk - if not, skip the event, and
	// read the next header if skip mode is on
	case infoChunk:
		if (_stream.peek()!='I')
		{
			skip();
			if (doSkip == on)
				return readHeader(mode, doSkip, checkInfo, infoCondition);
			else
				return false;
		}
		_stream.ignore(1);
		_status = infoPreBlock;
		break;

		// we want to read an event - if not, skip the info chunk (...), and
		// read the next header if skip mode is on
	case event:
		if (_stream.peek()!='E')
		{
			skip();
			if (doSkip == on)
				return readHeader(mode, doSkip, checkInfo, infoCondition);
			else
				return false;
		}
		_stream.ignore(1);
		_status = evPreBlock;
		break;

		// read what comes, and change the status accordinglys
	default:
		_status = preBlock;
		char nextId = nextBlockId();
		if (nextId=='E')
			_status = evPreBlock;
		else if (nextId=='I')
			_status = infoPreBlock;
	}

	// get size of info string
	pxl::int32_t infoSize = 0;
	_stream.read((char *)&infoSize, 4);

	//if info string is to be checked
	if (checkInfo==evaluate)
	{
		char* infoBuffer = new char[infoSize+1];
		infoBuffer[infoSize]=0;
		_stream.read(infoBuffer, infoSize);
		std::string info;
		info.assign(infoBuffer);
		delete infoBuffer;
		if (infoCondition!=info)
		{
			if (doSkip == on)
				return readHeader(mode, doSkip, checkInfo, infoCondition);
			else
				return false;
		}
	}
	else
		_stream.ignore(infoSize);

	if (_stream.eof() || _stream.bad() )
		return false;

	return true;
}

int ChunkReader::unzipEventData(int nBytes) throw(std::runtime_error)
{
	int ret, length = 0;
	pxl::int32_t have;


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
	do
	{
		int size = nBytes;
		if (size > iotl__iStreamer__lengthUnzipBuffer)
		{
			size = iotl__iStreamer__lengthUnzipBuffer;
		}

		strm.avail_in = _stream.read((char*)_inputBuffer, size).gcount();
		if (_stream.bad())
		{
			inflateEnd(&strm);
			return 0;
		}

		nBytes -= strm.avail_in;
		if (_stream.eof())
			nBytes = 0;

		if (strm.avail_in == 0)
			break;

		strm.next_in = _inputBuffer;

		// run inflate() on input until output buffer not full
		do
		{
			strm.avail_out = iotl__iStreamer__lengthUnzipBuffer;
			strm.next_out = _outputBuffer;

			ret = inflate(&strm, Z_NO_FLUSH);
			switch (ret)
			{
			case Z_STREAM_ERROR:
				throw std::runtime_error("pxl::ChunkReader::unzipEventData(): Internal inflate stream error.");
			case Z_NEED_DICT:
				ret = Z_DATA_ERROR; // fall through
			case Z_DATA_ERROR:
			case Z_MEM_ERROR:
				inflateEnd(&strm);
				return 0;
			default:
				break;
			}

			have = iotl__iStreamer__lengthUnzipBuffer - strm.avail_out;

			length += have;

			_buffer.resize(length);
			memcpy(_buffer.data()+(length-have), _outputBuffer, have);

		} while (strm.avail_out == 0);
	} while (nBytes > 0); // done when inflate() says it's done
	inflateEnd(&strm);

	return length;
}

	
/// Returns the size of the associated file.
size_t ChunkReader::getSize() const
{
      std::streampos pos = _stream.tellg();
      _stream.seekg (0, std::ios::end);

      size_t length = _stream.tellg();
      _stream.seekg (pos);

      return length;
}

}

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------



namespace pxl
{

bool ChunkWriter::newFileSection(const std::string& info, char cSectionMarker)
{
	_stream.write(&cSectionMarker, 1);
	_nBytes+=1;

	// info block:
	const char* cInfo = info.c_str();
	pxl::int32_t lengthInfo = info.length();
	_stream.write((char *) &lengthInfo, 4);
	_nBytes+=4;
	_stream.write(cInfo, lengthInfo);
	_nBytes+=lengthInfo;

	return true;
}

bool ChunkWriter::endEvent()
{
	// end event marker:
	writeFlag('e');
	_stream.write((char* ) &_nBytes, 4);
	_nBytes=0;
	return true;
}

/// Write char flag.
bool ChunkWriter::writeFlag(char cEvtMarker)
{
	_stream.write(&cEvtMarker, 1);
	_nBytes+=1;
	return true;
}
	

bool ChunkWriter::write(std::string info) throw(std::runtime_error)
{
	// write out block information
	const char* cInfo = info.c_str();
	pxl::int32_t lengthInfo = info.length();
	_stream.write((char *) &lengthInfo, 4);
	_nBytes+=4;
	_stream.write(cInfo, lengthInfo);
	_nBytes+=lengthInfo;

	// write out compression mode
	char compressed = 'Z';
	if (_compressionMode == ' ') compressed = ' ';
	_stream.write((char *) &compressed, 1);
	_nBytes+=1;

	// zip block:
	const char* cBuffer = _buffer.getData();
	pxl::int32_t lengthBuffer = _buffer.getSize();

	const char* cZip = cBuffer;
	pxl::int32_t lengthZip = lengthBuffer;

	char* cZipSpace = 0;
	unsigned long lengthZipSpace = 0;

	if (_compressionMode == ' ')
	{
		// no compression requires no action...
	}
	else if (_compressionMode >= '0' && _compressionMode <= '9')
	{
		// data compression a la Gero, i.e. compression level = 6:
		lengthZipSpace = long(double(lengthBuffer) * 1.05 + 16);
		cZipSpace = new char[lengthZipSpace];
		
		int status = compress2((Bytef*)cZipSpace, (uLongf*)&lengthZipSpace,
				(const Bytef*)cBuffer, lengthBuffer, _compressionMode - '0');
		switch (status)
		{
		case Z_MEM_ERROR:
			throw std::runtime_error("pxl::ChunkWriter::write(): zlib: not enough memory");
			break;
		case Z_BUF_ERROR:
			throw std::runtime_error("pxl::ChunkWriter::write(): zlib: buffer too small");
			break;
		case Z_STREAM_ERROR:
			throw std::runtime_error("pxl::ChunkWriter::write(): zlib: level parameter invalid");
			break;
		default:
			break;
		}

		cZip = cZipSpace;
		lengthZip = lengthZipSpace;
	}
	else
		throw std::runtime_error("pxl::FileChunkWriter::write(): Invalid compression mode.");

	_stream.write((char *) &lengthZip, 4);
	_nBytes+=4;
	_stream.write(cZip, lengthZip);
	_nBytes+=lengthZip;

	if (cZipSpace)
		delete[] cZipSpace;

	_buffer.clear();

	return true;
}

} //namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------



namespace pxl
{

/// This method fills the objects from the read-in block into the passed vector. The number of added objects is returned.
int InputHandler::readObjects(std::vector<Serializable*>& objects)
		throw(std::runtime_error)
{
	int nObjects = 0;
	while (getChunkReader().getInputStream().good())
	{
		_objectCount++;
		
		Id id(getChunkReader().getInputStream());

		Serializable* obj = ObjectFactory::instance().create(id);
		if (obj)
		{
			obj->deserialize(getChunkReader().getInputStream());
			objects.push_back(obj);
			nObjects++;
		}
		else
			throw std::runtime_error("InputHandler::readObjects(std::vector<Serializable*>& objects): unknown object in file!");
	}
	return nObjects;
}

/// This method fills the objects from the read-in block into the passed pxl::Event. The number of added objects is returned.
/// Caution: Only pxl::Relative derivatives will be restored, e.g. no pxl::Event.
int InputHandler::readObjects(Event* event) throw(std::runtime_error)
{
	int nObjects = 0;
	while (getChunkReader().getInputStream().good())
	{
		_objectCount++;
		
		Id id(getChunkReader().getInputStream());

		Relative *obj =
				dynamic_cast<Relative*>(ObjectFactory::instance().create(id));
		if (obj)
		{
			obj->deserialize(getChunkReader().getInputStream());
			event->setObject(obj);
			nObjects++;
		}
		else
			throw std::runtime_error("InputHandler::readObjects(pxl::Event* event): unknown object in file!");

	}
	return nObjects;
}

/// This method fills the objects from the read-in block into the passed pxl::BasicContainer. The number of added objects is returned.
int InputHandler::readObjects(BasicContainer* container)
		throw(std::runtime_error)
{
	int nObjects = 0;
	while (getChunkReader().getInputStream().good())
	{
		_objectCount++;
		
		Id id(getChunkReader().getInputStream());

		Serializable* obj = ObjectFactory::instance().create(id);
		if (obj)
		{
			obj->deserialize(getChunkReader().getInputStream());
			container->setObject(obj);
			nObjects++;
		}
		else
			throw std::runtime_error("InputHandler::readObjects(std::vector<pxl::Serializable*>& objects): unknown object in file!");
	}
	return nObjects;
}

Serializable* InputHandler::readNextObject() throw(std::runtime_error)
{
	while (!getChunkReader().getInputStream().good())
	{
		if (getChunkReader().getStatus()==ChunkReader::preHeader)
		{
			if (!getChunkReader().next())
			{
				return 0;
			}
		}

		getChunkReader().nextBlock();
		if (getChunkReader().eof())
			return 0;
	}
	
	_objectCount++;
	
	Id id(getChunkReader().getInputStream());
	Serializable* obj = ObjectFactory::instance().create(id);
	if (obj)
	{
		obj->deserialize(getChunkReader().getInputStream());
		return obj;
	}
	else
		throw std::runtime_error("InputHandler::readNextObject(): unknown object in file!");
	return obj;
}

Serializable* InputHandler::readPreviousObject() throw(std::runtime_error)
{
	return seekToObject(_objectCount-1);
}

Serializable* InputHandler::seekToObject(size_t index) throw(std::runtime_error)
{
	if (index <= _objectCount)
	{
///		bool success = 
		seekToFileSection(0);
		//DEACTIVATE for current VISPA version, should be fixed in IO 3.0
//		if (!success)
//			throw std::runtime_error("InputHandler::seekToObject: seekToFileSection(0) not successful");
		_objectCount = 0;
	}
		
	if (index > _objectCount)
	{
		Serializable* obj = 0;
		do
		{
			if (obj!=0)
				delete obj;
			
			obj = readNextObject();
			if (obj == 0)
				return 0;
		}
		while ( index  > _objectCount );
		return obj;
	}
	
	return 0;
}


bool InputHandler::readEvent(Event* event)
{
	if (getChunkReader().nextBlock() && getChunkReader().getInputStream().good())
	{
		Id id(getChunkReader().getInputStream());
		if (id==event->getStaticTypeId())
		{
			event->deserialize(getChunkReader().getInputStream());
			return true;
		}
	}
	return false;
}

bool InputHandler::readEventIf(Event* event, const std::string& blockInfo,
		skipMode doSkip)
{
	if (readBlockIf(blockInfo, doSkip))
	{
		Id id(getChunkReader().getInputStream());
		if (id==event->getStaticTypeId())
		{
			event->deserialize(getChunkReader().getInputStream());
			return true;
		}
	}
	return false;
}

bool InputHandler::readBasicContainer(BasicContainer* basicContainer)
{
	if (getChunkReader().nextBlock() && getChunkReader().getInputStream().good())
	{
		Id id(getChunkReader().getInputStream());
		if (id==basicContainer->getStaticTypeId())
		{
			basicContainer->deserialize(getChunkReader().getInputStream());
			return true;
		}
	}
	return false;
}

bool InputHandler::readBasicContainerIf(BasicContainer* basicContainer, const std::string& blockInfo,
			skipMode doSkip)
{
	if (readBlockIf(blockInfo, doSkip))
	{
		Id id(getChunkReader().getInputStream());
		if (id==basicContainer->getStaticTypeId())
		{
			basicContainer->deserialize(getChunkReader().getInputStream());
			return true;
		}
	}
	return false;
}

bool InputHandler::readInformationChunk(InformationChunk* chunk)
{
	if (getChunkReader().nextBlock() && getChunkReader().getInputStream().good())
	{
		Id id(getChunkReader().getInputStream());
		if (id==chunk->getStaticTypeId())
		{
			chunk->deserialize(getChunkReader().getInputStream());
			return true;
		}
	}
	return false;
}

bool InputHandler::readInformationChunkIf(InformationChunk* event,
		const std::string& blockInfo, skipMode doSkip)
{
	if (readBlockIf(blockInfo, doSkip))
	{
		Id id(getChunkReader().getInputStream());
		if (id==event->getStaticTypeId())
		{
			event->deserialize(getChunkReader().getInputStream());
			return true;
		}
	}
	return false;
}

int InputHandler::skipFileSections(int n)
{
	int skipped = 0;
	for (; n<0 && getChunkReader().previous(); ++n)
	--skipped;
	for (; n>0 && getChunkReader().skip(); --n)
	++skipped;
	return skipped;
}

bool InputHandler::seekToFileSection(int index)
{
	int count = index - getChunkReader().getSectionCount();
	return (skipFileSections(count) == count);
}

bool InputHandler::readBlockIf(const std::string& blockInfo,
		skipMode doSkip)
{
	bool success = false;
	while (!success && getChunkReader().getStatus()!=0)
	{
		success = getChunkReader().nextBlock(doSkip, ChunkReader::evaluate, blockInfo);
	}
	return success;
}

/// This method explicitly reads an object of type objecttype. Caution: This method should only be used if the type of the following object is known by hard.
template<class objecttype> bool InputHandler::readObject(objecttype* obj) throw(std::runtime_error)
{
	if (getChunkReader().getInputStream().good())
	{
		Id id(getChunkReader().getInputStream());
		if (id!=objecttype::getStaticTypeId())
		throw std::runtime_error("InputHandler::readObject(objecttype* obj): unexpected object in file!");
		obj->deserialize(getChunkReader().getInputStream());
		return true;
	}
	return false;
}

}
//namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------


namespace pxl {

ObjectFactory::ObjectFactory()
{
}
	
ObjectFactory& ObjectFactory::instance()
{
	static ObjectFactory f;
	return f;
}
		
Serializable *ObjectFactory::create (const Id& id)
{
	std::map<Id,const ObjectProducerInterface *>::iterator result;
	result = instance()._Producers.find (id);
	if (result == instance()._Producers.end())
		return 0;
	else
		return (*result).second->create();
}
	
void ObjectFactory::registerClass (const Id& id, const ObjectProducerInterface* producer)
{
	instance()._Producers[id] = producer;
}

} // namespace pxl

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------



namespace pxl {
OutputFile::OutputFile(const std::string& filename, size_t maxBlockSize, size_t maxNObjects) :
	OutputHandler(maxBlockSize, maxNObjects), _stream(filename.c_str(),
			std::ios::binary), _writer(_stream)
{
	if (_stream.good() == false)
		throw std::runtime_error("OutputFile: " + filename
				+ " could not be opened.");
}

OutputFile::~OutputFile()
{
	close();
}

void OutputFile::open(const std::string& filename)
{
	_stream.open(filename.c_str(), std::ios::binary);
}

void OutputFile::close()
{
	finish();
	_stream.close();
}

ChunkWriter& OutputFile::getChunkWriter()
{
	return _writer;
}

void OutputFile::setCompressionMode(char compressionMode)
{
	_writer.setCompressionMode(compressionMode);
}

void OutputFile::setCompressionMode(int compressionMode)
{
	_writer.setCompressionMode(compressionMode);
}

OutputFile::OutputFile(const OutputFile& original) :
	_stream(), _writer(_stream)
{
}

OutputFile& OutputFile::operator=(const OutputFile& other)
{
	return *this;
}

}
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------


namespace pxl
{

OutputHandler::OutputHandler(size_t maxSize, size_t maxNObjects) :
	_maxSize(maxSize), _newEvent(true), _maxNObjects(maxNObjects), _nObjects(0)
{
}

OutputHandler::~OutputHandler()
{
}


/// Use this method to write an information string describing the new event. Otherwise, this method need not necessarily be used.
bool OutputHandler::newFileSection(const std::string& info)
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
bool OutputHandler::writeStream(const std::string& info)
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
bool OutputHandler::writeFileSection(const std::string& info)
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

}//namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------




namespace pxl {

void AnalysisFork::beginJob(const ObjectOwner* input)
{
    for(ObjectOwnerTypeIterator<AnalysisFork> iter(&getObjectOwner());
       iter!=getObjectOwner().end(); iter++)
        (*iter)->beginJob(input);

    for(ObjectOwnerTypeIterator<AnalysisProcess> iter(&getObjectOwner());
    iter!=getObjectOwner().end(); iter++)
        (*iter)->beginJob(input);
}


void AnalysisFork::beginRun(const ObjectOwner* input)
{
    for(ObjectOwnerTypeIterator<AnalysisFork> iter(&getObjectOwner());
       iter!=getObjectOwner().end(); iter++)
        (*iter)->beginRun(input);

    for(ObjectOwnerTypeIterator<AnalysisProcess> iter(&getObjectOwner());
    iter!=getObjectOwner().end(); iter++)
        (*iter)->beginRun(input);
}

void AnalysisFork::analyseEvent(const ObjectOwner* input)
{
    for(ObjectOwnerTypeIterator<AnalysisFork> iter(&getObjectOwner());
       iter!=getObjectOwner().end(); iter++)
        (*iter)->analyseEvent(input);

    for(ObjectOwnerTypeIterator<AnalysisProcess> iter(&getObjectOwner());
    iter!=getObjectOwner().end(); iter++)
        (*iter)->analyseEvent(input);
}

void AnalysisFork::finishEvent(const ObjectOwner* input)
{
    for(ObjectOwnerTypeIterator<AnalysisFork> iter(&getObjectOwner());
       iter!=getObjectOwner().end(); iter++)
        (*iter)->finishEvent(input);

    for(ObjectOwnerTypeIterator<AnalysisProcess> iter(&getObjectOwner());
    iter!=getObjectOwner().end(); iter++)
        (*iter)->finishEvent(input);
}

void AnalysisFork::endRun(const ObjectOwner* input)
{
    for(ObjectOwnerTypeIterator<AnalysisFork> iter(&getObjectOwner());
       iter!=getObjectOwner().end(); iter++)
        (*iter)->endRun(input);

    for(ObjectOwnerTypeIterator<AnalysisProcess> iter(&getObjectOwner());
    iter!=getObjectOwner().end(); iter++)
        (*iter)->endRun(input);
}

void AnalysisFork::endJob(const ObjectOwner* input)
{
    for(ObjectOwnerTypeIterator<AnalysisFork> iter(&getObjectOwner());
       iter!=getObjectOwner().end(); iter++)
        (*iter)->endJob(input);

    for(ObjectOwnerTypeIterator<AnalysisProcess> iter(&getObjectOwner());
    iter!=getObjectOwner().end(); iter++)
        (*iter)->endJob(input);
}

Serializable* AnalysisFork::clone() const
{
    return new AnalysisFork(*this);
}

std::ostream& AnalysisFork::print(int level, std::ostream& os, int pan) const
{
    printPan1st(os, pan) << "AnalysisFork: " << getName() << std::endl;

    if (level>0) os << printContent(level, os, pan);

    for(ObjectOwner::const_iterator iter = getObjectOwner().begin();
        iter!=getObjectOwner().end(); ++iter)
    {

        if ((*iter)->getMotherRelations().size() == 0)
            (*iter)->printDecayTree(level, os, pan);
    }
    return os;
}

static ObjectProducerTemplate<AnalysisFork>
		_AnalysisForkProducer(AnalysisFork::getStaticTypeId());

} // namespace pxl

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------



namespace pxl
{

Serializable* AnalysisProcess::clone() const
{
	return new AnalysisProcess(*this);
}

std::ostream& AnalysisProcess::print(int level, std::ostream& os, int pan) const
{
	printPan1st(os, pan) << "AnalysisProcess: " << getName() << std::endl;

	if (level>0) os << printContent(level, os, pan);

	for (ObjectOwner::const_iterator iter = getObjectOwner().begin(); iter!=getObjectOwner().end(); ++iter)
	{
		if ((*iter)->getMotherRelations().size() == 0)
			(*iter)->printDecayTree(level, os, pan);
	}
	return os;
}

static ObjectProducerTemplate<AnalysisProcess>
		_AnalysisProcessProducer(AnalysisProcess::getStaticTypeId());

} // namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------


namespace pxl {

std::ostream& Collision::print(int level, std::ostream& os, int pan) const
{
    printPan1st(os, pan) << "Collision: " << getName() << std::endl;

    if (level>0) os << printContent(level, os, pan);

    return os;
}

static ObjectProducerTemplate<Collision>
		_CollisionProducer(Collision::getStaticTypeId());

} // namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------



namespace pxl {

std::ostream& EventView::print(int level, std::ostream& os, int pan) const
{
    printPan1st(os, pan) << "EventView: " << getName() << std::endl;
    
    if (level>0) printContent(level, os, pan);
    
    for(ObjectOwner::const_iterator iter = getObjectOwner().begin();
        iter!=getObjectOwner().end(); iter++) {

        if ((*iter)->getMotherRelations().size() == 0)
            (*iter)->printDecayTree(level, os, pan);
    }
    return os;
}

const Id& EventView::getStaticTypeId()
{
	static const Id id("c8db3cce-dc4b-421e-882a-83e213c9451f");
	return id;
}

static ObjectProducerTemplate<EventView>
		_EventViewProducer(EventView::getStaticTypeId());

} // namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------



namespace pxl
{

bool const operator==(const Particle& obj1, const Particle& obj2)
{
	return obj1.getVector() == obj2.getVector() && obj1.getCharge()
			== obj2.getCharge();
}

bool const operator!=(const Particle& obj1, const Particle& obj2)
{
	return obj1.getVector() != obj2.getVector() || obj1.getCharge()
			== obj2.getCharge();
}

std::ostream& Particle::print(int level, std::ostream& os, int pan) const
{
	printPan1st(os, pan);
	os << "Particle: '" << getName() << "', p = (" << getPt() << ", "
			<< getPz() << ") m = " << getMass() << std::endl;
	if (level>0)
		Object::printContent(level, os, pan);
	return os;
}

void Particle::setP4FromDaughters()
{
	setP4(0., 0., 0., 0.);
	for (Relations::const_iterator iter = getDaughterRelations().begin(); iter!=getDaughterRelations().end(); ++iter)
	{
		CommonParticle* daughter =
				dynamic_cast<CommonParticle*>(*iter);
		if (daughter)
			addP4(daughter->getPx(), daughter->getPy(), daughter->getPz(),
					daughter->getE());
	}
}

void Particle::setP4FromDaughtersRecursive()
{
	if (getDaughterRelations().size() > 0)
		setP4(0., 0., 0., 0.);

	for (Relations::const_iterator iter = getDaughterRelations().begin(); iter!=getDaughterRelations().end(); ++iter)
	{
		// update daughter particles
		Particle *particle = dynamic_cast<Particle *>(*iter);
		if (particle)
			particle->setP4FromDaughtersRecursive();
	}
	
	for (Relations::const_iterator iter = getDaughterRelations().begin(); iter!=getDaughterRelations().end(); ++iter)
	{
		// add all common particles
		CommonParticle* daughter =
				dynamic_cast<CommonParticle*>(*iter);
		if (daughter)
			addP4(daughter->getPx(), daughter->getPy(), daughter->getPz(),
					daughter->getE());

	}
}

static ObjectProducerTemplate<Particle>
		_ParticleProducer(Particle::getStaticTypeId());

} // namespace pxl

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------



namespace pxl
{

bool const operator==(const Vertex& obj1, const Vertex& obj2)
{
	return obj1.getVector() == obj2.getVector();
}

bool const operator!=(const Vertex& obj1, const Vertex& obj2)
{
	return obj1.getVector() != obj2.getVector();
}

std::ostream& Vertex::print(int level, std::ostream& os, int pan) const
{
	printPan1st(os, pan) << "Vertex: '" << getName() << "', x = (" << getX()
			<< ", " << getY() << ", " << getZ() << ")" << std::endl;
	if (level > 0) printContent(level, os, pan);
	return os;
}

static ObjectProducerTemplate<Vertex>
		_VertexProducer(Vertex::getStaticTypeId());

} // namespace pxl

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#include <cmath>
namespace pxl
{
	static ObjectProducerTemplate<AstroBasicObject> _AstroBasicObjectProducer(AstroBasicObject::getStaticTypeId());


	void AstroBasicObject::serialize(const OutputStream &out) const
	{
		Serializable::serialize(out);
		_vector.serialize(out);
		out.writeDouble(_dphi);
		out.writeDouble(_dtheta);
		out.writeLong((long)_time);
	}

	void AstroBasicObject::deserialize(const InputStream &in)
	{
		Serializable::deserialize(in);
		_vector.deserialize(in);
		in.readDouble(_dphi);
		in.readDouble(_dtheta);
		long t;
		in.readLong(t);
		_time = t;
	}


	void AstroBasicObject::setAzimuthZenith(double azimuth,double zenith)
	{
		_vector.setX(cos(azimuth) * sin(zenith));
		_vector.setY(sin(azimuth) * sin(zenith));
		_vector.setZ(cos(zenith));
	}

	void AstroBasicObject::setAzimuth(double azimuth)
	// Updates direction, assumes that zenith is set correctly
	{
		_vector.setX(cos(azimuth) * sin(getZenith()));
		_vector.setY(sin(azimuth) * sin(getZenith()));
	}

	void AstroBasicObject::setZenith(double zenith)
	// Updates direction, assumes that azimuth is set correctly
	{
		_vector.setX(cos(getAzimuth()) * sin(zenith));
		_vector.setY(sin(getAzimuth()) * sin(zenith));
		_vector.setZ(cos(zenith));
	}

	void AstroBasicObject::setLongitudeLatitude(double longitude,double latitude)
	{
		_vector.setX(cos(longitude) * cos(latitude));
		_vector.setY(sin(longitude) * cos(latitude));
		_vector.setZ(sin(latitude));
	}

	void AstroBasicObject::setLongitude(double longitude)
	// Updates direction, assumes that Longitude is set correctly
	{
		_vector.setX(cos(longitude) * cos(getLatitude()));
		_vector.setY(sin(longitude) * cos(getLatitude()));
	}

	void AstroBasicObject::setLatitude(double latitude)
	// Updates direction, assumes that Latitude is set correctly
	{
		_vector.setX(cos(getLongitude()) * cos(latitude));
		_vector.setY(sin(getLongitude()) * cos(latitude));
		_vector.setZ(sin(latitude));
	}


	double AstroBasicObject::getAzimuth() const
	{
		return atan2(_vector.getY(),_vector.getX());
	}

	double AstroBasicObject::getZenith() const
	{
		return acos(_vector.getZ());
	}

	double AstroBasicObject::getLongitude() const
	{
		return atan2(_vector.getY(),_vector.getX());
	}

	double AstroBasicObject::getLatitude() const
	{
		return asin(_vector.getZ());
	}

	double AstroBasicObject::angularDistanceTo(const AstroBasicObject &obj) const
	{
		double cosdistance = obj.getVector() * _vector;
		return acos(cosdistance);
	}
	double AstroBasicObject::angularDistanceTo(const AstroBasicObject *obj) const
	{
		double cosdistance = (obj->getVector()) * _vector;
		return acos(cosdistance);
	}

}
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------


namespace pxl
{

static ObjectProducerTemplate<AstroObject> _AstroObjectProducer(
		AstroObject::getStaticTypeId());

void AstroObject::serialize(const OutputStream &out) const
{
	AstroBasicObject::serialize(out);
	_userRecords.serialize(out);
	_softRelations.serialize(out);
}

void AstroObject::deserialize(const InputStream &in)
{
	AstroBasicObject::deserialize(in);
	_userRecords.deserialize(in);
	_softRelations.deserialize(in);
}

void AstroObject::linkSoft(AstroObject* astroobject, const std::string& type)
{
	_softRelations.set(astroobject, type);
	astroobject->_softRelations.set(this, type);
}

void AstroObject::unlinkSoft(AstroObject* astroobject, const std::string& type)
{
	if (astroobject)
	{
		_softRelations.remove(astroobject, type);
		astroobject->_softRelations.remove(this, type);
	}
}

void AstroObject::linkSoft(AstroObject &astroobject, const std::string& type)
{
	_softRelations.set(&astroobject, type);
	astroobject._softRelations.set(this, type);
}

void AstroObject::unlinkSoft(AstroObject &astroobject, const std::string& type)
{
	_softRelations.remove(&astroobject, type);
	astroobject._softRelations.remove(this, type);
}




}
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------


namespace pxl
{

static ObjectProducerTemplate<RegionOfInterest> _RegionOfInterestProducer(
		RegionOfInterest::getStaticTypeId());

std::ostream& RegionOfInterest::print(int level, std::ostream& os, int pan) const
{
	os << "Name: " << getName() << std::endl;
	os << "  Latitude : " << getLatitude() << std::endl;
	os << "  Longitude : " << getLongitude() << std::endl;
	os << "  ConeRadius: " << getConeRadius() << std::endl;
	return os;
}

void RegionOfInterest::serialize(const OutputStream &out) const
{
	AstroObject::serialize(out);
	out.writeDouble(_coneRadius);
	out.writeString(_name);
}

void RegionOfInterest::deserialize(const InputStream &in)
{
	AstroObject::deserialize(in);
	in.readDouble(_coneRadius);
	in.readString(_name);
}

bool RegionOfInterest::objectInCone(const AstroBasicObject& obj) const
{
	return (this->angularDistanceTo(obj) < _coneRadius);
}

bool RegionOfInterest::objectInCone(const AstroBasicObject* obj) const
{
	return (this->angularDistanceTo(obj) < _coneRadius);
}


bool RegionOfInterest::linkIfObjectInCone(AstroObject& obj,std::string name)
{
	if (this->angularDistanceTo(obj) < _coneRadius)
	{
		this->linkSoft(obj,name);
		return true;
	}
	return false;
}

bool RegionOfInterest::linkIfObjectInCone(AstroObject* obj,std::string name)
{
	if (this->angularDistanceTo(obj) < _coneRadius)
	{
		this->linkSoft(obj,name);
		return true;
	}
	return false;
}
}
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------


namespace pxl
{

static ObjectProducerTemplate<UHECR> _UHECRProducer(UHECR::getStaticTypeId());

std::ostream& UHECR::print(int level, std::ostream& os, int pan) const
{
	os << " Latitude : " << getLatitude() << std::endl;
	os << "Longitude : " << getLongitude() << std::endl;
	os << "   Energy : " << getEnergy() << std::endl;
	return os;
}

void UHECR::serialize(const OutputStream &out) const
{
	AstroObject::serialize(out);
	out.writeDouble(_energy);
	out.writeDouble(_denergy);
	out.writeDouble(_mass);
	out.writeDouble(_charge);
}

void UHECR::deserialize(const InputStream &in)
{
	AstroObject::deserialize(in);
	in.readDouble(_energy);
	in.readDouble(_denergy);
	in.readDouble(_mass);
	in.readDouble(_charge);
}

}
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2009 Martin Erdmann        -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------


namespace pxl
{

static ObjectProducerTemplate<UHECRSource> _UHECRSourceProducer(
		UHECRSource::getStaticTypeId());

std::ostream& UHECRSource::print(int level, std::ostream& os, int pan) const
{
	os << "Name: " << getName() << std::endl;
	os << "  Latitude  : " << getLatitude() << std::endl;
	os << "  Longitude : " << getLongitude() << std::endl;
	os << "  Distance  : " << getDistance() << std::endl;
	os << "  Rel. Lum. : " << getLuminosity() << std::endl;
	return os;
}

void UHECRSource::serialize(const OutputStream &out) const
{
	AstroObject::serialize(out);
	out.writeDouble(_distance);
	out.writeDouble(_luminosity);
	out.writeString(_name);
}

void UHECRSource::deserialize(const InputStream &in)
{
	AstroObject::deserialize(in);
	in.readDouble(_distance);
	in.readDouble(_luminosity);
	in.readString(_name);
}

}
