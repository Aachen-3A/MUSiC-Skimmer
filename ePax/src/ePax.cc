//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#include "ePaxPxl/ePax/interface/ePax.h"

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
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------


namespace pxl
{

MutableId::MutableId(const char* id)
{
	// reset
	for (int j = 0; j < 16; j++)
		bytes[j] = 0;

	// read 32 char = 16 bytes
	unsigned char first, second;
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

MutableId::MutableId(const std::string& id)
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


std::string MutableId::toString() const
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


} //namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#include <iostream>


namespace pxl
{

void ObjectOwner::init(const pxl::ObjectOwner& original)
{
	// copy objects: loop in STL style
	for (const_iterator iter = original._container.begin(); iter
			!= original._container.end(); iter++)
	{

		pxl::Relative* pOld = *iter;
		pxl::Relative* pNew = pOld->clone();

		set(pNew);

		pxl::CopyHistory::iterator insertPos =
				_copyHistory.lower_bound(pOld->id());
		if (insertPos == _copyHistory.end() || insertPos->first != pOld->id())
			_copyHistory.insert(insertPos,
					pxl::CopyHistory::iterator::value_type(pOld->id(), pNew));
		else
			insertPos->second = pNew;
	}

	// FIXME: possibly inefficient, might be done all in one loop
	// redirect relations: loop in PTL style
	for (const_iterator iter = original._container.begin(); iter
			!=original._container.end(); iter++)
	{
		pxl::Relative* pOld = *iter;
		pxl::Relative* pNew = 0;
		pxl::CopyHistory::const_iterator found = _copyHistory.find(pOld->id());
		if (found != _copyHistory.end())
			pNew = found->second;

		// mother relations
		for (pxl::Relations::const_iterator iter = pOld->getMotherRelations().begin(); iter!=pOld->getMotherRelations().end(); ++iter)
		{
			pxl::Relative* pOldRel = *iter;
			pxl::Relative* pNewRel = 0;
			pxl::CopyHistory::const_iterator foundRel =
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
							<< "pxl::ObjectOwner::ObjectOwner(...): WARNING: some original objects had relations to objects of other owners."
							<< std::endl;
			}
			else
				std::cerr
						<< "pxl::ObjectOwner::ObjectOwner(...): WARNING: some originally related objects no longer exist."
						<< std::endl;
		}

		// daughter relations
		// have been set automatically above
	}

	// redirect index:
	for (pxl::Index::const_iterator iter = original._index.begin(); iter
			!=original._index.end(); ++iter)
	{

		pxl::Relative* pOld = iter->second;

		pxl::Relative* pNew = 0;
		pxl::CopyHistory::const_iterator found = _copyHistory.find(pOld->id());
		if (found != _copyHistory.end())
			pNew = found->second;

		if (pNew)
			_index.insert(pxl::Index::const_iterator::value_type(iter->first,
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

void ObjectOwner::set(pxl::Relative* item)
{
	item->_refObjectOwner = this;
	_container.push_back(item);
	_uuidSearchMap.insert(std::pair<pxl::Id, pxl::Relative*>(item->getId(), item));
}

void ObjectOwner::remove(pxl::Relative* item)
{
	// search & remove possible indices (multiple occurrences possible!)
	for (pxl::Index::const_iterator iter = _index.begin(); iter != _index.end(); iter++)
	{
		if (item == iter->second)
			_index.erase(iter->first);
	}

	// search & remove possible copy history
	for (pxl::CopyHistory::const_iterator iter = _copyHistory.begin(); iter
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

bool ObjectOwner::has(const pxl::Relative* item) const
{
	return item->_refObjectOwner == this;
}

bool ObjectOwner::setIndex(const std::string& idx, pxl::Relative* obj,
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

	pxl::Index::iterator insertPos = _index.lower_bound(idx);
	if (insertPos == _index.end() || insertPos->first != idx)
	_index.insert(insertPos, pxl::Index::iterator::value_type(idx, obj));
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
	for (pxl::Index::const_iterator iter = _index.begin(); iter != _index.end(); ++iter)
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

	std::map<pxl::Id, pxl::Relative*> objIdMap;
	std::multimap<pxl::Relative*, pxl::Id> daughterRelationsMap;
	std::multimap<pxl::Relative*, pxl::Id> motherRelationsMap;
	std::multimap<pxl::Relative*, pxl::Id> flatRelationsMap;

	unsigned int size = 0;
	in.readUnsignedInt(size);
	for (unsigned int i=0; i<size; ++i)
	{
		pxl::Id typeId (in);
		// Contained object must be a Relative derivative
		pxl::Relative* object = dynamic_cast<pxl::Relative*>(pxl::ObjectFactory::instance().create(typeId));
		object->deserialize(in);
		set(object);
		objIdMap.insert(std::pair<pxl::Id, pxl::Relative*>(object->id(), object));

		int msize = 0;
		in.readInt(msize);
		for (int j=0; j<msize;++j)
		{
			pxl::Id id (in);
			motherRelationsMap.insert(std::pair<pxl::Relative*, pxl::Id>(object, id));
		}

		int dsize = 0;
		in.readInt(dsize);
		for (int j=0; j<dsize;++j)
		{
			pxl::Id id (in);
			daughterRelationsMap.insert(std::pair<pxl::Relative*, pxl::Id>(object, id));
		}
		int fsize = 0;
		in.readInt(fsize);
		for (int j=0; j<fsize;++j)
		{
			pxl::Id id (in);
			flatRelationsMap.insert(std::pair<pxl::Relative*, pxl::Id>(object, id));
		}

	}

	for (std::multimap<pxl::Relative*, pxl::Id>::const_iterator iter = daughterRelationsMap.begin();
	iter!=daughterRelationsMap.end(); ++iter)
	{
		pxl::Relative* target = objIdMap.find(iter->second)->second;
		iter->first->linkDaughter(target);
	}
	
	for (std::multimap<pxl::Relative*, pxl::Id>::const_iterator iter = flatRelationsMap.begin();
	iter!=flatRelationsMap.end(); ++iter)
	{
		pxl::Relative* target = objIdMap.find(iter->second)->second;
		iter->first->linkFlat(target);
	}

	

	in.readUnsignedInt(size);

	for (unsigned int i=0; i<size; ++i)
	{
		std::string name;
		in.readString(name);
		pxl::Id id (in);
		_index.insert(std::pair<std::string, pxl::Relative*>(name, objIdMap.find(id)->second));
	}
}

} // namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------


namespace pxl
{

void pxl::Relations::serialize(const OutputStream &out) const
{
	out.writeInt(size());
	for (const_iterator iter = begin(); iter!=end(); ++iter)
	{
		(*iter)->getId().serialize(out);
	}
}

}
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#include <stdexcept>


namespace pxl
{

void Relative::linkMother(pxl::Relative* target) throw(std::runtime_error)
{
	if (target->_refObjectOwner != this->_refObjectOwner)
		throw std::runtime_error("pxl::ObjectBase::linkDaughter(...): WARNING: mother and daughter have not the same object holder!");

	this->_motherRelations.set(target);
	target->_daughterRelations.set(this);
}

void Relative::linkDaughter(pxl::Relative* target) throw(std::runtime_error)
{
	if (target->_refObjectOwner != this->_refObjectOwner)
		throw std::runtime_error("pxl::ObjectBase::linkMother(...): WARNING: mother and daughter have not the same object holder!");

	this->_daughterRelations.set(target);
	target->_motherRelations.set(this);
}

void Relative::linkFlat(pxl::Relative* target) throw(std::runtime_error)
{
	if (target->_refObjectOwner != this->_refObjectOwner)
		throw std::runtime_error("pxl::ObjectBase::linkFlat(...): WARNING: potential relatives have not the same object holder!");

	this->_flatRelations.set(target);
	target->_flatRelations.set(this);
}

void Relative::unlinkMother(pxl::Relative* target)
{
	this->_motherRelations.erase(target);
	target->_daughterRelations.erase(this);
}

void Relative::unlinkDaughter(pxl::Relative* target)
{
	this->_daughterRelations.erase(target);
	target->_motherRelations.erase(this);
}

void Relative::unlinkFlat(pxl::Relative* target)
{
	this->_flatRelations.erase(target);
	target->_flatRelations.erase(this);
}

void Relative::unlinkMothers()
{
	for (pxl::Relations::const_iterator iter = _motherRelations.begin(); iter
			!=_motherRelations.end(); ++iter)
	{
		(*iter)->_daughterRelations.erase(this);
	}

	_motherRelations.clearContainer();
}

void Relative::unlinkDaughters()
{
	for (pxl::Relations::const_iterator iter = _daughterRelations.begin(); iter
			!=_daughterRelations.end(); ++iter)
	{

		(*iter)->_motherRelations.erase(this);
	}

	_daughterRelations.clearContainer();
}

void Relative::unlinkFlat()
{
	for (pxl::Relations::const_iterator iter = _flatRelations.begin(); iter
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

	for (pxl::Relations::const_iterator iter = _daughterRelations.begin(); iter
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
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------



namespace pxl
{
pxl::Relative* SoftRelations::getFirst(const pxl::ObjectOwner& owner,
const std::string& type) const
{
	if (type=="")
	{
		const_iterator found = _relationsMap.begin();
		if (found!=_relationsMap.end())
		{
			return owner.getById(found->second);
		}
		else
		return 0;
	}
	else
	{
		const_iterator found = _relationsMap.find(type);
		if (found!=_relationsMap.end())
		{
			return owner.getById(found->second);
		}
		else
		return 0;
	}

}

int SoftRelations::getSoftRelatives(std::vector<pxl::Relative*>& vec,
const pxl::ObjectOwner& owner, const std::string& type) const
{
	int size = vec.size();

	std::pair<const_iterator, const_iterator> iterators = _relationsMap.equal_range(type);
	for (const_iterator iter = iterators.first; iter!= iterators.second; ++iter)
	{
		pxl::Relative* relative = owner.getById(iter->second);
		if (relative!=0)
		vec.push_back(relative);
	}

	return vec.size()-size;
}

int SoftRelations::getSoftRelatives(std::vector<pxl::Relative*>& vec,
const pxl::ObjectOwner& owner) const
{
	int size = vec.size();
	for (const_iterator iter = _relationsMap.begin(); iter
	!= _relationsMap.end(); ++iter)
	{
		pxl::Relative* relative = owner.getById(iter->second);
		if (relative!=0)
		vec.push_back(relative);
	}
	return vec.size()-size;
}

template<class objecttype>
int SoftRelations::getSoftRelativesOfType(std::vector<objecttype*>& vec,
const pxl::ObjectOwner& owner, const std::string& type) const
{
	int size = vec.size();
	std::pair<const_iterator, const_iterator> iterators = _relationsMap.equal_range(type);
	for (const_iterator iter = iterators.first; iter!= iterators.second; ++iter)
	{
		pxl::Relative* relative = owner.getById(iter->second);
		if (relative!=0)
		{
			objecttype* object = dynamic_cast<objecttype*>(relative);
			if (object!=0)
			vec.push_back(object);
		}
	}
	return vec.size()-size;
}

template<class objecttype>
int SoftRelations::getSoftRelativesOfType(std::vector<objecttype*>& vec,
const pxl::ObjectOwner& owner) const
{
	int size = vec.size();
	for (const_iterator iter = _relationsMap.begin(); iter
	!= _relationsMap.end(); ++iter)
	{
		pxl::Relative* relative = owner.getById(iter->second);
		if (relative!=0)
		{
			objecttype* object = dynamic_cast<objecttype*>(relative);
			if (object!=0)
			vec.push_back(object);
		}
	}
	return vec.size()-size;
}

int SoftRelations::keepSoftRelatives(std::vector<pxl::Relative*>& vec) const
{
	std::vector<pxl::Relative*> keepItems;
	for (const_iterator iter = _relationsMap.begin(); iter
	!= _relationsMap.end(); ++iter)
	{
		for (std::vector<pxl::Relative*>::const_iterator itVec = vec.begin(); itVec!=vec.end(); ++itVec)
		{
			if ( (*itVec)->getId() == iter->second )
			keepItems.push_back(*itVec);
		}
	}
	vec.swap(keepItems);
	return vec.size();
}

int SoftRelations::keepSoftRelatives(std::vector<pxl::Relative*>& vec, const std::string& type) const
{
	std::vector<pxl::Relative*> keepItems;
	for (const_iterator iter = _relationsMap.begin(); iter
	!= _relationsMap.end(); ++iter)
	{
		if (iter->first==type)
		{
			for (std::vector<pxl::Relative*>::const_iterator itVec = vec.begin(); itVec!=vec.end(); ++itVec)
			{
				if ( (*itVec)->getId() == iter->second )
				keepItems.push_back(*itVec);
			}
		}
	}
	vec.swap(keepItems);
	return vec.size();
}

bool SoftRelations::has(const pxl::Relative* relative) const
{
	for (const_iterator iter = _relationsMap.begin(); iter
	!= _relationsMap.end(); ++iter)
	{
		if ( relative->getId() == iter->second )
		return true;
	}
	return false;
}

bool SoftRelations::has(const pxl::Relative* relative, const std::string& type) const
{
	std::pair<const_iterator, const_iterator> iterators = _relationsMap.equal_range(type);
	for (const_iterator iter = iterators.first; iter!= iterators.second; ++iter)
	{
		if ( relative->getId() == iter->second)
		return true;
	}
	return false;
}

bool SoftRelations::has(const pxl::Id& id) const
{
	for (const_iterator iter = _relationsMap.begin(); iter
	!= _relationsMap.end(); ++iter)
	{
		if ( id == iter->second )
		return true;
	}
	return false;
}

bool SoftRelations::has(const pxl::Id& id, const std::string& type) const
{
	std::pair<const_iterator, const_iterator> iterators = _relationsMap.equal_range(type);
	for (const_iterator iter = iterators.first; iter!= iterators.second; ++iter)
	{
		if ( id == iter->second )
		return true;
	}
	return false;
}

bool SoftRelations::hasType(const std::string& name) const
{
	if ( _relationsMap.count(name)>0 )
	return true;
	return false;
}

void SoftRelations::set(const pxl::Relative* relative, const std::string& type)
{
	_relationsMap.insert(std::pair<std::string, pxl::Id>(type, relative->getId()));
}

int SoftRelations::count(const std::string& name) const
{
	return _relationsMap.count(name);
}

void SoftRelations::remove(const pxl::Relative* relative)
{
	for (iterator iter = _relationsMap.begin(); iter != _relationsMap.end(); ++iter)
	{
		if (iter->second==relative->getId())
		{
			_relationsMap.erase(iter);
			break;
		}
	}
}

void SoftRelations::remove(const pxl::Relative* relative, const std::string& type)
{
	std::pair<iterator, iterator> iterators = _relationsMap.equal_range(type);
	for (iterator iter = iterators.first; iter!= iterators.second; ++iter)
	{
		if (iter->second==relative->getId())
		{
			_relationsMap.erase(iter);
			break;
		}
	}
}

std::ostream& SoftRelations::print(int level, std::ostream& os) const
{
	os << "SoftRelations of size " << _relationsMap.size() << "\n";
	for (const_iterator iter = _relationsMap.begin(); iter!=_relationsMap.end(); ++iter)
	{
		os << "--> ('" << iter->first << "', " << iter->second << ") \n";
	}
	return os;
}

} //namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------



namespace pxl
{

void UserRecord::serialize(const OutputStream &out) const
{
	out.writeUnsignedInt(getContainer()->size());
	for (const_iterator iter = getContainer()->begin(); iter!=getContainer()->end(); ++iter)
	{
		out.writeString(iter->first);
		pxl::Variant::Type type = iter->second.getType();

		char cType = ' ';
		switch (type)
		{
		case pxl::Variant::TYPE_BOOL:
			cType='b';
			out.writeChar(cType);
			out.writeBool(iter->second.get<bool>());
			break;
		case pxl::Variant::TYPE_CHAR:
			cType = 'c';
			out.writeChar(cType);
			out.writeChar(iter->second.get<char>());
			break;
		case pxl::Variant::TYPE_UCHAR:
			cType = 'C';
			out.writeChar(cType);
			out.writeUnsignedChar(iter->second.get<unsigned char>());
			break;
		case pxl::Variant::TYPE_INT:
			cType = 'i';
			out.writeChar(cType);
			out.writeInt(iter->second.get<int>());
			break;
		case pxl::Variant::TYPE_UINT:
			cType = 'I';
			out.writeChar(cType);
			out.writeUnsignedInt(iter->second.get<unsigned int>());
			break;
		case pxl::Variant::TYPE_SHORT:
			cType = 'o';
			out.writeChar(cType);
			out.writeShort(iter->second.get<short>());
			break;
		case pxl::Variant::TYPE_USHORT:
			cType = 'O';
			out.writeChar(cType);
			out.writeUnsignedShort(iter->second.get<unsigned short>());
			break;
		case pxl::Variant::TYPE_LONG:
			cType = 'l';
			out.writeChar(cType);
			out.writeLong(iter->second.get<long>());
			break;
		case pxl::Variant::TYPE_ULONG:
			cType = 'L';
			out.writeChar(cType);
			out.writeUnsignedLong(iter->second.get<unsigned long>());
			break;
		case pxl::Variant::TYPE_DOUBLE:
			cType = 'd';
			out.writeChar(cType);
			out.writeDouble(iter->second.get<double>());
			break;
		case pxl::Variant::TYPE_FLOAT:
			cType = 'f';
			out.writeChar(cType);
			out.writeFloat(iter->second.get<float>());
			break;
		case pxl::Variant::TYPE_STRING:
			cType = 's';
			out.writeChar(cType);
			out.writeString(iter->second.get<std::string>());
			break;
		case pxl::Variant::TYPE_USER:
			cType = 'u';
			out.writeChar(cType);
			break;
		case pxl::Variant::TYPE_PTR:
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
	for (unsigned int j=0; j<size; ++j)
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
			set<bool>(name, b);
		}
			break;
		case 'c':
		{
			char c;
			in.readChar(c);
			set<char>(name, c);
		}
			break;
		case 'C':
		{
			unsigned char c;
			in.readUnsignedChar(c);
			set<unsigned char>(name, c);
		}
			break;

		case 'i':
		{
			int ii;
			in.readInt(ii);
			set<int>(name, ii);
		}
			break;
		case 'I':
		{
			unsigned int ui;
			in.readUnsignedInt(ui);
			set<unsigned int>(name, ui);
		}
			break;
		case 'o':
		{
			short s;
			in.readShort(s);
			set<short>(name, s);
		}
			break;
		case 'O':
		{
			unsigned short us;
			in.readUnsignedShort(us);
			set<unsigned short>(name, us);
		}
			break;
		case 'l':
		{
			long l;
			in.readLong(l);
			set<long>(name, l);
		}
			break;
		case 'L':
		{
			unsigned long ul;
			in.readUnsignedLong(ul);
			set<unsigned long>(name, ul);
		}
			break;
		case 'd':
		{
			double d;
			in.readDouble(d);
			set<double>(name, d);
		}
			break;
		case 'f':
		{
			float f;
			in.readFloat(f);
			set<float>(name, f);
		}
			break;
		case 's':
		{
			std::string ss;
			in.readString(ss);
			set<std::string>(name, ss);
		}
			break;
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
	for (const_iterator iter = getContainer()->begin(); iter!=getContainer()->end(); ++iter)
	{
		os << "-->";
		pxl::Variant::Type type = iter->second.getType();

		os << " ('" << iter->first << "', ";
		switch (type)
		{
		case pxl::Variant::TYPE_BOOL:
			os << iter->second.get<bool>();
			break;
		case pxl::Variant::TYPE_CHAR:
			os << iter->second.get<char>();
			break;
		case pxl::Variant::TYPE_UCHAR:
			os << iter->second.get<unsigned char>();
			break;
		case pxl::Variant::TYPE_INT:
			os << iter->second.get<int>();
			break;
		case pxl::Variant::TYPE_UINT:
			os << iter->second.get<unsigned int>();
			break;
		case pxl::Variant::TYPE_SHORT:
			os << iter->second.get<short>();
			break;
		case pxl::Variant::TYPE_USHORT:
			os << iter->second.get<unsigned short>();
			break;
		case pxl::Variant::TYPE_LONG:
			os << iter->second.get<long>();
			break;
		case pxl::Variant::TYPE_ULONG:
			os << iter->second.get<unsigned long>();
			break;
		case pxl::Variant::TYPE_DOUBLE:
			os << iter->second.get<double>();
			break;
		case pxl::Variant::TYPE_FLOAT:
			os << iter->second.get<float>();
			break;
		case pxl::Variant::TYPE_STRING:
			os << iter->second.get<std::string>();
			break;
		case pxl::Variant::TYPE_USER:
			break;
		case pxl::Variant::TYPE_PTR:
			os << iter->second.get<void*>();
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
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------

#include <cstddef>
#include <sstream>
#include <string>
#include <vector>


namespace pxl
{

std::vector<pxl::VariantBase::TypeInfo> VariantBase::types;

/// This constant array serves the PXL variant data type. 
static const char *typeNames[] =
{ // Note: must match pxl::VariantBase::Type
		"null", "bool", "char", "uchar", "short", "ushort", "int", "uint",
				"long", "ulong", "float", "double", "string", "ptr" };

const pxl::VariantBase::TypeInfo& VariantBase::fallbackGetTypeInfo(Type t)  throw (std::runtime_error)
{
	if (PXL_UNLIKELY(types.begin() == types.end()))
	{
		for (unsigned int i = 0; i < sizeof typeNames / sizeof *typeNames; i++)
			types.push_back(typeNames[i]);
	}

	if ((std::size_t)t >= types.size())
	{
		std::ostringstream ss;
		ss <<"pxl::VariantBase::fallbackGetTypeInfo: " << "Variant type "
				<< (std::size_t)t << " undefined.";
		throw std::runtime_error(ss.str());
	}

	return types[t];
}

void VariantBase::wrongType(Type tShould, Type tIs) const  throw (std::runtime_error)
{
	const TypeInfo& is = getTypeInfo(tIs);
	const TypeInfo& should = getTypeInfo(tShould);

	std::ostringstream ss;
	ss << "pxl::VariantBase::wrongType: "
			<< "Trying to access pxl::Variant of type \"" << is.name
			<< "\" as \"" << should.name << "\".";

	throw std::runtime_error(ss.str());

}

} // namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------


namespace pxl {

void WkPtrBase::notifyDeleted()
{
    _objectRef = 0;
    if (_notifyChainOut)
        _notifyChainOut->notifyDeleted();
    _notifyChainIn = 0; 
    _notifyChainOut = 0;
}

void WkPtrBase::connect(pxl::Relative* pointer)
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
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
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
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------



namespace pxl
{

bool ChunkReader::skip()
{
	if (_stream.peek()==EOF)
		return false;

	//skip event header
	if (_status == preHeader)
	{
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
	if (_status != preHeader)
		return false;

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
			ss << "Unknown char identifier: " << id	<< "in pxl::ChunkReader::nextBlockIf.";
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
		if (skip == skipMode::on)
			return readBlock(skip, checkInfo, infoCondition);
		else
			return false;
	}
	else
	{
		// read chunk into buffer
		if (compressionMode==' ')
		{
			_buffer.resize(chunkSize);
			_stream.read(_buffer.data(), chunkSize);
			if (_stream.bad() || _stream.eof() )
				return false;
		}
		else if (compressionMode=='Z')
		{
			unzipEventData(chunkSize);
		}
		else
		{
			throw std::runtime_error("ChunkReader::readBlock(): Invalid compression mode.");
		}
	}

	return true;
}

bool ChunkReader::readHeader(readMode mode, skipMode doSkip,
		infoMode checkInfo, const std::string& infoCondition)
{
	endEvent();

	if (_stream.peek()==EOF || _stream.bad() )
		return false;

	switch (mode)
	{
	case infoChunk:
		_status = infoPreBlock;
		if (nextBlockId()!='I')
		{
			skip();
			if (doSkip == skipMode::on)
				return readHeader(mode, doSkip, checkInfo, infoCondition);
			else
				return false;
		}
		break;

	case event:
		_status = evPreBlock;
		if (nextBlockId()!='E')
		{
			skip();
			if (doSkip == skipMode::on)
				return readHeader(mode, doSkip, checkInfo, infoCondition);
			else
				return false;
		}
		break;

	default:
		_status = preBlock;
		_stream.ignore(1);
	}

	pxl::int32_t infoSize = 0;
	_stream.read((char *)&infoSize, 4);

	if (_stream.eof() || _stream.bad() )
		return false;

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
			if (doSkip == skipMode::on)
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

	unsigned char* _inputBuffer =
			new unsigned char[iotl__iStreamer__lengthUnzipBuffer];
	unsigned char* _outputBuffer =
			new unsigned char[iotl__iStreamer__lengthUnzipBuffer];

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
			size = iotl__iStreamer__lengthUnzipBuffer;

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

			_buffer.destroy();
			_buffer.resize(length);
			memcpy(_buffer.data(), _outputBuffer+(length-have), have);

		} while (strm.avail_out == 0);
	} while (nBytes > 0); // done when inflate() says it's done
	inflateEnd(&strm);

	delete[] _inputBuffer;
	delete[] _outputBuffer;

	return length;
}

}

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
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

bool ChunkWriter::write(std::string info, char compressionMode) throw(std::runtime_error)
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
	if (compressionMode == ' ') compressed = ' ';
	_stream.write((char *) &compressed, 1);
	_nBytes+=1;

	// zip block:
	const char* cBuffer = _buffer.getData();
	pxl::int32_t lengthBuffer = _buffer.getSize();

	const char* cZip = cBuffer;
	pxl::int32_t lengthZip = lengthBuffer;

	char* cZipSpace = 0;
	pxl::int32_t lengthZipSpace = 0;

	if (compressionMode == ' ')
	{
		// no compression requires no action...
	}
	else if (compressionMode >= '0' && compressionMode <= '9')
	{
		// data compression a la Gero, i.e. compression level = 6:
		lengthZipSpace = int(double(lengthBuffer) * 1.05 + 16);
		cZipSpace = new char[lengthZipSpace];

		int status = compress2((Bytef*)cZipSpace, (uLongf*)&lengthZipSpace,
				(const Bytef*)cBuffer, lengthBuffer, compressionMode - '0');
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

	_buffer.destroy();

	return true;
}

} //namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------



namespace pxl
{

/// This method explicitly reads an object of type objecttype. Caution: This method should only be used if the type of the following object is known by hard.
template<class objecttype> bool InputHandler::readObject(objecttype* obj)  throw(std::runtime_error)
{
	if (getChunkReader().getInputStream().good())
	{
		pxl::Id id(getChunkReader().getInputStream());
		if (id!=objecttype::getStaticTypeId())
			 throw std::runtime_error("InputFile::readObject(objecttype* obj): unexpected object in file!");
		obj->deserialize(getChunkReader().getInputStream());
		return true;
	}
	return false;
}

/// This method fills the objects from the read-in block into the passed vector. The number of added objects is returned.
int InputHandler::readObjects(std::vector<pxl::Serializable*>& objects)
{
	int nObjects = 0;
	while (getChunkReader().getInputStream().good())
	{
		pxl::Id id(getChunkReader().getInputStream());

		pxl::Serializable* obj = pxl::ObjectFactory::instance().create(id);
		if (obj)
		{
			obj->deserialize(getChunkReader().getInputStream());
			objects.push_back(obj);
			nObjects++;
		}
		else
			std::cerr
					<< "pxl::InputFile::readObjects: Warning: Unknown object Id."
					<< std::endl;
	}
	return nObjects;
}

/// This method fills the objects from the read-in block into the passed pxl::Event. The number of added objects is returned.
int InputHandler::readObjects(pxl::Event* event)
{
	int nObjects = 0;
	while (getChunkReader().getInputStream().good())
	{
		pxl::Id id(getChunkReader().getInputStream());
		pxl::Relative *obj =
				dynamic_cast<pxl::Relative*>(pxl::ObjectFactory::instance().create(id));
		if (obj)
		{
			obj->deserialize(getChunkReader().getInputStream());
			event->setObject(obj);
			nObjects++;
		}
		else
			std::cerr
					<< "pxl::InputFile::readObjects(pxl::Event*): Warning: Unknown object Id or no ObjectBase derivative."
					<< std::endl;
	}
	return nObjects;
}

} //namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------




namespace pxl {

void AnalysisFork::beginJob(const pxl::ObjectOwner* input)
{
    for(pxl::ObjectOwner::TypeIterator<pxl::AnalysisFork> iter(&getObjectOwner());
       iter!=getObjectOwner().end(); iter++)
        (*iter)->beginJob(input);

    for(pxl::ObjectOwner::TypeIterator<pxl::AnalysisProcess> iter(&getObjectOwner());
    iter!=getObjectOwner().end(); iter++)
        (*iter)->beginJob(input);
}


void pxl::AnalysisFork::beginRun(const pxl::ObjectOwner* input)
{
    for(pxl::ObjectOwner::TypeIterator<pxl::AnalysisFork> iter(&getObjectOwner());
       iter!=getObjectOwner().end(); iter++)
        (*iter)->beginRun(input);

    for(pxl::ObjectOwner::TypeIterator<pxl::AnalysisProcess> iter(&getObjectOwner());
    iter!=getObjectOwner().end(); iter++)
        (*iter)->beginRun(input);
}

void pxl::AnalysisFork::analyseEvent(const pxl::ObjectOwner* input)
{
    for(pxl::ObjectOwner::TypeIterator<pxl::AnalysisFork> iter(&getObjectOwner());
       iter!=getObjectOwner().end(); iter++)
        (*iter)->analyseEvent(input);

    for(pxl::ObjectOwner::TypeIterator<pxl::AnalysisProcess> iter(&getObjectOwner());
    iter!=getObjectOwner().end(); iter++)
        (*iter)->analyseEvent(input);
}

void pxl::AnalysisFork::finishEvent(const pxl::ObjectOwner* input)
{
    for(pxl::ObjectOwner::TypeIterator<pxl::AnalysisFork> iter(&getObjectOwner());
       iter!=getObjectOwner().end(); iter++)
        (*iter)->finishEvent(input);

    for(pxl::ObjectOwner::TypeIterator<pxl::AnalysisProcess> iter(&getObjectOwner());
    iter!=getObjectOwner().end(); iter++)
        (*iter)->finishEvent(input);
}

void pxl::AnalysisFork::endRun(const pxl::ObjectOwner* input)
{
    for(pxl::ObjectOwner::TypeIterator<pxl::AnalysisFork> iter(&getObjectOwner());
       iter!=getObjectOwner().end(); iter++)
        (*iter)->endRun(input);

    for(pxl::ObjectOwner::TypeIterator<pxl::AnalysisProcess> iter(&getObjectOwner());
    iter!=getObjectOwner().end(); iter++)
        (*iter)->endRun(input);
}

void pxl::AnalysisFork::endJob(const pxl::ObjectOwner* input)
{
    for(pxl::ObjectOwner::TypeIterator<pxl::AnalysisFork> iter(&getObjectOwner());
       iter!=getObjectOwner().end(); iter++)
        (*iter)->endJob(input);

    for(pxl::ObjectOwner::TypeIterator<pxl::AnalysisProcess> iter(&getObjectOwner());
    iter!=getObjectOwner().end(); iter++)
        (*iter)->endJob(input);
}

pxl::Relative* AnalysisFork::clone() const
{
    return new AnalysisFork(*this);
}

std::ostream& AnalysisFork::print(int level, std::ostream& os, int pan) const
{
    printPan1st(os, pan) << "AnalysisFork: " << getName() << std::endl;
    
    if (level>0) os << printContent(level, os, pan);
    
    for(pxl::ObjectOwner::const_iterator iter = getObjectOwner().begin();
        iter!=getObjectOwner().end(); ++iter) 
    {

        if ((*iter)->getMotherRelations().size() == 0)
            (*iter)->printDecayTree(level, os, pan);
    }
    return os;
}

} // namespace pxl

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------



namespace pxl
{

pxl::Relative* AnalysisProcess::clone() const
{
	return new AnalysisProcess(*this);
}

std::ostream& AnalysisProcess::print(int level, std::ostream& os, int pan) const
{
	printPan1st(os, pan) << "AnalysisProcess: " << getName() << std::endl;
	
	if (level>0) os << printContent(level, os, pan);
	
	for (pxl::ObjectOwner::const_iterator iter = getObjectOwner().begin(); iter!=getObjectOwner().end(); ++iter)
	{
		if ((*iter)->getMotherRelations().size() == 0)
			(*iter)->printDecayTree(level, os, pan);
	}
	return os;
}

} // namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------



namespace pxl
{

bool const operator==(const pxl::Basic3Vector& obj1, const pxl::Basic3Vector& obj2)
{
    return obj1.getX() == obj2.getX() && obj1.getY() == obj2.getY() && obj1.getZ() == obj2.getZ();
}

bool const operator!=(const pxl::Basic3Vector& obj1, const pxl::Basic3Vector& obj2)
{
    return obj1.getX() != obj2.getX() || obj1.getY() != obj2.getY() || obj1.getZ() != obj2.getZ();
}

} // namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------



namespace pxl
{

bool const operator==(const pxl::Basic4Vector& obj1, const pxl::Basic4Vector& obj2)
{
    return obj1.getX() == obj2.getX() && obj1.getY() == obj2.getY() && obj1.getZ() == obj2.getZ() && obj1.getE()
            == obj2.getE();
}

bool const operator!=(const pxl::Basic4Vector& obj1, const pxl::Basic4Vector& obj2)
{
    return obj1.getX() != obj2.getX() || obj1.getY() != obj2.getY() || obj1.getZ() != obj2.getZ() || obj1.getE()
            != obj2.getE();
}

} // namespace pxl

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------


namespace pxl {

std::ostream& pxl::Collision::print(int level, std::ostream& os, int pan) const
{
    printPan1st(os, pan) << "Collision: " << getName() << std::endl;

    if (level>0) os << printContent(level, os, pan);

    return os;
}

static ObjectFactory::ProducerTemplate<pxl::Collision>
		_CollisionProducer(pxl::Collision::getStaticTypeId());

} // namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
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

	for (pxl::ObjectOwner::const_iterator iter = _objects.begin(); iter
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

} // namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------



namespace pxl {

std::ostream& pxl::EventView::print(int level, std::ostream& os, int pan) const
{
    printPan1st(os, pan) << "EventView: " << getName() << std::endl;
    
    if (level>0) printContent(level, os, pan);
    
    for(pxl::ObjectOwner::const_iterator iter = getObjectOwner().begin();
        iter!=getObjectOwner().end(); iter++) {

        if ((*iter)->getMotherRelations().size() == 0)
            (*iter)->printDecayTree(level, os, pan);
    }
    return os;
}

static pxl::ObjectFactory::ProducerTemplate<pxl::EventView>
		_EventViewProducer(pxl::EventView::getStaticTypeId());

} // namespace pxl
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
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

static pxl::ObjectFactory::ProducerTemplate<pxl::Object> _ObjectProducer(pxl::Object::getStaticTypeId());
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------


static pxl::ObjectFactory::ProducerTemplate<pxl::ObjectManager>
		_ObjectManagerProducer(pxl::ObjectManager::getStaticTypeId());
//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------



namespace pxl
{

bool const operator==(const pxl::Particle& obj1, const pxl::Particle& obj2)
{
	return obj1.getVector() == obj2.getVector() && obj1.getCharge()
			== obj2.getCharge();
}

bool const operator!=(const pxl::Particle& obj1, const pxl::Particle& obj2)
{
	return obj1.getVector() != obj2.getVector() || obj1.getCharge()
			== obj2.getCharge();
}

std::ostream& pxl::Particle::print(int level, std::ostream& os, int pan) const
{
	printPan1st(os, pan);
	os << "Particle: '" << getName() << "', p = (" << getPt()
			<< ", " << getPz() << ") m = " << getMass() << std::endl;
	if (level>0)
		pxl::Object::printContent(level, os, pan);
	return os;
}

void pxl::Particle::setP4FromDaughters()
{
	setP4(0., 0., 0., 0.);
	for (pxl::Relations::const_iterator iter = getDaughterRelations().begin(); iter!=getDaughterRelations().end(); ++iter)
	{
		pxl::CommonParticle* daughter = dynamic_cast<pxl::CommonParticle*>(*iter);
		if (daughter)
			addP4(daughter->getPx(), daughter->getPy(), daughter->getPz(), daughter->getE());
	}
}

void pxl::Particle::setP4FromDaughtersRecursive()
{
	if (getDaughterRelations().size() > 0)
		setP4(0., 0., 0., 0.);

	for (pxl::Relations::const_iterator iter = getDaughterRelations().begin(); iter!=getDaughterRelations().end(); ++iter)
	{
		// update daughter particles
		pxl::Particle *particle = dynamic_cast<pxl::Particle *>(*iter);
		if (particle)
			particle->setP4FromDaughtersRecursive();

		// add all common particles
		pxl::CommonParticle* daughter = dynamic_cast<pxl::CommonParticle*>(*iter);
		if (daughter)
			addP4(daughter->getPx(), daughter->getPy(), daughter->getPz(), daughter->getE());
	}
}

static ObjectFactory::ProducerTemplate<pxl::Particle> _ParticleProducer(pxl::Particle::getStaticTypeId());

} // namespace pxl

//-------------------------------------------
// Project: Physics eXtension Library (PXL) -
//          http://pxl.sourceforge.net      -
// Copyright (C) 2006-2008                  -
//               RWTH Aachen, Germany       -
// Contact: pxl-users@lists.sourceforge.net -
//-------------------------------------------



namespace pxl
{

bool const operator==(const pxl::Vertex& obj1, const pxl::Vertex& obj2)
{
	return obj1.getVector() == obj2.getVector();
}

bool const operator!=(const pxl::Vertex& obj1, const pxl::Vertex& obj2)
{
	return obj1.getVector() != obj2.getVector();
}

std::ostream& pxl::Vertex::print(int level, std::ostream& os, int pan) const
{
	printPan1st(os, pan) << "Vertex: '" << getName() << "', x = (" << getX()
			<< ", " << getY() << ", " << getZ() << ")" << std::endl;
	if (level > 0) printContent(level, os, pan);
	return os;
}

static pxl::ObjectFactory::ProducerTemplate<pxl::Vertex>
		_VertexProducer(pxl::Vertex::getStaticTypeId());

} // namespace pxl

