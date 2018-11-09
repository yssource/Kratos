//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Pooyan Dadvand
//

#if !defined(KRATOS_MEMORY_H_INCLUDED )
#define  KRATOS_MEMORY_H_INCLUDED

/* System includes */
#include <utility>
#include <map>

/* External includes */
#include <memory>
#include <includes/kratos_export_api.h>

namespace Kratos {

class KRATOS_API(KRATOS_CORE) MemoryUsageInfo {
public:
	int mCreatedBytes;
	int mCurrentInstances;
	int mCreatedInstances;
	int mDeletedInstances;
	static MemoryUsageInfo& GetMemoryUsageInfoByClassName(std::string const& ClassName) {
		return AllClasesMemoryUsageInfoByName()[ClassName];
	}
	
	static std::map<std::string, MemoryUsageInfo>& AllClasesMemoryUsageInfoByName() {
		static std::map<std::string, MemoryUsageInfo> memory_usage;
		return memory_usage;
	}

	static std::string GetListOfAllAllocatedObjects() {
		auto& memory_usage_info_map = AllClasesMemoryUsageInfoByName();
		std::stringstream buffer;
		for (auto& i : memory_usage_info_map) {
			buffer << i.second.mCurrentInstances << " " << i.first << " (" << i.second.mCreatedInstances << " [" << i.second.mCreatedBytes << "Bytes] created and " << i.second.mDeletedInstances  << " deleted)" << std::endl;
		}
		return buffer.str();
	}
	void Add(int NumberOfInstances, int Size) {
		mCurrentInstances += NumberOfInstances;
		mCreatedInstances += NumberOfInstances;
		mCreatedBytes += Size;
	}
	void Delete(int NumberOfInstances) {
		mCurrentInstances -= NumberOfInstances;
		mDeletedInstances += NumberOfInstances;
	}
};



template<class T>
using shared_ptr = std::shared_ptr<T>; //std::shared_ptr<T>;

template<class T>
using weak_ptr = std::weak_ptr<T>; //std::weak_ptr<T>;

template<class T>
using unique_ptr = std::unique_ptr<T>;

template<typename C, typename...Args>
shared_ptr<C> make_shared(Args &&...args) {
	MemoryUsageInfo::GetMemoryUsageInfoByClassName("shared_ptr").Add(1, sizeof(shared_ptr<C>));
	return shared_ptr<C>(new C(std::forward<Args>(args)...));
	//return std::make_shared<C>(std::forward<Args>(args)...);

}

template<typename C, typename...Args>
unique_ptr<C> make_unique(Args &&...args) {
    // Note: std::make_unique is C++14, this can be updated once we upgrade from C++11
    return unique_ptr<C>(new C(std::forward<Args>(args)...));
}

template<typename C, typename...Args>
shared_ptr<C> static_pointer_cast(Args &&...args) {
    return std::static_pointer_cast<C>(std::forward<Args>(args)...);
}

template<typename C, typename...Args>
shared_ptr<C> dynamic_pointer_cast(Args &&...args) {
    return std::dynamic_pointer_cast<C>(std::forward<Args>(args)...);
}

template<typename C, typename...Args>
shared_ptr<C> const_pointer_cast(Args &&...args) {
    return std::const_pointer_cast<C>(std::forward<Args>(args)...);
}

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
	const MemoryUsageInfo &rThis) {

	return rOStream;
}
// template<typename C, typename...Args>
// shared_ptr<C> reinterpret_pointer_cast(Args &&...args) {
//     return std::reinterpret_pointer_cast<C>(std::forward<Args>(args)...);
// }
} // namespace Kratos




#define KRATOS_CLASS_POINTER_DEFINITION(a) typedef Kratos::shared_ptr<a > Pointer; \
static std::string GetClassName() {return #a;}\
static MemoryUsageInfo& GetClassMemoryUsageInfo(){\
	static MemoryUsageInfo& memory_usage_info = MemoryUsageInfo::GetMemoryUsageInfoByClassName(#a);\
	return memory_usage_info;\
}\
void* operator new(size_t sz)\
{\
	GetClassMemoryUsageInfo().Add(1, sz);\
	return ::operator new(sz);\
}\
\
void* operator new[](size_t sz)\
{\
	GetClassMemoryUsageInfo().Add(1, sz);\
	return ::operator new[](sz);\
}\
\
void operator delete(void* p)\
{\
	GetClassMemoryUsageInfo().Delete(1);\
	::operator delete(p);\
}\
\
void operator delete[](void* p)\
{\
	GetClassMemoryUsageInfo().Delete(1);\
	::operator delete[](p);\
}\
void* operator new(size_t sz, void* ptr)\
{\
	GetClassMemoryUsageInfo().Add(1, sz);\
	return ::operator new(sz,ptr);\
}\
\
void* operator new[](size_t sz, void* ptr)\
{\
	GetClassMemoryUsageInfo().Add(1, sz);\
	return ::operator new[](sz,ptr);\
}\
void operator delete(void* p, void* ptr)\
{\
	GetClassMemoryUsageInfo().Delete(1);\
	::operator delete(p,ptr);\
}\
\
void operator delete[](void* p, void* ptr)\
{\
	GetClassMemoryUsageInfo().Delete(1);\
	::operator delete[](p,ptr);\
}\
typedef Kratos::shared_ptr<a > SharedPointer; \
typedef Kratos::weak_ptr<a > WeakPointer; \
typedef Kratos::unique_ptr<a > UniquePointer\



#endif /* KRATOS_MEMORY_H_INCLUDED  defined */