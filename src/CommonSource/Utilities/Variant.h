// ***************************************************************************
// variant_t - An Improved Variant Type Based on Member Templates
// ---------------------------------------------------------------------------
// (c) 2000 Fernando Cacciola
// Dr. Dobb's (http://www.ddj.com/cpp/184401293)
// ***************************************************************************

#pragma once

#include <stdexcept>
#include <typeinfo>
#include <string>

using namespace std;

class variant_t {
public:
	variant_t() : data (NULL) {}
	variant_t(const variant_t & rhs) { 
		if(rhs.data != NULL) rhs.data->AddRef();
		data = rhs.data;
	}

	~variant_t() { 
		if(data != NULL) data->Release();
	}

	// NOTE: This code takes care of self-assignment.
	// DO NOT CHANGE THE ORDER of the statements.
	variant_t& operator=(const variant_t& rhs) {
		if(rhs.data != NULL) rhs.data->AddRef();
		if(data != NULL) data->Release();
		data = rhs.data;
		return * this;
	}

	// This member template constructor allows you to
	// instance a variant_t object with a value of any type.
	template<typename T>
	variant_t(T v)
		: data (new Impl<T>(v))
	{ 
		data->AddRef(); 
	}

	// This generic conversion operator let you retrieve
	// the value held. To avoid template specialization conflicts,
	// it returns an instance of type T, which will be a COPY
	// of the value contained.
	template<typename T> 
	operator T() const { 
		return CastFromBase<T>(data)->data;
	}

	// This forms returns a REFERENCE and not a COPY, which
	// will be significant in some cases.
	template<typename T> 
	const T& get() const { 
		return CastFromBase<T>(data)->data; 
	}

	template<typename T> 
	bool is_type() const { 
		return typeid(*data)==typeid(Impl<T>); 
	}

	template<typename T> 
	bool is_type(T v) const { 
		return typeid(*data)==typeid(v); 
	}

private:
	struct ImplBase {
		ImplBase() : refs (0) {}
		virtual ~ImplBase() {}
		void AddRef() { refs ++; }
		void Release() { 
			refs --;
			if(refs == 0) delete this;
		}
		size_t refs;
	};

	template<typename T>
	struct Impl : ImplBase {
		Impl (T v) : data (v) {}
		~Impl () {}
		T data;
	};

	// The following method is static because it doesn't
	// operate on variant_t instances.
	template<typename T> 
	static Impl<T>* CastFromBase(ImplBase* v) {
		// This upcast will fail if T is other than the T used
		// with the constructor of variant_t.
		Impl<T>* p = dynamic_cast<Impl<T>*> (v);
		if(p == NULL) throw invalid_argument(typeid(T).name()+string(" is not a valid type"));
		return p;
	}

	ImplBase* data;
};
