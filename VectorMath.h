
#ifndef __GSL_HELPER_H__
#define __GSL_HELPER_H__


#if defined _WIN32 || defined _WIN64
 #include <crtdbg.h>
 #define _gsl_asserte(expr,msg) _ASSERT_EXPR(expr,_CRT_WIDE(msg))
#else
 #include <assert.h>
 #define _gsl_asserte(expr,msg) assert(expr)
#endif
#define _gsl_assert(expr) _gsl_asserte(expr, #expr)


extern "C"
{
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix_double.h>
}
#include <limits>
#include <memory>
#include <functional>
#include <iostream>



namespace gsl
{

struct Precondition
{
	template<typename E, typename Condition>
	void check(Condition cond) {
		if (!cond) throw E();
	}

	template<typename E, typename Condition, typename... ConditionList>
	void check(Condition cond, ConditionList... conds) {
		if (!cond) throw E();
		else check<E>(conds...);
	}

	template<typename E, typename Condition, typename... ConditionList>
	void operator () (Condition cond, ConditionList... conds) {
		return check<E>(cond, conds...);
	}

	const char* file;
	int line;

	Precondition(const char* file, int line)
		: file(file), line(line) {}
};
#define precondition Precondition(__FILE__, __LINE__).operator()


struct out_of_range: public std::exception {};
struct bad_pointer: public std::exception {};
struct bad_size: public std::exception {};
struct zero_devide: public std::exception {};



/**
 * base_vector
 * объектно-ориентированная оболочка для gsl_vector*
 */
class base_const_vector
{
protected:
	base_const_vector() = default;
	base_const_vector(base_const_vector&&) = default;
	base_const_vector(const base_const_vector&) = default;

public:
	double at(size_t i) const throw(std::exception) {
		precondition<out_of_range>(i >= 0);
		precondition<out_of_range>(i < ptr()->size);
		return gsl_vector_get(ptr(), i);
	}

	double operator [] (size_t i) const throw (std::exception) {
		return at(i);
	}

	const double* data() const {
		return ptr()->data;
	}

	size_t size() const {
		return ptr()->size;
	}

	virtual ~base_const_vector() {}
	virtual const gsl_vector* ptr() const = 0;
};


class base_vector: public base_const_vector
{
protected:
	base_vector() = default;
	base_vector(base_vector&&) = default;
	base_vector(const base_vector&) = default;

public:
	double& at(size_t i) throw(std::exception) {
		precondition<out_of_range>(i >= 0);
		precondition<out_of_range>(i < ptr()->size);
		return *gsl_vector_ptr(ptr(), i);
	}

	double& operator [] (size_t i) throw (std::exception) {
		return at(i);
	}

	double* data() {
		return ptr()->data;
	}

	void set_zero() {
		gsl_vector_set_zero(ptr());
	}

	void set_all(double value) {
		gsl_vector_set_all(ptr(), value);
	}

	virtual ~base_vector() {}
	virtual gsl_vector* ptr() = 0;
};


class vector: public base_vector
{
public:
	vector(vector&&) = default;

	explicit vector(size_t n)
	{
		precondition<bad_size>(n > 0);
		myvec = gsl_vector_calloc(n);
		precondition<bad_pointer>(myvec);
	}

	vector(size_t n, double value)
	{
		precondition<bad_size>(n > 0);
		myvec = gsl_vector_alloc(n);
		gsl_vector_set_all(myvec, value);
		precondition<bad_pointer>(myvec);
	}

    vector(const double* array, size_t n)
	{
    	precondition<bad_pointer>(array);
    	precondition<bad_size>(n > 0);
		gsl_vector_const_view view = gsl_vector_const_view_array(array, n);
		myvec = gsl_vector_alloc(n);
		gsl_vector_memcpy(myvec, &view.vector);
		precondition<bad_pointer>(myvec);
	}

	explicit vector(const gsl_vector* vec)
	{
		precondition<bad_pointer>(vec);
		myvec = gsl_vector_alloc(vec->size);
		gsl_vector_memcpy(myvec, vec);
		precondition<bad_pointer>(myvec);
	}

	explicit vector(const gsl_vector_view& view)
	{
		myvec = gsl_vector_alloc(view.vector.size);
		gsl_vector_memcpy(myvec, &view.vector);
		precondition<bad_pointer>(myvec);
	}

	explicit vector(const gsl_vector_const_view& view)
	{
		myvec = gsl_vector_alloc(view.vector.size);
		gsl_vector_memcpy(myvec, &view.vector);
		precondition<bad_pointer>(myvec);
	}

	vector(const base_const_vector& other)
		: vector(other.ptr()) {}

	vector(const vector& other)
		: vector(other.myvec) {}

	~vector() {
		gsl_vector_free(myvec);
	}

	vector& operator = (const vector& other) {
		if( this != &other ) {
			precondition<bad_pointer>(other.myvec);
			precondition<bad_size>(myvec->size == other.myvec->size);
			gsl_vector_memcpy(myvec, other.myvec);
		}
		return *this;
	}

	gsl_vector* ptr() override {
		return myvec;
	}

	const gsl_vector* ptr() const override {
		return myvec;
	}

private:
	gsl_vector* myvec;
};



template<class V>
struct view_traits;


template<>
struct view_traits<gsl_vector_view> {

	using value_t = double;
	using vector_t = gsl_vector;
	using view_t = gsl_vector_view;

	inline static view_t view_array(value_t* arr, size_t n) {
		return gsl_vector_view_array(arr, n);
	}

	inline static view_t view_array_with_stride(value_t* arr, size_t n, size_t stride) {
		return gsl_vector_view_array_with_stride(arr, n, stride);
	}

	inline static view_t subvector(vector_t* v, size_t i, size_t n) {
		return gsl_vector_subvector(v, i, n);
	}

	inline static view_t subvector_with_stride(vector_t* v, size_t i, size_t stride, size_t n) {
		return gsl_vector_subvector_with_stride(v, i, stride, n);
	}
};


template<>
struct view_traits<gsl_vector_const_view> {

	using value_t = const double;
	using vector_t = const gsl_vector;
	using view_t = gsl_vector_const_view;

	inline static view_t view_array(value_t* arr, size_t n) {
		return gsl_vector_const_view_array(arr, n);
	}

	inline static view_t view_array_with_stride(value_t* arr, size_t n, size_t stride) {
		return gsl_vector_const_view_array_with_stride(arr, n, stride);
	}

	inline static view_t subvector(vector_t* v, size_t i, size_t n) {
		return gsl_vector_const_subvector(v, i, n);
	}

	inline static view_t subvector_with_stride(vector_t* v, size_t i, size_t stride, size_t n) {
		return gsl_vector_const_subvector_with_stride(v, i, stride, n);
	}
};


template<class View>
class base_view_t
{
protected:
	typedef view_traits<View> traits;
	typedef typename traits::value_t value_t;
	typedef typename traits::vector_t vector_t;
	typedef typename traits::view_t view_t;

	view_t make_view(vector_t* v, size_t i, size_t stride, size_t n) {
		precondition<bad_pointer>(v);
		precondition<out_of_range>(i >= 0 && i < v->size);
		precondition<out_of_range>(i + stride * (n - 1) < v->size);
		precondition<bad_size>(n > 0 && n <= v->size - i);

		return traits::subvector_with_stride(v, i, stride, n);
	}

public:
	base_view_t(const base_view_t& other)
		: view(other.view) {}

	base_view_t(const view_t& _view)
		: view(_view) {}

	base_view_t(vector_t* v)
		: view( make_view(v, 0, 1, v->size) ) {}

	base_view_t(vector_t* vector, size_t i, size_t n)
		: view( make_view(vector, i, 1, n) ) {}

	base_view_t(vector_t* vector, size_t i, size_t stride, size_t n)
		: view( make_view(vector, i, stride, n) ) {}

protected:
	view_t view;
};


class vector_const_view: public base_const_vector, public base_view_t<gsl_vector_const_view>
{
public:
	using base_view_t<gsl_vector_const_view>::base_view_t;

	vector_const_view(const base_const_vector& v)
		: base_view_t(v.ptr()) {}

	vector_const_view(const base_const_vector& v, size_t i, size_t n)
		: base_view_t(v.ptr(), i, n) {}

	vector_const_view(const base_const_vector& v, size_t i, size_t stride, size_t n)
		: base_view_t(v.ptr(), i, stride, n) {}

	const gsl_vector* ptr() const override {
		return &view.vector;
	}
};


class vector_view: public base_vector, public base_view_t<gsl_vector_view>
{
public:
	using base_view_t<gsl_vector_view>::base_view_t;

	vector_view(base_vector& v)
		: base_view_t(v.ptr()) {}

	vector_view(base_vector& v, size_t i, size_t n)
		: base_view_t(v.ptr(), i, n) {}

	vector_view(base_vector& v, size_t i, size_t stride, size_t n)
		: base_view_t(v.ptr(), i, stride, n) {}

	vector_view& operator = (const vector_view& other) {
		if (this != &other)
			view = other.view;
		return *this;
	}

	const gsl_vector* ptr() const override {
		return &view.vector;
	}

	gsl_vector* ptr() override {
		return &view.vector;
	}
};


inline std::ostream& operator << (std::ostream& output, const base_const_vector& v) {
	output << '{';
	for( size_t i=0; i<v.size(); ++i ) {
		if( i > 0 ) output << ", ";
		output << v[i];
	}
	output << '}';
	return output;
}

}//namespace gsl

#endif //__GSL_HELPER_H__
