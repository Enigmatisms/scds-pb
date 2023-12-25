#pragma once
#include <cmath>
#include <cstring>
#include <numeric>
#include <iostream>
#include <type_traits>

namespace scds {

#define POINT_TYPE_DEF(ndim) \
    using Point##ndim##f = Point<float, ndim>;  \
    using Point##ndim##d = Point<double, ndim>; \
    using Point##ndim##i = Point<int, ndim>;    \
    using Point##ndim##u = Point<unsigned int, ndim>;

/**
 * Operator definer
*/
#define POINT_BINARY_OPERATOR_OVERLOAD(operatorType, Ty, Ndim) \
template <typename PointType, typename = std::enable_if_t<!std::is_arithmetic_v<PointType>, int>> \
constexpr Point operator operatorType (PointType&& pt) const { \
    static_assert(std::is_same_v<std::decay_t<PointType>, Point<std::decay_t<decltype(pt._data[0])>, Ndim>>, \
                      "PointType must have the same Ndim as the current instance."); \
    Point<Ty, Ndim> res{}; \
    for (size_t i = 0; i < Ndim; i++) \
        res[i] = this->_data[i] operatorType pt[i]; \
    return res; \
}

// '##' in macro is for concatenation, it can not only concatenate operators, but variable names, like var + 1 -> var1
#define POINT_BINARY_INPLACE_OVERLOAD(operatorType, Ty, Ndim) \
template <typename PointType, typename = std::enable_if_t<!std::is_arithmetic_v<PointType>, int>> \
Point& operator operatorType##= (PointType&& pt) { \
    static_assert(std::is_same_v<std::decay_t<PointType>, Point<std::decay_t<decltype(pt._data[0])>, Ndim>>, \
                      "PointType must have the same Ndim as the current instance."); \
    for (size_t i = 0; i < Ndim; i++) \
        this->_data[i] operatorType##= pt[i]; \
    return *this; \
}

// Since we are using the universal reference T&&, for compilation, the compiler will first try to call the unary operator
// even if the input val is of a PointType (so the binary operator is hidden), which perplexes the compiler
// so we will enable unary operators only when T is of a arithmetic type
#define POINT_UNARY_OPERATOR_OVERLOAD(operatorType, Ty, Ndim) \
template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0> \
constexpr Point operator operatorType (T&& val) const { \
    Point<Ty, Ndim> res{}; \
    for (size_t i = 0; i < Ndim; i++) \
        res[i] = this->_data[i] operatorType static_cast<Ty>(val); \
    return res; \
}

#define POINT_UNARY_INPLACE_OVERLOAD(operatorType, Ty, Ndim) \
template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0> \
Point& operator operatorType##= (T&& val) { \
    for (size_t i = 0; i < Ndim; i++) \
        _data[i] operatorType##= static_cast<Ty>(val); \
    return *this; \
}

#define POINT_STD_REDUCE_OPERATOR_OVERLOAD(operatorType, Ty, Ndim) \
constexpr Ty operatorType () const { \
    Ty ret = _data[0]; \
    for (size_t i = 1; i < Ndim; i++) \
        ret = std::operatorType(ret, _data[i]); \
    return ret; \
}

#define POINT_REDUCE_OPERATOR_OVERLOAD(function, operatorType, Ty, Ndim, init_val) \
constexpr Ty function () const { \
    Ty ret = init_val; \
    for (size_t i = 0; i < Ndim; i++) \
        ret operatorType##= _data[i]; \
    return ret; \
}

// Declaring this function as constexpr is useless in C++17, since the rules for constexpr are too restrictive
// We can not have more than one return position in the function
#define POINT_ANY_OPERATOR_OVERLOAD(operatorType, Ty, Ndim) \
bool any_##operatorType () const { \
    for (size_t i = 0; i < Ndim; i++) { \
        bool flag = std::is##operatorType(_data[i]); \
        if (flag) return true; \
    } \
    return false; \
}

#define POINT_ACCESS_OVERLOAD(dimensionName, index) \
constexpr Ty& dimensionName() {\
    static_assert(Ndim >= (index + 1), "Invalid dimension for " #dimensionName "()");\
    return _data[index];\
}\
constexpr const Ty& dimensionName() const {\
    static_assert(Ndim >= (index + 1), "Invalid dimension for " #dimensionName "()");\
    return _data[index];\
}

#define POINT_BINARY_COMPARE_OVERLOAD(operatorType, Ty, Ndim) \
template <typename PointType, typename = std::enable_if_t<!std::is_arithmetic_v<PointType>, int>> \
constexpr Point<bool, Ndim> operator operatorType (PointType&& pt) const { \
    static_assert(std::is_same_v<std::decay_t<PointType>, Point<std::decay_t<decltype(pt._data[0])>, Ndim>>, \
                      "PointType must have the same Ndim as the current instance."); \
    Point<bool, Ndim> res{}; \
    for (size_t i = 0; i < Ndim; i++) { \
        res[i] = _data[i] operatorType static_cast<Ty>(pt[i]); \
    } \
    return res; \
}

#define POINT_UNARY_COMPARE_OVERLOAD(operatorType, Ty, Ndim) \
template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0> \
constexpr Point<bool, Ndim> operator operatorType (T&& val) const { \
    Point<bool, Ndim> res{}; \
    for (size_t i = 0; i < Ndim; i++) { \
        res[i] = _data[i] operatorType static_cast<Ty>(val); \
    } \
    return res; \
}

#define POINT_TYPE_CONVERT_OVERLOAD(typeName, type, Ndim) \
constexpr Point<type, Ndim> to_##typeName() const { \
    Point<type, Ndim> res{}; \
    for (size_t i = 0; i < Ndim; i++) { \
        res[i] = static_cast<type>(this->_data[i]); \
    } \
    return res; \
}

#define POINT_EWISE_OVERLOAD(operationName, type, Ndim) \
template <typename PointType> \
constexpr Point<type, Ndim> operationName(PointType&& pt) const { \
    static_assert(std::is_same_v<std::decay_t<PointType>, Point<std::decay_t<decltype(pt._data[0])>, Ndim>>, \
                      "PointType must have the same Ndim as the current instance."); \
    Point<type, Ndim> res{}; \
    for (size_t i = 0; i < Ndim; i++) { \
        res[i] = std::operationName(this->_data[i], static_cast<type>(pt[i])); \
    } \
    return res; \
} \
template <typename PointType> \
Point<type, Ndim>& operationName##_inplace(PointType&& pt) { \
    static_assert(std::is_same_v<std::decay_t<PointType>, Point<std::decay_t<decltype(pt._data[0])>, Ndim>>, \
                      "PointType must have the same Ndim as the current instance."); \
    for (size_t i = 0; i < Ndim; i++) { \
        this->_data[i] = std::operationName(this->_data[i], static_cast<type>(pt[i])); \
    } \
    return *this; \
}

/**
 * 2D-4D points
*/
template<typename Ty, size_t Ndim, typename = std::enable_if_t<Ndim >= 2 && Ndim <=4>>
class Point {
public:
    constexpr Point() { memset(_data, static_cast<Ty>(0), Ndim * sizeof(Ty)); }

    // TODO: What if we take {1, 2}, 3 as its input?
    template <typename... Args>
    constexpr Point(Args&&... args) : _data{static_cast<Ty>(std::forward<Args>(args))...} {
        static_assert(sizeof...(Args) <= Ndim, "Too many arguments.");
        
    }

    template <typename PointType>
    Point(PointType&& pt) {
        for (size_t i = 0; i < Ndim; i++)
            _data[i] = static_cast<Ty>(pt[i]);
    }
public:
    POINT_BINARY_OPERATOR_OVERLOAD(+, Ty, Ndim)
    POINT_BINARY_OPERATOR_OVERLOAD(-, Ty, Ndim)
    POINT_BINARY_OPERATOR_OVERLOAD(*, Ty, Ndim)
    POINT_BINARY_OPERATOR_OVERLOAD(/, Ty, Ndim)   

    POINT_BINARY_INPLACE_OVERLOAD(+, Ty, Ndim)
    POINT_BINARY_INPLACE_OVERLOAD(-, Ty, Ndim)
    POINT_BINARY_INPLACE_OVERLOAD(*, Ty, Ndim)
    POINT_BINARY_INPLACE_OVERLOAD(/, Ty, Ndim)  

    POINT_UNARY_OPERATOR_OVERLOAD(+, Ty, Ndim)
    POINT_UNARY_OPERATOR_OVERLOAD(-, Ty, Ndim)
    POINT_UNARY_OPERATOR_OVERLOAD(*, Ty, Ndim)
    POINT_UNARY_OPERATOR_OVERLOAD(/, Ty, Ndim)

    POINT_UNARY_INPLACE_OVERLOAD(+, Ty, Ndim)
    POINT_UNARY_INPLACE_OVERLOAD(-, Ty, Ndim)
    POINT_UNARY_INPLACE_OVERLOAD(*, Ty, Ndim)
    POINT_UNARY_INPLACE_OVERLOAD(/, Ty, Ndim)

    POINT_BINARY_COMPARE_OVERLOAD(>,  Ty, Ndim)
    POINT_BINARY_COMPARE_OVERLOAD(>=, Ty, Ndim)
    POINT_BINARY_COMPARE_OVERLOAD(<,  Ty, Ndim)
    POINT_BINARY_COMPARE_OVERLOAD(<=, Ty, Ndim)

    POINT_UNARY_COMPARE_OVERLOAD(>,   Ty, Ndim)
    POINT_UNARY_COMPARE_OVERLOAD(>=,  Ty, Ndim)
    POINT_UNARY_COMPARE_OVERLOAD(<,   Ty, Ndim)
    POINT_UNARY_COMPARE_OVERLOAD(<=,  Ty, Ndim)

    POINT_ACCESS_OVERLOAD(x, 0)
    POINT_ACCESS_OVERLOAD(y, 1)
    POINT_ACCESS_OVERLOAD(z, 2)
    POINT_ACCESS_OVERLOAD(w, 3)

    // explicit type conversion
    POINT_TYPE_CONVERT_OVERLOAD(int,    int,        Ndim)
    POINT_TYPE_CONVERT_OVERLOAD(bool,   bool,       Ndim)
    POINT_TYPE_CONVERT_OVERLOAD(float,  float,      Ndim)
    POINT_TYPE_CONVERT_OVERLOAD(short,  short,      Ndim)
    POINT_TYPE_CONVERT_OVERLOAD(double, double,     Ndim)
    POINT_TYPE_CONVERT_OVERLOAD(i64,    long,       Ndim)
    POINT_TYPE_CONVERT_OVERLOAD(u64,    uint64_t,   Ndim)
    POINT_TYPE_CONVERT_OVERLOAD(u32,    uint32_t,   Ndim)
    POINT_TYPE_CONVERT_OVERLOAD(u16,    unsigned short, Ndim)

    POINT_STD_REDUCE_OPERATOR_OVERLOAD(max, Ty, Ndim)
    POINT_STD_REDUCE_OPERATOR_OVERLOAD(min, Ty, Ndim)

    POINT_REDUCE_OPERATOR_OVERLOAD(sum,  +, Ty, Ndim, 0)
    POINT_REDUCE_OPERATOR_OVERLOAD(prod, *, Ty, Ndim, 1)
    POINT_REDUCE_OPERATOR_OVERLOAD(any,  |, Ty, Ndim, false)
    POINT_REDUCE_OPERATOR_OVERLOAD(all,  &, Ty, Ndim, true)

    POINT_ANY_OPERATOR_OVERLOAD(inf, Ty, Ndim)
    POINT_ANY_OPERATOR_OVERLOAD(nan, Ty, Ndim)

    POINT_EWISE_OVERLOAD(max, Ty, Ndim)
    POINT_EWISE_OVERLOAD(min, Ty, Ndim)

    

    template <typename PointType>
    constexpr Ty dot(PointType&& pt) const {
        static_assert(std::is_same_v<std::decay_t<PointType>, Point<std::decay_t<decltype(pt._data[0])>, Ndim>>,
                        "PointType must have the same Ndim as the current instance.");
        return ((*this) * pt).sum();
    }

    constexpr Ty length2() const {
        Ty sum = 0;
        for (size_t i = 0; i < Ndim; i++)
            sum += _data[i] * _data[i];
        return sum;
    }

    constexpr Point abs() const {
        Point<Ty, Ndim> res{};
        for (size_t i = 0; i < Ndim; i++)
            res[i] = std::abs(this->_data[i]);
        return res;
    }

    constexpr Ty length() const {
        return std::sqrt(this->length2());
    }

    Point& normalize() {
        (*this) /= this->length();
        return *this;
    }

    constexpr Point normalized() const {
        return (*this) / this->length();
    }

    friend constexpr std::ostream& operator<<(std::ostream& os, const Point<Ty, Ndim>& point) {
        os << "(";
        for (size_t i = 0; i < Ndim; ++i) {
            os << point._data[i];
            if (i < Ndim - 1)
                os << ", ";
        }
        os << ")";
        return os;
    }


    // =========== expand dimension by one (homogeneous coordinates) ==========
    template<typename PointType>
    void assign(PointType&& pt) {
        size_t incoming_ndim = std::min(pt.ndim(), Ndim);
        for (size_t i = 0; i < incoming_ndim; ++i) {
            _data[i] = pt[i];
        }
    }

    auto expand(Ty val) const {
        if constexpr (Ndim == 4) return *this;
        else {
            Point<Ty, Ndim + 1> result;
            result.assign(*this);
            result[Ndim] = val;
            return result;
        }
    }

    template<typename T>
    void fill(T val) {
        for (size_t i = 0; i < Ndim; ++i)
            _data[i] = static_cast<Ty>(val);
    }

    template<typename T>
    static Point from_pointer(const T* const ptr) {
        Point res{};
        for (size_t i = 0; i < Ndim; ++i)
            res[i] = static_cast<Ty>(ptr[i]);
        return res;
    }

    // ============ indexing ================
    constexpr Ty& operator[](int idx) {
        return this->_data[idx];
    }

    constexpr const Ty& operator[](int idx) const {
        return this->_data[idx];
    }

    constexpr const Ty* const_data() const {
        return _data;
    }

    constexpr Ty* data() const {
        return _data;
    }

    constexpr size_t ndim() const noexcept {
        return Ndim;
    }
private:
    Ty _data[Ndim];
};

template<typename T>
using Point2 = Point<T, 2>;
template<typename T>
using Point3 = Point<T, 3>;
template<typename T>
using Point4 = Point<T, 4>;

template <typename T>
struct __is_point_type : std::false_type {};
template <typename Ty, std::size_t Ndim>
struct __is_point_type<Point<Ty, Ndim>> : std::true_type {};

// For point type checking, will be true if T is of type Point<Ty, Ndim>
#define STATIC_ASSERT_POINT_TYPE(T) \
    static_assert(__is_point_type<std::decay_t<T>>::value, "Input type is not a Point type.")

/**
 * Directly using std::enable_if_t<std::is_arithmetic_v<T>> might fail. Compiler might not be able to
 * deduct type T due to SFINAE (make type deduction more difficult)
 * This is actually a overload for class member operator \ 
*/
#define POINT_UNARY_OPERATOR_INVERSE(operatorType) \
template <typename T = float, typename Ty, size_t Ndim, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0> \
constexpr Point<Ty, Ndim> operator operatorType (T val, const Point<Ty, Ndim>& pt) { \
    Point<Ty, Ndim> res{}; \
    for (size_t i = 0; i < Ndim; i++) \
        res[i] = static_cast<Ty>(val) operatorType pt[i]; \
    return res; \
}

POINT_UNARY_OPERATOR_INVERSE(+)
POINT_UNARY_OPERATOR_INVERSE(-)
POINT_UNARY_OPERATOR_INVERSE(*)
POINT_UNARY_OPERATOR_INVERSE(/)

POINT_TYPE_DEF(2)
POINT_TYPE_DEF(3)
POINT_TYPE_DEF(4)

#undef POINT_TYPE_DEF
#undef POINT_EWISE_OVERLOAD
#undef POINT_ASSESS_OVERLOAD
#undef POINT_ANY_OPERATOR_OVERLOAD
#undef POINT_TYPE_CONVERT_OVERLOAD
#undef POINT_UNARY_OPERATOR_INVERSE
#undef POINT_UNARY_INPLACE_OVERLOAD
#undef POINT_UNARY_COMPARE_OVERLOAD
#undef POINT_BINARY_COMPARE_OVERLOAD
#undef POINT_BINARY_INPLACE_OVERLOAD
#undef POINT_UNARY_OPERATOR_OVERLOAD
#undef POINT_BINARY_OPERATOR_OVERLOAD
#undef POINT_REDUCE_OPERATOR_OVERLOAD
#undef POINT_STD_REDUCE_OPERATOR_OVERLOAD

}   // end of namespace scds  