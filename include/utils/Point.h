#pragma once
#include <cstring>
#include <numeric>
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

#define POINT_REDUCE_OPERATOR_OVERLOAD(operatorType, Ty, Ndim) \
constexpr Ty operatorType () const { \
    Ty ret = _data[0]; \
    for (size_t i = 1; i < Ndim; i++) \
        ret = std::operatorType(ret, _data[i]); \
    return ret; \
}

// Declaring this function as constexpr is useless in C++17, since the rules for constexpr are too restrictive
// We can not have more than one return position in the function
#define POINT_ANY_OPERATOR_OVERLOAD(operatorType, Ty, Ndim) \
bool any_##operatorType () const { \
    for (size_t i = 1; i < Ndim; i++) { \
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

    POINT_ACCESS_OVERLOAD(x, 0)
    POINT_ACCESS_OVERLOAD(y, 1)
    POINT_ACCESS_OVERLOAD(z, 2)
    POINT_ACCESS_OVERLOAD(w, 3)

    POINT_REDUCE_OPERATOR_OVERLOAD(max, Ty, Ndim)
    POINT_REDUCE_OPERATOR_OVERLOAD(min, Ty, Ndim)
    POINT_ANY_OPERATOR_OVERLOAD(inf, Ty, Ndim)
    POINT_ANY_OPERATOR_OVERLOAD(nan, Ty, Ndim)

    template <typename PointType>
    constexpr Ty dot(PointType&& pt) const {
        static_assert(std::is_same_v<std::decay_t<PointType>, Point<std::decay_t<decltype(pt._data[0])>, Ndim>>,
                        "PointType must have the same Ndim as the current instance.");
        return ((*this) * pt).sum();
    }

    constexpr Ty sum() const {
        Ty sum{0};
        for (size_t i = 0; i < Ndim; i++)
            sum += _data[i];
        return sum;
    }

    constexpr Ty prod() const {
        Ty prod{1};
        for (size_t i = 0; i < Ndim; i++)
            prod *= _data[i];
        return prod;
    }

    constexpr Ty length2() const {
        Ty sum = 0;
        for (size_t i = 0; i < Ndim; i++)
            sum += _data[i] * _data[i];
        return sum;
    }

    constexpr Ty length() const {
        return std::sqrt(this->length2());
    }

    Point& normalize() {
        (*this) /= this->length();
        return *this;
    }

    constexpr Point& normalized() const {
        Ty len = this->length();
        return (*this) / len;
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
    constexpr void assign(PointType&& pt) {
        constexpr size_t incoming_ndim = std::min(pt.ndim(), Ndim);
        for (size_t i = 0; i < incoming_ndim; ++i) {
            _data[i] = pt[i];
        }
    }

    constexpr auto expand(Ty val) const {
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

    constexpr size_t ndim() const {
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

#undef POINT_ASSESS_OVERLOAD
#undef POINT_BINARY_OPERATOR_OVERLOAD
#undef POINT_BINARY_INPLACE_OVERLOAD
#undef POINT_UNARY_OPERATOR_OVERLOAD
#undef POINT_UNARY_INPLACE_OVERLOAD
#undef POINT_UNARY_OPERATOR_INVERSE
#undef POINT_REDUCE_OPERATOR_OVERLOAD
#undef POINT_ANY_OPERATOR_OVERLOAD
#undef POINT_TYPE_DEF



}   // end of namespace scds  