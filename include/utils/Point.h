#pragma once
#include <cstring>
#include <type_traits>

namespace scds {
using u32 = unsigned int;
using u64 = unsigned long;

/**
 * 2D-4D points
*/
template<typename Ty, u64 Ndim, typename = std::enable_if_t<Ndim >= 2 && Ndim <=4>>
class Point {
public:
    Point() { memset(data, static_cast<Ty>(0), Ndim * sizeof(Ty)); }

    // TODO: What if we take {1, 2}, 3 as its input?
    template <typename... Args>
    Point(Args&&... args) : data{std::forward<Args>(args)...} {
        static_assert(sizeof...(Args) <= Ndim, "Too many arguments.");
    }
public:
    template <typename PointType>
    Point operator+(PointType&& pt) const {
        static_assert(std::is_same_v<std::decay_t<PointType>, Point<std::decay_t<decltype(p.data[0])>, Ndim>>,
                      "PointType must have the same Ndim as the current instance.");
        Point res{};
        for (u64 i = 0; i < Ndim; i++)
            res[i] = this->data[i] + pt[i];
        return res;
    }

    Point operator-(PointType&& pt) const {
        static_assert(std::is_same_v<std::decay_t<PointType>, Point<std::decay_t<decltype(p.data[0])>, Ndim>>,
                      "PointType must have the same Ndim as the current instance.");
        Point res{};
        for (u64 i = 0; i < Ndim; i++)
            res[i] = this->data[i] - pt[i];
        return res;
    }


    // ============ indexing ================
    Ty& operator[](int idx) {
        return this->data[idx];
    }

    // TODO: can this get simplified?
    const Ty& operator[](int idx) const {
        return this->data[idx];
    }

    Ty operator[](int idx) const {
        return this->data[idx];
    }

private:
    Ty data[Ndim];
};


}   // end of namespace scds  