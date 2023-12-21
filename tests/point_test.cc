#include <iostream>
#include "utils/Point.h"
#include "utils/utils.h"

using namespace scds;

int main()
{   
    Point3f p1{1, 2, 3};
    Point3f p2{2, 2, 4};
    Point3f p3{0, 3, 1};
    std::cout << "Floating point test:\n";
    PRINT_OPS(p1)
    PRINT_OPS(p2)
    PRINT_OPS((p1 - p2).abs())
    PRINT_OPS((p1 - p2).abs())
    PRINT_OPS((p1 >= 2))
    PRINT_OPS(p1.length2())
    PRINT_OPS(p2.length())
    PRINT_OPS(1 / p2)
    PRINT_OPS(p2.max())
    PRINT_OPS(p2.any_inf())
    PRINT_OPS(p2.any_nan())
    PRINT_OPS(p2.expand(3))
    PRINT_OPS(p2.max())
    PRINT_OPS(p2.min())
    PRINT_OPS(p1.sum())
    PRINT_OPS(p1.prod())
    PRINT_OPS(((p1 - p2) >= 2))
    PRINT_OPS(((p1 - p2) >= 2).any())
    PRINT_OPS(((p1 - p2) >= 2).all())
    PRINT_OPS(p2.dot(p1));
    PRINT_OPS(p1.normalized());
    PRINT_OPS(p1 / p2);
    PRINT_OPS(p1.to_int());
    PRINT_OPS(p1.to_double());
    std::cout << "non-inplace operation:\n";
    auto temp = p1.max(p3);
    temp -= int(3);
    PRINT_OPS(temp);
    PRINT_OPS(p1);
    std::cout << "inplace operation:\n";
    p1.max_inplace(p3);
    PRINT_OPS(p1);

    std::cout << "Integral and constexpr test:\n";

    constexpr Point3i p4{1, 2, 3};
    constexpr Point3i p5{-1, -2, -3};
    PRINT_OPS(p4 + p5);
    PRINT_OPS(p4 * p5);
    PRINT_OPS(p4.length2());
    PRINT_OPS(p4.expand(4));
    PRINT_OPS(p5.to_u32());
    return 0;
}