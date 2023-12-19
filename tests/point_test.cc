#include <iostream>
#include "utils/Point.h"

#define PRINT_OPS(ops) \
    std::cout << #ops << ":\t\t" << ops << std::endl;

int main()
{   
    scds::Point3f p1{1, 2, 3};
    scds::Point3f p2{2, 2, 4};
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
    PRINT_OPS(p1.to_bool());

    std::cout << "Integral and constexpr test:\n";

    constexpr scds::Point3i p4{1, 2, 3};
    constexpr scds::Point3i p5{-1, -2, -3};
    PRINT_OPS(p4 + p5);
    PRINT_OPS(p4 * p5);
    PRINT_OPS(p4.length2());
    PRINT_OPS(p4.expand(4));
    PRINT_OPS(p5.to_u32());
    return 0;
}