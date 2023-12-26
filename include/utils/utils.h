#pragma once
#include <iostream>
#include <fmt/core.h>

namespace scds {

#define PRINT_OPS(ops) std::cout << #ops << ":\t" << ops << std::endl;

#define SCDS_RUNTIME_ERROR(...) throw std::runtime_error(fmt::format(__VA_ARGS__))
    
} // namespace scds
