/**
 * This executable file is for CMake build testing. Will be deleted in the future.
 * @date: 2023.12.19
*/
#include <iostream>

int main(int argc, char** argv) {
    std::cout << "Hello world: ";
    if (argc > 1)
        std::cout << argv[1] << std::endl;
    else
        std::cout << "Nothing." << std::endl;
    return 0;
}
