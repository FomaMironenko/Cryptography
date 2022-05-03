#include <iostream>
#include <string>

#include <gmpxx.h>

int main() {
    mpz_class a, b;
    std::cout << "Enter a := ";
    std::cin  >> a;
    std::cout << "Emter b := ";
    std::cin  >> b;

    std::cout << a << " * " << b << " = " << a*b << std::endl;

    return EXIT_SUCCESS;
}