 #include <iostream>

extern "C" {
  // remember that in fortran values are passed by reference, so we need to receive them
  // using pointers in C/C++
  void print_double_(double * num1, double * num2) {
    std::cout << "num1:" << *num1 << " and num2:" << *num2 << std::endl;
  }
  void print_int_(int * num1, int * num2) {
    std::cout << "num1:" << *num1 << " and num2:" << *num2 << std::endl;
  }
}
