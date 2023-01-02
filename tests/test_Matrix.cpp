#include "Matrix.hpp"
#include<string>

void UnitTests_Matrix() {
  Matrix<double, 1> md1({4});
  std::cout << md1;
  
  Matrix<double, 2> md2({4,5});
  std::cout << md2;

  Matrix<double, 3> md3({4,5,2});
  std::cout << md3;
  
  Matrix<int, 1> mi1({4});
  std::cout << mi1;
  
  Matrix<int, 2> mi2({4,5});
  std::cout << mi2;

  Matrix<int, 3> mi3({4,5,2});
  std::cout << mi3;
  
  Matrix<std::string, 3> ms3({4,5,2});
  std::cout << ms3;
}

int main() {
  UnitTests_Matrix();
}
