#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <limits>
#include "../src/MatrixManip.hpp"

/// Support routines /////////////////////////////////////////////////////////

template<class T>
T abs(T a) {
  return (a>0)?a:-a;
}

/// TESTS ////////////////////////////////////////////////////////////////////

TEST_CASE("mateye [n=1..10]") {
  for(int n=1; n<10; n++) {
    M2 I = MatrixManip::mateye(n);
    for(int i=0; i<n; i++) {
      for(int j=0; j<n; j++) {
        if(i==j)
          CHECK(doctest::Approx(I[i][j])==1);
        else
          CHECK(doctest::Approx(I[i][j])==0);
      }
    }
  }
}

TEST_CASE("mateye [n=1..10,v=-100:0.1:100]") {
  for(int vk=-100; vk<=100; vk++) {
    datatype v = (datatype)vk/(datatype)10.0;
    for(int n=1; n<=10; n++) {
      M2 I = MatrixManip::mateye(n, v);
      for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
          if(i==j)
            CHECK(doctest::Approx(I[i][j])==v);
          else
            CHECK(doctest::Approx(I[i][j])==0);
        }
      }
    }
  }
}

TEST_CASE("mattranspose") {
  M2 M = { {1, 2, 3},
           {4, 5, 6},
           {7, 8, 9} };
  M2 T = MatrixManip::allocMatrix(3,3);
  MatrixManip::mattranspose(M, T);
  M = { {1, 4, 7},
        {2, 5, 8},
        {3, 6, 9} };
  CHECK(MatrixManip::matCompare2D(M, T));
}

TEST_CASE("cholesky") { 
  int n = 3;
  M2 M = MatrixManip::allocMatrix(n, n);
  M = { {4,12,-16},
        {12,37,-43},
        {-16,-43,98} };
  M2 L = MatrixManip::allocMatrix(n, n);
  MatrixManip::cholesky(M, L);
  M2 Lp = MatrixManip::allocMatrix(n, n);
  Lp={  {2,0,0},
        {6,1,0},
        {-8,5,3} };
  CHECK(MatrixManip::matCompare2D(L, Lp));
}
