#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <limits>
#include "../src/MatrixManip.hpp"

template<class T>
T abs(T a) {
  return (a>0)?a:-a;
}

bool matCompare2D(M2 Ma, M2 Mb) {
  datatype tolerance = 1e-12;
  int dim1=Ma.size(), dim2=Ma[0].size();
  for(int i=0; i<dim1; i++) {
    for(int j=0; j<dim2; j++) {
      if(abs(Ma[i][j]-Mb[i][j])>tolerance) {
        return false;
      }
    }
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////

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
  int n = 3;
  M2 M = MatrixManip::allocMatrix(n, n);
  M[0][0] = 1; M[0][1] = 2; M[0][2] = 3;
  M[1][0] = 4; M[1][1] = 5; M[1][2] = 6;
  M[2][0] = 7; M[2][1] = 8; M[2][2] = 9;
  M2 T = MatrixManip::allocMatrix(n, n);
  MatrixManip::mattranspose(M, T);
  M[0][0] = 1; M[0][1] = 4; M[0][2] = 7;
  M[1][0] = 2; M[1][1] = 5; M[1][2] = 8;
  M[2][0] = 3; M[2][1] = 6; M[2][2] = 9;
  CHECK(matCompare2D(M, T));
}

TEST_CASE("cholesky") { 
  int n = 3;
  M2 M = MatrixManip::allocMatrix(n, n);
  M[0][0] = 4;   M[0][1] = 12;  M[0][2] = -16;
  M[1][0] = 12;  M[1][1] = 37;  M[1][2] = -43;
  M[2][0] = -16; M[2][1] = -43; M[2][2] = 98;
  M2 L = MatrixManip::allocMatrix(n, n);
  MatrixManip::cholesky(M, L);
  M2 Lp = MatrixManip::allocMatrix(n, n);
  Lp[0][0] = 2;  Lp[0][1] = 0; Lp[0][2] = 0;
  Lp[1][0] = 6;  Lp[1][1] = 1; Lp[1][2] = 0;
  Lp[2][0] = -8; Lp[2][1] = 5; Lp[2][2] = 3;
  CHECK(matCompare2D(L, Lp));
}
