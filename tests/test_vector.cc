#include <iostream>

#include <vector.h>
#include <matrix.h>
#include <lapack_interface.h>

namespace bla = ASC_bla;

int main()
{
  size_t n = 10;
  bla::Vector<double> x(n), y(n);

  bla::Matrix<double, bla::ORDERING::ColMajor> A(n, n);
  bla::Matrix<double, bla::ORDERING::ColMajor> B(n, n);
  bla::Matrix<double, bla::ORDERING::ColMajor> X(n, n);

  for (size_t i = 0; i < x.Size(); i++)
  {
    x(i) = i;
    y(i) = 10;
  }

  for (size_t i = 0; i < A.Size_Rows(); i++)
    for (size_t j = 0; j < A.Size_Cols(); j++)
      A(i, j) = i + j;

  for (size_t i = 0; i < B.Size_Rows(); i++)
    for (size_t j = 0; j < B.Size_Cols(); j++)
      B(i, j) = 2*i + j;

  for (size_t i = 0; i < X.Size_Rows(); i++)
    for (size_t j = 0; j < X.Size_Cols(); j++)
      X(i, j) = i;

  bla::Vector<double> z = x + y;

  bla::Matrix<double, bla::ORDERING::ColMajor> C = A*B;
  bla::Matrix<double, bla::ORDERING::ColMajor> D = A*B + A*B;
  bla::Matrix<double, bla::ORDERING::ColMajor> E = B*A;

  bla::Matrix<double, bla::ORDERING::ColMajor> Identity = bla::Matrix<double, bla::ORDERING::ColMajor>(n, n);
  Identity.Diag() = 1;
  //bla::Matrix<double, bla::ORDERING::ColMajor> E = InverseLapack(A);
  //bla::Matrix<double, bla::ORDERING::ColMajor> I = A*E;

  std::cout << "A = " << A << std::endl;
  std::cout << "B = " << B << std::endl;
  std::cout << "C = " << C << std::endl;
  std::cout << "D = " << D << std::endl;
  std::cout << "E = " << E << std::endl;
  std::cout << "I = " << Identity << std::endl;
  
  std::cout << "A*I = " << A*Identity << std::endl;
  std::cout << "I*A = " << Identity*A << std::endl;
  std::cout << "X = " << X << std::endl;
  std::cout << "A*X = " << A*X << std::endl;
  std::cout << "X*A = " << X*A << std::endl;

  std::cout << "A = " << A << std::endl;
  std::cout << "X = " << X << std::endl;

  std::cout << "X.Row(0) = " << X.Row(0) << std::endl;
  std::cout << "X.Row(1) = " << X.Row(1) << std::endl;
  std::cout << "X.Row(2) = " << X.Row(2) << std::endl;
  std::cout << "X.Row(3) = " << X.Row(3) << std::endl;
  std::cout << "X.Row(4) = " << X.Row(4) << std::endl;

  std::cout << "X.Col(0) = " << X.Col(0) << std::endl;
  std::cout << "X.Col(1) = " << X.Col(1) << std::endl;
  std::cout << "X.Col(2) = " << X.Col(2) << std::endl;
  std::cout << "X.Col(3) = " << X.Col(3) << std::endl;
  std::cout << "X.Col(4) = " << X.Col(4) << std::endl;
  std::cout << "X.Row(4).Data() = " << X.Row(4).Data()[10] << std::endl;
  std::cout << "X.Flatten() = " << X.Flatten() << std::endl;

  size_t i = 1;
  std::cout << "X.Row(i) = " << X.Row(i) << std::endl;
  std::cout << "A.Col(i) = " << A.Col(i) << std::endl;

  std::cout << "InnerProduct(X.Row(i), A.Col(i)) = " << ASC_bla::InnerProduct<8>((size_t)n, X.Row(i), A.Col(i)) << std::endl;


  std::cout << "A*x = " << A*x << std::endl;
  std::cout << "y*A = " << y*A << std::endl;

  std::cout << "x+y = " << z << std::endl;

  std::cout << "type of (x+3*y) is  " << typeid(x + 3 * y).name() << std::endl;

  std::cout << "x+3*y = " << x + 3 * y << std::endl;

  std::cout << "sizeof(x+3*y) = " << sizeof(x + 3 * y) << std::endl;

  x.Range(2, 9) = 3;
  x.Slice(1, 5) = 10;

  std::cout << "x = " << x << std::endl;
}
