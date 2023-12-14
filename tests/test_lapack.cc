#include <iostream>

#include <vector.h>
#include <lapack_interface.h>


using namespace ASC_bla;
using namespace std;


int main()
{
  Vector<double> x(5);
  Vector<double> y(5);

  Matrix<double, ORDERING::ColMajor> A(5, 5);

  for (int i = 0; i < x.Size(); i++)
    {
      x(i) = i;
      y(i) = 2;
    }

  for (int i = 0; i < A.Size_Rows(); i++)
  {
    for (int j = 0; j < A.Size_Cols(); j++)
    {
      if (i == j)
        A(i, j) = 1;
      else
        A(i, j) = 0;
    }
  }

  cout << "x = " << x << endl;
  cout << "y = " << y << endl;
  
  AddVectorLapack (2, x, y);  
  cout << "y+2*x = " << y << endl;

  cout << "A = " << endl << A << endl;
  cout << "A^-1 = "<< endl << Matrix(InverseLapack(A)) << endl;
}

  
