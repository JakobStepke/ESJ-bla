#include <iostream>

#include <vector.h>
#include <matrix.h>
#include <lapack_interface.h>

namespace bla = ASC_bla;

int main()
{
    size_t n = 5;
    size_t m = 3;
    bla::Vector<double> x(n), y(n);

    bla::Matrix<double, bla::ORDERING::ColMajor> A(n, m);
    bla::Matrix<double, bla::ORDERING::ColMajor> B(m, n);

    for (size_t i = 0; i < x.Size(); i++)
    {
        x(i) = i;
        y(i) = 10;
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

    for (int i = 0; i < B.Size_Rows(); i++)
    {
        for (int j = 0; j < B.Size_Cols(); j++)
        {
            if (i == 0)
                B(i, j) = 1;
            else
                if (i == 1)
                    B(i, j) = 2;
                else
                    B(i, j) = 0;
        }
    }

    bla::Vector<double> z = x + y;

    bla::Matrix<double, bla::ORDERING::ColMajor> C = A * B;
    bla::Matrix<double, bla::ORDERING::ColMajor> D = A * B + A * B;
    bla::Matrix<double, bla::ORDERING::ColMajor> E = InverseLapack(A);
    bla::Matrix<double, bla::ORDERING::ColMajor> I = A * E;

    std::cout << "A = " << A << std::endl;
    std::cout << "B = " << B << std::endl;
    std::cout << "C = " << C << std::endl;
    std::cout << "D = " << D << std::endl;
    std::cout << "E = " << E << std::endl;
    std::cout << "I = " << I << std::endl;

    std::cout << "A*x = " << A * x << std::endl;
    std::cout << "y*A = " << y * A << std::endl;

    std::cout << "x+y = " << z << std::endl;

    std::cout << "type of (x+3*y) is  " << typeid(x + 3 * y).name() << std::endl;

    std::cout << "x+3*y = " << x + 3 * y << std::endl;

    std::cout << "sizeof(x+3*y) = " << sizeof(x + 3 * y) << std::endl;

    x.Range(2, 9) = 3;
    x.Slice(1, 5) = 10;

    std::cout << "x = " << x << std::endl;
}
