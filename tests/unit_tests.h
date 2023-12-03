#include <cassert>
#include <iostream>

#include <vector.h>
#include <matrix.h>
#include <lapack_interface.h>

void test_pass() {
    assert(0 == 0);
}

void test_vector() {
	size_t n = 5;
	ASC_bla::Vector<double> x(n), y(n);

	for (size_t i = 0; i < x.Size(); i++)
	{
		x(i) = i;
		y(i) = 10;
	}

	ASC_bla::Vector<double> z = x + y;

	assert(z(0) == 10);
	assert(z(1) == 11);
	assert(z(2) == 12);
	assert(z(3) == 13);
	assert(z(4) == 14);
}

void test_matrix() {
	size_t n = 5;
	size_t m = 3;
	ASC_bla::Matrix<double, ASC_bla::ORDERING::ColMajor> A(n, m);
	ASC_bla::Matrix<double, ASC_bla::ORDERING::ColMajor> B(m, n);

	for (size_t i = 0; i < A.Size_Rows(); i++)
	{
		for (size_t j = 0; j < A.Size_Cols(); j++)
		{
			if (i == j)
				A(i, j) = 1;
			else
				A(i, j) = 0;
		}
	}

	for (size_t i = 0; i < B.Size_Rows(); i++)
	{
		for (size_t j = 0; j < B.Size_Cols(); j++)
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

	std::cout << "A = " << A << std::endl;
	std::cout << "B = " << B << std::endl;

	ASC_bla::Matrix<double, ASC_bla::ORDERING::ColMajor> C = B*A;
	ASC_bla::Matrix<double, ASC_bla::ORDERING::ColMajor> D = A * B + A * B;
	ASC_bla::Matrix<double, ASC_bla::ORDERING::ColMajor> E = ASC_bla::InverseLapack(A);
	ASC_bla::Matrix<double, ASC_bla::ORDERING::ColMajor> I = A * E;
	std::cout << "C = " << C << std::endl;
	std::cout << "D = " << D << std::endl;
	std::cout << "E = " << E << std::endl;
	std::cout << "I = " << I << std::endl;

	assert(C(0, 0) == 1);
	assert(C(0, 1) == 1);
	assert(C(1, 0) == 2);
	assert(C(1, 1) == 2);
	assert(C(2, 0) == 0);
}
