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

	ASC_bla::Matrix<double, ASC_bla::ORDERING::ColMajor> C = A*B;
	ASC_bla::Matrix<double, ASC_bla::ORDERING::ColMajor> D = A * B + A * B;
	ASC_bla::Matrix<double, ASC_bla::ORDERING::ColMajor> E = ASC_bla::InverseLapack(A);
	ASC_bla::Matrix<double, ASC_bla::ORDERING::ColMajor> I = A * E;
	auto RowA1 = A.Row(1);
	auto ColB1 = B.Col(1);
	std::cout << "C = " << C << std::endl;
	std::cout << "D = " << D << std::endl;
	std::cout << "E = " << E << std::endl;
	std::cout << "I = " << I << std::endl;
	//std::cout << "RowA1 = " << RowA1 << std::endl;
	//std::cout << "ColB1 = " << ColB1 << std::endl;
	//std::cout << ASC_bla::InnerProduct<1>((size_t)3, RowA1, ColB1) << std::endl;

	assert(C(0, 0) == 1);
	assert(C(0, 1) == 1);
	assert(C(1, 0) == 2);
	assert(C(1, 1) == 2);
	assert(C(2, 0) == 0);
	assert(C(2, 1) == 0);
	assert(C(3, 0) == 0);
	assert(C(3, 1) == 0);
	assert(C(4, 0) == 0);
	assert(C(4, 1) == 0);
}

void test_matrix_row_col()
{
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

	
	auto RowA0 = A.Row(0);
	auto ColB0 = B.Col(0);
	
	std::cout << "RowA0 = " << RowA0 << std::endl;
	std::cout << "ColB0 = " << ColB0 << std::endl;

	assert(RowA0(0) == 1);
	assert(RowA0(1) == 0);
	assert(RowA0(2) == 0);

	assert(ColB0(0) == 1);
	assert(ColB0(1) == 2);
	assert(ColB0(2) == 0);
}

void test_matrix_vector()
{
	size_t n = 5;
	size_t m = 3;
	ASC_bla::Matrix<double, ASC_bla::ORDERING::ColMajor> A(n, m);
	ASC_bla::Vector<double> x(m);

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

	for (size_t i = 0; i < x.Size(); i++)
	{
		x(i) = i;
	}

	std::cout << "A = " << A << std::endl;
	std::cout << "x = " << x << std::endl;

	auto Ax = A * x;

	std::cout << "Ax = " << Ax << std::endl;

	assert(Ax(0) == 0);
	assert(Ax(1) == 1);
	assert(Ax(2) == 2);
	assert(Ax(3) == 0);
	assert(Ax(4) == 0);
}

void test_100_100_matmul()
{
	size_t n = 100;
	size_t m = 100;
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
			B(i, j) = i + 1;
		}
	}

	// std::cout << "A = " << A << std::endl;
	// std::cout << "B = " << B << std::endl;

	ASC_bla::Matrix<double, ASC_bla::ORDERING::ColMajor> C = A*B;
	
	std::cout << "C = " << C << std::endl;

	assert(C(0, 0) == B(0, 0));
	assert(C(0, 1) == B(0, 1));
	assert(C(0, 2) == B(0, 2));
	assert(C(0, 3) == B(0, 3));
	assert(C(0, 4) == B(0, 4));
	assert(C(0, 5) == B(0, 5));

	assert(C(1, 0) == B(1, 0));
	assert(C(1, 1) == B(1, 1));
	assert(C(1, 2) == B(1, 2));


	assert(C(20, 20) == B(20, 20));
	assert(C(20, 21) == B(20, 21));
}