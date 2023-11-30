#include <iostream>

#include <vector.h>
#include <matrix.h>
#include <lapack_interface.h>

#include <chrono>

using namespace std;
using namespace ASC_bla;

int main()
{
    int arr[] = {10, 100, 1000};

    for (int n : arr)
    {
        cout << n << "\n";
        cout << "Own multiplication\n";
        size_t flops = n * n * n;
        size_t runs = size_t(1e9 / flops) + 1;

        Matrix<double, ORDERING::ColMajor> A(n, n);
        Matrix<double, ORDERING::ColMajor> B(n, n);
        Matrix<double, ORDERING::ColMajor> C(n, n);

        for (size_t i = 0; i < A.Size_Rows(); i++)
            for (size_t j = 0; j < A.Size_Cols(); j++)
                A(i, j) = i + j;

        for (size_t i = 0; i < B.Size_Rows(); i++)
            for (size_t j = 0; j < B.Size_Cols(); j++)
                B(i, j) = 2 * i + j;

        cout << "Starting run\n";
        auto start = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < runs; i++)
        {
            C = A * B;
        }
        auto end = std::chrono::high_resolution_clock::now();
        double time = std::chrono::duration<double>(end - start).count();

        cout << "n = " << n << ", time = " << time << " s, GFlops = "
             << (flops * runs) / time * 1e-9 << endl;

        cout << "Lapack multiplication\n";
        start = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < runs; i++)
        {
            MultMatMatLapack(A, B, C);
        }
        end = std::chrono::high_resolution_clock::now();

        time = std::chrono::duration<double>(end - start).count();

        cout << "n = " << n << ", time = " << time << " s, GFlops = "
             << (flops * runs) / time * 1e-9 << endl;
    }
}