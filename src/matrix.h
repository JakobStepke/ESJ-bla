#ifndef FILE_MATRIX_H
#define FILE_MATRIX_H

#include <iostream>

#include "expression.h"

#include "vector.h"
#include "ordering.h"

#include "simd.h"

using namespace ASC_HPC;

namespace ASC_bla
{

    template <typename T, ORDERING ORD>
    class MatrixView : public MatExpr<MatrixView<T, ORD>, ORD>
    {
    protected:
        size_t rows_;
        size_t cols_;
        size_t dist_;
        T *data_;

    public:
        MatrixView(size_t rows, size_t cols, T *data)
            : rows_(rows), cols_(cols), data_(data)
            , dist_((ORD == ORDERING::ColMajor) ? rows : cols)
            {
         }

        template <typename TB, ORDERING ORDB>
        MatrixView &operator=(const MatExpr<TB, ORDB> &m2)
        {
            rows_ = m2.Size_Cols();
            cols_ = m2.Size_Cols();
            data_ = m2.Data();

            return *this;
        }

        MatrixView &operator=(T scal)
        {
            rows_ = scal;
            cols_ = scal;
            data_ = 0;
            return *this;
        }

        auto View() const { return MatrixView(rows_, cols_, data_); }
        size_t Size_Cols() const { return cols_; }
        size_t Size_Rows() const { return rows_; }
        size_t Dist() const { return dist_; }

        T &operator()(size_t i, size_t j)
        {
            if constexpr (ORD == ORDERING::ColMajor)
            {
                return data_[i+j * rows_];
            }
            else
            {
                return data_[j+i * cols_];
            }
        }
        const T &operator()(size_t i, size_t j) const
        {
            if constexpr (ORD == ORDERING::ColMajor)
            {
                return data_[i+j * rows_];
            }
            else
            {
                return data_[j+i * cols_];
            }
        }

        auto Transpose() const
        {
            if constexpr (ORD == ORDERING::ColMajor)
            {
                return MatrixView<T, ORDERING::RowMajor>(cols_, rows_, data_);
            }
            else
            {
                return MatrixView<T, ORDERING::ColMajor>(cols_, rows_, data_);
            }
        }

        auto Data() const
        {
            return data_;
        }

        auto Row(size_t i) const
        {
            if constexpr (ORD == ORDERING::ColMajor)
            {
                auto data_tmp = data_ + i * cols_;

                return VectorView(cols_, 1, data_tmp);
            }
            else
            {
                return Transpose().Row(i);
            }
        }

        auto Col(size_t i) const
        {
            if constexpr (ORD == ORDERING::ColMajor)
            {
                return Row(i);
            }
            else
            {
                return Transpose().Col(i);
            }
        }

        auto Rows(size_t first, size_t next) const
        {
            if constexpr (ORD == ORDERING::ColMajor)
            {
                return MatrixView(next - first, cols_, data_ + first);
            }
            else
            {
                return MatrixView(next - first, cols_, data_ + first * cols_);
            }
        }

        auto Cols(size_t first, size_t next) const
        {
            if constexpr (ORD == ORDERING::ColMajor)
            {
                return MatrixView(rows_, next - first, data_ + first * rows_);
            }
            else
            {
                return MatrixView(rows_, next - first, data_ + first);
            }
        }

        auto Size() const
        {
            return rows_ * cols_;
        }
    };

    template <typename T, ORDERING ORD>
    class Matrix : public MatrixView<T, ORD>
    {

        using MatrixView<T, ORD>::rows_;
        using MatrixView<T, ORD>::cols_;
        using MatrixView<T, ORD>::data_;
        using MatrixView<T, ORD>::dist_;
        using MatrixView<T, ORD>::Size;
        using MatrixView<T, ORD>::Row;
        using MatrixView<T, ORD>::Col;

     public:
        using MatrixView<T, ORD>::operator();


    public:
        // Constructors
        Matrix(size_t rows, size_t cols) : MatrixView<T, ORD>(rows, cols, new T[rows*cols]) {}

        // Copy constructor
        Matrix(const Matrix &other) : Matrix (other.rows_, other.cols_)
        {
            for (size_t i = 0; i < rows_; ++i)
            {
                for (size_t j = 0; j < cols_; ++j)
                {
                    (*this)(i, j) = other(i, j);
                }
            }
        }

        // Move constructor
        Matrix(Matrix &&other) noexcept : MatrixView<T, ORD>(other.rows_, other.cols_, other.data_) {}

        // Destructor
        ~Matrix() = default;

        // Assignment operator
        Matrix &operator=(const Matrix &other)
        {
            if (this != &other)
            {
                rows_ = other.rows_;
                cols_ = other.cols_;
                data_ = other.data_;
            }
            return *this;
        }

        // Output stream operator for easy printing
        friend std::ostream &operator<<(std::ostream &os, const Matrix &matrix)
        {
            for (size_t i = 0; i < matrix.rows_; ++i)
            {
                if (i > 0)
                {
                    os << "\n";
                }
                for (size_t j = 0; j < matrix.cols_; ++j)
                {
                    if (j > 0)
                    {
                        os << " ";
                    }
                    if constexpr (ORD == ORDERING::ColMajor)
                    {
                        os << matrix.Data()[i + j * matrix.rows_];
                    }
                    else
                    {
                        os << matrix.Data()[j + i * matrix.cols_];
                    }
                }
            }
            return os;
        }
        // Matrix-Matrix Multiplication
        Matrix operator*(const Matrix &other) const
        {
            if (Size() != other.Size())
            {
                // Invalid multiplication, return an empty matrix or throw an exception
                throw std::invalid_argument("Invalid multiplication");
            }

            Matrix result(rows_, other.cols_);

            for (size_t i = 0; i < rows_; ++i)
            {
                auto row = std::remove_const_t<decltype(Row(0))>(Col(i));
                for (size_t j = 0; j < other.cols_; ++j)
                {
                    const int dy = 1;
                    const size_t SW = 16;

                    if constexpr (std::is_same<double, T>::value)
                    {
                        /*
                        for (size_t k = 0; k < cols_; ++k)
                        {
                            sum = FMA(SIMD<double,SW>(px[k]), SIMD<double,SW>(py+k*dy), sum);
                        }*/
                        auto col = std::remove_const_t<decltype(Col(0))>(other.Col(j));

                        

                        result(i, j) = InnerProduct<SW>(cols_, row, 1, col, 1);
                    }
                    else{
                        T sum = 0;
                        for (size_t k = 0; k < cols_; ++k)
                        {
                            if constexpr (ORD == ORDERING::ColMajor)
                            {
                                sum += (*this)(i, k) * other(k, j);
                            }
                            else
                            {
                                sum += (*this)(i, k) * other(k, j);
                            }
                        }
                        result(i, j) = sum;
                    }
                }
            }
            return result;
        }

        Matrix operator*(const T &scal) const
        {
            Matrix result(rows_, cols_);

            for (size_t i = 0; i < rows_; ++i)
            {
                for (size_t j = 0; j < cols_; ++j)
                {
                    if constexpr (ORD == ORDERING::ColMajor)
                    {
                        result(i, j) = (*this)(i, j) * scal;
                    }
                    else
                    {
                        result(i, j) = (*this)(i, j) * scal;
                    }
                }
            }
            return result;
        }

        // Matrix-Matrix Addition
        Matrix operator+(const Matrix &other) const
        {
            if (rows_ != other.rows_ || cols_ != other.cols_)
            {
                // Invalid addition, return an empty matrix or throw an exception
                throw std::invalid_argument("Invalid addition");
            }

            Matrix result(rows_, cols_);

            for (size_t i = 0; i < rows_; ++i)
            {
                for (size_t j = 0; j < cols_; ++j)
                {
                    if constexpr (ORD == ORDERING::ColMajor)
                    {
                        result(i, j) = (*this)(i, j) + other(i, j);
                    }
                    else
                    {
                        result(i, j) = (*this)(i, j) + other(i, j);
                    }
                }
            }
            return result;
        }


        // Transpose method
        Matrix<T, ORD> transpose() const
        {
            Matrix<T, ORD> result(cols_, rows_);

            for (size_t i = 0; i < rows_; ++i)
            {
                for (size_t j = 0; j < cols_; ++j)
                {
                    if constexpr (ORD == ORDERING::ColMajor)
                    {
                        result(j, i) = (*this)(i, j);
                    }
                    else
                    {
                        result(j, i) = (*this)(i, j);
                    }
                }
            }
            return result;
        }
        // Matrix-Vector Multiplication
        Vector<T> operator*(const Vector<T> &vector) const
        {
            if constexpr (ORD == ORDERING::ColMajor)
            {
                if (cols_ != vector.Size())
                {
                    // Invalid multiplication, return an empty vector or throw an exception
                    return Vector<T>(0);
                }

                Vector<T> result(rows_);
                for (size_t i = 0; i < rows_; ++i)
                {
                    T sum = 0;
                    for (size_t j = 0; j < cols_; ++j)
                    {
                        sum += (*this)(i, j) * vector(j);
                    }
                    result(i) = sum;
                }
                return result;
            }
            else
            {
                if (rows_ != vector.Size())
                {
                    // Invalid multiplication, return an empty vector or throw an exception
                    return Vector<T>(0);
                }

                Vector<T> result(cols_);
                for (size_t j = 0; j < cols_; ++j)
                {
                    T sum = 0;
                    for (size_t i = 0; i < rows_; ++i)
                    {
                        sum += (*this)(i, j) * vector(i);
                    }
                    result(j) = sum;
                }
                return result;
            }
        }

        Matrix operator-(const Matrix &other) const
        {
            if (rows_ != other.rows_ || cols_ != other.cols_)
            {
                // Invalid addition, return an empty matrix or throw an exception
                throw std::invalid_argument("Invalid addition");
            }

            Matrix result(rows_, cols_);

            for (size_t i = 0; i < rows_; ++i)
            {
                for (size_t j = 0; j < cols_; ++j)
                {
                    if constexpr (ORD == ORDERING::ColMajor)
                    {
                        result(i, j) = (*this)(i, j) - other(i, j);
                    }
                    else
                    {
                        result(i, j) = (*this)(i, j) - other(i, j);
                    }
                }
            }
            return result;
        }

        Matrix operator/(const T &scal) const
        {
            Matrix result(rows_, cols_);

            for (size_t i = 0; i < rows_; ++i)
            {
                for (size_t j = 0; j < cols_; ++j)
                {
                    if constexpr (ORD == ORDERING::ColMajor)
                    {
                        result(i, j) = (*this)(i, j) / scal;
                    }
                    else
                    {
                        result(i, j) = (*this)(i, j) / scal;
                    }
                }
            }
            return result;
        }

        Matrix operator-() const
        {
            Matrix result(rows_, cols_);

            for (size_t i = 0; i < rows_; ++i)
            {
                for (size_t j = 0; j < cols_; ++j)
                {
                    if constexpr (ORD == ORDERING::ColMajor)
                    {
                        result(i, j) = -(*this)(i, j);
                    }
                    else
                    {
                        result(i, j) = -(*this)(i, j);
                    }
                }
            }
            return result;
        }

        Matrix operator+=(const Matrix &other)
        {
            if (rows_ != other.rows_ || cols_ != other.cols_)
            {
                // Invalid addition, return an empty matrix or throw an exception
                throw std::invalid_argument("Invalid addition");
            }

            for (size_t i = 0; i < rows_; ++i)
            {
                for (size_t j = 0; j < cols_; ++j)
                {
                    if constexpr (ORD == ORDERING::ColMajor)
                    {
                        (*this)(i, j) += other(i, j);
                    }
                    else
                    {
                        (*this)(i, j) += other(i, j);
                    }
                }
            }
            return *this;
        }

        Matrix operator-=(const Matrix &other)
        {
            if (rows_ != other.rows_ || cols_ != other.cols_)
            {
                // Invalid addition, return an empty matrix or throw an exception
                throw std::invalid_argument("Invalid addition");
            }

            for (size_t i = 0; i < rows_; ++i)
            {
                for (size_t j = 0; j < cols_; ++j)
                {
                    if constexpr (ORD == ORDERING::ColMajor)
                    {
                        (*this)(i, j) -= other(i, j);
                    }
                    else
                    {
                        (*this)(i, j) -= other(i, j);
                    }
                }
            }
            return *this;
        }

        Matrix operator*=(const T &scal)
        {
            for (size_t i = 0; i < rows_; ++i)
            {
                for (size_t j = 0; j < cols_; ++j)
                {
                    if constexpr (ORD == ORDERING::ColMajor)
                    {
                        (*this)(i, j) *= scal;
                    }
                    else
                    {
                        (*this)(i, j) *= scal;
                    }
                }
            }
            return *this;
        }

        Matrix operator/=(const T &scal)
        {
            for (size_t i = 0; i < rows_; ++i)
            {
                for (size_t j = 0; j < cols_; ++j)
                {
                    if constexpr (ORD == ORDERING::ColMajor)
                    {
                        (*this)(i, j) /= scal;
                    }
                    else
                    {
                        (*this)(i, j) /= scal;
                    }
                }
            }
            return *this;
        }

        Matrix operator*=(const Matrix &other)
        {
            if (Size() != other.Size())
            {
                // Invalid multiplication, return an empty matrix or throw an exception
                throw std::invalid_argument("Invalid multiplication");
            }

            Matrix result(rows_, other.cols_);

            for (size_t i = 0; i < rows_; ++i)
            {
                for (size_t j = 0; j < other.cols_; ++j)
                {
                    T sum = 0;
                    for (size_t k = 0; k < cols_; ++k)
                    {
                        if constexpr (ORD == ORDERING::ColMajor)
                        {
                            sum += (*this)(i, k) * other(k, j);
                        }
                        else
                        {
                            sum += (*this)(i, k) * other(k, j);
                        }
                    }
                    result(i, j) = sum;
                }
            }
            *this = result;
            return *this;
        }

        T Determinant() const
        {
            if (rows_ != cols_)
            {
                // Invalid determinant, return an empty matrix or throw an exception
                throw std::invalid_argument("Invalid determinant");
            }

            if (rows_ == 1)
            {
                return (*this)(0, 0);
            }
            else if (rows_ == 2)
            {
                return (*this)(0, 0) * (*this)(1, 1) - (*this)(0, 1) * (*this)(1, 0);
            }
            else
            {
                T det = 0;
                for (size_t i = 0; i < rows_; ++i)
                {
                    Matrix<T, ORD> sub_matrix(rows_ - 1, cols_ - 1);
                    for (size_t j = 1; j < rows_; ++j)
                    {
                        for (size_t k = 0; k < cols_; ++k)
                        {
                            if (k < i)
                            {
                                sub_matrix(j - 1, k) = (*this)(j, k);
                            }
                            else if (k > i)
                            {
                                sub_matrix(j - 1, k - 1) = (*this)(j, k);
                            }
                        }
                    }
                    if constexpr (ORD == ORDERING::ColMajor)
                    {
                        det += (*this)(0, i) * sub_matrix.Determinant();
                    }
                    else
                    {
                        det += (*this)(0, i) * sub_matrix.Determinant();
                    }
                }
                return det;
            }
        }

        // Obsolete
        [[deprecated("This does not work as intended")]]
        Matrix Inverse() const
        {
            if (rows_ != cols_)
            {
                // Invalid determinant, return an empty matrix or throw an exception
                throw std::invalid_argument("Invalid determinant");
            }

            Matrix<T, ORD> result(rows_, cols_);

            if (rows_ == 1)
            {
                result(0, 0) = 1 / (*this)(0, 0);
            }
            else if (rows_ == 2)
            {
                T det = (*this)(0, 0) * (*this)(1, 1) - (*this)(0, 1) * (*this)(1, 0);
                result(0, 0) = (*this)(1, 1) / det;
                result(0, 1) = -(*this)(0, 1) / det;
                result(1, 0) = -(*this)(1, 0) / det;
                result(1, 1) = (*this)(0, 0) / det;
            }
            else
            {
                T det = Determinant();
                for (size_t i = 0; i < rows_; ++i)
                {
                    for (size_t j = 0; j < cols_; ++j)
                    {
                        Matrix<T, ORD> sub_matrix(rows_ - 1, cols_ - 1);
                        for (size_t k = 0; k < rows_; ++k)
                        {
                            for (size_t l = 0; l < cols_; ++l)
                            {
                                if (k < i && l < j)
                                {
                                    sub_matrix(k, l) = (*this)(k, l);
                                }
                                else if (k < i && l > j)
                                {
                                    sub_matrix(k, l - 1) = (*this)(k, l);
                                }
                                else if (k > i && l < j)
                                {
                                    sub_matrix(k - 1, l) = (*this)(k, l);
                                }
                                else if (k > i && l > j)
                                {
                                    sub_matrix(k - 1, l - 1) = (*this)(k, l);
                                }
                            }
                        }
                        if constexpr (ORD == ORDERING::ColMajor)
                        {
                            result(j, i) = sub_matrix.Determinant() / det;
                        }
                        else
                        {
                            result(i, j) = sub_matrix.Determinant() / det;
                        }
                    }
                }
            }
            return result;
        }
    };

    template <typename T, ORDERING ORD>
    Vector<T> operator*(const Vector<T> &vector, const Matrix<T, ORD> &matrix)
    {
        if constexpr (ORD == ORDERING::ColMajor)
        {
            if (matrix.Size_Rows() != vector.Size())
            {
                // Invalid multiplication, return an empty vector or throw an exception
                return Vector<T>(0);
            }

            Vector<T> result(matrix.Size_Cols());
            for (size_t j = 0; j < matrix.Size_Cols(); ++j)
            {
                T sum = 0;
                for (size_t i = 0; i < matrix.Size_Rows(); ++i)
                {
                    sum += matrix(i, j) * vector(i);
                }
                result(j) = sum;
            }
            return result;
        }
        else
        {
            if (matrix.Size_Cols() != vector.Size())
            {
                // Invalid multiplication, return an empty vector or throw an exception
                return Vector<T>(0);
            }

            Vector<T> result(matrix.Size_Rows());
            for (size_t i = 0; i < matrix.Size_Rows(); ++i)
            {
                T sum = 0;
                for (size_t j = 0; j < matrix.Size_Cols(); ++j)
                {
                    sum += matrix(i, j) * vector(j);
                }
                result(i) = sum;
            }
            return result;
        }
    }
}
#endif