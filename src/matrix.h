#ifndef FILE_MATRIX_H
#define FILE_MATRIX_H

#include <iostream>

#include "expression.h"
#include "vector.h"
#include "ordering.h"

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
        MatrixView(size_t rows, size_t cols, T *data, size_t dist)
            : rows_(rows), cols_(cols), data_(data), dist_(dist) {}

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

        auto Row(size_t i)
        {
            if constexpr (ORD == ORDERING::ColMajor)
            {
                auto data_tmp = data_ + i * cols_;

                return VectorView(cols_, 1, data_tmp);
            }
            else
            {
                return Transpose(this).Row(i);
            }
        }

        auto Col(size_t i)
        {
            if constexpr (ORD == ORDERING::ColMajor)
            {
                return this.Row(i);
            }
            else
            {
                return Transpose(this).Col(i);
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
    };

    template <typename T, ORDERING ORD>
    class Matrix : protected MatrixView<T, ORD>
    {

    public:
        // Constructors
        Matrix(size_t rows, size_t cols) : MatrixView<T, ORD>(rows, cols, new T[rows*cols], if constexpr (ORD == ORDERING::ColMajor) {rows} else {cols}) {}

        // Copy constructor
        Matrix(const Matrix &other) : Matrix (other.rows_, other.cols_)
        {
            for (size_t i = 0; i < rows_; ++i)
            {
                for (size_t j = 0; j < Size_Cols(); ++j)
                {
                    (*this)(i, j) = other(i, j);
                }
            }
        }

        // Move constructor
        Matrix(Matrix &&other) noexcept : MatrixView<T, ORD>(other.Size_Rows(), other.Size_Cols(), other.Data(), other.Dist()) {}

        // Destructor
        ~Matrix() = default;

        // Access operator for element (i, j)
        T &operator()(size_t i, size_t j)
        {
            if constexpr (ORD == ORDERING::ColMajor)
            {
                return Data()[i + j * rows_];
            }
            else
            {
                return Data()[j + i * cols_];
            }
        }

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
            if (Size_Cols() != other.Size_Rows())
            {
                // Invalid multiplication, return an empty matrix or throw an exception
                return nullptr;
            }

            Matrix result(Size_Rows(), other.Size_Cols());

            for (size_t i = 0; i < rows_; ++i)
            {
                for (size_t j = 0; j < other.Size_Cols(); ++j)
                {
                    T sum = 0;
                    for (size_t k = 0; k < Size_Cols(); ++k)
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
            return result;
        }

        Matrix operator*(const T &scal) const
        {
            Matrix result(Size_Rows(), Size_Cols());

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
            if (Size_Rows != other.Size_Rows() || Size_Cols() != other.Size_Cols())
            {
                // Invalid addition, return an empty matrix or throw an exception
                return nullptr;
            }

            Matrix result(Size_Rows(), Size_Cols());

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
            Matrix<T, ORD> result(Size_Cols(), Size_Rows());

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

                Vector<T> result(Size_Rows());
                for (size_t i = 0; i < Size_Rows(); ++i)
                {
                    T sum = 0;
                    for (size_t j = 0; j < Size_Cols(); ++j)
                    {
                        sum += (*this)(i, j) * vector(j);
                    }
                    result(i) = sum;
                }
                return result;
            }
            else
            {
                if (Size_Rows() != vector.Size())
                {
                    // Invalid multiplication, return an empty vector or throw an exception
                    return Vector<T>(0);
                }

                Vector<T> result(Size_Cols());
                for (size_t j = 0; j < Size_Cols(); ++j)
                {
                    T sum = 0;
                    for (size_t i = 0; i < Size_Rows(); ++i)
                    {
                        sum += (*this)(i, j) * vector(i);
                    }
                    result(j) = sum;
                }
                return result;
            }
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