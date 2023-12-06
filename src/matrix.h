#ifndef FILE_MATRIX_H
#define FILE_MATRIX_H

#include <iostream>

#include "expression.h"

#include "vector.h"
#include "ordering.h"

#include "simd.h"

#include <optional>

using namespace ASC_HPC;

namespace ASC_bla
{

    template <typename T = double, ORDERING ORD = ORDERING::ColMajor>
    class MatrixView : public MatExpr<MatrixView<T, ORD>, ORD>
    {
    protected:
        size_t rows_;
        size_t cols_;
        size_t dist_;
        T *data_;

        auto MultiplyEight(const MatrixView<T, ORD> &other) const
        {
            MatrixView<T, ORD> result(rows_, other.cols_);
            
            for (size_t i = 0; i < rows_; i+=8)
            {
                auto row = Row(i);
                auto col_0 = other.Col(0);
                auto col_1 = other.Col(1);
                auto col_2 = other.Col(2);
                auto col_3 = other.Col(3);
                auto col_4 = other.Col(4);
                auto col_5 = other.Col(5);
                auto col_6 = other.Col(6);
                auto col_7 = other.Col(7);

                result(i, 0) = InnerProduct<8>(8, row, 1, col_0, 1);
                result(i, 1) = InnerProduct<8>(8, row, 1, col_1, 1);
                result(i, 2) = InnerProduct<8>(8, row, 1, col_2, 1);
                result(i, 3) = InnerProduct<8>(8, row, 1, col_3, 1);
                result(i, 4) = InnerProduct<8>(8, row, 1, col_4, 1);
                result(i, 5) = InnerProduct<8>(8, row, 1, col_5, 1);
                result(i, 6) = InnerProduct<8>(8, row, 1, col_6, 1);
                result(i, 7) = InnerProduct<8>(8, row, 1, col_7, 1);
            }
            return result;
        }

    public:
        MatrixView(size_t rows, size_t cols, T *data)
            : rows_(rows), cols_(cols), data_(data), dist_((ORD == ORDERING::ColMajor) ? rows : cols)
        {
        }

        MatrixView(size_t rows, size_t cols)
            : rows_(rows), cols_(cols), data_(new T[rows * cols]), dist_((ORD == ORDERING::ColMajor) ? rows : cols)
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

        template <typename TB, ORDERING ORDB>
        void CopyToNew(const MatrixView<TB, ORDB> &m2)
        {
            rows_ = m2.Size_Cols();
            cols_ = m2.Size_Rows();
            data_ = new T[rows_ * cols_];

            for (size_t i = 0; i < rows_; ++i)
            {
                for (size_t j = 0; j < cols_; ++j)
                {
                    (*this)(i, j) = m2(i, j);
                }
            }
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
                return data_[i + j * rows_];
            }
            else
            {
                return data_[j + i * cols_];
            }
        }
        const T &operator()(size_t i, size_t j) const
        {
            if constexpr (ORD == ORDERING::ColMajor)
            {
                return data_[i + j * rows_];
            }
            else
            {
                return data_[j + i * cols_];
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
                return VectorView(cols_, rows_, data_ + i);
            }
            else
            {
                return VectorView(cols_, (size_t)1, data_ + i  * cols_);
            }
        }

        auto Col(size_t i) const
        {
            if constexpr (ORD == ORDERING::ColMajor)
            {
                return VectorView(rows_, (size_t)1, data_ + i * rows_);
            }
            else
            {
                return VectorView(rows_, cols_, data_ + i);
            }
        }

        auto Flatten() const
        {
            return VectorView(rows_ * cols_, (size_t)1, data_);
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

        auto SubMatrix(size_t first_row, size_t row_length, size_t first_col, size_t col_length) const
        {
            MatrixView<T, ORD> result(row_length, col_length);
            for (size_t i = 0; i < row_length; ++i)
            {
                for (size_t j = 0; j < col_length; ++j)
                {
                    result(i, j) = (*this)(first_row + i, first_col + j);
                }
            }

            return result;
        }

        auto Diag() const
        {
            return VectorView<T, size_t>(std::min(rows_, cols_), dist_ + 1, data_);
        }

        auto Size() const
        {
            return rows_ * cols_;
        }

        auto operator*=(const T &scal)
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

        auto operator/=(const T &scal)
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

        // Matrix-Matrix Multiplication, using SIMD (InnerProduct), if T is double
        MatrixView<T, ORD> operator*(const MatrixView &other) const
        {
            if (Size() != other.Size())
            {
                // Invalid multiplication, return an empty matrix or throw an exception
                throw std::invalid_argument("Invalid multiplication");
            }

            MatrixView result(rows_, other.cols_);

            if constexpr (std::is_same<double, T>::value)
            {
                MatrixView result_to_be_merged(rows_, other.cols_);
                
                std::optional<MatrixView<double, ORDERING::ColMajor>> sub_result_list = std::nullopt;

                // Check if dimension >= 96, if yes, use caching
                if (rows_ >= 96 && other.cols_ >= 96)
                {
                    sub_result_list = MatrixView<double, ORDERING::ColMajor>(96 * 96, size_t(rows_ / 96) * size_t(other.cols_ / 96));
                    //std::cout << "sub_result_list init\n" << sub_result_list.value() << "\n";
                }

                for (size_t i = 0; i < rows_; i+=96)
                {
                    for (size_t j = 0; j < other.cols_; j+=96)
                    {
                        // Check if sub matrix can be created
                        // +1 to avoid overflow
                        if (i + 96 + 1 > rows_  || j + 96 + 1 > other.cols_)
                        {
                            // If not, use the normal multiplication

                            // std::cout << "Normal\n";
                            for (size_t l = 0; l + j < other.cols_; ++l)
                            {
                                for (size_t k = 0; k + i < rows_; ++k)
                                {
                                    if(i + k + 7 < rows_)
                                    {
                                        auto tmp = InnerProduct8<8>(cols_, Row(i + k), Row(i + k + 1), Row(i + k + 2), Row(i + k + 3), Row(i + k + 4), Row(i + k + 5), Row(i + k + 6), Row(i + k + 7), other.Col(j + l));
                                        result_to_be_merged(i + k, j + l) = std::get<0>(tmp);
                                        k++;
                                        result_to_be_merged(i + k, j + l) = std::get<1>(tmp);
                                        k++;
                                        result_to_be_merged(i + k, j + l) = std::get<2>(tmp);
                                        k++;
                                        result_to_be_merged(i + k, j + l) = std::get<3>(tmp);
                                        k++;
                                        result_to_be_merged(i + k, j + l) = std::get<4>(tmp);
                                        k++;
                                        result_to_be_merged(i + k, j + l) = std::get<5>(tmp);
                                        k++;
                                        result_to_be_merged(i + k, j + l) = std::get<6>(tmp);
                                        k++;
                                        result_to_be_merged(i + k, j + l) = std::get<7>(tmp);
                                    }
                                    else
                                    {
                                        result_to_be_merged(i + k, j + l) = InnerProduct<8>(cols_, Row(i + k), other.Col(j + l));
                                    }
                                }
                            }
                        }
                        else
                        {
                            // If yes, use caching

                            // std::cout << "Caching\n";

                            auto sub_matrix = SubMatrix(i, 96, j, 96);
                            auto sub_matrix_other = other.SubMatrix(i, 96, j, 96);

                            // std::cout << "Sub matrix\n" << sub_matrix << "\n";
                            // std::cout << "Sub matrix other\n" << sub_matrix_other << "\n";

                            auto sub_matrix_result = sub_matrix * sub_matrix_other;

                            // std::cout << "Sub matrix result\n" << sub_matrix_result << "\n";

                            auto flat = sub_matrix_result.Flatten();
                            auto col_var = sub_result_list.value().Col(size_t((i + j) / 96));

                            for (size_t k = 0; k < flat.Size(); ++k)
                            {
                                col_var(k) = flat(k);
                            }

                            // std::cout << "Sub matrix result flatten\n" << sub_matrix_result.Flatten() << "\n";

                            // std::cout << "Sub result list\n" << sub_result_list.value().Col(size_t((i + j) / 96)) << "\n";
                        }
                    }
                }
                
                // std::cout << "Rows: " << rows_ << "\n";
                // std::cout << "Cols: " << other.cols_ << "\n";

                // Check if dimension < 96, if yes, copy results from result_to_be_merged to result
                if (rows_ <= 96 || other.cols_ <= 96)
                {
                    // std::cout << "result_to_be_merged\n" << result_to_be_merged << "\n";

                    return result_to_be_merged;
                }

                // std::cout << "sub_result_list\n" << sub_result_list.value() << "\n";

                // Merge the result from sub_result_list to result by using the sub_result_list in the following way:
                // result(i, j) = sum(sub_result_list.Row(i * j)(cols//96 * k)) where k is an index from 0 to cols//96 * rows//96

                for (size_t i = 0; i < rows_; ++i)
                {
                    for (size_t j = 0; j < other.cols_; ++j)
                    {
                        // std::cout << "i: " << i << "\n";

                        auto row_var = sub_result_list.value().Row(((i + j)/96) * size_t(rows_ / 96) * size_t(other.cols_ / 96));

                        // std::cout << "row_var\n" << row_var << "\n";

                        auto sum = 0.0;

                        for (size_t k = 0; k < size_t(cols_ / 96) * size_t(rows_ / 96); ++k)
                        {
                            sum += row_var(k + j/96);
                        }

                        result(i, j) = sum;
                    }
                }

                return result;
            }
            else
            {

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
            }
            return result;
        }

        auto operator*(const T &scal) const
        {
            MatrixView result(rows_, cols_);

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
        auto operator+(const MatrixView &other) const
        {
            if (rows_ != other.rows_ || cols_ != other.cols_)
            {
                // Invalid addition, return an empty matrix or throw an exception
                throw std::invalid_argument("Invalid addition");
            }

            MatrixView result(rows_, cols_);

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
        auto transpose() const
        {
            MatrixView result(cols_, rows_);

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

        auto operator-(const MatrixView &other) const
        {
            if (rows_ != other.rows_ || cols_ != other.cols_)
            {
                // Invalid addition, return an empty matrix or throw an exception
                throw std::invalid_argument("Invalid addition");
            }

            MatrixView result(rows_, cols_);

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

        auto operator/(const T &scal) const
        {
            MatrixView result(rows_, cols_);

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

        auto operator-() const
        {
            MatrixView result(rows_, cols_);

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

        auto operator+=(const MatrixView &other)
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

        auto operator-=(const MatrixView &other)
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

        auto operator*=(const MatrixView &other)
        {
            if (Size() != other.Size())
            {
                // Invalid multiplication, return an empty matrix or throw an exception
                throw std::invalid_argument("Invalid multiplication");
            }

            MatrixView result(rows_, other.cols_);

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

        friend std::ostream &operator<<(std::ostream &os, const MatrixView &matrix)
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
                        os << matrix(i, j);
                    }
                    else
                    {
                        os << matrix(i, j);
                    }
                }
            }
            return os;
        }
    };

    template <typename T = double, ORDERING ORD = ORDERING::ColMajor>
    class Matrix : public MatrixView<T, ORD>
    {

        using MatrixView<T, ORD>::rows_;
        using MatrixView<T, ORD>::cols_;
        using MatrixView<T, ORD>::data_;
        using MatrixView<T, ORD>::dist_;

    public:
        using MatrixView<T, ORD>::operator();
        using MatrixView<T, ORD>::Size;
        using MatrixView<T, ORD>::Row;
        using MatrixView<T, ORD>::Col;

    public:
        // Constructors
        Matrix(size_t rows, size_t cols) : MatrixView<T, ORD>(rows, cols, new T[rows * cols]) {}
        Matrix(const MatrixView<T, ORD> &other) : MatrixView<T, ORD>(other.Size_Rows(), other.Size_Cols(), other.Data()) {}

        // Copy constructor
        Matrix(const Matrix &other) : Matrix(other.rows_, other.cols_)
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