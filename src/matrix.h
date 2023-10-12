#ifndef FILE_MATRIX_H
#define FILE_MATRIX_H

#include <iostream>

#include "vector.h" // Assuming your Vector class is in "vector.h"

namespace ASC_bla {
    enum class ORDERING { ColMajor, RowMajor };

    template <typename T, ORDERING ORD>
    class Matrix {
        size_t rows_;
        size_t cols_;
        Vector<T> data_;

    public:
        // Constructors
        Matrix(size_t rows, size_t cols) : rows_(rows), cols_(cols), data_(rows * cols) {}

        // Copy constructor
        Matrix(const Matrix& other) : rows_(other.rows_), cols_(other.cols_), data_(other.data_) {}

        // Move constructor
        Matrix(Matrix&& other) noexcept : rows_(other.rows_), cols_(other.cols_), data_(std::move(other.data_)) {}

        // Destructor
        ~Matrix() = default;

        // Access operator for element (i, j)
T& operator()(size_t i, size_t j) {
    if constexpr (ORD == ORDERING::ColMajor) {
        return data_(i + j * rows_);
    } else {
        return data_(j + i * cols_);
    }
}

// Const version of the access operator
const T& operator()(size_t i, size_t j) const {
    if constexpr (ORD == ORDERING::ColMajor) {
        return data_(i + j * rows_);
    } else {
        return data_(j + i * cols_);
    }
}


        // Assignment operator
        Matrix& operator=(const Matrix& other) {
            if (this != &other) {
                rows_ = other.rows_;
                cols_ = other.cols_;
                data_ = other.data_;
            }
            return *this;
        }

        // Output stream operator for easy printing
        friend std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
            for (size_t i = 0; i < matrix.rows_; ++i) {
                if (i > 0) {
                    os << "\n";
                }
                for (size_t j = 0; j < matrix.cols_; ++j) {
                    if (j > 0) {
                        os << " ";
                    }
                    if constexpr (ORD == ORDERING::ColMajor) {
                        os << matrix.data_(i + j * matrix.rows_);
                    } else {
                        os << matrix.data_(j + i * matrix.cols_);
                    }
                }
            }
            return os;
        }
        // Matrix-Matrix Multiplication
Matrix operator*(const Matrix& other) const {
    if (cols_ != other.rows_) {
        // Invalid multiplication, return an empty matrix or throw an exception
        return Matrix(0, 0);
    }

    Matrix result(rows_, other.cols_);

    for (size_t i = 0; i < rows_; ++i) {
        for (size_t j = 0; j < other.cols_; ++j) {
            T sum = 0;
            for (size_t k = 0; k < cols_; ++k) {
                if constexpr (ORD == ORDERING::ColMajor) {
                    sum += (*this)(i, k) * other(k, j);
                } else {
                    sum += (*this)(i, k) * other(k, j);
                }
            }
            result(i, j) = sum;
        }
    }
    return result;
}
// Transpose method
        Matrix<T, ORD> transpose() const {
            Matrix<T, ORD> result(cols_, rows_);
            for (size_t i = 0; i < rows_; ++i) {
                for (size_t j = 0; j < cols_; ++j) {
                    result(j, i) = (*this)(i, j);
                }
            }
            return result;
        }
        // Matrix-Vector Multiplication
        Vector<T> operator*(const Vector<T>& vector) const {
            if constexpr (ORD == ORDERING::ColMajor) {
                if (cols_ != vector.Size()) {
                    // Invalid multiplication, return an empty vector or throw an exception
                    return Vector<T>(0);
                }

                Vector<T> result(rows_);
                for (size_t i = 0; i < rows_; ++i) {
                    T sum = 0;
                    for (size_t j = 0; j < cols_; ++j) {
                        sum += (*this)(i, j) * vector(j);
                    }
                    result(i) = sum;
                }
                return result;
            } else {
                if (rows_ != vector.Size()) {
                    // Invalid multiplication, return an empty vector or throw an exception
                    return Vector<T>(0);
                }

                Vector<T> result(cols_);
                for (size_t j = 0; j < cols_; ++j) {
                    T sum = 0;
                    for (size_t i = 0; i < rows_; ++i) {
                        sum += (*this)(i, j) * vector(i);
                    }
                    result(j) = sum;
                }
                return result;
            }
        }
    };
}
#endif