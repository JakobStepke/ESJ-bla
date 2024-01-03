#ifndef FILE_VECTOR_H
#define FILE_VECTOR_H

#include <iostream>
#include <cmath>

#include "forward_decl.h"

#include "expression.h"
#include "simd.h"

using namespace ASC_HPC;

namespace ASC_bla
{

  template <typename T, typename TDIST>
  class VectorView : public VecExpr<VectorView<T, TDIST>>
  {
  protected:
    T *data_;
    size_t size_;
    TDIST dist_;

  public:
    VectorView(size_t size, T *data)
        : data_(data), size_(size) {}

    VectorView(size_t size, TDIST dist, T *data)
        : data_(data), size_(size), dist_(dist) {}

    VectorView(const VectorView &v)
        : data_(v.data_), size_(v.size_), dist_(v.dist_) {}

    template <typename TB>
    VectorView &operator=(const VecExpr<TB> &v2)
    {
      for (size_t i = 0; i < size_; i++)
        data_[dist_ * i] = v2(i);
      return *this;
    }

    VectorView &operator=(const VectorView &v2)
    {
      for (size_t i = 0; i < size_; i++)
        data_[dist_ * i] = v2(i);
      return *this;
    }

    VectorView &operator=(T scal)
    {
      for (size_t i = 0; i < size_; i++)
        data_[dist_ * i] = scal;
      return *this;
    }

    VectorView &operator+=(T scal)
    {
      for (size_t i = 0; i < size_; i++)
        data_[dist_ * i] += scal;
      return *this;
    }

    VectorView &operator-=(T scal)
    {
      for (size_t i = 0; i < size_; i++)
        data_[dist_ * i] -= scal;
      return *this;
    }

    template <typename TB>
    VectorView &operator+=(const VecExpr<TB> &v2)
    {
      for (size_t i = 0; i < size_; i++)
        data_[dist_ * i] += v2(i);
      return *this;
    }

    template <typename TB>
    VectorView &operator-=(const VecExpr<TB> &v2)
    {
      for (size_t i = 0; i < size_; i++)
        data_[dist_ * i] -= v2(i);
      return *this;
    }

    VectorView &operator*=(T scal)
    {
      for (size_t i = 0; i < size_; i++)
        data_[dist_ * i] *= scal;
      return *this;
    }

    VectorView &operator/=(T scal)
    {
      for (size_t i = 0; i < size_; i++)
        data_[dist_ * i] /= scal;
      return *this;
    }

    auto View() const { return VectorView(size_, dist_, data_); }
    size_t Size() const { return size_; }
    auto Dist() const { return dist_; }
    T &operator()(size_t i) { return data_[dist_ * i]; }
    const T &operator()(size_t i) const { return data_[dist_ * i]; }

    auto Data() const { return data_; }

    auto Range(size_t first, size_t next) const
    {
      return VectorView(next - first, dist_, data_ + first * dist_);
    }

    auto Slice(size_t first, size_t slice) const
    {
      return VectorView<T, size_t>(size_ / slice, dist_ * slice, data_ + first * dist_);
    }

    auto Norm2() const
    {
      double sum = 0;
      for (size_t i = 0; i < size_; i++)
        sum += pow(data_[dist_ * i], 2);
      return sqrt(sum);
    }

    auto AsMatrix(size_t rows, size_t cols) const
    {
      static_assert(std::is_same<TDIST, std::integral_constant<size_t, 1>>::value == true, "gapped vectors cannot be converted to matrices");
      
      
      auto mat = MatrixView<T>(rows, cols, data_);
      // std::cout << "mat = " << mat << std::endl;
      return mat;
    }
  };

  template <typename T = double>
  class Vector : public VectorView<T>
  {
    typedef VectorView<T> BASE;
    using BASE::data_;
    using BASE::size_;

  public:
    Vector(size_t size)
        : VectorView<T>(size, new T[size]) { ; }
    // Vector from list

    Vector(std::initializer_list<T> list)
        : VectorView<T>(list.size(), new T[list.size()])
    {
      size_t i = 0;
      for (auto it = list.begin(); it != list.end(); ++it)
        data_[i++] = *it;
    }

    Vector(const Vector &v)
        : Vector(v.Size())
    {
      *this = v;
    }

    Vector(Vector &&v)
        : VectorView<T>(0, nullptr)
    {
      std::swap(size_, v.size_);
      std::swap(data_, v.data_);
    }

    template <typename TB>
    Vector(const VecExpr<TB> &v)
        : Vector(v.Size())
    {
      *this = v;
    }

    ~Vector() { delete[] data_; }

    using BASE::operator=;
    Vector &operator=(const Vector &v2)
    {
      for (size_t i = 0; i < size_; i++)
        data_[i] = v2(i);
      return *this;
    }

    Vector &operator=(Vector &&v2)
    {
      for (size_t i = 0; i < size_; i++)
        data_[i] = v2(i);
      return *this;
    }

    Vector &operator=(T scal)
    {
      for (size_t i = 0; i < size_; i++)
        data_[i] = scal;
      return *this;
    }

    Vector &operator*(Vector &&v2)
    {
      Vector ret = Vector(size_);
      for (size_t i = 0; i < size_; i++)
        ret(i) = this(i) * v2(i);
      return ret;
    }
  };

  template <typename... Args>
  std::ostream &operator<<(std::ostream &ost, const VectorView<Args...> &v)
  {
    if (v.Size() > 0)
      ost << v(0);
    for (size_t i = 1; i < v.Size(); i++)
      ost << ", " << v(i);
    return ost;
  }

  template <size_t SW>
  auto InnerProduct(size_t n, const VectorView<double, size_t> x,
                    const VectorView<double, size_t> y)
  {
    SIMD<double, SW> sum{0.0};
    for (size_t i = 0; i < n; i++)
    {
      // sum += px[i] * SIMD<double,SW>(py+i*dy);
      sum = FMA(SIMD<double, SW>(x(i)), SIMD<double, SW>(y(i)), sum);
    }

    /*
    std::cout << "sum = " << sum << std::endl;
    std::cout << "x = " << x << std::endl;
    std::cout << "y = " << y << std::endl;
    */
    double first_val = sum.GetFirst();

    return first_val;
  }

  template <>
  auto InnerProduct<1>(size_t n, const VectorView<double, size_t> x,
                       const VectorView<double, size_t> y)
  {
    double sum = 0.0;
    for (size_t i = 0; i < n; i++)
    {
      // std::cout << "x(i) = " << x(i) << std::endl;
      // std::cout << "y(i) = " << y(i) << std::endl;
      // std::cout << "x(i)*y(i) = " << x(i) * y(i) << std::endl;
      // std::cout << "sum = " << sum << std::endl;
      sum += x(i) * y(i);
    }
    return sum;
  }

  template <size_t SW>
  auto InnerProduct8(size_t n, const VectorView<double, size_t> x0,
                     const VectorView<double, size_t> x1,
                     const VectorView<double, size_t> x2,
                     const VectorView<double, size_t> x3,
                     const VectorView<double, size_t> x4,
                     const VectorView<double, size_t> x5,
                     const VectorView<double, size_t> x6,
                     const VectorView<double, size_t> x7,
                     const VectorView<double, size_t> y)
  {
    SIMD<double, SW> sum0{0.0};
    SIMD<double, SW> sum1{0.0};
    SIMD<double, SW> sum2{0.0};
    SIMD<double, SW> sum3{0.0};
    SIMD<double, SW> sum4{0.0};
    SIMD<double, SW> sum5{0.0};
    SIMD<double, SW> sum6{0.0};
    SIMD<double, SW> sum7{0.0};

    for (size_t i = 0; i < n; i++)
    {
      // sum += px[i] * SIMD<double,SW>(py+i*dy);
      sum0 = FMA(SIMD<double, SW>(x0(i)), SIMD<double, SW>(y(i)), sum0);
      sum1 = FMA(SIMD<double, SW>(x1(i)), SIMD<double, SW>(y(i)), sum1);
      sum2 = FMA(SIMD<double, SW>(x2(i)), SIMD<double, SW>(y(i)), sum2);
      sum3 = FMA(SIMD<double, SW>(x3(i)), SIMD<double, SW>(y(i)), sum3);
      sum4 = FMA(SIMD<double, SW>(x4(i)), SIMD<double, SW>(y(i)), sum4);
      sum5 = FMA(SIMD<double, SW>(x5(i)), SIMD<double, SW>(y(i)), sum5);
      sum6 = FMA(SIMD<double, SW>(x6(i)), SIMD<double, SW>(y(i)), sum6);
      sum7 = FMA(SIMD<double, SW>(x7(i)), SIMD<double, SW>(y(i)), sum7);
    }
    return std::tuple(sum0.GetFirst(), sum1.GetFirst(), sum2.GetFirst(), sum3.GetFirst(),
                      sum4.GetFirst(), sum5.GetFirst(), sum6.GetFirst(), sum7.GetFirst());
  }

  template <size_t S, typename T = double>
  class Vec : public VecExpr<Vec<S, T>>
  {
    T data[S];

  public:
    Vec() {}

    Vec(const Vec<S, T> &v2)
    {
      for (int i = 0; i < S; i++)
      {
        data[i] = v2(i);
      }
    }

    Vec(T all)
    {
      for (size_t i = 0; i < S; i++)
      {
        data[i] = all;
      }
    };

    // initializer list constructor
    Vec(std::initializer_list<T> list)
    {
      if (list.size() != S)
        throw std::invalid_argument("initializer list with wrong size");
      // copy list
      for (size_t i = 0; i < S; i++)
      {
        data[i] = list.begin()[i];
      }
    }

    // constructor from VectorView
    template <typename T2, typename TDIST>
    Vec(VectorView<T2, TDIST> v2)
    {
      if (v2.Size() != S)
        throw std::invalid_argument("VectorView has wrong size");
      // copy list
      for (size_t i = 0; i < S; i++)
      {
        data[i] = v2(i);
      }
    }

    // constructor from VecExpr
    template <typename T2>
    Vec(const VecExpr<T2> &v2)
    {
      if (v2.Size() != S)
        throw std::invalid_argument("VecExpr has wrong size");

      for (int i = 0; i < S; i++)
      {
        data[i] = v2(i);
      }
    }

    // VectorView operator=
    template <typename T2, typename TDIST>
    Vec &operator=(VectorView<T2, TDIST> v2)
    {
      if (v2.Size() != S)
        throw std::invalid_argument("VectorView has wrong size");

        // std::cout << "Using this operator=" << std::endl;
      // copy list
      for (size_t i = 0; i < S; i++)
      {
        // std::cout << "v2(i) = " << v2(i) << std::endl;
        data[i] = v2(i);
        // std::cout << "data[i] = " << data[i] << std::endl;
      }
      return *this;
    }

    Vec &operator=(const Vec<S, T> &v2)
    {
      // std::cout << "Using the other operator=" << std::endl;
      for (int i = 0; i < S; i++)
      {
        data[i] = v2(i);
      }
      return *this;
    }

    template <typename T2>
    Vec &operator=(const VecExpr<T2> &v2)
    {
      // std::cout << "Using the other other operator=" << std::endl;
      if (v2.Size() != S)
        throw std::invalid_argument("VecExpr has wrong size");

      for (int i = 0; i < S; i++)
      {
        data[i] = v2(i);
      }
      return *this;
    }

    T operator()(size_t i) const
    {
      // std::cout << "Vec::operator()(" << i << ")" << std::endl;
      // std::cout << "data[i] = " << data[i] << std::endl;
      return data[i];
    }

    T &operator()(size_t i)
    {
      return data[i];
    }

    size_t Size() const { return S; }

    auto CopyToVectorView()
    {
      auto copyData = new T[S];
      for (size_t i = 0; i < S; i++)
      {
        // std::cout << "data[i] = " << data[i] << std::endl;
        copyData[i] = data[i];
      }
      return VectorView<T>(S, copyData);
    }
  };
}

#endif
