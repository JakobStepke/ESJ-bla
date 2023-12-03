#ifndef FILE_VECTOR_H
#define FILE_VECTOR_H

#include <iostream>
#include <cmath>

#include "expression.h"
#include "simd.h"

using namespace ASC_HPC;

namespace ASC_bla
{


 
   
  template <typename T = double, typename TDIST = std::integral_constant<size_t,1> >
  class VectorView : public VecExpr<VectorView<T,TDIST>>
  {
  protected:
    T * data_;
    size_t size_;
    TDIST dist_;
  public:
    VectorView (size_t size, T * data)
      : data_(data), size_(size) { }
    
    VectorView (size_t size, TDIST dist, T * data)
      : data_(data), size_(size), dist_(dist) { }
    
    template <typename TB>
    VectorView & operator= (const VecExpr<TB> & v2)
    {
      for (size_t i = 0; i < size_; i++)
        data_[dist_*i] = v2(i);
      return *this;
    }

    VectorView & operator= (T scal)
    {
      for (size_t i = 0; i < size_; i++)
        data_[dist_*i] = scal;
      return *this;
    }

    VectorView & operator+= (T scal)
    {
      for (size_t i = 0; i < size_; i++)
        data_[dist_*i] += scal;
      return *this;
    }

    VectorView & operator-= (T scal)
    {
      for (size_t i = 0; i < size_; i++)
        data_[dist_*i] -= scal;
      return *this;
    }

    VectorView & operator+= (const VectorView & v2)
    {
      for (size_t i = 0; i < size_; i++)
        data_[dist_*i] += v2(i);
      return *this;
    }

    VectorView & operator-= (const VectorView & v2)
    {
      for (size_t i = 0; i < size_; i++)
        data_[dist_*i] -= v2(i);
      return *this;
    }

    VectorView & operator*= (T scal)
    {
      for (size_t i = 0; i < size_; i++)
        data_[dist_*i] *= scal;
      return *this;
    }

    VectorView & operator/= (T scal)
    {
      for (size_t i = 0; i < size_; i++)
        data_[dist_*i] /= scal;
      return *this;
    }
    
    auto View() const { return VectorView(size_, dist_, data_); }
    size_t Size() const { return size_; }
    auto Dist() const { return dist_; }    
    T & operator()(size_t i) { return data_[dist_*i]; }
    const T & operator()(size_t i) const { return data_[dist_*i]; }
    
    auto Data() const { return data_; }

    auto Range(size_t first, size_t next) const {
      return VectorView(next-first, dist_, data_+first*dist_);
    }

    auto Slice(size_t first, size_t slice) const {
      return VectorView<T,size_t> (size_/slice, dist_*slice, data_+first*dist_);
    }

    auto Norm2() const {
      double sum = 0;
      for (size_t i = 0; i < size_; i++)
        sum += data_[dist_*i]*data_[dist_*i];
      return sqrt(sum);
    }
      
  };
  
  template <typename T = double>
  class Vector : public VectorView<T>
  {
    typedef VectorView<T> BASE;
    using BASE::size_;
    using BASE::data_;
  public:
    Vector (size_t size) 
      : VectorView<T> (size, new T[size]) { ; }
      // Vector from list

      Vector (std::initializer_list<T> list)
      : VectorView<T> (list.size(), new T[list.size()]) {
        size_t i = 0;
        for (auto it = list.begin(); it != list.end(); ++it)
          data_[i++] = *it;}
    
    Vector (const Vector & v)
      : Vector(v.Size())
    {
      *this = v;
    }

    Vector (Vector && v)
      : VectorView<T> (0, nullptr)
    {
      std::swap(size_, v.size_);
      std::swap(data_, v.data_);
    }

    template <typename TB>
    Vector (const VecExpr<TB> & v)
      : Vector(v.Size())
    {
      *this = v;
    }
    
    
    ~Vector () { delete [] data_; }

    using BASE::operator=;
    Vector & operator=(const Vector & v2)
    {
      for (size_t i = 0; i < size_; i++)
        data_[i] = v2(i);
      return *this;
    }

    Vector & operator= (Vector && v2)
    {
      for (size_t i = 0; i < size_; i++)
        data_[i] = v2(i);
      return *this;
    }

    Vector & operator* (Vector && v2)
    {
      Vector ret = Vector(size_);
      for (size_t i = 0; i < size_; i++)
        ret(i) = this(i)*v2(i);
      return ret;
    }
    
  };


  template <typename ...Args>
  std::ostream & operator<< (std::ostream & ost, const VectorView<Args...> & v)
  {
    if (v.Size() > 0)
      ost << v(0);
    for (size_t i = 1; i < v.Size(); i++)
      ost << ", " << v(i);
    return ost;
  }

  template <size_t SW>
auto InnerProduct (size_t n, const VectorView<double, size_t> x, size_t dx,
                   const VectorView<double, size_t> y, size_t dy)
{
  SIMD<double,SW> sum{0.0};
  for (size_t i = 0; i < n; i++)
    {
      // sum += px[i] * SIMD<double,SW>(py+i*dy);
      sum = FMA(SIMD<double,SW>(x.Data()[i]), SIMD<double,SW>(y.Data()[i*dy]), sum);
    }
    /*
    std::cout << "sum = " << sum << std::endl;
    std::cout << "x = " << x << std::endl;
    std::cout << "y = " << y << std::endl;
    */

  double first_val = sum.GetFirst();

  return first_val;
}

template <size_t SW>
auto InnerProduct8 (size_t n, const VectorView<double, size_t> x0, 
                    const VectorView<double, size_t> x1, 
                    const VectorView<double, size_t> x2, 
                    const VectorView<double, size_t> x3, 
                    const VectorView<double, size_t> x4, 
                    const VectorView<double, size_t> x5, 
                    const VectorView<double, size_t> x6, 
                    const VectorView<double, size_t> x7, 
                    const VectorView<double, size_t> y, size_t dy)
{
  SIMD<double,SW> sum0{0.0};
  SIMD<double,SW> sum1{0.0};
  SIMD<double,SW> sum2{0.0};
  SIMD<double,SW> sum3{0.0};
  SIMD<double,SW> sum4{0.0};
  SIMD<double,SW> sum5{0.0};
  SIMD<double,SW> sum6{0.0};
  SIMD<double,SW> sum7{0.0};

  for (size_t i = 0; i < n; i++)
    {
      // sum += px[i] * SIMD<double,SW>(py+i*dy);
      sum0 = FMA(SIMD<double,SW>(x0.Data()[i]), SIMD<double,SW>(y.Data()[i*dy]), sum0);
      sum1 = FMA(SIMD<double,SW>(x1.Data()[i]), SIMD<double,SW>(y.Data()[i*dy]), sum1);
      sum2 = FMA(SIMD<double,SW>(x2.Data()[i]), SIMD<double,SW>(y.Data()[i*dy]), sum2);
      sum3 = FMA(SIMD<double,SW>(x3.Data()[i]), SIMD<double,SW>(y.Data()[i*dy]), sum3);
      sum4 = FMA(SIMD<double,SW>(x4.Data()[i]), SIMD<double,SW>(y.Data()[i*dy]), sum4);
      sum5 = FMA(SIMD<double,SW>(x5.Data()[i]), SIMD<double,SW>(y.Data()[i*dy]), sum5);
      sum6 = FMA(SIMD<double,SW>(x6.Data()[i]), SIMD<double,SW>(y.Data()[i*dy]), sum6);
      sum7 = FMA(SIMD<double,SW>(x7.Data()[i]), SIMD<double,SW>(y.Data()[i*dy]), sum7);
    }
  return std::tuple(sum0.GetFirst(), sum1.GetFirst(), sum2.GetFirst(), sum3.GetFirst(),
                    sum4.GetFirst(), sum5.GetFirst(), sum6.GetFirst(), sum7.GetFirst());
}
  
}

#endif
