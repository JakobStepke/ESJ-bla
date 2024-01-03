#ifndef FILE_EXPRESSION_H
#define FILE_EXPRESSION_H

#include <ostream>
#include "forward_decl.h"
namespace ASC_bla
{

  template <typename T>
  class VecExpr
  {
  public:
    auto Upcast() const { return static_cast<const T&> (*this); }
    size_t Size() const { return Upcast().Size(); }
    auto operator() (size_t i) const { return Upcast()(i); }
  };

  template <typename T, ORDERING ORD>
  class MatExpr
  {
  public:
    auto Upcast() const { return static_cast<const T&> (*this); }
    size_t Rows() const { return Upcast().Size_Rows(); }
    size_t Cols() const { return Upcast().Size_Cols(); }
    auto operator() (size_t row, size_t col) const { return Upcast()(row, col); }
  };
  
  template <typename TA, typename TB, ORDERING ORD>
  class SumMatExpr : public MatExpr<SumMatExpr<TA,TB, ORD>, ORD>
  {
    TA a_;
    TB b_;
  public:
    SumMatExpr (TA a, TB b) : a_(a), b_(b) { }

    auto operator() (size_t i, size_t j) const { return a_(i, j)+b_(i, j); }
    size_t Size_Rows() const { return a_.Size_Rows(); }      
    size_t Size_Cols() const { return a_.Size_Cols(); }      
  };
  
  template <typename TA, typename TB, ORDERING ORD>
  auto operator+ (const MatExpr<TA, ORD> & a, const MatExpr<TB, ORD> & b)
  {
    return SumMatExpr(a.Upcast(), b.Upcast());
  }

  template <typename TA, typename TB, ORDERING ORD>
  class ProdMatExpr : public MatExpr<ProdMatExpr<TA,TB, ORD>, ORD>
  {
    TA a_;
    TB b_;
  public:
    ProdMatExpr (TA a, TB b) : a_(a), b_(b) { }

    auto operator() (size_t i, size_t j) const {
      a_.Row(i);
      b_.Col(j);
      return a_.Row(i)*b_.Col(j); 
      }
    size_t Size_Rows() const { return a_.Size_Rows(); }      
    size_t Size_Cols() const { return a_.Size_Cols(); }      
  };
  
  template <typename TA, typename TB, ORDERING ORD>
  auto operator* (const MatExpr<TA, ORD> & a, const MatExpr<TB, ORD> & b)
  {
    return ProdMatExpr(a.Upcast(), b.Upcast());
  }

  template <typename TA, typename TB, ORDERING ORD>
  class ProdMatVecExpr : public MatExpr<ProdMatVecExpr<TA,TB, ORD>, ORD>
  {
    TA a_;
    TB b_;
  public:
    ProdMatVecExpr (TA a, TB b) : a_(a), b_(b) { }

    auto operator() (size_t i) const {
      a_.Row(i);
      return a_.Row(i)*b_; 
      }
    size_t Size() const { return b_.Size(); }    
  };

  template <typename TA, typename TB, ORDERING ORD>
  auto operator* (const MatExpr<TA, ORD> & a, const VecExpr<TB> & b)
  {
    return ProdMatVecExpr(a.Upcast(), b.Upcast());
  }
 
  template <typename TA, typename TB>
  class SumVecExpr : public VecExpr<SumVecExpr<TA,TB>>
  {
    TA a_;
    TB b_;
  public:
    SumVecExpr (TA a, TB b) : a_(a), b_(b) { }

    auto operator() (size_t i) const { return a_(i)+b_(i); }
    size_t Size() const { return a_.Size(); }      
  };
  
  template <typename TA, typename TB>
  auto operator+ (const VecExpr<TA> & a, const VecExpr<TB> & b)
  {
    return SumVecExpr(a.Upcast(), b.Upcast());
  }

  template <typename TA, typename TB>
  auto operator- (const VecExpr<TA> & a, const VecExpr<TB> & b)
  {
    return SumVecExpr(a.Upcast(), -b.Upcast());
  }



  
  template <typename TSCAL, typename TV>
  class ScaleVecExpr : public VecExpr<ScaleVecExpr<TSCAL,TV>>
  {
    TSCAL scal_;
    TV vec_;
  public:
    ScaleVecExpr (TSCAL scal, TV vec) : scal_(scal), vec_(vec) { }

    auto operator() (size_t i) const { return scal_*vec_(i); }
    size_t Size() const { return vec_.Size(); }      
  };
  
  template <typename T>
  auto operator* (double scal, const VecExpr<T> & v)
  {
    return ScaleVecExpr(scal, v.Upcast());
  }

  template <typename T>
  auto operator* (const VecExpr<T> & v, double scal)
  {
    return ScaleVecExpr(scal, v.Upcast());
  }

  template <typename T>
  auto operator- (const VecExpr<T> & v)
  {
    return ScaleVecExpr(-1.0, v.Upcast());
  }



  template <typename T>
  std::ostream & operator<< (std::ostream & ost, const VecExpr<T> & v)
  {
    if (v.Size() > 0)
      ost << v(0);
    for (size_t i = 1; i < v.Size(); i++)
      ost << ", " << v(i);
    return ost;
  }
  
}
 
#endif
