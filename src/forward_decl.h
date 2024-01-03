#ifndef  FORWARD_DECL_H
#define  FORWARD_DECL_H

namespace ASC_bla
{

    enum class ORDERING { ColMajor, RowMajor };

    template <typename T = double, typename TDIST = std::integral_constant<size_t, 1>>
    class VectorView;

    template <typename T = double, ORDERING ORD = ORDERING::ColMajor>
    class MatrixView;
}

#endif