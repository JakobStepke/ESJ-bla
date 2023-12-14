#include <sstream>
#include <pybind11/pybind11.h>

#include "vector.h"
#include "ordering.h"
#include "matrix.h"

using namespace ASC_bla;
namespace py = pybind11;

PYBIND11_MODULE(bla, m)
{
     m.doc() = "Basic linear algebra module"; // optional module docstring

     py::class_<Vector<double>>(m, "Vector")
         .def(py::init<size_t>(),
              py::arg("size"), "create vector of given size")
         .def("__len__", &Vector<double>::Size,
              "return size of vector")

         .def("__setitem__", [](Vector<double> &self, int i, double v)
              {
        if (i < 0) i += self.Size();
        if (i < 0 || i >= self.Size()) throw py::index_error("vector index out of range");
        self(i) = v; })
         .def("__getitem__", [](Vector<double> &self, int i)
              { return self(i); })

         .def("__setitem__", [](Vector<double> &self, py::slice inds, double val)
              {
        size_t start, stop, step, n;
        if (!inds.compute(self.Size(), &start, &stop, &step, &n))
          throw py::error_already_set();
        self.Range(start, stop).Slice(0,step) = val; })

         .def("__add__", [](Vector<double> &self, Vector<double> &other)
              { return Vector<double>(self + other); })

         .def("__rmul__", [](Vector<double> &self, double scal)
              { return Vector<double>(scal * self); })

         .def("__str__", [](const Vector<double> &self)
              {
        std::stringstream str;
        str << self;
        return str.str(); })

         .def(py::pickle(
             [](Vector<double> &self) { // __getstate__
                  /* return a tuple that fully encodes the state of the object */
                  return py::make_tuple(self.Size(),
                                        py::bytes((char *)(void *)&self(0), self.Size() * sizeof(double)));
             },
             [](py::tuple t) { // __setstate__
                  if (t.size() != 2)
                       throw std::runtime_error("should be a 2-tuple!");

                  Vector<double> v(t[0].cast<size_t>());
                  py::bytes mem = t[1].cast<py::bytes>();
                  std::memcpy(&v(0), PYBIND11_BYTES_AS_STRING(mem.ptr()), v.Size() * sizeof(double));
                  return v;
             }));

     py::class_<Matrix<double, ORDERING::RowMajor>>(m, "Matrix")
         .def(py::init<size_t, size_t>(),
              py::arg("width"), py::arg("length"), "create matrix of given size")
         .def("width", &Matrix<double, ORDERING::RowMajor>::Size_Rows,
              "return width of matrix")
         .def("height", &Matrix<double, ORDERING::RowMajor>::Size_Cols,
              "return height of matrix")
         .def("transpose", &Matrix<double, ORDERING::RowMajor>::transpose,
              "return transpose of matrix")
         .def("T", &Matrix<double, ORDERING::RowMajor>::transpose,
              "return transpose of matrix")
         //.def("determinant", &Matrix<double, ORDERING::RowMajor>::Determinant, "return determinant of matrix")
         //.def("inverse", &Matrix<double, ORDERING::RowMajor>::Inverse, "return inverse of matrix")

         .def("__setitem__", [](Matrix<double, ORDERING::RowMajor> &self, std::tuple<int, int> ind, double v)
              {
        if (std::get<0>(ind) < 0) std::get<0>(ind) += self.Size_Rows();
        if (std::get<1>(ind) < 0) std::get<1>(ind) += self.Size_Cols();
        if (std::get<0>(ind) < 0 || std::get<0>(ind) >= self.Size_Rows() ||
            std::get<1>(ind) < 0 || std::get<1>(ind) >= self.Size_Cols())
          throw py::index_error("matrix index out of range");
        self(std::get<0>(ind), std::get<1>(ind)) = v; })

         .def("__getitem__",
              [](Matrix<double, ORDERING::RowMajor> self, std::tuple<int, int> ind)
              {
                   return self(std::get<0>(ind), std::get<1>(ind));
              })

         .def("__add__", [](Matrix<double, ORDERING::RowMajor> &self, Matrix<double, ORDERING::RowMajor> &other)
              { return Matrix<double, ORDERING::RowMajor>(self + other); })

         .def("__mul__", [](Matrix<double, ORDERING::RowMajor> &self, Matrix<double, ORDERING::RowMajor> &other)
              { return Matrix<double, ORDERING::RowMajor>(self * other); })

         .def("__mul__", [](Matrix<double, ORDERING::RowMajor> &self, double &scal)
              { return Matrix<double, ORDERING::RowMajor>(self * scal); })

         .def("__rmul__", [](Matrix<double, ORDERING::RowMajor> &self, double &scal)
              { return Matrix<double, ORDERING::RowMajor>(self * scal); })

         .def("__mul__", [](Matrix<double, ORDERING::RowMajor> &self, Vector<double> &vec)
              { return Vector<double>(self * vec); })

         .def("__rmul__", [](Matrix<double, ORDERING::RowMajor> &self, Vector<double> &vec)
              { return Vector<double>(vec * self); })

         .def("__str__", [](const Matrix<double, ORDERING::RowMajor> &self)
              {
        std::stringstream str;
        str << self;
        return str.str(); })

         .def(py::pickle(
             [](Matrix<double, ORDERING::RowMajor> &self) { // __getstate__
                  /* return a tuple that fully encodes the state of the object */
                  return py::make_tuple(self.Size_Rows(), self.Size_Cols(),
                                        py::bytes((char *)(void *)&self(0, 0), self.Size_Rows() * self.Size_Cols() * sizeof(double)));
             },
             [](py::tuple t) { // __setstate__
                  if (t.size() != 3)
                       throw std::runtime_error("should be a 3-tuple!");

                  Matrix<double, ORDERING::RowMajor> m(t[0].cast<size_t>(), t[1].cast<size_t>());
                  py::bytes mem = t[2].cast<py::bytes>();
                  std::memcpy(&m(0, 0), PYBIND11_BYTES_AS_STRING(mem.ptr()), m.Size_Rows() * m.Size_Cols() * sizeof(double));
                  return m;
             }));
}
