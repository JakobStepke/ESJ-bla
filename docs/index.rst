.. ASC-bla documentation master file, created by
   sphinx-quickstart on Tue Aug 29 06:39:02 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ASC-bla's documentation!
===================================

ASC-bla is a C++ library for basic linear algebra operations.
The library provides template classes **Vector** and **Matrix**.

Installation is via git-clone:

..  code-block::
    
    git clone https://github.com/TUWien-ASC/ASC-bla.git


To configure and build some tests do

..  code-block::

    cd ASC-bla
    mkdir build
    cd build
    cmake ..
    make
    

To use ASC-bla in your code, set the compiler include path properly, and include the header files

..  code-block::

    #include <vector.h>
    #include <matrix.h>

All objects are implemented in the namespace ASC_bla. To use them with less typing, you can set

..  code-block::
    
    namespace bla = ASC_bla;

or even

..  code-block::
    
    using namespace ASC_bla;

    

You can create vectors and compute with vectors like:

..  code-block:: cpp
                 
   Vector<double> x(5), y(5), z(5);
   for (int i = 0; i < x.Size(); i++)
      x(i) = i;
   y = 5.0
   z = x+3*y;
   cout << "z = " << z << endl;


For matrices you can choose between row-major (`RowMajor`) or column-major (`ColMajor`) storage,
default is row-major.

..  code-block:: cpp

   Matrix<double,RowMajor> m1(5,3), m2(3,3);
   for (int i = 0; i < m1.Height(); i++)
     for (int j = 0; j < m1.Width(); j++)
       m1(i,j) = i+j;
   m2 = 3.7;
   Matrix product = m1 * m2;
   
You can extract a rows or a columns from a matrix:

..  code-block:: cpp

   Vector col1 = product.Col(1);

You can perform different operations with matrices. As an example, we declare and output some vectors and matrices:

.. code-block:: cpp

  size_t n = 3;
  bla::Vector<double> x(n), y(n);

  bla::Matrix<double, bla::ORDERING::ColMajor> A(n, n);
  bla::Matrix<double, bla::ORDERING::ColMajor> B(n, n);

  for (size_t i = 0; i < x.Size(); i++)
  {
    x(i) = i;
    y(i) = 10;
  }

  for (size_t i = 0; i < A.Size_Rows(); i++)
    for (size_t j = 0; j < A.Size_Cols(); j++)
      A(i, j) = i + j;

  for (size_t i = 0; i < B.Size_Rows(); i++)
    for (size_t j = 0; j < B.Size_Cols(); j++)
      B(i, j) = 2*i + j;

  std::cout << "A:" << std::endl << A << std::endl;
  std::cout << "B:" << std::endl << B << std::endl;
  std::cout << "x:" << std::endl << x << std::endl;
  std::cout << "y:" << std::endl << y << std::endl;

This will result in the following output:

.. code-block::

  A:
  0 1 2
  1 2 3
  2 3 4
  B:
  0 1 2
  2 3 4
  4 5 6
  x:
  0, 1, 2
  y:
  10, 10, 10

You can multiply two matrices:

.. code-block:: cpp

  bla::Matrix<double, bla::ORDERING::ColMajor> C = A*B;
  std::cout << "C: " << std::endl << C << std::endl;

.. code-block::

  C:
  10 13 16
  16 22 28
  22 31 40
   
.. toctree::
   :maxdepth: 2
   :caption: Contents:





Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
