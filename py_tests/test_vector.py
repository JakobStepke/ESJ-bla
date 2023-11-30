# search for libraray like bla.cpython-312-darwin.so in the build directory:
# import sys
# sys.path.append('/Users/joachim/texjs/lva/ws2324/ScientificComputing/ASC-bla/build')
# from bla import Vector

# import from the installed ASCsoft package:
from ASCsoft.bla import Vector
from ASCsoft.bla import Matrix

x = Vector(3)
y = Vector(3)

for i in range(len(x)):
    x[i] = i
y[:] = 2    

print ("x =", x)
print ("y =", y)
print ("x+3*y =", x+3*y)


x = Vector(10)
x[0:] = 1
print (x)

x[3:7] = 2
print (x)

x[0:10:2] = 3
print (x)

A = Matrix(3,3)
B = Matrix(3,3)

for i in range(3):
    for j in range(3):
        A[i,j] = i+j
        B[i,j] = i-j

print ("A =", A)
print ("B =", B)
print ("A*B =", A*B)
print ("A+3*B =", A+3*B)
print ("A*B*A =", A*B*A)

# Vector-matrix multiplication:
print ("A*x =", A*y)
print ("x*A =", y*A)

# Width of A:
print ("width of A =", A.width())

# Height of A:
print ("height of A =", A.height())

# Transpose of A:
print ("transpose of A =", A.transpose())
print ("transpose of A =", A.T())
print ("det of A =", A.determinant())
print ("inverse of A =", A.inverse())

print ("A*Inverse(A) =", A.inverse()*A)

print ("A[0,0] =", A[0,0])