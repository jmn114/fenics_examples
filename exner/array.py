import math
import numpy as np

a = np.array([1., 0., 0., 1., 2., 1])
b = np.array([ 1.5, 0.0, -0.9, 0.0, 0.8, 0.0])

print b

print b.shape
print b.ndim

print a * b

print np.absolute(b)

print 2*b
print b**3    