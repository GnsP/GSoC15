#!/usr/bin/python
from numpy import array
from scipy.linalg import lu
a = array([[1,1,1,10],[1,1,-1,0],[1,-1,1,4],[-1,1,1,6]])

pl, u = lu(a, permute_l=True)

print u
