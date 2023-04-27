import numpy as np

print("Solving x for matrix eq: Ax=b ...") 


# Define matrix A using Numpy arrays 
# linalg.solve is the function of NumPy to solve a system of linear scalar equations 
# Sample example:
A = np.array([[7,0,-2,-1,0], [0,13,-4,0,-1], [-1,-1,3,0,-1], [-1,0,0,6,-1], [0,-1,-4,-2,15]])
b = np.array([4,8,0,0,0])

print("A=",A)
print("b=",b)
print()
#np.allclose(np.dot(A, x), b)

x=np.linalg.solve(A,b)
print ("Solutions: x=\n",x)
I=2-x[0]-x[1]
print('I=',I)
print("R_eff=",1/I)


