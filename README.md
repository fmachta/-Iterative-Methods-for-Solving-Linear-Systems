# Iterative-Methods-for-Solving-Linear-Systems
In the realm of numerical solutions for linear systems of equations, iterative methods provide a systematic process to approximate solutions. This report describes the implementation of two iterative methods: the Jacobi method and the Gauss-Seidel method, in C++.


## Introduction

In the realm of numerical solutions for linear systems of equations, iterative methods provide a systematic process to approximate solutions. This report describes the implementation of two iterative methods: the Jacobi method and the Gauss-Seidel method, in C++.


## Methods Used


### Jacobi Iterative Method

- It decomposes the matrix A into a diagonal component D and a remainder R.
- The iterative formula: x = (b - R*x0) / D
- The iterations continue until the L2 norm (Euclidean distance) of the difference between subsequent solutions becomes less than a specified error.


### Gauss-Seidel Iterative Method

- This method enhances the Jacobi method by using the updated values as soon as they are known.
- The matrix A is decomposed into a lower triangular matrix L and an upper triangular matrix U.
- The iterative formula: x = inv(L) * (b - U*x0)
- Like the Jacobi method, iterations continue until the L2 norm of the difference is below the specified error.


## Implementation Overview
Diagonal Dominance Check: Before beginning iterations, the system checks if the matrix is strictly diagonally dominant because for these methods to guarantee convergence, diagonal dominance is essential.
User Interaction:
- The user specifies the number of equations.
- Coefficients can be input either through the command line or from a file.
- The user provides the desired error for convergence and an initial guess for the solution.

  
Iteration and Convergence: The program iteratively computes solutions. If the desired accuracy isn't reached within 50 iterations, the method terminates, reporting the result of the 50th iteration.


## Results & Output


For each method, the solution after each iteration is printed. The final output is the solution vector [x1, x2, ...]T for each method. The T indicates that it's a column vector.


### Command Line
```
Enter the number of equations (n <= 10): 3
Enter coefficients through:
1. Command Line
2. File
Choose (1/2): 1
Enter the augmented coefficient matrix (including b values):
5 -1 0 7

-1 3 -1 4

0 -1 2 5
Enter the desired stopping error: 0 0 0
Enter the starting solution:
0 0 0
Solution using Jacobi: [2.0870 3.4348 4.2174 ]T
Solution using Gauss-Seidel: [2.0870 3.4348 4.2174 ]T
```


### File (using Equations.txt):

```
Enter the number of equations (n <= 10): 3
Enter coefficients through:
1. Command Line
2. File
Choose (1/2): 2
Enter the filename: equations.txt
Enter the desired stopping error: 0 0 0
Enter the starting solution:
0 0 0
Solution using Jacobi: [2.0870 3.4348 4.2174 ]T
Solution using Gauss-Seidel: [2.0870 3.4348 4.2174 ]T
```

## Conclusion

Iterative methods, such as Jacobi and Gauss-Seidel, provide practical tools for approximating solutions to linear systems. While their convergence for arbitrary systems isn't guaranteed, they can efficiently solve diagonally dominant systems. This C++ program offers a robust interface for users to employ these methods and understand their workings.

