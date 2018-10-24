# Polar Decomposition

We provide an algorithm for computing the polar decomposition $$A=UH$$.
We initially iterate using the Newton method and switch to the matrix multiplication rich
Newton-Schulz method once it is guaranteed to converge.
The motivation for this is that matrix multiplication is very fast
on high performance computers.

We also provide a routine to compute the square root of a symmetric positive
definite matrix.
