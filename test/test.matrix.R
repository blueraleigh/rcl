library(rcl)


A = matrix(runif(25), 5, 5)
B = matrix(runif(25), 5, 5)
C = .Call("rcl_test_matrix_product", A, B, PACKAGE="rcl")
stopifnot(all.equal(A %*% B, C))
C = .Call("rcl_test_matrix_cum", A, B, PACKAGE="rcl")
stopifnot(all.equal(A %*% B %*% B %*% B %*% B, C))

A = matrix(runif(24), 6, 4)
B = matrix(runif(48), 4, 8)
C = .Call("rcl_test_matrix", A, B, PACKAGE="rcl")
stopifnot(all.equal(A %*% B, C))




R = matrix(0, 3, 3)
R[1, 2] = 0.10
R[1, 3] = 0.05
R[2, 1] = 0.01
R[2, 3] = 0.002
R[3, 1] = 0.2
R[3, 2] = 0.04
diag(R) = -rowSums(R)

I = matrix(0, nrow(R), nrow(R))
diag(I) = 1
t = 100
k = 10L
all.equal(as.matrix(Matrix::expm(R*t)), .Call("rcl_test_matrix_exponential", R, t, I, k, PACKAGE="rcl"))
system.time(Matrix::expm(R*t))
system.time(.Call("rcl_test_matrix_exponential", R, t, PACKAGE="rcl"))
.Call("rcl_test_matrix_exponential", R, t, PACKAGE="rcl")


A = matrix(runif(15), 3, 5)
x = runif(5)
.Call("rcl_test_matrix_vector_product", A, x, PACKAGE="rcl")

A %*% x
