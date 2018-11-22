library(Matrix)
A = matrix(c(1,1,0,0,0,0,1,1,0,1,0,0,0,0,1,0,1,0,1,0,0,1,0,1,1,0,0,0,0,0,0,1), 
           nrow=8, ncol=4,byrow = TRUE);

b = matrix(c(20,20,10,10,20,20,10,10), 
           nrow=8, ncol=1,byrow = TRUE);

p = matrix(c(18,18,8,8,18,18,8,8), 
           nrow=8, ncol=1,byrow = TRUE);

AA = t(A)%*%A
expand(lu(AA))
solve(AA)
L = expand(lu(AA))$L
paste( c(matrix(t(L))), collapse = ", ")

L = expand(lu(AA))$L
paste( c(matrix(t(L))), collapse = ", ")

U = expand(lu(AA))$U
paste( c(matrix(t(U))), collapse = ", ")


INV =  solve(AA)
paste( c(matrix(t(INV))), collapse = ", ")


AA = AA[c(4,3,2,1),]
paste( c(t(AA)), collapse = ", ")
paste( c(t(solve(AA))), collapse = ", ")

AA = AA[c(1,3,2,4),]
paste( c(t(AA)), collapse = ", ")
paste( c(t(solve(AA))), collapse = ", ")



