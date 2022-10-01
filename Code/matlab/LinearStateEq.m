function g = LinearStateEq(u,xLag, A, B)

g = A*xLag + B*u;