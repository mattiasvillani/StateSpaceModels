function u = randmn(my,COVAR)
u=my+chol(COVAR)'*randn(size(COVAR,1),1);

