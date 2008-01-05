    N=10 ;
    A = zeros(10,10)
    A = ones(10,10)
    A = eyes(10,10)
    A(N/2,:) = -1
    A(1:3,N/2) = -2
    A'
  
