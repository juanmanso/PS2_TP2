syms ax ay w T

A = [0 0 1 0 0  0  0  0;
     0 0 0 1 0  0  0  0;
     0 0 0 0 ax ay 0  0;
     0 0 0 0 0  0  ax ay;
     0 0 0 0  0 w 0 0;
     0 0 0 0 -w 0 0 0;
     0 0 0 0  0 0  0 w;
     0 0 0 0  0 0 -w 0];
 
 Ad = eye(size(A)) + A * T + A * A * T^2 / 2