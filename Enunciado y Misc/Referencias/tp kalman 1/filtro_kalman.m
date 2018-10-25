function x = filtro_kalman(A,B,C,R,Q)

I = eye(3);
O = zeros(3);

xhat00 = 0;
P00 = [I