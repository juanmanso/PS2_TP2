function [x, P, g] = kalman_qr(x_prev, P_prev, A, B, C, Q, R, y)
    
    n = length(x_prev);
    p = length(y);
    q = size(Q, 1);
    
	if (y ~= 0)
        %Es una actualización.
        
        %Arreglo del algoritmo QR.
        Arreglo = [chol(R)',           C * chol(P_prev)',           zeros(p, q) ;
                   zeros(n, p),        A * chol(P_prev)',           B * chol(Q)';
                   -y * inv(chol(R)),  x_prev' * inv(chol(P_prev)), zeros(1, q)];

        %Factorización QR
        [Q_qr, R_qr] = qr(Arreglo');

        %Separo X, Z, W1 y W2
        R_qr = R_qr';
        X = R_qr(1 : p, 1 : p);
        Z = R_qr(p + 1 : p + n, p + 1 : p + n);
        W1 = R_qr(end, 1 : p);
        W2 = R_qr(end, p + 1 : p + n);

        %Computo la estimación, la P y las innovaciones.
        P =  Z * Z';
        x =  Z * W2';
        g = -X * W1';
    else
        %Solo propago.
        x = A * x_prev; 
        P = A * P_prev * A' + B * Q * B';
        g = y' - C * x_prev;
    end
end