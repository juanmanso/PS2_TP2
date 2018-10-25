function xhat = kalman2(A,B,C,Q,R,xhat0,P0,y)

n = size(A,2); % Dimensión del vector a estimar
q = size(R,2); % Dimensión del ruido de medición
p = size(Q,2); % Dimensión del ruido del proceso
N = size(y,2); % Cantidad de muestras

if (length(size(A))==2)
    Ak = A;
else
    Ak = squeeze(A(1,:,:)); % A[0]
end

if (length(size(B))==2)
    Bk = B;
else
    Bk = squeeze(B(1,:,:)); % B[0]
end

if (length(size(C))==2)
    Ck = C;
else
    Ck = squeeze(C(2,:,:)); % C[1]
end

if (length(size(R))==2)
    Rk = R;
else
    Rk = squeeze(R(2,:,:)); % R[1]
end

if (length(size(Q))==2)
    Qk = Q;
else
    Qk = squeeze(Q(2,:,:)); % Q[1]
end

xhat = zeros(n,N);
xhat(:,1) = xhat0; % xhat[0|0]

% Predicción inicial:
xhat(:,2) = Ak*xhat(:,1); % xhat[1|0]
Pk = Ak*P0*Ak' + Bk*Qk*Bk'; % P[1|0]
Z = chol(Pk); % Z[1] = P[1|0]^(1/2)
R_sqrt = chol(Rk); % R[1]^(1/2)
Q_sqrt = chol(Qk); % Q[1]^(1/2)

if (length(size(A))==3)
    Ak = squeeze(A(2,:,:)); % A[1]
end

if (length(size(B))==3)
    Bk = squeeze(B(2,:,:)); % B[1]
end

M = [R_sqrt     Ck*Z zeros(q,p);
     zeros(p,q) Ak*Z  Bk*Q_sqrt];
 
[Q_fact, R_fact] = qr(M',0);
for i=1:length(R_fact)
    if(R_fact(i,i)<0)
        R_fact(i,i) = -R_fact(i,i);
        Q_fact(:,i) = -Q_fact(:,i);
    end
end

% Para actualización necesito:
X_inv = inv(R_fact(1:q,1:q)'); % X_inv[1] = (R[1]+C[1]*P[1|0]*C[1]')^(-1/2)
K = (Z'*Z)*Ck'*(X_inv'*X_inv); % K[1] = P[1|0]*C[1]'*(R[1]+C[1]*P[1|0]*C[1]')^(-1)

% Para predicción necesito:
fQ = [-y(:,2)'*inv(R_sqrt') xhat(:,2)'*inv(Z') zeros(1,p)]*Q_fact; % fQ = [W1 W2 W3] 
W = fQ(q+1:q+n); % W = W2

% Actualización:
xhat(:,2) = xhat(:,2) + K*(y(:,2)-Ck*xhat(:,2)); % x[1|1] = x[1|0] + K[1]*(y[1]-C[1]*xhat[1|0])

% Predicción:
Z = R_fact(q+1:end,q+1:end)'; % Z[2] = P[2|1]^(1/2)
xhat(:,3) = Z*W'; % x[2|1]

% Predicción iterativa:
for k=3:N-1
    if (length(size(A))==3)
        Ak = squeeze(A(k,:,:)); % A[k]
    end
    
    if (length(size(B))==3)
        Bk = squeeze(B(k,:,:)); % B[k]
    end
    
    if (length(size(C))==3)
        Ck = squeeze(C(k,:,:)); % C[k]
    end
    
    if (length(size(R))==3)
        Rk = squeeze(R(k,:,:)); % R[k]
    end
    
    if (length(size(Q))==3)
        Qk = squeeze(Q(k,:,:)); % Q[k]
    end
    
    R_sqrt = chol(Rk); % R[k]^(1/2)
    Q_sqrt = chol(Qk); % Q[k]^(1/2)
    M = [R_sqrt     Ck*Z zeros(q,p);
        zeros(p,q) Ak*Z  Bk*Q_sqrt];
    
    [Q_fact, R_fact] = qr(M',0);
    for i=1:length(R_fact)
        if(R_fact(i,i)<0)
            R_fact(i,i) = -R_fact(i,i);
            Q_fact(:,i) = -Q_fact(:,i);
        end
    end
    
    
    % Actualización:
    X_inv = inv(R_fact(1:q,1:q)'); % X_inv[k] = (R[k]+C[k]*P[k|k-1]*C[k]')^(-1/2)
    K = (Z'*Z)*Ck'*(X_inv'*X_inv); % K[k] = P[k|k-1]*C[k]'*(R[k]+C[k]*P[k|k-1]*C[k]')^(-1)
    xhat(:,k) = xhat(:,k) + K*(y(:,k)-Ck*xhat(:,k)); % x[k|k] = x[k|k-1] + K[k]*(y[k]-C[k]*xhat[k|k-1])

    % Predicción:
    Z = R_fact(q+1:end,q+1:end)'; % Z[k+1] = P[k+1|k]^(1/2)
    fQ = [-y(:,k)'*inv(R_sqrt') W zeros(1,p)]*Q_fact; % fQ = [W1 W2 W3] 
    W = fQ(q+1:q+n); % W[k+1] = W2
    xhat(:,k+1) = Z*W'; % xhat[k+1|k] = Z[k+1]*W[k+1]'
end
 
    
 