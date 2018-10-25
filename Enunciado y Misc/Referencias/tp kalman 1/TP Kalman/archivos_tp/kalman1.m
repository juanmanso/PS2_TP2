function xhat = kalman1(A,B,C,D,E,u,Q,R,xhat0,P0,y)

if (length(size(A))==2)
    Ak = A;
end

if (length(size(B))==2)
    Bk = B;
end

if (length(size(C))==2)
    Ck = C;
end

if (length(size(R))==2)
    Rk = R;
end

if (length(size(Q))==2)
    Qk = Q;
end

if (length(size(D))==2)
    Dk = D;
end

if (length(size(E))==2)
    Ek = E;
end

xhat = zeros(size(A,2),size(y,2));
xhat(:,1) = xhat0;
P = P0;

if (length(u)==1 && u==0)
    
    for k=1:length(xhat)-1
        
        if (length(size(A))==3)
            Ak = squeeze(A(k,:,:));
        end
        
        if (length(size(B))==3)
            Bk = squeeze(B(k,:,:));
        end
        
        if (length(size(C))==3)
            Ck = squeeze(C(k,:,:));
        end
        
        if (length(size(R))==3)
            Ck = squeeze(R(k,:,:));
        end
        
        if (length(size(Q))==3)
            Qk = squeeze(Q(k,:,:));
        end
               
        % Predicción:
        xhat(:,k+1) = Ak*xhat(:,k);
        P = Ak*P*Ak' + Bk*Qk*Bk';
        % Actualización:
        K = P*Ck'*(inv(Rk+Ck*P*Ck'));
        xhat(:,k+1) = xhat(:,k+1) + K*(y(:,k+1) -Ck*xhat(:,k+1));
        P = (eye(length(Ak)) - K*Ck)*P;       
        
    end

else
    
    for k=1:length(xhat)-1
        if (length(size(A))==3)
            Ak = squeeze(A(k,:,:));
        end
        
        if (length(size(B))==3)
            Bk = squeeze(B(k,:,:));
        end
        
        if (length(size(C))==3)
            Ck = squeeze(C(k,:,:));
        end
        
        if (length(size(R))==3)
            Ck = squeeze(R(k,:,:));
        end
        
        if (length(size(Q))==3)
            Qk = squeeze(Q(k,:,:));
        end
        
        if (length(size(D))==3)
            Dk = squeeze(D(k,:,:));
        end
        
        if (length(size(E))==3)
            Ek = squeeze(E(k,:,:));
        end
        
        % Predicción
        xhat(:,k+1) = Ak*xhat(:,k) + Dk*u(:,k);
        P = Ak*P*Ak' + Bk*Qk*Bk';
        % Actualización:
        K = P*Ck'*(inv(Rk+Ck*P*Ck'));
        xhat(:,k+1) = xhat(:,k+1) + K*(y(:,k+1)-Ck*xhat(:,k+1)-Ek*u(:,k+1));
        P = (eye(length(Ak)) - K*Ck)*P;
    end
end










