function eta = ruido_blanco_gaussiano(R,n)

[U, S, V] = svd(R);
eta = zeros(length(R),n);
for i=1:n
    eta(:,i) = U*sqrt(S)*randn(length(R),1);
end
