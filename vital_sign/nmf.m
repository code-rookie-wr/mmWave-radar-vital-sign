function [W, H] = nmf(V, K, MAXITER)
F = size(V,1);
T = size(V,2);
% rand('seed¡¯,0)
W = 1+rand(F, K);
H = 1+rand(K, T);
ONES = ones(F,T);

Ds=[];
for i=1:MAXITER
    % update activations
    H = H .* (W'*( V./(W*H+eps))) ./ (W'*ONES);
    % update dictionaries
    W = W .* ((V./(W*H+eps))*H') ./(ONES*H');
    
    D = sum(sum(V.*(log(V)-log(W*H))-V+(W*H)));
    Ds = [Ds; D];
end
% normalize W to sum to 1
sumW = sum(W);
W = W*diag(1./sumW);
H = diag(sumW)*H;