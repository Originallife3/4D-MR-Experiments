function [X,tnn] = prox_tnn(Y,lambda1)

   [n1,n2,n3,n4] = size(Y);
for j=1:n4
    n12 = min(n1,n2);
    % get core matrix
    [U, S, V] = tsvd(Y(:,:,:,j));
    reshapeS = zeros(n12, n3);
    for i = 1:n12
        reshapeS(i,:) = S(i,i,:);
    end
    [A,D,B] = svd(reshapeS, 'econ');
    VT=B';
%   Dt = diag(D);
    D=diag(SoftThresh(diag(D),D(1)*lambda1));
    L = A* D * VT;
    T = zeros(n12,n12,n3,n4);
    for i = 1:n12
      Tj(i,i,:) = L(i,:);
    end
    T(:,:,:,j) = tprod(tprod(U,Tj), tran(V));  
 
    Xj = zeros(n1,n2,n3);
    X = zeros(n1,n2,n3,n4);
%  Y2 = fft2(T,[],3);
   Y2 = fft2(T);
   tnn = 0;

[U,S,V] = svd(Y2(:,:,1,j),'econ');
St = diag(S);
St=diag(SoftThresh(diag(St),St(1)*lambda1));
Xj(:,:,1) = U*diag(St)*V';

tnn = tnn+sum(St);
halfn3 = round(n3/2);
for i = 2 : halfn3
 [U,S,V] = svd(Y2(:,:,i,j),'econ');
 St = diag(S);
 St=diag(SoftThresh(diag(St),St(1)*lambda1));
 Xj(:,:,i) = U*diag(St)*V';
 tnn = tnn+sum(St)*2;
end
    Xj(:,:,n3+2-i) = conj(Xj(:,:,i));
% if n3 is even
if mod(n3,2) == 0
    i = halfn3+1;
    [U,S,V] = svd(Y2(:,:,i,j),'econ');
    St = diag(S);
    St=diag(SoftThresh(diag(St),St(1)*lambda1));
    Xj(:,:,i) = U*diag(St)*V';
    tnn = tnn+sum(St);
end
X(:,:,:,j) = Xj;
tnn = tnn/n3;
end
X = ifft2(X);
end


function y=SoftThresh(x,p)
y=(abs(x)-p).*x./abs(x).*(abs(x)>p);
y(isnan(y))=0;
end    
