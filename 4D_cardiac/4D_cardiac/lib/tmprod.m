function X = tmprod(X,U,mode)
%TMPROD mode-n tensor-matrix product.
sX=size(X); sU=size(U);
sX(mode)=sU(1);
X=mat2tens(U*tens2mat(X,mode),mode,sX);
end