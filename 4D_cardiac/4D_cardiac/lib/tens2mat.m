function X = tens2mat(X,mode)
%TENS2MAT Matricize a tensor.
X=reshape(permute(X,[mode 1:mode-1 mode+1:ndims(X)]),size(X,mode),[]);
end