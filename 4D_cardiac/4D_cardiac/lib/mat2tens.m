function X = mat2tens(X,mode,size_tens)
%MAT2TENS Tensorize a matrix.
X=reshape(X,size_tens([mode 1:mode-1 mode+1:end])); 
X=permute(X,[2:mode 1 mode+1:numel(size_tens)]);
end