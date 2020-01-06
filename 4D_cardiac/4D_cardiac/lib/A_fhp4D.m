function A = A_fhp4D(z, S)
    z=double(z); S = single(S);
    [n1,n2,n3,n4] = size(z);
    A = zeros(n1,n2,n3,n4);
    p=1/sqrt(n1*n2)*fft2(z) ;
    A(S) = p(S);
end
