function [M] =trpca_4D(A,At,x_init,b,opts);
opts.nite=60;
M=x_init;
[nx,ny,nz,nt]=size(M);nxy=max(nx,ny);
S=zeros(nx,ny,nz,nt);L=S;
ite=0;
mu=4;
% rho=0.02;alpha=1.1;
% lambda1=1/sqrt(max(nxy,nt));
lambda1=0.2;
lambda2=2/(sqrt(nt*max(nx,ny)));
fprintf('\n ********** L+S Tensor reconstruction **********\n')

while(1)
    ite=ite+1;
    M0=M;
   [L,tnnL] = prox_tnn(M-S,lambda1);
   S = opts.T'*(SoftThresh(opts.T*(M-L),lambda2/mu));
   % data consistency
    resk=A(L+S) - b;
	M=L+S-At(resk);
% resident
%   Lpre=L;
%   tmp2=opts.T*S;
%cost2= norm(resk(:),2)^2+lambda1*tnnL+lambda2*norm(tmp2(:),1);
 if (ite > opts.nite) || (norm(M(:)-M0(:))<opts.tol*norm(M0(:))), break;end
 fprintf(' ite: %d ,update: %f3\n', ite,norm(M(:)-M0(:))/norm(M0(:))); 
end
end
% soft-thresholding function
function y=SoftThresh(x,p)
y=(abs(x)-p).*x./abs(x).*(abs(x)>p);
y(isnan(y))=0;
end    