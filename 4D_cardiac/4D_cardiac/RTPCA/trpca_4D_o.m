function [M] =trpca_4D_o(A,At,x_init,b,opts);
opts.nite=60;
M=x_init;
[nx,ny,nz,nt]=size(M);
S=zeros(nx,ny,nz,nt); 
ite=0;
rho=1.1;
mu=0.02;
max_mu=1e10;
lambda=1/sqrt(nt*max(nx,ny));
fprintf('\n ********** L+S Tensor reconstruction **********\n')

while(1),
    ite=ite+1;
    M0=M;
    % low-rank update
    [L,tnnL,trank] = prox_tnn_o(M-S,1/mu);
     % sparse update
    S = opts.T'*(SoftThresh(opts.T*(M-L),lambda/mu));
    %  data consistency
    resk=A(L+S) - b;
	M=L+S-At(resk);
    mu=min(rho*mu,max_mu);
    
if (ite > opts.nite) || (norm(M(:)-M0(:))<opts.tol*norm(M0(:))), break;end
fprintf(' ite: %d ,update: %f3\n', ite,norm(M(:)-M0(:))/norm(M0(:))); 
end
end

% soft-thresholding function
function y=SoftThresh(x,p)
y=(abs(x)-p).*x./abs(x).*(abs(x)>p);
y(isnan(y))=0;
end    