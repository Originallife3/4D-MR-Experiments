function [Mk] = lps_2D(A,At,Mj,b,opts);
opts.nite=50;
opts.lambda_L=0.1;opts.lambda_S=0.01;
[nx,ny,nt,n4]=size(Mj);Mk=Mj;Lk=zeros(nx,ny,nt,n4);Sk=Lk;
ite=0;
fprintf('\n ********** L+S Matrix reconstruction **********\n')

while(1),
ite=ite+1;
Mpre=Mk;
for j=1:n4
        M=reshape(Mpre(:,:,:,j),[nx*ny,nt]);
        Lpre=M;
        S=zeros(nx*ny,nt);
 % low-rank update
            [Ut,St,Vt]=svd(M-S,0);
            St=diag(SoftThresh(diag(St),St(1)*opts.lambda_L));
            L=Ut*St*Vt';
             % sparse update
            S=reshape(opts.T'*(SoftThresh(opts.T*reshape(M-L,[nx,ny,nt]),opts.lambda_S)),[nx*ny,nt]);
            % data consistency
            Lk(:,:,:,j)=reshape(L,nx,ny,nt);
            Sk(:,:,:,j)=reshape(S,nx,ny,nt);                 
end  
    resk=A(Lk+Sk)-b;
    Mk=Lk+Sk-At(resk);
    fprintf(' ite: %d ,update: %f3\n', ite,norm(Mk(:)-Mpre(:))/norm(Mpre(:)));                
    if (ite > opts.nite) || (norm(Mk(:)-Mpre(:))<opts.tol*norm(Mpre(:))), break;end
end
end

% soft-thresholding function
function y=SoftThresh(x,p)
y=(abs(x)-p).*x./abs(x).*(abs(x)>p);
y(isnan(y))=0;
end    