function [U,cost,opts] = minSN_4D(A,At,x_init,b,T3D,opts)
%%% optimisation for the reconstruction


    U=x_init; [m n d t] = size(U); % initial estimate
    Lam1 = zeros(m,n,d,t); % parametre for sparsity regularisation



    o=0;cost=[];
    for out = single(1:opts.outer_iter) % outer iteration
        o=o+1
        for in = single(1:opts.inner_iter) % inner iteration       
            [U1,S1,temp]=svd(reshape(U+Lam1/opts.beta1,m,d*n*t),'econ');
            [U2,S2,temp]=svd(reshape(permute(U+Lam1/opts.beta1,[2 3 4 1]),n,m*d*t),'econ');
            [U3,S3,temp]=svd(reshape(permute(U+Lam1/opts.beta1,[3 4 1 2]),d,n*m*t),'econ');
            [U4,S4,temp]=svd(reshape(permute(U+Lam1/opts.beta1,[4 1 2 3]),t,n*m*d),'econ');
            s=tmprod(tmprod(tmprod(tmprod(U+Lam1/opts.beta1,U1',1),U2',2),U3',3),U4',4);
            thres=(1/opts.beta1).*(abs(s).^(opts.p-1));
            s(abs(s)<=abs(thres))=0;
            idx=find(abs(s)>abs(thres));
            s(idx)=s(idx)-abs(thres(idx))./abs(s(idx)).*s(idx);
            Lambda=tmprod(tmprod(tmprod(tmprod(s,U1,1),U2,2),U3,3),U4,4);

            % ----------------
            %  Solve for the recon: Conjugate gradient update
            % ----------------


            [U,earray1] = CG_solver_4D(b,A, At,Lambda,Lam1,opts, U, 1e-7,10);
            %% plot a frame of the recon while it iterates
%             figure(3); 
%             imagesc(abs(double((U(:,:,4,2))))); title('A frame of the reconstruction'); colormap(gray);

            %% cost calculations
            e = A(U) - b;     
            [U1,S1,temp]=svd(reshape(U,m,d*n*t),'econ');
            [U2,S2,temp]=svd(reshape(permute(U,[2 3 4 1]),n,m*d*t),'econ');
            [U3,S3,temp]=svd(reshape(permute(U,[3 4 1 2]),d,n*m*t),'econ');
            [U4,S4,temp]=svd(reshape(permute(U,[4 1 2 3]),t,n*m*d),'econ');
            sigmaq=tmprod(tmprod(tmprod(tmprod(U,U1',1),U2',2),U3',3),U4',4);


            cost = [cost, sum(abs(e(:)).^2)  +  sum(abs(sigmaq(:)).^(opts.p)./(opts.p))*opts.mu1 ];
%             figure(30);plot(double(cost)); hold on; pause(0.1); hold off; title('Cost');



            if in>1
                if abs(cost(end) - cost(end-1))/abs(cost(end-1)) < 1e-3
                    break;
                end
            end


            %%        Update rules for the Lagrange multipliers
            Lam1 = Lam1 - 1.618*opts.beta1*(Lambda - U);
        end     

        %% increment beta 1
        opts.beta1=opts.beta1*opts.beta1rate;

    end
end 