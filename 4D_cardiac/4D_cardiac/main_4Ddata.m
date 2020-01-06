clear; clc;
load data4D.mat;
%%%CS-HOSVD on 4D data
% x=(double(abs(data4D(:,:,:,1:10))));
x=(double(abs(data4D(:,:,:,1:10))));
x=x./max(x(:));
[n1,n2,n3,n4]=size(x);
line =35; 
% generate sampling mask/define how many radial lines to sample in k-space
for i=1:size(x,4)
    [mask(:,:,:,i)] = strucrand(n1,n2,n3,line);
     mask(:,:,:,i) = fftshift(fftshift(mask(:,:,:,i),1),2);
end
clear i;
S=find(mask~=0);
%% Define the forward and backward Fourier operators A and Atranspose (At)
A = @(z)A_fhp4D(z,S);
At=@(z)At_fhp4D(z,S,n1,n2,n3,n4);
%% First guess, direct IFFT
b=A(x); 
x_init = At(b);
%% tensor robust priciple component analysis
%% based on TENSOR_L+S
opts.tol=1e-4;
opts.T=TempFFT(4);

tic;
[Recon1] = trpca_4D( A,At,x_init,b,opts);
time1=toc;
R1=abs(Recon1);
%% original scheme
tic;
[Recon2] = trpca_4D_o( A,At,x_init,b,opts);
time2=toc;
R2=abs(Recon2);
%% RPCA
 tic;
 [Recon3] = lps_2D( A,At,x_init,b,opts);
 time3=toc;
 R3=abs(Recon3);
%% HOSVD
tic;
 opts.mu1 = 3e-1; % Regularization parameter for sparsity
 opts.p=0.1; % The value of p in sparsity regularization
 opts.beta1=1;% The penelty for sparsity
 opts.beta1rate = 10;% The increment for sparsity
 opts.outer_iter =20; % # of outer iterations
 opts.inner_iter = 100; % # of inner iterations
[Recon4,cost,opts] = minSN_4D(A,At,x_init,b,1,opts);
time4=toc;
R4=abs(Recon4); 
%% imshow
figure(1);
subplot(1,5,5);imagesc(abs(R1(:,:,8,4)));axis off;
subplot(1,5,4);imagesc(abs(R2(:,:,8,4)));axis off;
subplot(1,5,3);imagesc(abs(R3(:,:,8,4)));axis off;
subplot(1,5,2);imagesc(abs(R4(:,:,8,4)));axis off;
subplot(1,5,1);imagesc(abs(x(:,:,8,4)));axis off;

%% method1
R1=abs(Recon1);[psnr_T,ssim_T] = MQAI(x(:,:,:,4),R1(:,:,:,4));
R2=abs(Recon2);[psnr_o,ssim_o] = MQAI(x(:,:,:,1),R2(:,:,:,4));
R3=abs(Recon3);[psnr_m,ssim_m] = MQAI(x(:,:,:,4),R3(:,:,:,4));
R4=abs(Recon4);[psnr_hosvd,ssim_hosvd] = MQAI(x(:,:,:,4),R4(:,:,:,4));
% %% method2
% [n1,n2,n3,n4]=size(R1); psnrT=0;ssimT=0;
% psnr_T=zeros(1,n4);ssim_T=zeros(1,n4);
% for i=1:n4
%     [psnr_t,ssim_t] = MQAI(x(:,:,:,i),R1(:,:,:,i));
%     psnr_T(1,i)=psnr_t;
%     ssim_T(1,i)=ssim_t;
% end
% for j=1:n4
%     psnrT=psnrT+psnr_T(1,i);
%     ssimT=ssimT+ssim_T(1,i);
% end
% psnrT=psnrT/n4;
% ssimT=ssimT/n4; 
% % figure(2);
% cost_t=x-R1;cost_o=x-R2; cost_m=x-R3; cost_h=x-R4;
% imagesc(abs(cost_t(:,:,8,4)));axis off;
% imagesc(abs(cost_o(:,:,8,4)));axis off;
% imagesc(abs(cost_m(:,:,8,4)));axis off;
% imagesc(abs(cost_h(:,:,8,4)));axis off;
% imagesc(abs(x(:,:,8,4)));

