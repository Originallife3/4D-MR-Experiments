function [psnr_r,ssim_r,fsim_r,mse_r] = MQAI(A,Recon)
psnr_r=psnr(A,Recon);
ssim_r=ssim_self(A,Recon);
end