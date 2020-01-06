function [ssim] =ssim_sel(imagery1, imagery2)

[m, n, k,l] = size(imagery1);
[mm, nn, kk,ll] = size(imagery2);
m = min(m, mm);
n = min(n, nn);
k = min(k, kk);
imagery1 = imagery1(1:m, 1:n, 1:k);
imagery2 = imagery2(1:m, 1:n, 1:k);
ssim =0;

        for i = 1:k
           ssim= ssim + ssim_index(imagery1(:, :, i), imagery2(:, :, i));
        end
        ssim= ssim/k;
end
