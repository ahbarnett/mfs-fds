load('scat_k300_Nsds.mat'); 
imagesc(Ns, ds, log10(err')); colorbar; caxis([-10 0]); colorbar; set(gca, 'FontSize', 16); axis xy;
title('err vs N,d when k = 300'); 

figure; 
load('scat_k30.mat'); 
semilogy(Ns, err); set(gca, 'FontSize', 16); 
title('err vs N when k = 30'); 
figure; 
plot(Ns, ts'); set(gca, 'FontSize', 16); 
title('t(s) vs N when k = 30'); 

figure; 
load('scat_k300.mat'); 
semilogy(Ns, err); set(gca, 'FontSize', 16); 
title('err vs N when k = 300'); 
figure; 
plot(Ns, ts'); set(gca, 'FontSize', 16); 
title('t(s) vs N when k = 300'); 

figure; 
load('scat_k3000.mat');
semilogy(Ns, err); set(gca, 'FontSize', 16); 
title('err vs N when k = 3000'); 
figure; 
plot(Ns, ts'); set(gca, 'FontSize', 16); 
title('t(s) vs N when k = 3000'); 

figure; 
load('scat_N1e5_ks.mat'); 
semilogy(ks, err); set(gca, 'FontSize', 16); 
title('err vs k when N = 1e5'); 
figure; 
plot(ks, ts'); set(gca, 'FontSize', 16); 
title('t(s) vs N when N = 1e5'); 