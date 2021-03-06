% chin jan 23 2013

% top script for montecarlo sim comparing EMDCS and vanilla

clear
close all
clc

addpath emd_flow/
addpath Utils/

%%%%% Load signal

n = 100; w = 40;
k = 10; B = 120;
tmp = load('sigsbee_part_small1.mat');
X = tmp.tmp;
norm_X_fro = norm(X, 'fro');

N = n*w; K = k*w;
E = 1;
Mvec = 100:10:110;
err_cosamp = zeros(E,length(Mvec));
err_emdcs_cosamp = zeros(E,length(Mvec));
err_iht = zeros(E,length(Mvec));
err_emdcs_iht = zeros(E,length(Mvec));

for ee=1:E
    
    for mm=1:length(Mvec)
        
        disp([ee,mm])
        
        tic;
        
        %%%%% acquire compressive measurements
        M = Mvec(mm);
        
        Phi = 1/sqrt(M)*randn(M,N); % can use random Fourier etc.
        y = Phi*X(:);
        
        %%%%% reconstruct using CoSaMP
        opt.iter = 20;
        Xhat_cosamp = cosamp(y,Phi,K,opt.iter);
        Xhat_cosamp = reshape(Xhat_cosamp,n,w);
        err_cosamp(ee,mm) = norm(Xhat_cosamp - X, 'fro')/norm_X_fro;
        
        %%%%% reconstruct using EMDCS-CoSaMP
        opt.tol = 1e-3; opt.K = K; opt.B = B;
        opt.w = w; opt.k = k;
        opt.verbose = true; opt.pause = 0;
        Xhat_emdcs_cosamp = emdcs(y,Phi,opt);
        Xhat_emdcs_cosamp = reshape(Xhat_emdcs_cosamp,n,w);
        err_emdcs_cosamp(ee,mm) = norm(Xhat_emdcs_cosamp - X, 'fro')/norm_X_fro;
        
        %%%%% reconstruct using IHT
        opt.iter = 50;
        opt.stepsize = 0.5;
        Xhat_iht = iht(y, Phi, K, opt.iter, opt.stepsize);
        Xhat_iht = reshape(Xhat_iht, n, w);
        err_iht(ee,mm) = norm(Xhat_iht - X, 'fro') / norm_X_fro;
        
        %%%%% reconstruct using EMDCS-IHT
        opt.B = B;
        opt.w = w;
        opt.k = k;
        opt.verbose = true;
        Xhat_emdcs_iht = emdcs_iht(y, Phi, opt);
        Xhat_emdcs_iht = reshape(Xhat_emdcs_iht, n, w);
        err_emdcs_iht(ee,mm) = norm(Xhat_emdcs_iht - X, 'fro') / norm_X_fro;
        
        toc
    end

end

% treshold error because we are interested in the probability of recovery
% err_threshold = 0.05;
% recovery_cosamp = error_to_recovery_indicator(err_cosamp, err_threshold);
% recovery_emdcs_cosamp = error_to_recovery_indicator(err_emdcs_cosamp, err_threshold);
% recovery_iht = error_to_recovery_indicator(err_iht, err_threshold);
% recovery_emdcs_iht = error_to_recovery_indicator(err_emdcs_iht, err_threshold);

figure(1), clf
hold on
box on
plot(Mvec,mean(recovery_emdcs_cosamp),'-+r','LineWidth',2,'MarkerSize',15);
plot(Mvec,mean(recovery_emdcs_iht),'-sb','LineWidth',2,'MarkerSize',15);
plot(Mvec,mean(recovery_cosamp),'-*g','LineWidth',2,'MarkerSize',15)
plot(Mvec,mean(recovery_iht),'-dm','LineWidth',2,'MarkerSize',15);
axisfortex('','Number of measurements','Mean error')
legend('EMD-CoSaMP', 'EMD-IHT', 'CoSaMP', 'IHT', 'Location','SouthEast')
axis tight
grid on