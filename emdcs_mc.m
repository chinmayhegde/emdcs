% chin jan 23 2013

% top script for montecarlo sim comparing EMDCS and vanilla

clear
close all
clc

addpath emd_flow/
addpath Utils/
addpath /Users/Chin/Documents/Chin/Acads/Utils/spgl1-1.7/

%%%%% Construct synthetic signal

n = 100; w = 10;
k = 2; B = 20;
X = zeros(n,w);
X(1:(2*w),:) = [eye(w); eye(w)];

N = n*w; K = k*w;
E = 10;
Mvec = 60:10:150;
err_cosamp = zeros(E,length(Mvec));

for ee=1:E
    
    for mm=1:length(Mvec)
        
        disp([ee,mm])
        
        %%%%% acquire compressive measurements
        M = Mvec(mm);
        
        Phi = 1/sqrt(M)*randn(M,N); % can use random Fourier etc.
        y = Phi*X(:);
        
        %%%%% reconstruct using CoSaMP
        opt.iter = 20;
        Xhat_cosamp = cosamp(y,Phi,K,opt.iter);
        Xhat_cosamp = reshape(Xhat_cosamp,n,w);
        err_cosamp(ee,mm) = norm(Xhat_cosamp - X, 'fro')/norm(X);
        
        %%%%% reconstruct using EMDCS
        opt.tol = 1e-3; opt.K = K; opt.B = B;
        opt.w = w; opt.k = k;
        opt.verbose = 0; opt.pause = 0;
        Xhat_emdcs = emdcs(y,Phi,opt);
        Xhat_emdcs = reshape(Xhat_emdcs,n,w);
        err_emdcs(ee,mm) = norm(Xhat_emdcs - X, 'fro')/norm(X);
        
    end

end

figure(1), clf
%subplot(1,3,1), imagesc(X), axisfortex('','Original',''), rmaxis
%subplot(1,3,2), imagesc(Xhat_cosamp), axisfortex('','CoSaMP',''), rmaxis
%subplot(1,3,3), imagesc(Xhat_emdcs), axisfortex('','EMD-CS',''), rmaxis
hold on
box on
plot(Mvec,mean(err_cosamp),'r-.*','LineWidth',2)
plot(Mvec,mean(err_emdcs),'-d','LineWidth',2);
axisfortex('','Number of measurements','Normalized recovery error')