% ludwig jan 27 2013

% top script for calling emdcs

clear
close all
clc

addpath emd_flow/
addpath Utils/

%%%%% Construct random walk

n = 100; w = 10;
p_step = 0.8;
k = 3; B = k * w * p_step;
X = gen_random_walk(w, n, k, p_step, 0.5, 1);

%%%%% acquire compressive measurements
N = n*w; K = k*w;
M = 100;

Phi = 1/sqrt(M)*randn(M,N); % can use random Fourier etc.
y = Phi*X(:);

%%%%% reconstruct using CoSaMP
opt.iter = 50;
Xhat_cosamp = cosamp(y,Phi,K,opt.iter);
Xhat_cosamp = reshape(Xhat_cosamp,n,w);
err_cosamp = norm(Xhat_cosamp - X, 'fro')/norm(X);

%%%%% reconstruct using EMDCS
opt.tol = 1e-3; opt.K = K; opt.B = B;
opt.w = w; opt.k = k;
opt.verbose = true; opt.pause = 0;
Xhat_emdcs = emdcs(y,Phi,opt);
Xhat_emdcs = reshape(Xhat_emdcs,n,w);
err_emdcs = norm(Xhat_emdcs - X, 'fro')/norm(X);

figure(1), clf
colormap gray
subplot(1,3,1), imagesc(X), caxis([-1 1]), axisfortex('','Original',''), rmaxis
subplot(1,3,2), imagesc(Xhat_cosamp), caxis([-1 1]), axisfortex('','CoSaMP',''), rmaxis
subplot(1,3,3), imagesc(Xhat_emdcs), caxis([-1 1]), axisfortex('','EMD-CS',''), rmaxis

