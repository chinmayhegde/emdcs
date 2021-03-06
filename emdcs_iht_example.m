% ludwig jan 25 2013

% top script for calling emdcs_iht

clear
close all
clc

addpath emd_flow/
addpath Utils/

%%%%% Construct synthetic signal

n = 100; w = 10;
k = 2; B = 20;
X = zeros(n,w);
X(1:(2*w),:) = [eye(w); eye(w)];

%%%%% acquire compressive measurements
N = n*w; K = k*w;
M = 80;

Phi = 1/sqrt(M)*randn(M,N); % can use random Fourier etc.
y = Phi*X(:);

%%%%% reconstruct using IHT
opt.iter = 50;
opt.stepsize = 0.5;
Xhat_iht = iht(y, Phi, K, opt.iter, opt.stepsize);
Xhat_iht = reshape(Xhat_iht, n, w);
err_cosamp = norm(Xhat_iht - X, 'fro') / norm(X);

%%%%% reconstruct using EMDCS IHT
opt.B = B;
opt.w = w;
opt.k = k;
opt.verbose = true;
Xhat_emdcs = emdcs_iht(y, Phi, opt);
Xhat_emdcs = reshape(Xhat_emdcs,n,w);
err_emdcs = norm(Xhat_emdcs - X, 'fro')/norm(X);

figure(1), clf
subplot(1,3,1), imagesc(X), axisfortex('','Original',''), rmaxis
subplot(1,3,2), imagesc(Xhat_iht), axisfortex('','IHT',''), rmaxis
subplot(1,3,3), imagesc(Xhat_emdcs), axisfortex('','EMD-IHT',''), rmaxis

