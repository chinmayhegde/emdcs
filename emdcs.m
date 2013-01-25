% chin jan 23 2013

% EMDCS: function for CS recon using EMD model

function [xhat,xprev] = emdcs(yy,Phi,opt)

tol = opt.tol;
Its = opt.iter;
B = opt.B;
k = opt.k;
w = opt.w;
yy = yy(:); % 
[~,N] = size(Phi);

xprev = zeros(N,Its);
kk=1; 
maxiter= 1000;
xhat = zeros(N,1);

while le(kk,Its),
    
    %----- EMD Model approximation---%
    r = yy - Phi*xhat;
    proxy = Phi'*r;
    
    tmp = reshape(proxy,[],w);
    mags = tmp.^2;
    supp = emd_flow(mags,2*k,2*B,opt.verbose);
    supp = double(supp(:));
    tt= union(find(ne(xhat,0)),find(ne(supp,0)));
    
    %------Least Squares------%
    coeff = cgsolve(Phi(:,tt)'*Phi(:,tt), Phi(:,tt)'*yy,...
                                        tol,maxiter, 0);
  
    bb2= zeros(N,1);
    bb2(tt)= coeff;
    
    %---Pruning----%
    kk = kk+1;   
    tmp = reshape(bb2,[],w);
    mags = tmp.^2;
    supp = emd_flow(mags,k,B,opt.verbose);
    xhat= tmp.*double(supp);
    xhat = xhat(:);
    xprev(:,kk) = xhat; % current signal estimate
    if (norm(xprev(:,kk)-xprev(:,kk-1)) < 1e-3*norm(xprev(:,kk)))
       break;
    end
    
    if opt.verbose
        disp(kk)
    end
end
xprev(:,kk+1:end)=[];