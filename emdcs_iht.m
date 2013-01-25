function xhat = emdcs_iht(y, phi, opts)
    xhat = zeros(size(phi, 2), 1);
    
    k = opts.k;
    B = opts.B;
    width = opts.w;
    num_iter = opts.iter;
    stepsize = opts.stepsize;
    if ~isfield(opts, 'verbose')
        opts.verbose = false;
    end
    
    for i = 1:num_iter
        prevx = xhat;
        a = xhat + stepsize * (phi' * (y - phi * xhat));
        
        % find support with EMD flow
        tmp = reshape(a, [], width);
        tmp = tmp.^2;
        supp = emd_flow(tmp, k, B, opts.verbose);
        supp = reshape(supp, [], 1);
        
        xhat = a .* double(supp);
        
        if (norm(prevx - xhat) < 1e-3 * norm(xhat))
            break;
        end
        
        if opts.verbose
            disp(i)
        end
    end
end

