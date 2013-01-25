function xhat = iht(y, phi, k, num_iter, stepsize)
    xhat = zeros(size(phi, 2), 1);

    for i = 1:num_iter
        prevx = xhat;
        a = xhat + stepsize * (phi' * (y - phi * xhat));
        xhat = threshold(a, k);
        
        if (norm(prevx - xhat) < 1e-3 * norm(xhat))
            break;
        end
    end
end

