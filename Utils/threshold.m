function x = threshold(x, k)
    [~, indices] = sort(abs(x), 'descend');
    x(indices(k+1 : end)) = 0;
end