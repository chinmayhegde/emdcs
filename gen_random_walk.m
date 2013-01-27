function x = gen_random_walk(w, h, num_walks, p_step, p_up, stepsize)
    x = zeros(h, w);
    for iwalk = 1:num_walks
        cury = randi(h);
        x(cury, 1) = 1;
        for ix = 2:w
            if rand() <= p_step
                if rand() <= p_up
                    cury = max(1, cury - stepsize);
                else
                    cury = min(h, cury + stepsize);
                end
            end
            x(cury, ix) = 1;
        end
    end
end

