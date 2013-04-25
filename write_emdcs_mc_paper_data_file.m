function write_emdcs_mc_paper_data_file(filename, Mvec, recovery_emdcs_cosamp, ...
        recovery_emdcs_iht, recovery_cosamp, recovery_iht)
    file = fopen(filename, 'w');
    fprintf(file, 'm recovery_emdcs_cosamp recovery_emdcs_iht recovery_cosamp recovery_iht\n');
    for ii = 1:numel(Mvec)
        fprintf(file, '%d %f %f %f %f\n', Mvec(ii), recovery_emdcs_cosamp(ii), ...
            recovery_emdcs_iht(ii), recovery_cosamp(ii), recovery_iht(ii));
    end
    fclose(file);
end

