function run_csn_parallel(csn_weighted, IO_path, n_cores)
    
    addpath(genpath('External_tools/CSN/'));
    
    data = load(fullfile(IO_path, 'gene_expr.mat'));
    gene_expr = data.gene_expr;
    sample_ids = data.sample_ids;
    gene_ids = data.gene_ids;
    num_samples = size(gene_expr, 2);
    results_dir = fullfile(IO_path, 'CSN_results_tmp');
    
    
    p = gcp('nocreate');
    if isempty(p)
        parpool('local', n_cores);
    elseif p.NumWorkers ~= n_cores
        delete(p);
        parpool('local', n_cores);
    end
    
    
    parfor idx = 1:num_samples
        sample_select = sample_ids{idx};
        out_file = fullfile(results_dir, ['CSN_result_', sample_select, '.parquet']);
        log_file = fullfile(results_dir, ['log_', sample_select, '.txt']);
    
        try
            csn_output = csnet(gene_expr, idx, 0.05, 0.1, csn_weighted);

            if csn_weighted
                csn_results = csn_output{idx};
            else
                csn_results = double(full(csn_output{idx}));
            end

            save_helper(out_file, csn_results, gene_ids, sample_select);
    
            fid = fopen(log_file, 'w');
            fprintf(fid, 'Success: %s completed at %s', sample_select, datestr(now));
            fclose(fid);
    
        catch ME
            fid = fopen(log_file, 'w');
            fprintf(fid, 'Error on %s: %s', sample_select, ME.message);
            fclose(fid);
        end
    end
end


function save_helper(fname, csn_results, gene_ids, sample_name)
    T = array2table(csn_results, 'VariableNames', gene_ids, 'RowNames', gene_ids);
    sample_id = repmat(string(sample_name), size(T, 1), 1);
    T = addvars(T, sample_id, 'Before', 1);
    parquetwrite(fname, T);
end