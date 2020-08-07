clear all;
close all;
warning off;
addpath .\genetic;
addpath .\functions;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PARAMETERS defined by users.
%  the input/output path should be pre-defined.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
project_name = "example";
path_of_occupancy = '.\data\occupancy_data.xlsx';
path_of_gene_expression = '.\data\gene_expression.xlsx';
path_of_output = '.\output\';
name_of_output = strcat(project_name, "-", datestr(datetime("today")), "-", randstring(8));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[binding_data, binding_txt] = xlsread(path_of_occupancy);
binding_txt = upper(binding_txt(2:end, :));

exp_datafiles = path_of_gene_expression;

time_intervales = [1,1,2,2,4,4,2];
KO_sheets = {'TF_exp'};
    
for i = 1: 1
    [orig_gene_exp_data{i}, orig_gene_txt{i}] = xlsread(exp_datafiles, 'Gene_exp');
    a = isnan(orig_gene_exp_data{i}(:,1));
    if sum(a) > 0.5*size(orig_gene_exp_data{i}, 2)
        orig_gene_exp_data{i} = orig_gene_exp_data{i}(:, 2:end);
    end
    
    for j = 1: length(KO_sheets)    
        [orig_TF_exp_data{i}{j}, orig_TF_txt{i}{j}] = xlsread(exp_datafiles, KO_sheets{j});
    end
    
    orig_genes{i} = orig_gene_txt{i}(2:end, 1);
    time_points(i) = size(orig_gene_exp_data{i}, 2);
    
    sample_hours(i) = time_points(i) * time_intervales(i);
end

TGs = orig_genes{1};

TFs = orig_TF_txt{1}{1}(2:end,1);

num_TGs = length(TGs);
num_TFs = length(TFs);

for i = 1: 1% experiments
    TG_exp_data_all{i} = orig_gene_exp_data{i}(:, :);
    for j = 1: length(KO_sheets)
        TF_exp_data_all{i}{j} = orig_TF_exp_data{i}{j};
    end
end

code = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'};
operator_set = {'|', '&', '>'};
delay = 1;

replicates = [1, 1, 1, 2, 3, 3, 1];

TFs = upper(TFs);
TGs = upper(TGs);
binding_txt(:, 2) = upper(binding_txt(:, 2));

ind_experiment = 1;
for ind_KO = 1: length(KO_sheets)
    TG_gene_data = TG_exp_data_all{ind_experiment};
    TF_exp_data = TF_exp_data_all{ind_experiment}{ind_KO};
    num_TG = length(TGs);
    num_TF = length(TFs);
    
    TF_exp_data = TF_exp_data / max(max(TF_exp_data));
    TF_exp_data_1 = reshape(TF_exp_data, size(TF_exp_data, 1) * size(TF_exp_data, 2), 1);
    
    D = quantile(TF_exp_data_1, [0.5, 0.95]);
    Gene_expr_median = D(1);
    Gene_expr_max95 = D(2);
    k1 = - log(0.5) / Gene_expr_median;
    k2 = - log(0.05) / Gene_expr_max95;
    k = min(k1, k2);
    
    Gene_TF_score = zeros(num_TG, num_TF);
    all_score = [];
    
    disp('start calculate Gene TF score');
    
    for TG_i = 1: num_TG 
        Target = TGs{TG_i};
        TF_bind_on_Target_i = strcmp(Target, binding_txt(:, 2));
        
        TF_ind_bind_on_Target = find(TF_bind_on_Target_i);
        TF_bind_on_Target = binding_txt(TF_ind_bind_on_Target, 1);
        
        for j = 1 : length(TF_ind_bind_on_Target)
            ind = strcmp(upper(TF_bind_on_Target{j}), TFs);
            ind1 = find(ind);
            Gene_TF_score(TG_i, ind1) = Gene_TF_score(TG_i, ind1) + binding_data(TF_ind_bind_on_Target(j), 1);
        end
    end
    
    disp('Finish calculate Gene TF score')
    
    tmp = reshape(Gene_TF_score, 1, num_TG*num_TF);
    tmp_i = find(tmp >0);
    Gene_TF_score = Gene_TF_score / std(tmp(tmp_i));
    
    TF_score_mean = sum(sum(Gene_TF_score)) / (num_TG * num_TF);
    
    disp('start calculate TF occupancy')
    Gene_TF_binding_strength = Gene_TF_score;
    occu_mtx = [];
    for TG_i = 1: num_TG
        for j = 1: num_TF
            occupancy{TG_i, j} = 1 - exp(- Gene_TF_binding_strength(TG_i,j) * TF_exp_data(j,:));
            occu_mtx = [occu_mtx; occupancy{TG_i,j}];
        end
    end
    
    disp('Finish calculate TF occupancy');
    disp('start calculate the logics');
    
    X = cell(num_TG,1);
    beta_low = 2;
    beta_high = 5;
    mZdif = zeros(num_TG, size(TG_gene_data,2)-1);
    
    for t = 2: size(TG_gene_data, 2)
        mZdif(:,t-1) = (TG_gene_data(:,t) - TG_gene_data(:,t - delay));% / time_intervals(t-1);
    end
   
    for TG_i = 1: num_TG
        X{TG_i} = occu_mtx((TG_i-1)*num_TF + 1: TG_i*num_TF, 1 : size(TG_gene_data,2)-1);
        
        for j = 1:num_TF-1
            for k = j+1: num_TF
                X{TG_i} = [X{TG_i}; occu_mtx((TG_i-1) * num_TF + j, 1 :size(TG_gene_data,2)-1) .* occu_mtx((TG_i-1) * num_TF + k, 1: size(TG_gene_data,2)-1)]; % ab
            end
        end
    end
    
    [logic_gate, comp_var, map] = fun_2R_logic_generator(num_TF, 1);
    map = abs(map);
    lgn = size(map, 1);
    vn = size(map, 2);
    
    logic_num = num_TF + sum((1:(num_TF-1)) * 4);
    
    p0 = zeros(num_TG, logic_num);
    
    p = zeros(num_TG, lgn);
    
    predicted = cell(1, num_TG);
    
    for TG_i = 1: num_TG
        y = mZdif(TG_i, :);
        xarr = (X{TG_i}(:, 1:end-delay))';
        yarr = (y(delay+1:end))';
        
        yarr = yarr - mean(yarr);
        
        [n, ~] = size(xarr);
        xarr = xarr - ones(n, 1) * mean(xarr);
        
        xarr_std = sqrt(sum(xarr.^2));
        
        for j = 1: size(xarr, 2)
            if xarr_std(j) == 0
                xarr(:, j) = xarr(:, j);
            else
                xarr(:, j) = xarr(:, j) / xarr_std(j);
            end
        end
        beta = lars(xarr, yarr);
        
        if beta_high > size(beta, 1)
            beta_high = size(beta, 1);
        end
        
        if sum(sum(xarr_std)) == 0
            beta = zeros(size(beta,1), size(beta,2));
            p0(TG_i,1:logic_num) = zeros(1, logic_num);
        else
            for row = beta_low: beta_high
                if sum(isnan(beta(row,:))) >= 1
                    continue;
                end
                
                sbeta = beta(row, :);
                sbeta_std = std(sbeta(find(sbeta)));
                
                if sum(abs(sbeta)) == 0
                    sigm_sbeta = sbeta;
                elseif sbeta_std == 0
                    sigm_sbeta = sbeta / max(abs(sbeta));
                else
                    sigm_sbeta = (1- exp(-abs(sbeta/sbeta_std)))./(1+ exp(-abs(sbeta/sbeta_std)));
                end
                
                thresh = 0.9;
                p_elmt = zeros(lgn, vn);
                
                for lg = 1: lgn
                    for v = 1: vn
                        if map(lg, v) == 0
                            p_elmt(lg, v) = 1 - abs(sigm_sbeta(v));
                        else
                            p_elmt(lg, v) = abs(sigm_sbeta(v));
                        end
                    end
                end
                
                for lg = 1: lgn
                    p(TG_i, lg) = 1;
                    features = length(find(map(lg,:)));
                    for v = 1: vn
                        if p_elmt(lg, v) >= thresh
                            p(TG_i, lg) = p(TG_i, lg) * p_elmt(lg, v);
                        else
                            p(TG_i, lg) = p(TG_i, lg) * (sqrt(features))/(sqrt(features) + 1);
                        end
                    end
                end
                
                if sum(p(TG_i,:)) == 0
                else
                    p1 = p(TG_i, :)/sum(p(TG_i, :));
                end
                [~, b(TG_i)] = max(p1);
                
                p0(TG_i, :) = p0(TG_i, :) + p1;
                
            end
        end
    end
    
    p0max = max(p0, [], 2);
    mmax = [];
    
    for TG_i = 1: num_TG
        if p0max(TG_i) == 0
        else
            p0(TG_i,:) = p0(TG_i,:) / p0max(TG_i);
        end
    end
    
    all_p = reshape(p0, 1, size(p0,1) * size(p0,2));
    [sort_p, sort_ind] = sort(all_p, 'descend');
    
    
    for TG_i = 1: num_TG
        target_gene(TG_i) = TG_i;
        
        cutting_thresh = 0.8;
        
        logic_ind = find(p0(TG_i,:) > cutting_thresh );
        all_logic_string{TG_i}='';
        
        for j = 1: length(logic_ind)
            logic_code = logic_gate{logic_ind(j)};
            
            if j > 1
                all_logic_string{TG_i} = strcat(all_logic_string{TG_i}, '.OR.');
            end
            
            if length(logic_code) == 1
                regulator_code = logic_code;
                [c, ia, ib] = intersect(regulator_code, code);
                all_logic_string{TG_i} = strcat(all_logic_string{TG_i}, '(', TFs(ib), ')');
            else
                regulator1_code = logic_code(1);
                regulator2_code = logic_code(3);
                operator_code = logic_code(2);
                
                [c, ia1, ib1] = intersect(regulator1_code, code);
                [c, ia2, ib2] = intersect(regulator2_code, code);
                all_logic_string{TG_i} = strcat(all_logic_string{TG_i}, '(', TFs(ib1), '.', operator_code,'.', TFs(ib2), ')');
            end
        end
    end
    
    for TG_i = 1: num_TG
        cutting_thresh = 0.8;
        logic_ind = find(p0(TG_i,:) > cutting_thresh);
        
        logic_label{TG_i} = [];
        tmp_label = 0;
        
        logic_numberCode{TG_i} = [];
        if isempty(logic_ind)
            para{TG_i} = [];
        else
            total_logics(TG_i) = length(logic_ind);
            regulator_record = [];  
            for j = 1: total_logics(TG_i)
                logic_code = logic_gate{logic_ind(j)};
                if length(logic_code) == 1                
                    regulator_code = logic_code; 
                    [c, ia, ib1] = intersect(regulator_code, code); 
                    logic_numberCode{TG_i} = [logic_numberCode{TG_i}; [ib1 0 0 p0(TG_i,logic_ind(j))]];
                    number_of_regulators{TG_i}(j) = 1;   
                    ib2 = 0;
                    regulator_record_i = [ib1, ib2];
                    
                else 
                    regulator1_code = logic_code(1);
                    regulator2_code = logic_code(3);
                    operator_code = logic_code(2);
                    
                    [c, ia1, ib1] = intersect(regulator1_code, code);
                    [c, ia2, ib2] = intersect(regulator2_code, code);
                    [c, ia3, operator_number] = intersect(operator_code, operator_set);
                    
                    number_of_regulators{TG_i}(j) = 2;
                    logic_numberCode{TG_i} = [logic_numberCode{TG_i}; [ib1, ib2, operator_number, p0(TG_i,logic_ind(j))]];
                    if operator_number == 2
                        regulator_record_i = [ib1 0];
                    elseif operator_number == 1
                        regulator_record_i = [ib1, ib2];
                    else
                        regulator_record_i = [ib1, 0];
                    end
                end
                regulator_record = [regulator_record, regulator_record_i];
            end
        end
    end
    
    for i = 1: length(all_logic_string)
        disp(all_logic_string{i})
        result_logic{i,1} = TGs{i};
        result_logic{i,2} = char(all_logic_string{i});
    end
    disp('Finish calculate logics')
    all_p0{ind_experiment}{ind_KO} = p0;
end

weights = time_points .* replicates;
weights = weights / sum(weights);

p1 = cell(1, length(KO_sheets));

for ind_KO = 1: length(KO_sheets)
    
    tmp_p1 =  all_p0{ind_experiment}{ind_KO};
    
    tmp_p1 = tmp_p1 .* (tmp_p1 >= 0.5);
    
    p1{ind_KO} = tmp_p1  * weights(ind_experiment);
    
end

for ind_KO = 1: length(KO_sheets)
    for TG_i = 1: num_TG
        cutting_thresh = 0.8;
        logic_ind = find(p1{ind_KO}(TG_i,:) > cutting_thresh);
        logic_label{TG_i} = [];
        tmp_label = 0;
        logic_numberCode{TG_i} = [];
        
        if isempty(logic_ind)
            para{TG_i} = [];
        else
            total_logics(TG_i) = length(logic_ind);
            regulator_record = [];
            
            for j = 1: total_logics(TG_i)
                logic_code = logic_gate{logic_ind(j)};
                if length(logic_code) == 1
                    
                    regulator_code = logic_code;
                    [c, ia, ib1] = intersect(regulator_code, code);
                    logic_numberCode{TG_i} = [logic_numberCode{TG_i}; [ib1 0 0 p0(TG_i,logic_ind(j))]];
                    
                    number_of_regulators{TG_i}(j) = 1;
                    
                    ib2 = 0;
                    regulator_record_i = [ib1, ib2];
                    
                else
                    
                    regulator1_code = logic_code(1);
                    regulator2_code = logic_code(3);
                    operator_code = logic_code(2);
                    
                    [c, ia1, ib1] = intersect(regulator1_code, code);
                    [c, ia2, ib2] = intersect(regulator2_code, code);
                    [c, ia3, operator_number] = intersect(operator_code, operator_set);
                    
                    number_of_regulators{TG_i}(j) = 2;
                    logic_numberCode{TG_i} = [logic_numberCode{TG_i}; [ib1, ib2, operator_number, p0(TG_i,logic_ind(j))]];
                    if operator_number == 2
                        regulator_record_i = [ib1 0];
                    elseif operator_number == 1
                        regulator_record_i = [ib1, ib2];
                    else
                        regulator_record_i = [ib1, 0];
                    end
                end
                regulator_record = [regulator_record, regulator_record_i];
            end
        end
    end
    
    for i = 1: length(all_logic_string)
        disp(all_logic_string{i})
        result_logic{i,1} = TGs{i};
        result_logic{i,2} = char(all_logic_string{i});
    end
    save_logic_file = strcat(path_of_output, name_of_output);
    xlswrite(save_logic_file, result_logic);
    
end

