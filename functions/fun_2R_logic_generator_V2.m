%% construct the logic gate and composite variables
%clear all
%close all

function [logic_gate, comp_var, fact_map, logic_mtx] = fun_2R_logic_generator_V2(TF_num, taylor_degree)
% taylor_degree = 1, means: a, b, ab
% taylor_degree = 2, means: a, b, ab, a2, b2
% taylor_degree = 3, means: a, b, ab, a2, b2, a2b, ab2
% taylor_degree = 4, means: a, b, ab, a2, b2, a2b, ab2, a2b2

%% construct logic_gate
code = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'};

for i = 1: TF_num
    logic_mtx(i, :) = [i 0 0];
end

TF_code= code(1:TF_num);
%logic_operator = {'&','|','>'}; % '&' is 1; '|' is 2; '>' is 3

logic_gate = TF_code;
for i = 1: TF_num    
    for j = i + 1: TF_num
        logic_gate(length(logic_gate) + 1) = strcat(TF_code(i), '&', TF_code(j));
        logic_gate(length(logic_gate) + 1) = strcat(TF_code(i), '|', TF_code(j));
        logic_gate(length(logic_gate) + 1) = strcat(TF_code(i), '>', TF_code(j));
        logic_gate(length(logic_gate) + 1) = strcat(TF_code(j), '>', TF_code(i));
        
        logic_mtx = [logic_mtx; [i j 1]];
        logic_mtx = [logic_mtx; [i j 2]];
        logic_mtx = [logic_mtx; [i j 3]];
        logic_mtx = [logic_mtx; [j i 3]];
    end
end

%% construct composite_variable
if taylor_degree >= 1
    comp_var = TF_code;
    
    for i = 1: TF_num -1 
        for j = i+1: TF_num
            comp_var(length(comp_var) + 1) = strcat(TF_code(i), TF_code(j));
        end
    end
    
    %{
    if taylor_degree >= 2
      
        for i = 1: TF_num            
            comp_var(length(comp_var) + 1) = strcat(TF_code(i), '2');            
        end
    
        if taylor_degree >= 3
            for j = 1: TF_num - 1 
                for k = 1: TF_num
                    if k ~= j
                        comp_var(length(comp_var) + 1) = strcat(TF_code(j), '2', );            
                        
                        X{i} = [X{i} ; odata((i-1)* num_TF + j, 1 :size2_mZdif).^2 .* odata((i-1)* num_TF + k, 1 :size2_mZdif)]; % a2b, b2a
                    end
                end
            end
        end
        %}    
end

%% construct the fact_map
fact_map= zeros(length(logic_gate), length(comp_var));

for i = 1: length(logic_gate)
    str_a = logic_gate(i);
    if length(str_a{1}) == 1
        [str1, I0, I1] = intersect(str_a, comp_var);
        fact_map(i, I1) = 1;
    else
        if strcmp(logic_gate{i}(2), '&') == 1
            str_ab = strcat(logic_gate{i}(1), logic_gate{i}(3));
            [str_ab_match, I0, I1] = intersect(str_ab, comp_var);
            fact_map(i, I1) = 1;
        elseif strcmp(logic_gate{i}(2), '|') == 1
            str_a = logic_gate{i}(1);
            str_b = logic_gate{i}(3);
            
            [str_a_match, I0, I1] = intersect(str_a, comp_var);
            [str_b_match, I0, I2] = intersect(str_b, comp_var);
            
            fact_map(i, I1) = 1;
            fact_map(i, I2) = 1;
        elseif strcmp(logic_gate{i}(2), '>') == 1
            str_a = logic_gate{i}(1);
            str_b = logic_gate{i}(3);
            str_ab = strcat(logic_gate{i}(1), logic_gate{i}(3));
            
            [str_a_match, I0, I1] = intersect(str_a, comp_var);
            [str_ab_match, I0, I2] = intersect(str_ab, comp_var);
            
            
            
            if isempty(str_ab_match) == 0
                fact_map(i, I1) = 1;
                fact_map(i, I2) = -1;
            else
                str_a = logic_gate{i}(3);
                str_b = logic_gate{i}(1);
                str_ba = strcat(logic_gate{i}(3), logic_gate{i}(1));
                [str_b_match, I0, I3] = intersect(str_b, comp_var);
                [str_ba_match, I0, I4] = intersect(str_ba, comp_var);
                
                fact_map(i, I3) = 1;
                fact_map(i, I4) = -1;
            end
        end
    end
end

endstr = fun_x_range_calc(length(comp_var) + 1);

%{
xlswrite('model_file.xlsx', logic_gate(:), 'map23', ['A2:', 'A', num2str(length(logic_gate)+1)]);
xlswrite('model_file.xlsx', comp_var(1,:), 'map23', ['B1:', strcat(endstr{1}, '1')]);
xlswrite('model_file.xlsx', fact_map, 'map23', ['B2:' strcat(endstr{1}, num2str(length(logic_gate)+1))]);
%}


