
function X_est = fun_X_estimation2(X_1, logic_numberCode, number_of_logics, number_of_regulators, YA, p0)

len = length(YA{1}(1,:));

%omega = p0(:, 1:number_of_logics);
Imax_mtx = rep(p0(:, number_of_logics + number_of_logics * 2 + 1), [1 len]);
Imax_arr = p0(:, number_of_logics + number_of_logics * 2 + 1);
%kd_arr = p0(:, number_of_logics + number_of_logics * 2 + 2);
kd_arr = p0(:, end);

X_est_0 = 0;
X_est = 0;

for i = 1: number_of_logics        

    if number_of_regulators(i) == 1 %% one logic with one regulator
        %ka, Imax, kd        
        omega_ith = rep(p0(:, i), [1 len]);
        kb_1_arr = p0(:, number_of_logics + (i-1) * 2 + 1);
        X_est_ith_logic{i} = omega_ith .* Imax_mtx .* (1 - exp(- kb_1_arr ./ Imax_arr * (YA{1}(1,:))));
    else   %% one logic with two regulators
         if logic_numberCode(i, 3) == 2 % logic '&'
             %kab, Imax, kd
             omega_ith = rep(p0(:, i), [1 len]);
             kb_1_arr = p0(:, number_of_logics + (i-1) * 2 + 1);
             X_est_ith_logic{i} = omega_ith .* Imax_mtx .* (1 - exp(- kb_1_arr./ Imax_arr * (YA{i}(1,:) .* YA{i}(2,:))));
         elseif logic_numberCode(i, 3) == 1 % logic '|'
             %ka, kb, Imax, kd
             omega_ith = rep(p0(:, i), [1 len]);
             kb_1_arr = p0(:, number_of_logics + (i-1) * 2 + 1);
             kb_2_arr = p0(:, number_of_logics + (i-1) * 2 + 2);
             X_est_ith_logic{i} = omega_ith .* (Imax_mtx /2) .* (1 - exp(- kb_1_arr./ Imax_arr * (YA{i}(1,:))) + 1 - exp(- kb_2_arr./ Imax_arr * (YA{i}(2,:))));
         else % logic '>'
             %ka, Imax, kd
             omega_ith = rep(p0(:, i), [1 len]);
             kb_1_arr = p0(:, number_of_logics + (i-1) * 2 + 1);
             X_est_ith_logic{i} = omega_ith .* Imax_mtx .* (1 - exp(- kb_1_arr./ Imax_arr * (YA{i}(1,:) .* (1-YA{i}(2,:)))));
         end   
    end    

    X_est_0 = X_est_0 + X_est_ith_logic{i};
end

X_est = X_est_0 + (1- kd_arr) * X_1;

