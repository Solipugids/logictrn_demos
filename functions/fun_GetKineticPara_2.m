
%function [para, err, objv1, objv2] = fun_GetKineticPara_1(X, X_1, n_reg, YA, YB, rg, operator_number)

function [para, err, objv1] = fun_GetKineticPara_2(X, X_1, YA, logic_numberCode, number_of_logics, number_of_regulators, range, regulator_record)

NVAR = size(range, 2);    

NIND = 40; % Number of individuals
MAXGEN = 500; % Maximum no. of generations

PRECI = 20; % Precision of variables
GGAP = 0.9; % Generation gap
% Build field descriptor


FieldD = [rep([PRECI],[1,NVAR]); range; rep([1;0;1;1],[1,NVAR])];
% Initialise population
Chrom = crtbp(NIND, NVAR*PRECI);
gen = 1; % Counter
% Evaluate initial population
p0 = bs2rv(Chrom,FieldD); 

sum_p0_omega = sum(p0(:, 1: number_of_logics), 2);

sum_p0_omega_mtx = rep(sum_p0_omega, [1 number_of_logics]);

p0(:, 1: number_of_logics) = p0(:, 1: number_of_logics) ./ sum_p0_omega_mtx;

    %{
    for i = 1: length(regulator_record)
        if regulator_record(i) ~= 0
            p0(:, number_of_logics + i) = p0(:, number_of_logics+regulator_record(i));
            %p0(:, i) = (p0(:, i) + p0(:, regulator_record(i))) /2;
        end
    end
    %}
for i = 1: max(regulator_record)
    ind = find(regulator_record == i);
    for j = 2: length(ind)
        p0(:, number_of_logics + ind(j)) = p0(:, number_of_logics + ind(1));
    end
end

X_est = fun_X_estimation2(X_1, logic_numberCode, number_of_logics, number_of_regulators, YA, p0);


Xmtr = rep(X, [size(X_est, 1), 1]);
%ObjV = objfun1(X_est - Xmtr);
ObjV = objfun6(X_est - Xmtr);


% Generational loop
while gen < MAXGEN,
    % Assign fitness values to entire population
    FitnV = ranking(ObjV);
    % Select individuals for breeding
    SelCh = select('sus', Chrom, FitnV, GGAP);
    % Recombine individuals (crossover)
    SelCh = recombin('xovsp',SelCh,0.7);
    % Apply mutation
    SelCh = mut(SelCh);
    % Evaluate offspring, call objective function
    p0 = bs2rv(SelCh,FieldD);
    
    %{
    for i = 1: length(regulator_record)
        if regulator_record(i) ~= 0
            %p0(:, i) = p0(:, regulator_record(i));
            p0(:, number_of_logics + i) = p0(:, number_of_logics+regulator_record(i));
            %p0(:, i) = (p0(:, i) + p0(:, regulator_record(i))) /2;
        end
    end
    %}
    for i = 1: max(regulator_record)
       ind = find(regulator_record == i);
       for j = 2: length(ind)
          p0(:, number_of_logics + ind(j)) = p0(:, number_of_logics + ind(1));
       end
    end
    
    sum_p0_omega = sum(p0(:, 1: number_of_logics), 2);

    sum_p0_omega_mtx = rep(sum_p0_omega, [1 number_of_logics]);

    p0(:, 1: number_of_logics) = p0(:, 1: number_of_logics) ./ sum_p0_omega_mtx;
    
    X_est = fun_X_estimation2(X_1, logic_numberCode, number_of_logics, number_of_regulators, YA, p0);

        
    Xmtr = rep(X, [size(X_est, 1), 1]);
    %ObjVSel = objfun1(X_est - Xmtr);
     
    ObjVSel = objfun1(X_est - Xmtr);
    %ObjVSel = objfun6(X_est - Xmtr);
    
    %ObjVSel = log10(objfun1(X_est - Xmtr) + 1) + (fun_entropy(p(:,1:6)).^2)';
    
    objv1(gen) = min(ObjVSel);
    %objv1(gen) = min(log10(objfun1(X_est - Xmtr) + 1));
    %objv2(gen) = min(fun_entropy(p(:,1:6)).^2);

    % Reinsert offspring into population
    [Chrom ObjV]=reins(Chrom,SelCh,1,1,ObjV,ObjVSel);
    % Increment counter

    % Increment counter
    gen = gen+1;
    gen;
end

p0 = bs2rv(Chrom, FieldD);

sum_p0_omega = sum(p0(:, 1: number_of_logics), 2);

sum_p0_omega_mtx = rep(sum_p0_omega, [1 number_of_logics]);

p0(:, 1: number_of_logics) = p0(:, 1: number_of_logics) ./ sum_p0_omega_mtx;

%{
    for i = 1: length(regulator_record)
        if regulator_record(i) ~= 0
            %p0(:, i) = p0(:, regulator_record(i));
            p0(:, number_of_logics + i) = p0(:, number_of_logics+regulator_record(i));
            %p0(:, i) = (p0(:, i) + p0(:, regulator_record(i))) /2;
        end
    end
%}
for i = 1: max(regulator_record)
    ind = find(regulator_record == i);
    for j = 2: length(ind)
        p0(:, number_of_logics + ind(j)) = p0(:, number_of_logics + ind(1));
    end
end

X_est = fun_X_estimation2(X_1, logic_numberCode, number_of_logics, number_of_regulators, YA, p0);


Xmtr = rep(X, [size(X_est, 1), 1]);
 
ObjV1 = objfun1(X_est - Xmtr);
FitnV1 = ranking(ObjV1);
[Mv, MI] = max(FitnV1);

para = p0(MI, :);
err = sum((X_est(MI,:)-X).^2);


end


