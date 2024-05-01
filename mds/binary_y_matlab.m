cd /n/janson_lab/lab/sma/CompDA_paper/R/logistic_lasso;
i_job = string(i_job);
x = readtable(append("/n/janson_lab/lab/sma/CompDA_paper/results/simulation/binary_Y/data/x_", i_job, ".csv"));
x = x{:,:};
y = readtable(append("/n/janson_lab/lab/sma/CompDA_paper/results/simulation/binary_Y/data/y_", i_job, ".csv"));
y = y{:,:};

w = x>0;
ind_d = sum(w);

%pick main bacteria
[n_old p_old] = size(x);
pick = find(ind_d>= p_old*0.1);
x_or = x(: , pick);

x_min=min( x_or(x_or>0));
% eliminate 0
x_or = x_or + 0.5*x_min-min(x_or, 0.5*x_min);

%convert to compositional data
[n p] = size(x_or);
x = x_or./(ones(p,1)*sum(x_or'))';
x = log(x);
x = [ones(n,1) x];

%% onec variable selection
[n p] = size(x);
constr =  [0;ones(p-1,1)];
[constr2, S, V] = svd(constr);
constr= constr2(: , 1:size(S,2));
Pc=constr*constr';
x = x - x*Pc;

for j=1:p
    norm_x(j) = 1/10*norm(x(:,j),2);
end
b2 = 1./norm_x;
c =diag(b2);
x_one =x*c;
constr3 = constr'*c;
[constr2, S, V] = svd(constr3');
constr3= constr2(: , 1:size(S,2));
constr3= constr3';

mu=1; level =0.95; penalize =0; length1 = 15; length2 =25; d=1.5;
[res.beta_n,res.lambda_best, res.EBIC] = biased_estimate_BIC(penalize,y,x_one,constr3, mu, length1, length2);
[res.beta_u,res.CI_l, res.CI_u, res.CI_M] = debiased_cvx(y,x_one, res.beta_n, constr3, res.lambda_best, level,c,d);

CI = zeros(p_old,3);
CI(pick,1) = res.beta_u(2:p);
CI(pick,2) = res.CI_l(2:p);
CI(pick,3) = res.CI_u(2:p);

filename = append("/n/janson_lab/lab/sma/CompDA_paper/results/simulation/binary_Y/fit/matlab_", i_job, ".csv");
csvwrite(filename, CI);
