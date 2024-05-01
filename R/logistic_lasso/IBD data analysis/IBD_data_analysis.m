%set path to correct directory
%check path
pwd

%import data 
load([pwd,'/IBD data analysis/IBD.mat'])
load([pwd,'/IBD data analysi/IBD_bac_name.mat'])

y = [zeros(26,1); ones(85,1)];
w = IBD1>0;
ind_d = sum(w);

%pick main bacteria
pick = find(ind_d>= 111*0.2);
x_or = IBD1(: , pick);

x_min=min( x_or(x_or>0));
% eliminate 0
x_or = x_or + 0.5*x_min-min(x_or, 0.5*x_min);

%convert to compositional data
[n p] = size(x_or);
x = x_or./(ones(p,1)*sum(x_or'))';
x = log(x);
x = [ones(n,1) x];

train_set_y = [y(1:16);y(27:82)];
test_set_y = [y(17:24);y(83:110)];

train_set_x = [x(1:16,:);x(27:82,:)];
test_set_x = [x(17:24,:);x(83:110,:)];


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

mu=1; level =0.95; penalize =0; length1 = 15; length2 =25;
d=1.5;
[res.beta_n,res.lambda_best, res.EBIC] = biased_estimate_BIC(penalize,y,x_one,constr3, mu, length1, length2);
[res.beta_u,res.CI_l, res.CI_u, res.CI_M] = debiased_cvx(y,x_one, res.beta_n, constr3, res.lambda_best, level,c,d);
res.beta_n = c*res.beta_n;
sum(res.beta_n(2:78))+sum(res.beta_u(2:78))


width = (res.CI_u-res.CI_l)/2;
width = width/norminv(0.975)*norminv(1-(1-level)/2);
% find the correct bacteria;
variable = ((res.CI_u < 0)|(res.CI_l > 0))*1;
variable_id = find(variable==1);
variable_id = variable_id(variable_id>1)-1;
pick_id = pick(variable_id);
bacteria = IBD_bac_name(1,pick_id);

CI = zeros(78*2,2);
CI(1:78,1) = res.beta_u;
CI(79:(78*2),1) = res.beta_n;

CI(1:78,2) = (res.CI_u-res.CI_l)/2;

filename2 = [pwd,'\IBD data analysis\results\IBD_result_vs.mat'];
save(filename2,'res');

filename3 = [pwd,'\IBD data analysis\results\IBD_result_CI.csv'];
csvwrite(filename3, CI);

filename4 = [pwd,'\IBD data analysis\results\IBD_result_id_vs_new.csv'];
csvwrite(filename4, pick_id);

filename5 = [pwd,'\IBD data analysis\results\IBD_result_bac_vs_new.csv'];
export(bacteria,'File',filename5,'Delimiter',',')

%% onec train/test
id_control = [1:1:26];
id_case = [27:1:111];
%resampling
nrep=50;
rng(1111);
for i=1:nrep
    % train data set
    sample_case= sort(randsample(id_case, 56));
    sample_control = sort(randsample(id_control, 16));
    sample = [sample_control,sample_case];
    % test set
    test_case = setdiff(id_case,sample_case);
    test_case = sort(randsample(test_case, 28));
    test_control = setdiff(id_control, sample_control);
    test_control = sort(randsample(test_control, 8));
    test = [test_control,test_case];
    
    rep_set{i}.x = x(sample,:);
    rep_set{i}.y = y(sample);
    rep_set{i}.xtest = x(test, :);
    rep_set{i}.ytest = y(test, :);
end
for i=1:nrep
x = rep_set{i}.x;
y = rep_set{i}.y;

[n p] = size(x);
constr =  [0;ones(p-1,1)];
[constr2, S, V] = svd(constr);
constr= constr2(: , 1:size(S,2));
Pc=constr*constr';
x = x - x*Pc;

[n p] = size(x);

    for j=1:p
        norm_x(j) = 1/10*norm(x(:,j),2);
    end
    b2 = 1./norm_x;
    c =diag(b2);
    x2 =x*c;
    constr3 = c*constr;
    %constr3 = constr3';
    
    if (sum(abs(constr3))==0)
        %constr = zeros(p,1);
    else
        [constr4,S,V] = svd(constr3);
        constr3 = constr4(:,1:size(S,2));
    end
    constr3 = constr3';
    
    %useset{i}.x = x;
    useset{i}.x2 = x2;
    useset{i}.y = y;
    useset{i}.constr3 = constr3;
    useset{i}.c = c;
end

penalize =0;
result = cell(nrep);
result_select = cell(nrep);
parpool('local',4);
pick = find(ind_d>= 111*0.2);
parfor i=1:nrep
    y = useset{i}.y;
    x2 = useset{i}.x2;
    xtest = rep_set{i}.xtest;
    ytest = rep_set{i}.ytest;
    
    constr3 = useset{i}.constr3;
    c = useset{i}.c;
    mu=1; level =0.95; penalize =0; tolerance2 = 1e-6; step =1e5;length1 = 15; length2 =25;d=1.5;
    %one constraint
    [beta_n_1,lambda_best_1,EBIC_1] = biased_estimate_BIC(penalize,y,x2,constr3, mu, length1, length2);
    [beta_u_1,CI_l_1,CI_u_1,CI_M_1] = debiased_cvx(y,x2, beta_n_1, constr3,lambda_best_1, level,c,d);
    beta_n_1 = c*beta_n_1;

    %no constraint
    [beta_n,lambda_best, EBIC] = biased_estimate_BIC(penalize,y,x2,0, mu,length1, length2);
    [beta_u,CI_l, CI_u, CI_M] = debiased_cvx(y,x2, beta_n, 0, lambda_best, level,c,d);
    beta_n = c*beta_n;

    variable = ((CI_u_1 < 0)|(CI_l_1 > 0))*1;
    variable_id = find(variable==1);
    beta_u_1_s = beta_u_1([variable_id]);
    x_1_s = xtest(:,variable_id);
    
    
    variable = ((CI_u < 0)|(CI_l > 0))*1;
    variable_id = find(variable==1);
    beta_u_s = beta_u([variable_id]);
    x_s = xtest(:,variable_id);
    
    %ROC
    %full variable
    %no constraint, biased
    prob_noc_n = xtest * beta_n;
    %no constraint, debiased
    prob_noc_u = xtest * beta_u;
    % one constraint, no debiased
    prob_onec_n = xtest * beta_n_1;
    % one constraint, debiased
    prob_onec_u = xtest * beta_u_1;
    
    %selected variable
    %no constraint, biased
    prob_noc_n_s = xtest * beta_n;
    %no constraint, debiased
    prob_noc_u_s = x_s * beta_u_s;
    % one constraint, no debiased
    prob_onec_n_s = xtest * beta_n_1;
    % one constraint, debiased
    prob_onec_u_s = x_1_s * beta_u_1_s;
    
    result{i} = [prob_noc_n;prob_noc_u; prob_onec_n; prob_onec_u];
    result_select{i} = [prob_noc_n_s;prob_noc_u_s; prob_onec_n_s; prob_onec_u_s];
end

delete(gcp)

save([pwd,'\IBD data analysis\ROC_result.mat'],'result')
save([pwd,'\IBD data analysis\ROC_result_select.mat'],'result_select')

load([pwd,'\IBD data analysis\ROC_result.mat'])
load([pwd,'\IBD data analysis\ROC_result_select.mat'])

ROC_result = zeros(nrep, 144);
for i=1:nrep
    prob_temp = result{i};
    % no constraint, lasso
    prob = prob_temp(1:36,:);  %36 = 28 + 8
    for j =1: length(prob)
        if prob(j)>=10
            prob1(j) =1;
        elseif prob(j) <=-10
            prob1(j) =0;
        else
            prob1(j) = exp(prob(j))/(exp(prob(j))+1);
        end
    end
    prob_noc_n = prob1;
    
    %no constraint, de-biased
    prob = prob_temp(37:72,:);  %36 = 28 + 8
    for j =1: length(prob)
        if prob(j)>=10
            prob1(j) =1;
        elseif prob(j) <=-10
            prob1(j) =0;
        else
            prob1(j) = exp(prob(j))/(exp(prob(j))+1);
        end
    end
    prob_noc_u = prob1;
    
    prob = prob_temp(73:108,:);  %36 = 28 + 8
    for j =1: length(prob)
        if prob(j)>=10
            prob1(j) =1;
        elseif prob(j) <=-10
            prob1(j) =0;
        else
            prob1(j) = exp(prob(j))/(exp(prob(j))+1);
        end
    end
    prob_onec_n = prob1;
    
    prob = prob_temp(109:144,:);  %36 = 28 + 8
    for j =1: length(prob)
        if prob(j)>=10
            prob1(j) =1;
        elseif prob(j) <=-10
            prob1(j) =0;
        else
            prob1(j) = exp(prob(j))/(exp(prob(j))+1);
        end
    end
    prob_onec_u = prob1;
    
    ROC_result(i,:) = [prob_noc_n,prob_noc_u,prob_onec_n,prob_onec_u];
end

ROC_result_select = zeros(nrep, 144);
for i=1:nrep
    prob_temp = result_select{i};
    % no constraint, lasso
    prob = prob_temp(1:36,:);  %36 = 28 + 8
    for j =1: length(prob)
        if prob(j)>=10
            prob1(j) =1;
        elseif prob(j) <=-10
            prob1(j) =0;
        else
            prob1(j) = exp(prob(j))/(exp(prob(j))+1);
        end
    end
    prob_noc_n = prob1;
    
    %no constraint, de-biased
    prob = prob_temp(37:72,:);  %36 = 28 + 8
    for j =1: length(prob)
        if prob(j)>=10
            prob1(j) =1;
        elseif prob(j) <=-10
            prob1(j) =0;
        else
            prob1(j) = exp(prob(j))/(exp(prob(j))+1);
        end
    end
    prob_noc_u = prob1;
    
    prob = prob_temp(73:108,:);  %36 = 28 + 8
    for j =1: length(prob)
        if prob(j)>=10
            prob1(j) =1;
        elseif prob(j) <=-10
            prob1(j) =0;
        else
            prob1(j) = exp(prob(j))/(exp(prob(j))+1);
        end
    end
    prob_onec_n = prob1;
    
    prob = prob_temp(109:144,:);  %36 = 28 + 8
    for j =1: length(prob)
        if prob(j)>=10
            prob1(j) =1;
        elseif prob(j) <=-10
            prob1(j) =0;
        else
            prob1(j) = exp(prob(j))/(exp(prob(j))+1);
        end
    end
    prob_onec_u = prob1;
    
    ROC_result_select(i,:) = [prob_noc_n,prob_noc_u,prob_onec_n,prob_onec_u];
end

filename = [pwd,'\IBD data analysis\ROC_result.csv'];
csvwrite(filename, ROC_result);
filename2 = [pwd,'\IBD data analysis\ROC_result_select.csv'];
csvwrite(filename2, ROC_result_select);

% random forest
ROC_result_random_forest = zeros(nrep, 36);
for i=1:nrep
    x = rep_set{i}.x;
    y = rep_set{i}.y;
    xtest = rep_set{i}.xtest;
    ytest = rep_set{i}.ytest;
    %run random forest
    NumTrees = 50;
    Mdl = TreeBagger(NumTrees,x,y);
    [Yfit,scores] = predict(Mdl,xtest);
    ROC_result_random_forest(i,:) = scores(:,2);
end
filename = [pwd,'\IBD data analysis\ROC_result_rf.csv'];
csvwrite(filename, ROC_result_random_forest);


