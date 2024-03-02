function output = UKB_DCM_dem_mediation(x,m,y,covs,weights, nperm)

%Conduct a simple regression analysis with X predicting Y to
%test for path c alone
mdl = fitglm(zscore([x, covs]),y,'weights', weights,  'CategoricalVars', [3],'Distribution', 'binomial');
path_c = mdl.Coefficients(2,:);
path_c.Properties.RowNames = {'c'};

%Conduct a simple regression analysis with X predicting M to 
% test for path a
mdl = fitglm(zscore([x,covs]),m, 'weights', weights);
path_a = mdl.Coefficients(2,:);
path_a.Properties.RowNames = {'a'};

%Conduct a simple regression analysis with M predicting Y to
%test the significance of path b alone
mdl = fitglm(zscore([m,covs]),y, 'weights', weights,  'CategoricalVars', [3],'Distribution', 'binomial');
path_b = mdl.Coefficients(2,:);
path_b.Properties.RowNames = {'b'};

%Conduct a multiple regression analysis with X and M 
% predicting Y
mdl = fitglm(zscore([x, m,covs]),y, 'weights', weights, 'CategoricalVars', [3], 'Distribution', 'binomial');
path_b2 = mdl.Coefficients(3,:);
path_b2.Properties.RowNames = {'b2'};
path_cprime = mdl.Coefficients(2,:);
path_cprime.Properties.RowNames = {'cprime'};

%formally calculate indirect effect by multiplying coefficients
path_ab_estimate = path_a.Estimate.*path_b2.Estimate;


%Run permutation test to create a null distribution for path_ab
for iperm = 1:nperm
    shuffle = randperm(length(y));
    yperm = y(shuffle);
    weightsperm = weights(shuffle);

    mdl1 = fitglm(zscore([x,covs]),m, 'weights', weightsperm,'CategoricalVars', [3]);
    mdl2 = fitglm(zscore([x, m,covs]),yperm, 'weights', weightsperm, 'CategoricalVars', [3],'Distribution', 'binomial');
    path_ab_null(iperm) = mdl1.Coefficients.Estimate(2)*mdl2.Coefficients.Estimate(3);
end

%generate p-value for path_ab
ind = find(path_ab_estimate == sort([path_ab_null, path_ab_estimate]));
p = 1- ind/(length(path_ab_null)+1);
p_lim = 1/nperm;
p = max([p, p_lim]);

path_ab = table(path_ab_estimate, nan, nan, p,'Rownames', {'ab'}, 'VariableNames', path_a.Properties.VariableNames);
output = [path_a;path_b;path_c;path_cprime;path_ab];


end