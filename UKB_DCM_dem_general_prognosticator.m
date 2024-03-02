% Code written by Sam Ereira (s.ereira@qmul.ac.uk) - last updated 2/3/24
% Fits an elastic-net linear regression model with K-fold cross-validation


function [metrics] = UKB_DCM_dem_general_prognosticator(X,Y,K,k)

metrics = [];
for iter = 1:9 %iterate because the kfold is sensitive to the way the data is randomised

    % Split into K folds
    K_folds =[];

    tmp_list = [1:length(Y)];
    cases_outer = tmp_list(randperm(length(tmp_list)));

    %work out how many cases will be in each subset
    ncases_outer = round(length(cases_outer)/(K));

    for i = 1:K
        if length(cases_outer) < 2*ncases_outer %at the end just use all remaining cases for a larger fold if necessary
            ncases_outer = length(cases_outer);
        end
        K_folds{i} = [cases_outer(1:ncases_outer)];
        cases_outer(1:ncases_outer) = [];
    end


    %start OUTER FOLDS of cross-validation now.
    alpha_list = 0:0.1:1; %list of alpha (mixing) parameters for elastic net
    Yhat =[]; Yreal=[];
    for i_outer = 1:K

        test_set = K_folds{i_outer};
        Xtest_outer = X(test_set,:);
        Ytest_outer = Y(test_set);

        train_set = cell2mat(K_folds(setdiff([1:K], i_outer)))';
        Xtrain_outer = X(train_set,:);
        Ytrain_outer = Y(train_set);

        pars = struct('fdev',-1); %ensures that we use the full sequence of 100 lambdas each time
        glmnetControl(pars);

        %At this point split the training data into k (5) subsets and run an inner
        %loop of cross validation with k (5) folds

        lambdas = nan(numel(alpha_list), 100,k);

        %final column is class labels
        cases_inner = train_set;
        ncases_inner = floor(size(cases_inner,1)/k);

        k_folds =[];
        for i = 1:k
            if length(cases_inner) < 2*ncases_inner %at the end just use all remaining cases for a larger fold
                ncases_inner = length(cases_inner);
            end

            k_folds{i} = cases_inner(1:ncases_inner);
            cases_inner(1:ncases_inner) = [];
        end


        %reset to default options to use full lambda sequence
        options = glmnetSet;

        %start inner folds of cross-validation
        MSqE = [];lambdas = []; SqE = [];
        for i_inner = 1:k
            %Define test-set
            Xtest_inner = X(k_folds{i_inner},:);
            Ytest_inner = Y(k_folds{i_inner});

            %Define train-set
            Xtrain_inner = X(cell2mat(k_folds(setdiff([1:k], i_inner))'),:);
            Ytrain_inner = Y(cell2mat(k_folds(setdiff([1:k], i_inner))'));

            for alpha_ind = 1:numel(alpha_list) %Loop over hyperparameters to optimise
                options.alpha = alpha_list(alpha_ind);

                %Train linear regression model
                fit = glmnet(Xtrain_inner, Ytrain_inner, 'gaussian', options);

                %add padding if glmnet stopped short of 100 lambdas
                offset = 100 - size(fit.beta,2);
                if offset >0
                    tmp = fit.beta(:,1); pad = repmat(tmp, 1, offset); fit.beta = [pad, fit.beta];
                    tmp = fit.lambda(1); pad = repmat(tmp, offset, 1); fit.lambda = [pad;fit.lambda];
                    tmp = fit.a0(1); pad = repmat(tmp, offset,1); fit.a0 = [pad; fit.a0];
                end

                %Test model predictions and compute squared error
                SqE{i_inner}(:,:,alpha_ind) = (glmnetPredict(fit, Xtest_inner) - Ytest_inner).^2;

                %save lambda sequences as they are always slightly different
                lambdas(alpha_ind,:,i_inner) = fit.lambda;
            end

        end

        %Find hyperparameters that minimise median squared error
        tmp = cellfun(@(x) squeeze(median(x, 1)), SqE, 'UniformOutput', false);
        MSqE = mean(cat(3, tmp{:}),3); [~, ind] = min(MSqE(:));
        [opt_lam_ind, opt_alpha_ind] = ind2sub([size(MSqE,1), size(MSqE,2)],ind);
        options.alpha = alpha_list(opt_alpha_ind);
        options.lambda = mean(lambdas(opt_alpha_ind,opt_lam_ind,:),3);

        fit = glmnet(Xtrain_outer, Ytrain_outer, 'gaussian', options);

        %Model prediction for test set (lef out subject)
        Yhat{i_outer} = glmnetPredict(fit, Xtest_outer);
        Yreal{i_outer} = Ytest_outer;
    end

    Yhat = cat(1, Yhat{:}); %concatenated folds
    Yreal = cat(1,Yreal{:});

    metrics(iter).Yhat = Yhat;
    metrics(iter).Yreal = Yreal
    [metrics(iter).PearsonR, metrics(iter).PearsonP] = corr(Yhat,Yreal)
    [metrics(iter).SpearmanR, metrics(iter).SpearmanP] = corr(Yhat,Yreal, 'type','Spearman');
end

end