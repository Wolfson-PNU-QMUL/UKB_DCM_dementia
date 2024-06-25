% Code written by Sam Ereira (s.ereira@qmul.ac.uk) - last updated 2/3/24
% Fits an elastic-net logistic regression model with K-fold cross-validation.

function [metrics] = UKB_DCM_dem_general_classifier(X_full,Y_full,K,k)


metrics = [];
for iter = 1:9 %iterate because the k-fold is sensitive to the way the data is randomised

    K_folds =[];

    %subset 1 is the prevalent cases and their matched controls
    tmp_list = GCM_sj_list;
    K_folds{1} = tmp_list(end-241:end,:);
    tmp_list(end-241:end,:)=[];

    %now do subsets 2-10
    % split cases and controls and randomly shuffle their order
    ctrls_outer = tmp_list(tmp_list(:,2)==0,:); ctrls_outer = ctrls_outer(randperm(size(ctrls_outer,1)),:);
    cases_outer = tmp_list(tmp_list(:,2)==1,:); cases_outer = cases_outer(randperm(size(cases_outer,1)),:);

    %work out how many cases will be in each subset
    ncases_outer = round(size(cases_outer,1)/(K-1));

    for i = 2:K
        if size(cases_outer,1) < 2*ncases_outer %at the end just use all remaining cases for a larger fold if necessary
            ncases_outer = size(cases_outer,1);
        end

        nctrls_outer = 10*ncases_outer;

        K_folds{i} = [cases_outer(1:ncases_outer,:);ctrls_outer(1:nctrls_outer,:)];

        cases_outer(1:ncases_outer,:) = [];
        ctrls_outer(1:nctrls_outer,:) = [];
    end


    %start OUTER FOLDS of cross-validation now. subset 1 (prevalents) is never
    %used as a test set so start with subset 2
    opt_lam = []; opt_alpha = [];
    alpha_list = 0:0.1:1; %list of alpha (mixing) parameters for elastic net
    for i_outer = 2:K

        test_set = K_folds{i_outer};
        Xtest_outer = X_full(test_set(:,3),:);
        Ytest_outer = Y_full(test_set(:,3),:);

        train_set = cell2mat(K_folds(setdiff([1:K], i_outer))');
        Xtrain_outer = X_full(train_set(:,3),:);
        Ytrain_outer = Y_full(train_set(:,3),:);

        pars = struct('fdev',-1); %ensures that we use the full sequence of 100 lambdas each time
        glmnetControl(pars);

        %At this point split the training data into k (5) subsets and run an inner
        %loop of cross validation with k (5) folds

        lambdas = nan(numel(alpha_list), 100,k);
        x_tmp = cell(100, numel(alpha_list),k); y_tmp = cell(100, numel(alpha_list),k);
        T_tmp = cell(100, numel(alpha_list),k); AUC_tmp = zeros(100, numel(alpha_list),k);

        %final column is class label
        ctrls_inner = Xtrain_outer(Ytrain_outer==0,:); ctrls_inner = [ctrls_inner,zeros(size(ctrls_inner,1),1)];
        cases_inner = Xtrain_outer(Ytrain_outer==1,:); cases_inner = [cases_inner,ones(size(cases_inner,1),1)];
        ncases_inner = floor(size(cases_inner,1)/k);

        k_folds =[];
        for i = 1:k
            if size(cases_inner,1) < 2*ncases_inner %at the end just use all remaining cases for a larger fold
                ncases_inner = size(cases_inner,1);
            end
            nctrls_inner = 10*ncases_inner;
            k_folds{i} = [cases_inner(1:ncases_inner,:);ctrls_inner(1:nctrls_inner,:)];
            cases_inner(1:ncases_inner,:) = [];
            ctrls_inner(1:nctrls_inner,:) = [];
        end


        %reset to default options to use full lambda sequence
        options = glmnetSet;

        %start inner folds of cross-validation
        for i_inner = 1:k

            %Define test-set
            Xtest_inner = k_folds{i_inner}(:,1:end-1);
            Ytest_inner = k_folds{i_inner}(:,end);

            %Define train-set
            Xtrain_inner = cell2mat(k_folds(setdiff([1:k], i_inner))');
            Ytrain_inner = Xtrain_inner(:,end);
            Xtrain_inner(:,end) = [];

            weights = ones(length(Ytrain_inner),1);
            weights(Ytrain_inner==1) = 1-sum(Ytrain_inner==1)./numel(Ytrain_inner); weights(Ytrain_inner==0) = 1-sum(Ytrain_inner==0)./numel(Ytrain_inner);

            prob_tmp = nan(size(Xtest_inner,1), 100, numel(alpha_list));
            for alpha_ind = 1:numel(alpha_list) %Loop over hyperparameters to optimise
                options.alpha = alpha_list(alpha_ind);

                %Train logistic regression model
                options.weights = weights;
                fit = glmnet(Xtrain_inner, Ytrain_inner, 'binomial', options);

                %add padding if glmnet stopped short of 100 lambdas
                offset = 100 - size(fit.beta,2);
                if offset >0
                    tmp = fit.beta(:,1); pad = repmat(tmp, 1, offset); fit.beta = [pad, fit.beta];
                    tmp = fit.lambda(1); pad = repmat(tmp, offset, 1); fit.lambda = [pad;fit.lambda];
                    tmp = fit.a0(1); pad = repmat(tmp, 1, offset); fit.a0 = [pad, fit.a0];
                end

                %Output probabilities from regression models...P(case)
                prob_tmp(:,:, alpha_ind) = glmnetPredict(fit, Xtest_inner,[],'response');

                %save lambda sequences as they are always slightly different
                lambdas(alpha_ind,:,i_inner) = fit.lambda;
            end


            % Fit a ROC curve for each combination of hyperparameters for each
            % inner fold
            for l = 1:size(prob_tmp, 2)
                for alpha = 1:numel(alpha_list)
                    [x_tmp{l,alpha,k},y_tmp{l,alpha,k},T_tmp{l,alpha,k},AUC_tmp(l,alpha,k)] = perfcurve(Ytest_inner,squeeze(prob_tmp(:,l,alpha)),1);
                end
            end
        end

        %find optimal hyperparameters by maximising average AUC-ROC
        AUC_tmp = mean(AUC_tmp,3); %average over 5 inner folds
        [~,ind] = max(AUC_tmp(:)); %find the maximum AUC
        [opt_lam_ind, opt_alpha_ind] = ind2sub([size(AUC_tmp,1), size(AUC_tmp,2)],ind);
        opt_lam(i_outer) = mean(lambdas(opt_alpha_ind,opt_lam_ind,:));
        opt_alpha(i_outer) = alpha_list(opt_alpha_ind);

        %train classifier on full training set using these optimised parameters
        options.lambda = opt_lam(i_outer); options.alpha = opt_alpha(i_outer);
        weights = ones(numel(Ytrain_outer),1);
        weights(Ytrain_outer==1) = 1-sum(Ytrain_outer==1)./numel(Ytrain_outer);
        weights(Ytrain_outer==0) = 1-sum(Ytrain_outer==0)./numel(Ytrain_outer);
        options.weights = weights;
        fit = glmnet(Xtrain_outer, Ytrain_outer, 'binomial', options);

        %test classifier on left out test set and save p(case) and binary prediction
        prob = glmnetPredict(fit, Xtest_outer,[],'response');
        [metrics(iter).x(:,i_outer-1), metrics(iter).y(:,i_outer-1), metrics(iter).T(:,i_outer-1), metrics(iter).AUC(i_outer-1)] = perfcurve(Ytest_outer,prob,1);

        % Point biserial Correlation
        [metrics(iter).corr.r(i_outer-1),~,metrics(iter).corr.p(i_outer-1),metrics(iter).corr.ci(:,i_outer-1)] = pointbiserial(Ytest_outer,prob);
    end
end

end
