function model = train_PCA_reg_model(train_data, train_times, nmodes)
    [m, n] = size(train_data);

    mean_train_data = mean(train_data);

    % calculate PCA coefficients
    train_data = train_data - repmat(mean_train_data, m, 1);
    [U, S, V_train] = svds(train_data, nmodes);

    % fit model
    mdl = fitlm(train_data*V_train, train_times, 'linear');

    % store data
    model.mu = mean_train_data;
    model.V = V_train;
    model.reg_mdl = mdl;
end
