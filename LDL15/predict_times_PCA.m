function varargout = predict_times_PCA(train_data, train_times, test_data, nmodes)

[m, n] = size(train_data);
[m2, n2] = size(test_data);

mean_train_data = mean(train_data);

train_data = train_data - repmat(mean_train_data, m, 1);
test_data = test_data - repmat(mean_train_data, m2, 1);

[U, S, V] = svds(train_data, nmodes);

mdl = fitlm(train_data*V,train_times,'linear');

[pred_time,pred_time_int] = predict(mdl ,test_data*V);

varargout{1} = pred_time;
if nargout > 1
    varargout{2} = pred_time_int;
end