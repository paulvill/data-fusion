function [t, tmin, tmax] = predict_time_image(model, image, symmetrize)
% this function takes an image and returns the predicted time, along with
% the 95% confidence intervals for the time
% image is a single-channel image of nuclei, oriented with the dorsal side
% at the top
% t is the predicted time
% tmin is the minimum predicted time in the 95% confidence interval
% tmax is the maximum predicted time in the 95% confidence interval

if nargin < 3
    symmetrize = false;
end

rot_angle = 0;
image = image_normalize(image, rot_angle, symmetrize);
data = (create_PCA_data(image) - model.mu)*model.V;

[ypred, yci] = predict(model.reg_mdl, data, 'prediction', 'observation', ...
    'alpha', 1-0.682);

t = ypred;
tmin = yci(1);
tmax = yci(2);


