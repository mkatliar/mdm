%% Initialize
clear;
model = ModelF();

% model.SersicIndex = 0.89;
% model.MoffatBeta = 3.6;
model.PSF.ShapeParameters = [[2 1.5/2], 0.3 * pi];
model.Galaxy.ShapeParameters = [[4 3/4], 0.4 * pi];
model.ScaleS = 2000;
model.ScaleG = 17000;
model.NoiseS.Sigma = 1e2;
model.NoiseG.Sigma = 1e2;
model.OffsetS = [2.3 1];
model.OffsetG = [-2.3 0.5];
% Compute and show image
img = model.GenerateImage(true);

figure(1);
img.Show();

[l, grad] = model.LogLikelihood(img)

%% Estimate parameters
tic

gs = GlobalSearch('Display', 'iter', 'StartPointsToRun', 'bounds-ineqs', ...
    'NumStageOnePoints', 5, ...
    'NumTrialPoints', 100);

ms = MultiStart('Display', 'iter', 'StartPointsToRun', 'bounds-ineqs', ...
    'UseParallel', 'never');

est = Estimate(img, ...
    optimset('Algorithm', 'interior-point', ...
        'Display', 'notify', ...
        'GradObj', 'on', 'UseParallel', 'never'), ...
    ms, 5);

disp('----------------------------------------------------')
est.FunctionValue
[model.ParameterVector; est.Value; est.Delta; est.Gradient]

% Compute ellipticity
disp('True vs estimated ellipticity:');
disp([model.Ellipticity(); est.Model.Ellipticity]);
fprintf('Ellipticity estimation error: %f\n', sqrt(mean((model.Ellipticity - est.Model.Ellipticity).^2)));
toc

figure(2);
est.Model.GenerateImage(true).Show();