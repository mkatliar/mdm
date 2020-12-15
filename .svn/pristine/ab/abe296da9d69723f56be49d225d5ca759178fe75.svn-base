%% Initialize
mode = 'test';

switch mode
    case 'training'
        N = 40000;
    case 'test'
        N = 60000;
    otherwise
        error('Unknown mode!');
end

N1 = 1;
N2 = N1 + 1000 - 1;

X = zeros(N2, ModelF.ParameterVectorLength);
S = zeros(N2, ModelF.ParameterVectorLength);
e = zeros(N2, 2);

% Set global optimization options.
ms = MultiStart('Display', 'off', 'StartPointsToRun', 'bounds-ineqs', ...
    'UseParallel', 'never');

%% Process objects
load '../data/sol.mat';

tic;
parfor i = N1 : N2
    fprintf('%s\t%d\t', mode, i);
    imp = ImagePair.Load(mode, i);
    
    est = Estimate(imp, optimset('GradObj', 'on', 'UseParallel', 'always'), ms, 3);
    
    X(i, :) = est.Value;
    e(i, :) = est.Model.Ellipticity();
    H(:, :, i) = est.Hessian;
    %Delta(i, :) = est.Delta;
    
%     figure(1);
%     imp.Show();
%     
%     figure(2);
%     model.GenerateImage(true).Show();
%     drawnow;
    
%     if strcmp(mode, 'training')
%         fprintf('ellipticity error = %f\tMSE so far = %f', ...
%             sqrt(mean((sol(i, 2 : 3) - e(i, :)).^2)), ...
%             sqrt(mean(mean((sol(1 : i, 2 : 3) - e(1 : i, :)).^2))));
%     end

    fprintf('ellipticity error = %f', ...
        sqrt(mean((sol.(mode)(i, 2 : 3) - e(i, :)).^2)));
    
    fprintf('\n');
end

estimate.(mode).X(N1 : N2, :) = X(N1 : N2, :);
estimate.(mode).H(:, :, N1 : N2) = H(:, :, N1 : N2);
estimate.(mode).e(N1 : N2, :) = e(N1 : N2, :);
save('estimate', 'estimate');

fprintf('Ellipticity MSE = %f\n', sqrt(mean(mean((sol.(mode)(N1 : N2, 2 : 3) ...
    - estimate.(mode).e(N1 : N2, :)).^2))));

toc
