%%
imp = ImagePair.Load('training', 1);
x = X(1, :);

i1 = ModelF.Index.GalaxyShape(2);
i2 = ModelF.Index.GalaxyShape(3);

[lb, ub] = ModelF.Bounds(true);

% lb(i1) = 0.6;
% ub(i1) = 0.8;
lb([i1, i2]) = [0.6, 1.5];
ub([i1, i2]) = [0.8, 2.0];

n1 = 100;
n2 = 100;

x1 = lb(i1) + (ub(i1) - lb(i1)) / n1 * (0 : n1);
x2 = lb(i2) + (ub(i2) - lb(i2)) / n2 * (0 : n2);
L = zeros(n1 + 1, n2 + 1);

parfor k1 = 1 : n1 + 1
    m = ModelF();
    for k2 = 1 : n2 + 1
        m.ParameterVector = x;
        m.ParameterVector([i1, i2]) = [x1(k1), x2(k2)];
        L(k2, k1) = m.LogLikelihood(imp);
    end
end

%% Plot
imagesc(x1, x2, exp(L - max(L(:)))); axis xy