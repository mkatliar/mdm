function ShowPrior(x, H, k)
%%
[lb, ub] = ModelF.Bounds();

is_angle = any([ModelF.Index.PSFShape(3), ModelF.Index.GalaxyShape(3)] == k);
if is_angle
    lb = 0;
    ub = pi;
else
    lb = lb(k);
    ub = ub(k);
    
    if ~isfinite(lb)
        lb = min(x(:, k));
    end
    
    if ~isfinite(ub)
        ub = max(x(:, k));
    end
end

mu = x(:, k); 
a = squeeze(H(k, k, :));
ind = a > 0;

if is_angle
    f = @(x) mean((2*pi)^(-1/2) * a(ind).^(1/2) .* exp(-1/2 * circ_d(x - mu(ind), pi).^2.*a(ind)));
else
    f = @(x) mean((2*pi)^(-1/2) * a(ind).^(1/2) .* exp(-1/2 * (x - mu(ind)).^2.*a(ind)));
end

x = lb : (ub - lb) / 100 : ub;
y = zeros(size(x));

for i = 1 : length(x)
    y(i) = f(x(i));
end

if is_angle
    polar([x, pi + x], [y, y]);
else
    plot(x, y);
end