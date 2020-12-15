function [x1, e] = refine1(x0, H)
%
N = size(x0, 1);
x1 = x0;
e = zeros(N, 2);

ind = 1 : ModelF.ParameterVectorLength;
ind([ModelF.Index.PSFShape(3), ModelF.Index.GalaxyShape(3)]) = [];

C = inv(cov(x0(:, ind)));
mu = mean(x0(:, ind));

for i = 1 : N
    Hi = H(ind, ind, i);
    xi = x0(i, ind);
    
    x1(i, ind) = (xi * Hi + mu * C) / (Hi + C);

    r = x1(i, ModelF.Index.GalaxyShape(2));
    theta = x1(i, ModelF.Index.GalaxyShape(3));
            
    e(i, :) = (1 - r) / (1 + r) * [cos(2 * theta), -sin(2 * theta)];
    
    disp(strcat('i = ', mat2str(i)));
end