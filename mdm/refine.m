function x1 = refine(x0, H, mask)
%
N = size(x0, 1);
x1 = zeros(size(x0));

for i = 1 : N
    c = zeros(size(x0));
    log_w = zeros(1, N);

    Hi = H(:, :, i);
    xi = x0(i, :);
    
    for j = 1 : N
        Hj = H(:, :, j);
        Hw = ones(size(Hj));
        Hw(~mask, :) = 0;
        Hw(:, ~mask) = 0;
        Hj = Hj .* Hw;
        
        xj = x0(j, :);
        
        c(j, :) = (xi * Hi + xj * Hj) / (Hi + Hj);
        b = xi * Hi * xi.' + xj * Hj * xj.' - ((xi * Hi + xj * Hj) / (Hi + Hj)) * (xi * Hi + xj * Hj).';
        log_w(j) = (log(det(Hj(mask, mask))) - log(det(Hi + Hj)) - b) / 2;
    end
    
    w = exp(log_w - max(log_w));
    x1(i, :) = w * x0 / sum(w);
    ws = sort(w, 'descend');
    
    disp(strcat('i = ', mat2str(i), ', w = ', mat2str(ws(1 : 3))));
end