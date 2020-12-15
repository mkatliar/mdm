function [A, dA_xi, dA_theta] = SersicF(xi, ~)
    % FT of space-domain galaxy intensity profile.
    % A = AMPLITUDE, NOT SQUARED!
    % xi = spatial frequency MAGNITUDE, NOT SQUARED!

    A = 1 ./ (1 + 4 * pi^2 * xi.^2).^(3 / 2);

    if nargout > 1
        dA_xi = -12 * pi^2 * xi ./ (1 + 4 * pi^2 * xi.^2).^(5 / 2);
    end                    

    if nargout > 2
        % Производная по параметрам = [], т.к. сейчас
        % SersicIndex фиксирован и = 1.
        dA_theta = [];
    end 
end