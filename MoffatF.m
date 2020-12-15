function [A, dA_xi, dA_beta] = MoffatF(xi, beta)
    % FT of space-domain star intensity (PSF) profile
    % A = AMPLITUDE, NOT SQUARED!
    % xi = spatial frequency MAGNITUDE, NOT SQUARED!

    A = 2 * pi^(beta-1) * (beta-1) * xi.^(beta-1) .* besselk(1 - beta, 2 * pi * xi) / gamma(beta);

    if nargout > 1
        % Посчитать производную по радиусу
        dA_xi = -4 * pi^beta * xi.^(beta-1) .* besselk(-2 + beta, 2 * pi * xi) / gamma(beta-1);
    end

    if nargout > 2
        % Посчитать производную по beta
        dA_beta = pi^(beta-1) * xi.^(beta-1) / gamma(beta) ...
            .* (2 * besselk(beta-1, 2 * pi * xi) .* (1 + (beta-1) * log(pi * xi) - (beta-1) * psi(beta)) ...
            - 2 * (beta-1) * dbesselk1(1-beta, 2 * pi * xi));
    end
end