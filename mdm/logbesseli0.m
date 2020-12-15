function y = logbesseli0(x)
% y = log(besseli(0, x));
large = x > 10;
y(~large) = log(besseli(0, x(~large)));
y(large) = x(large) - log(2 * pi) / 2 - log(x(large)) / 2;
y = reshape(y, size(x));