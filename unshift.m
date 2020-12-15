function s = unshift(s)
%
F = fft2(s);
W = F .* conj(F);
A = abs(F);

s = ifft2(A);
s(1, 1) = NaN;
s = ifftshift(s);