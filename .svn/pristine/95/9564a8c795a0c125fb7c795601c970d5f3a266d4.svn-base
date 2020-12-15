function C = fftconv2(A, B)
%
[M1 N1 K1] = size(A);
[M2 N2 K2] = size(B);

if isa(A, 'parallel.gpu.GPUArray') || isa(B, 'parallel.gpu.GPUArray')
    Aex = parallel.gpu.GPUArray.zeros([M1+M2 N1+N2 K1]);
    Bex = parallel.gpu.GPUArray.zeros([M1+M2 N1+N2 K2]);
else
    Aex = zeros([M1+M2 N1+N2 K1]);
    Bex = zeros([M1+M2 N1+N2 K2]);
end

Aex(1:M1, 1:N1, :) = A;
fftA = fft2(Aex);

Bex(1:M2, 1:N2, :) = B;
fftB = fft2(Bex);

if K2 == 1
    C = ifft2(fftA .* repmat(fftB, [1 1 K1]));
else
    C = ifft2(repmat(fftA, [1 1 K2]) .* fftB);
end

s = ceil(([M2 N2] + 1) / 2);
C = real(C(s(1) : s(1) + M1 - 1, s(2) : s(2) + N1 - 1, :));