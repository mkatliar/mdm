function c = dconv2(dx, y, x, dy)
%
c = fftconv2(dx, y) + fftconv2(x, dy);