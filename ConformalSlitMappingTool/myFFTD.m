function [deri_Y] = myFFTD(X)
Y = fft(X);
N = length(Y);
% 计算数值导数
deri_Y = ifft(Y.*[1i*(0:N/2), 1i*(-N/2+1:-1)]);
