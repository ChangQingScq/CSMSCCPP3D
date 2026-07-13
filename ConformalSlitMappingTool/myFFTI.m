function [integ_Y] = myFFTI(X)
    Y = fft(X);
    N = length(Y);
    Y1 = Y./[1i*(0:N/2), 1i*(-N/2+1:-1)];
    Y1(1) = 0;
    integ_Y = (ifft(Y1));
    k = mean(X);
    Timelist = 0:pi*2/N:pi*2-pi*2/N;
    integ_Y = integ_Y + k*Timelist;
end
