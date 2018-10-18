function fftabs = plotfft(x, Rs, slope, name, ifplot)


c = 3e8;
fftx = fft(x);
[L, ~] = size(x);
dist = Rs*(0:(L/2))/L * c/slope/2;
P1 = abs(fftx(1:L/2+1,:)/L);
P1(2:end-1,:) = 2*P1(2:end-1,:);

if ifplot
    figure
    plot(dist,P1(:,ifplot));
    title(sprintf('Single-Sided Amplitude Spectrum of %s', name));
    xlabel('d (m)')
    ylabel('|P1(f)|')
end

fftabs = P1;

end