clear all
[X,Fs] = audioread('70 alvin.wav');
amax=0.5;
amin=15;
L=length(X);
wp=2*pi*400/Fs;
ws=2*pi*600/Fs;
[n,wp]=cheb1ord(wp/pi,ws/pi,amax,amin);
[b,a]=cheby1(n,amax,wp);
[y,W] = filter(b,a,X);
nVec=1:L;
f=fft(y);
f2=fft(X);

figure(1)
plot(nVec/Fs,X)

figure(2)
plot(nVec/Fs,y)

P2 = abs(f/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Z = Fs*(0:(L/2))/L;
figure(3)
plot(Z,P1)
xlim([0 1000])
title('Single-Sided Amplitude Spectrum of Y(t)')
xlabel('f (Hz)')
ylabel('|P1(f2)|')

P2 = abs(f2/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Z = Fs*(0:(L/2))/L;
figure(4)
plot(Z,P1)
xlim([0 1000])
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f2)|')

figure(5)
plot(nVec/Fs,(X-y))
ylim([-0.1 0.1])
W


