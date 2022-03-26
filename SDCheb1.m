clear all

Fs=48000;
amax=0.5;
amin=15;
wp=2*pi*400/Fs;
ws=2*pi*600/Fs;
[n,wp]=cheb1ord(wp/pi,ws/pi,amax,amin);
[b,a]=cheby1(n,amax,wp);
%max atten passband: 0.5dB
%min atten rejectband: 15dB

disp('Order of H(z):')
disp(n)
disp('Numerator of H(z)')
disp(b);
disp('Denominator of H(z)')
disp(a);

%IMPORTANT: Verify that design satisfy your specifications
[HVec,wVec] = freqz(b,a,1000);  %<= freqz evaluates H(z) at z=exp(jw) at 100 uniformly spaced w's

figure(2);
subplot(2,1,1);plot(wVec/(pi),abs(HVec),'-r','linewidth',2)
xlabel('w/pi','fontsize',20)
ylabel('|H(exp(jw))|','fontsize',20);
axis([0 0.7 0 1])
grid on
subplot(2,1,2);plot(wVec/(pi),20*log10(abs(HVec)),'-r','linewidth',2)
axis([0 0.7 -20 1])
grid on
xlabel('w/pi','fontsize',20)
ylabel('|H(exp(jw))| (dB)','fontsize',20);
set(gca,'fontsize',20)


%IMPORTANT: Verify that design satisfy your specifications at wp and ws
[HVec,wVec] = freqz(b,a,[wp, ws]);  %<= you can specify which frequencies to use in 
disp(sprintf('|H(exp(jwp)|=%f',20*log10(abs(HVec(1)))));
disp(sprintf('|H(exp(jws)|=%f',20*log10(abs(HVec(2)))));

