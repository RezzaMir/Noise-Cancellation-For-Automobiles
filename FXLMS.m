%**************************************************************************
%  FXLMS.m - Matlab program simulation of single-channel ANC system       *
%            using the FXLMS algorithm                                    *
%**************************************************************************
%  System configuration:
%
%  Initialization -> Off-line secondary path modeling -> On-line ANC
%--------------------------------------------------------------------------
%  Off-line modeling:
%                           _______
%                          |       |       d(n)
%  x(n) ------------------>|  S(z) |--------|
%              |           |_______|        |
%              |                 ^          |
%              |            ____/__         |
%              |           |       | y(n)   v +
%              |---------->| S^(z) |----->(sum)--------> e(n)
%              |           |_______|     -        |
%              |            ___|___               |
%              |           |       |              |
%              |---------->|  LMS  |<-------------|
%                          |_______|
%
%     x(n) is an internally-generated white noise
%     S(z) is the secondary path represented by an IIR filter
%     d(n) is the output of S(z)
%     S^(z) is an adaptive filter that models the secondary path
%     y(n) is the output of S^(z)
%     e(n) = d(n) - y(n) is used by the LMS algorithm to update S^(z)
%--------------------------------------------------------------------------
%  On-line ANC:
%                   ________
%                  |        |              d(n)
%   x(n) --------->|  P(z)  |-----------------------
%            |     |________|                      |
%            |           ^                         |
%            |      ____/___         ______        |
%            |     |        | y(n)  |      |y'(n)  v +      e(n)
%            |---->|  W(z)  |------>| S(z) |---->(sum)----------->
%          __|__   |________|       |______|    -         |
%         |     |      |                                  |
%         |S^(z)|   ___|___                               |
%         |_____|  |       |                              |
%            |---->|  LMS  |<-----------------------------|
%           x'(n)  |_______|
%
%     x(n) is the reference signal
%     P(z) is the primary path represented by an IIR filter
%     d(n) is the output of P(z)
%     W(z) is an adaptive filter updated by the FXLMS algorithm
%     y(n) is the output of W(z)
%     S(z) is the secondary path represented by an IIR filter
%     y'(n) is the output of S(z)
%     e(n) = d(n) - y'(n) is the error signal
%     S^(z) is a fixed FIR filter from the off-line modeling mode
%     x'(n) is a filtered reference signal used by the FXLMS algorithm
%--------------------------------------------------------------------------
%  Synopsis
%
%      These parameter can be defined by user.
%      LS  - order of S^(z) (integer)
%      MUS - step size of S^(z) (real)
%      LW  - order of W(z) (integer)
%      MUW - step size of W(z) (real)
%
%--------------------------------------------------------------------------
%  Files
%    
%    "TF.mat" contains transfer function of P(z) and S(z)
%
%    P_z - coefficients of numerator of P(z)
%    P_p - coefficients of denominator of P(z)
%              (The leading coefficient is assumed to be equal to 1)
%    S_z - coefficients of numerator of S(z)
%    S_p - coefficients of denominator of S(z)
%              (The leading coefficient is assumed to be equal to 1)
%
%**************************************************************************
%*    Clear workspace and command window                                  *
%**************************************************************************
clear all 
close all
clc

LS=128;        % order of S^(z) 
MUS=0.005;     % step size of S^(z) 
LW=256;        % order of W(z) 
MUW=e-10;     % step size of W(z) 
% Import Transfer function of primary path P(z), secondary path S(z),
load('TF');
x=importdata('SEC13R.mat');    % Input x(n) from data file
Fs=4800;     % Fs for 'SEC13R.mat'
% Fs=10480;  % Fs for 'SEC18R.mat'

Spec=fft(x,1024);
figure 
plot(Fs/1024:Fs/1024:512*Fs/1024,20*log10(abs(Spec(1:512))));
xlabel('freq (Hz)');
ylabel('dB');
title('Power Spectra of X(z)');

[Hp wp]=freqz(P_z,P_p);
[Hs ws]=freqz(S_z,S_p);
figure
subplot(211)
plot(pi/512:pi/512:pi,abs(Hp));
xlabel('freq');
ylabel('Magnitude');
title('Spectrum of P(z)');
subplot(212)
plot(pi/512:pi/512:pi,abs(Hs));
xlabel('freq');
ylabel('Magnitude');
title('Spectrum of S(z)');
%%
%**************************************************************************
%*    Off-line modeling of secondary path S^(z)                           *
%**************************************************************************
% 1. Generate sampled zero-mean white noise signal wh_n 
count = 30000;
wh_n = rand(1,count)-0.5;

% 2. Obtaining desired signal d(n) 
d_off= filter(S_z,S_p,wh_n); % FIR filtering of white noise by S(z)

% 3. LMS algorithm 
S_hat = zeros(1,LS); % Define coefficient vector of secondary path S^(z) 
wh_n_bf=zeros(1,LS); % Define white noise buffer
a=1/LS;
P0=std(d_off(1))^2;
for n=1:count
    % Calculate off-line modeling error 
    wh_n_bf=[wh_n(n) wh_n_bf(1,1:LS-1)]; % update white noise buffer
    y_off(n) = wh_n_bf*S_hat'; % filtering wh_n(n) by S^(z) to get y(n)
    e_off(n)=d_off(n)-y_off(n); % e(n) = d(n) - y(n)

    % Use off-line modeling stepsize MUS to update S^(z)
    S_hat = S_hat+MUS*e_off(n)*wh_n_bf;
    
    if n==1
        Pe(n)=(1-a)*P0+a*e_off(n)^2;
    else
        Pe(n)=(1-a)*Pe(n-1)+a*e_off(n)^2;
    end
end

% 4. Off-line modeling result
figure
subplot(211);
plot(1/Fs:1/Fs:length(Pe)/Fs,10*log10(Pe));
xlabel('time (s)');
ylabel('dB');
title('Power of e_off(n)');

subplot(212);
plot(d_off);
hold on;
plot(e_off,'r');
title('Off-line error plot');
xlabel('time (s)');
ylabel('Amplitude');
legend('white noise','error');

%%
%**************************************************************************
%*    On-line ANC                                                         *
%**************************************************************************
% 1.IIR filtering of x(n) by P(z) to get d(n)
d=filter(P_z,P_p,x)';    

% 2. Filted-X LMS algorithm 
W=zeros(1,LW); % Define coefficient vector of adaptive filter W(z)
xW_bf=zeros(1,LW); % Define buffer to filter x(n) by W(z)
yz_bf=zeros(1,length(S_z)); % Define buffer to FIR filtering y(n) by S_z
yp_bf=zeros(1,length(S_p)-1); % Define buffer to FIR filtering yp(n-1) by S_p
xS_bf=zeros(1,LS); % Define buffer to FIR filtering x(n) by S^(z) to get x'(n)
xp_bf=zeros(1,LW); % Define buffer to update W vector

for n = 1 : length(x);
    xW_bf=[x(n) xW_bf(1,1:LW-1)];
    y(n) = xW_bf * W'; % FIR filtering x(n) by W(z) to get y(n)
    
    yz_bf=[y(n) yz_bf(1:length(S_z)-1)];
    yp(n)=yz_bf*S_z-yp_bf*S_p(2:end); % y'(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
                                      %         - a(2)*y'(n-1) - ... - a(na+1)*y'(n-na)
    yp_bf=[yp(n) yp_bf(1:length(yp_bf)-1)]; % update y'(n) buffer
    
    e(n)=d(n)-yp(n);
    
    % Filtered x(n) by S^(z) to get x'(n)
    xS_bf=[x(n) xS_bf(1,1:LS-1)];
    xp(n)=xS_bf*S_hat';
    
    % Use on-line modeling stepsize MUW to update W(z)
    xp_bf=[xp(n),xp_bf(1,1:LW-1)];
    W=W+MUW*e(n)*xp_bf;
end

figure
subplot(211)
plot(d)
hold on    
plot(e,'r')   
title('ANC result');    
xlabel('time (s)');
ylabel('Amplitude');
legend('d','e'); 
hold off

Spec_e=fft(e(end-1024:end),1024);
Spec_d=fft(d(end-1024:end),1024);

subplot(212)
plot(Fs/1024:Fs/1024:length(Spec_d)/2*Fs/1024,20*log10(abs(Spec_d(1:512))));
hold on
plot(Fs/1024:Fs/1024:length(Spec_e)/2*Fs/1024,20*log10(abs(Spec_e(1:512))),'r');
title('Power Spectrum of last 1024 samples')
xlabel('Freq (Hz)');
ylabel('dB');
legend('d','e'); 

disp(' Comput average cancelation dB for last 1024 samples')
mean(20*log10(abs(Spec_d(1:512)))-20*log10(abs(Spec_e(1:512))))

% 
H=Hp/Hs;
figure
subplot(211)
plot(abs(H))
title('P(z)/S(z)');

[Hw ww]=freqz(W,1);
subplot(212)
plot(abs(Hw));
title('W(z)');
