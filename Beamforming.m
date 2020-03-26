% Performance evaluation of different beamforming algorithms

clc;close all;clear all;

%% Model 

N=2; % # of antennas
L=5000; % # of symbols
EbN0_dB=10;
alpha=1; % signal level of x2 relative to x1
a1=steer(N,0.5,2); a2=steer(N,0.5,45); % Steering Vectors
s1=((randi([0 1],1,L)*2-1) +1j*(randi([0 1],1,L)*2-1))/sqrt(2); %Signal 1
s2=((randi([0 1],1,L)*2-1) +1j*(randi([0 1],1,L)*2-1))/sqrt(2); %Signal 2

x1 = a1*s1; x2 = a2*s2 * alpha; % Received signal vectors
n = (randn(N,L)+1j*randn(N,L))/sqrt(2)* 10^(-EbN0_dB/20); % noise
x = x1+x2+n;

sigma_d=s1*s1';
Rxx=x*x';
Rjj=x2*x2';
sigma_n=n*n';
Rjn=((a2*a2')*Rjj)+sigma_n*eye(N);

%% Zero Forcing 

w = [eye(N) - inv(a2'*a2)*a2*a2']'*a1; % Weight Vector

y1 = w'*x1; y2 = w'*x2; yn = w'*n;
P1 = sum(y1*y1')/L; P2 = sum(y2*y2')/L;
Pn = sum(yn*yn')/L; SINR = P1/(P2+Pn);
fprintf('Output SINR (ZF): %8.4f dB\n',10*log10(SINR));
angle_range=[-90:90];
for ii=1:length(angle_range)
a = steer(N,0.5,angle_range(ii));
p(ii) = abs(w'*a);
end
figure(1)
plot(angle_range,p/max(p))
xlabel('Angle of arrival (degree)')
ylabel('Normalized array pattern')
title('Zero Forcing');
ylim([0 1])

%% MMSE

% w = inv(x*x'/L)*(x'*s1'); % MMSE(unnormalized)
w = (sigma_d*inv(Rxx))*a1; %MMSE

y1 = w'*x1; y2 = w'*x2; yn = w'*n;
P1 = sum(y1*y1')/L; P2 = sum(y2*y2')/L;
Pn = sum(yn*yn')/L; SINR = P1/(P2+Pn);
fprintf('Output SINR (MMSE): %8.4f dB\n',10*log10(SINR));
angle_range=[-90:90];
for ii=1:length(angle_range)
a = steer(N,0.5,angle_range(ii));
p(ii) = abs(w'*a);
end
figure(2)
plot(angle_range,p/max(p))
xlabel('Angle of arrival (degree)')
ylabel('Normalized array pattern')
title('MMSE');
ylim([0 1])

%% MPDR

w = (inv(Rxx)*a1)/(a1'*(inv(Rxx)*a1)); % Weight Vector

y1 = w'*x1; y2 = w'*x2; yn = w'*n;
P1 = sum(y1*y1')/L; P2 = sum(y2*y2')/L;
Pn = sum(yn*yn')/L; SINR = P1/(P2+Pn);
fprintf('Output SINR (MPDR): %8.4f dB\n',10*log10(SINR));
angle_range=[-90:90];
for ii=1:length(angle_range)
a = steer(N,0.5,angle_range(ii));
p(ii) = abs(w'*a);
end
figure(3)
plot(angle_range,p/max(p))
xlabel('Angle of arrival (degree)')
ylabel('Normalized array pattern')
title('MPDR');
ylim([0 1])

%% MVDR

w = (inv(Rjn)*a1)/(a1'*(inv(Rjn)*a1)); % Weight Vector

y1 = w'*x1; y2 = w'*x2; yn = w'*n;
P1 = sum(y1*y1')/L; P2 = sum(y2*y2')/L;
Pn = sum(yn*yn')/L; SINR = P1/(P2+Pn);
fprintf('Output SINR (MVDR): %8.4f dB\n',10*log10(SINR));
angle_range=[-90:90];
for ii=1:length(angle_range)
a = steer(N,0.5,angle_range(ii));
p(ii) = abs(w'*a);
end
figure(4)
plot(angle_range,p/max(p))
xlabel('Angle of arrival (degree)')
ylabel('Normalized array pattern')
title('MVDR');
ylim([0 1])

