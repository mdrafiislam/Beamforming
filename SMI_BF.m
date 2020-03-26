%SMI beamformer with diagonal loading

clc;close all;clear all;

%% Model

N=6; % # of antennas
L=5000; % # of symbols
EbN0_dB=10;
alpha=1; % signal level of x2 relative to x1
a1=steer(N,0.5,0); a2=steer(N,0.5,45); % Steering Vectors
s1=((randi([0 1],1,L)*2-1) +1j*(randi([0 1],1,L)*2-1))/sqrt(2); %Signal 1
s2=((randi([0 1],1,L)*2-1) +1j*(randi([0 1],1,L)*2-1))/sqrt(2); %Signal 2

x1 = a1*s1; x2 = a2*s2 * alpha; % Received signal vectors
n = (randn(N,L)+1j*randn(N,L))/sqrt(2)* 10^(-EbN0_dB/20); % noise
x = x1+x2+n;

Rxx_avg = (x*x')/N; %sample covariance matrix

%% SMI beamformer without diagonal loading

w = (inv(Rxx_avg)*a1)/(a1'*(inv(Rxx_avg)*a1)); 
y1 = w'*x1; y2 = w'*x2; yn = w'*n;
P1 = sum(y1*y1')/L; P2 = sum(y2*y2')/L;
Pn = sum(yn*yn')/L; SINR = P1/(P2+Pn);
fprintf('Output SINR (without DL): %8.4f dB\n',(SINR));
angle_range=[-90:90];
for ii=1:length(angle_range)
a = steer(N,0.5,angle_range(ii));
p(ii) = abs(w'*a);
end
figure(1)
plot(angle_range,p/max(p))
title('Array Pattern')
xlabel('Angle of arrival (degree)')
ylabel('Normalized array pattern')
title('SMI beamformer without diagonal loading')
ylim([0 1])
%% SMI beamformer with diagonal loading

DL_factor =10*Pn; %Diagonal loading factor
Rxx_DL = Rxx_avg + (DL_factor*eye(N));%Diagonally loaded(DL)sample covariance matrix
w = (inv(Rxx_DL)*a1)/(a1'*(inv(Rxx_DL)*a1)); %SMI with DL
y1 = w'*x1; y2 = w'*x2; yn = w'*n;
P1 = sum(y1*y1')/L; P2 = sum(y2*y2')/L;
Pn = sum(yn*yn')/L; SINR = P1/(P2+Pn);
fprintf('Output SINR (with DL): %8.4f dB\n',(SINR));

angle_range=[-90:90];
for ii=1:length(angle_range)
a = steer(N,0.5,angle_range(ii));
p(ii) = abs(w'*a);
end
figure(2)
plot(angle_range,p/max(p))
title('Array Pattern')
xlabel('Angle of arrival (degree)')
ylabel('Normalized array pattern')
title('SMI beamformer with diagonal loading')
ylim([0 1])
