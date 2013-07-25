clc;clear all;close all;
%% Initialization
tmax=0.8; fmax=10;
order=7;smoothing=61;

%% Measured data processing
load measured_data.mat
m_f0=m_f;m_p0=m_p;%original value
% Pick partial data points. if interval is 1, we pick all of them.
m_f=m_f(1:1:end);m_p=m_p(1:1:end); 
% Savitzky-Golay filtering. 
%second arg: polynomial order
%third arg: number of data points used for smoothing
m_f=sgolayfilt(m_f,order,smoothing); %3,21
m_p=sgolayfilt(m_p,order,smoothing);

N=length(m_f);
t=linspace(0,0.8,N); ts=t(2)-t(1);
fs=1/ts;
f=(0:N-1)*fs/N; jw=j*2*pi*f;
M_F=fft(m_f);M_P=fft(m_p);IMP_measured=M_P./M_F;

% f_LP_ind=find(f<fmax);%get index below 15Hz, basically 1~13 index
% f_LP=f(f_LP_ind);
% M_F_LP=M_F(f_LP_ind);M_P_LP=M_P(f_LP_ind); IMP_measured_LP = M_P_LP./M_F_LP;
% fcn_IMP=polyfit(f_LP,IMP_measured_LP,7);
% ff=linspace(f_LP(1),f_LP(end),30);
% IMP_1=polyval(fcn_IMP,ff);f_LP=ff;
% Substitute to lower values
%IMP_measured=IMP_1;f=f_LP;


figure(1)
subplot(2,2,1),plot(t,m_f),title('flow'),xlabel('t'),xlim([0, tmax])
subplot(2,2,2),plot(t,m_p),title('Pressure'),xlabel('t'),xlim([0, tmax])
subplot(2,2,3),plot(f,20*log10(abs(IMP_measured))),xlim([0,fmax])
%subplot(2,2,3),plot(f_LP,20*log10(abs(IMP_1)),'.'),xlim([0,fmax])
title('magnitude'),xlabel('f'),ylabel('Magnitude Modulus(dB)')
subplot(2,2,4),plot(f,angle(IMP_measured)*180/pi),xlim([0,fmax])
%subplot(2,2,4),plot(f_LP,angle(IMP_1)*180/pi,'.'),xlim([0,fmax])
title('phase'),xlabel('f'),ylabel('Phase(degree)'),ylim([-90 90])

%% Predicted data processing
Rp_W4P=0.63;Rc_W4P=0.045;C_W4P=2.53;L_W4P=0.0054; %W4P
Rp_W4S=0.79;Rc_W4S=0.0288;C_W4S= 3.525;L_W4S= 9.64e-4; %W4S
% R=mmHg/ml*s, C=ml/mmHg, L=mmHg*s^2/ml

%% FFT
%Measured
M_F=fft(m_f); M_P=fft(m_p); M_R=M_P./M_F;
% W4P
W4P_R=1/C_W4P./(jw+1/(C_W4P*Rp_W4P))+jw*Rc_W4P./(Rc_W4P/L_W4P+jw);
%C_W4P2=0.001;
%W4P_R=1/C_W4P./(jw+1/(C_W4P*Rp_W4P))+1./(1/Rc_W4P+jw*C_W4P2+1./(jw*L_W4P)
%);%test
W4P_P= M_F.*W4P_R;
% W4S
W4S_R=1/C_W4S./(jw+1/(C_W4S*Rp_W4S))+Rc_W4S+L_W4S*jw;
W4S_P=M_F.*W4S_R;

figure(2)
subplot(2,2,1),plot(f,10*log10(abs(M_R)),...
    f,10*log10(abs(W4P_R)),'g',f,10*log10(abs(W4S_R)),'r'),
xlim([0 fmax]),xlabel('f(Hz)'),title('Magnitude in dB(Impedance)')
legend('Measured','W4P','W4S')
subplot(2,2,2),plot(f,angle(M_R)*180/pi,...
    f,angle(W4P_R)*180/pi,'g',f,angle(W4S_R)*180/pi,'r'),
xlim([0 fmax]),xlabel('f(Hz)'),title('Phase(Impedance)'),%ylim([-90 90])
legend('Measured','W4P','W4S')
subplot(2,2,3),plot(f,10*log10(abs(M_P)),...
    f,10*log10(abs(W4P_P)),'g',f,10*log10(abs(W4S_P)),'r'),
xlim([0 fmax]),xlabel('f(Hz)'),title('Magnitude in dB(Pressure)')
legend('Measured','W4P','W4S')
subplot(2,2,4),plot(f,angle(M_P)*180/pi,...
    f,angle(W4P_P)*180/pi,'g',f,angle(W4S_P)*180/pi,'r'),
xlim([0 fmax]),xlabel('f(Hz)'),title('Phase(Pressure)'),%ylim([-90 90])
legend('Measured','W4P','W4S')