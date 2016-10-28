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

fcn=polyfit(t,m_f,9);
%% ODE45: RC part of W4P
% This part intentionally ignore flow of inertance
initial_P=65+(-5:0.02:5); %Initial Pao
for a=1:length(initial_P)
    init_y1=initial_P(a) - Rc_W4P*m_f(1);
    [t1 P_3WK]=ode45(@(t,P) wk3(t,P,m_f,Rp_W4P,C_W4P),t,init_y1);
 %[t1 P_3WK]=ode45(@(t,P) wk3(t,P,fcn,Rp_W4P,C_W4P),t,init_y1);
    Pao = P_3WK' + Rc_W4P*m_f;
    SC_3WK(a) = log(sum((m_p-Pao).^2))+log(N)*1;
    [min_SC_3WK ind_SC_3WK] = min(SC_3WK);
end
figure(3)
subplot(1,2,1),plot(initial_P, SC_3WK),xlabel('initial P_a_o(mmHg)')
title('SC with P_a_o(W4P)')


%% ODE45: Rest part of W4P
initial_I= 180+(-10:0.02:10);
for a=1:length(initial_I)
    init_y2=initial_P(ind_SC_3WK)-Rc_W4P*(m_f(1)-initial_I(a));%V1=Pao-Rc*(Iao-I2)
    [t1 P_W4P]=ode45(@(t,P) wk3(t,P,m_f,Rp_W4P,C_W4P),t,init_y2); %get V1 
    %[t1 P_W4P]=ode45(@(t,P) wk3(t,P,fcn,Rp_W4P,C_W4P),t,init_y2); %get V1 
    init_y2_inductor=m_f(1)-(m_p(1)-init_y2)/Rc_W4P;%I2=Iao-(Vao-V1)/Rc
    [t1 I2]= ode45(@(t,I) wk4_inductor(t,I,m_f,Rc_W4P,L_W4P),t,init_y2_inductor); %get I2
    %[t1 I2]= ode45(@(t,I) wk4_inductor(t,I,fcn,Rc_W4P,L_W4P),t,init_y2_inductor); %get I2
    Pao2=P_W4P'+(m_f-I2')*Rc_W4P;
    SC_W4P(a) = log(sum((m_p-Pao2).^2))+log(N)*1;
    [min_SC_W4P ind_SC_W4P] = min(SC_W4P);    
end
subplot(1,2,2),plot(initial_I,SC_W4P),xlabel('initial flow in inductor(ml/s)')
title('SC with current flow in L(W4P)')
init_y2=initial_P(ind_SC_3WK)-Rc_W4P*(m_f(1)-initial_I(ind_SC_W4P));
[t1 P_W4P]=ode45(@(t,P) wk3(t,P,m_f,Rp_W4P,C_W4P),t,init_y2);
%[t1 P_W4P]=ode45(@(t,P) wk3(t,P,fcn,Rp_W4P,C_W4P),t,init_y2);
init_y2_inductor=m_f(1)-(m_p(1)-init_y2)/Rc_W4P;%I2=Iao-(Vao-V1)/Rc
[t1 I2]= ode45(@(t,I) wk4_inductor(t,I,m_f,Rc_W4P,L_W4P),t,init_y2_inductor); %get I2
%[t1 I2]= ode45(@(t,I) wk4_inductor(t,I,fcn,Rc_W4P,L_W4P),t,init_y2_inductor); %get I2
Pao2=P_W4P'+(m_f-I2')*Rc_W4P;
figure(4)
plot(t1,Pao2,t,m_p),
title('pressure comparison between W4P and measured pressure')
xlabel('time(s)'),ylabel('Pressure(mmHg)'),legend('calculated pressure','measured pressure')

%% ODE45: W4S
initial_P=75+(-5:0.02:5); %Initial Pao
% extrapolation
m_temp=interp1(t,m_f,[t(1)-ts t t(end)+ts],'spline','extrap');
%m_temp=interp1(t,m_f,[t(1) t t(end)+ts+ts],'spline','extrap');
for a=2:length(m_temp)-1
   m_df(a-1)=(m_temp(a+1)-m_temp(a-1))/(2*ts); %dIao/dt
end
%m_df=diff(m_f)/ts; m_df=[m_df m_df(end)];%dIao/dt ->not perfect
%L_W4S= 9.04e-4;
for a=1:length(initial_P)
    init_y3=initial_P(a)-Rc_W4S*m_f(1)-L_W4S*m_df(1); %Pao-Rc*Iao - L*(dIao/dt)
    [t1 P_W4S] = ode45(@(t,P) wk3(t,P,m_f,Rp_W4S,C_W4S),t,init_y3);
   %[t1 P_W4S] = ode45(@(t,P) wk3(t,P,fcn,Rp_W4S,C_W4S),t,init_y3);
    Pao3=P_W4S'+Rc_W4S*m_f+ L_W4S*m_df;% Pao=P1+ Rc*f_ao
    SC_W4S(a)=log(sum((m_p-Pao3).^2))+log(N)*1;
end
[min_SC_W4S ind_SC_W4S]=min(SC_W4S);
init_y3=initial_P(ind_SC_W4S)-Rc_W4S*m_f(1)-L_W4S*m_df(1);
[t3 P_W4S] = ode45(@(t,P) wk3(t,P,m_f,Rp_W4S,C_W4S),t,init_y3);
%[t3 P_W4S] = ode45(@(t,P) wk3(t,P,fcn,Rp_W4S,C_W4S),t,init_y3);
Pao3=P_W4S'+Rc_W4S*m_f+ L_W4S*m_df;

figure(5)
subplot(1,2,1),plot(initial_P,SC_W4S),xlabel('initial P_a_o(mmHg)')
title('SC in W4S')
subplot(1,2,2),plot(t,Pao3,t,m_p)
title('pressure comparison between W4S and measured pressure')
xlabel('time(s)'),ylabel('Pressure(mmHg)'),legend('Calculated pressure','Measured pressure')

figure(6)
plot(t,m_p,t,Pao2,t,Pao3),legend('Measured','W4P','W4S')
title('Pressure comparison among 4-WK models and measured value')
xlabel('time(s)'),ylabel('Pressure(mmHg)')

