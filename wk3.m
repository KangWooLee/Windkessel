function Pdot=wk3(t,P,flow,Rp,C)
N=length(flow);
t_array=linspace(0,0.8,N);
diff_t=abs(t-t_array); 
[min_t ind_t]=min(diff_t);
Pdot=flow(ind_t)/C - P/(Rp*C);
end