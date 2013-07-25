function Idot=wk4_inductor(t,I,flow,Rc,L)
N=length(flow);
t_array=linspace(0,0.8,N);
diff_t=abs(t-t_array); 
[min_t ind_t]=min(diff_t);
Idot=Rc/L*(flow(ind_t)-I);
end


