% Plot of Mz for IR Question.

t = [0:0.01:1];
Mz = -1.63*exp(-t/0.5)+1;
Mz(51:101) = 1-exp(-(t(51:101)-.5)/.5);
plot(t,Mz)
lplot('Time (s)','M_z / M_0','Simple IR Signal Example',[0 1 -1 1]);
setprops
