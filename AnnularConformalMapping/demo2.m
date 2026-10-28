clear
close all
clc
[v, f] = generatingMeshForSlitMapping();

inner_bd_id = 2; %inner_bd_id should not be outer_bd_id
[tau] = annularSlitMapping(v, f, inner_bd_id);
figure()
hold on
trisurf(f, real(exp(tau)),imag(exp(tau)), 0*imag(tau))
figure(1)
bd = meshboundaries(f);
hold on
myPlot(v(bd{1},:),'r')
myPlot(v(bd{2},:),'g')
myPlot(v(bd{3},:),'b')
myPlot(v(bd{4},:),'y')

figure(2)
tau_v = [real(exp(tau)),imag(exp(tau))];
bd = meshboundaries(f);
hold on
myPlot(tau_v(bd{1},:),'r')
myPlot(tau_v(bd{2},:),'g')
myPlot(tau_v(bd{3},:),'b')
myPlot(tau_v(bd{4},:),'y')