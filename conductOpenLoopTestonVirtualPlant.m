clear
close all
clc


modelName = 'openLoopTestBed';

Ts=0.01; %sampling period, in seconds (rather, time units)
tfinal = 100;
timeInput = [0:Ts:tfinal]';


%%---
% construct the input structure
uValues = 20*sin(3*timeInput);

u = [];
u.time = timeInput;
u.signals.values = uValues;
u.signals.dimensions = 1;

%%- run simulation
Out = sim(modelName,'StopTime',num2str(tfinal));
%%--


%% plot the results

timeOut = Out.y.time;
yValues = Out.y.signals.values;

figure
plot(timeInput,uValues,'b*--',timeOut,yValues,'r*-');
xlabel('time');ylabel('y and u-command');
title('Beware: u is the u commanded, it might be getting saturated');
legend('u','y');











