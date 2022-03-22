clear all
close all
clc
%% Plant Data
modelName = 'closedLoop_DesignTest';
A = [-0.5355 -25.1019 -6.2361;1 0 0; 0 1 0];
B = [1 0 0]';
C = [0.1336 -2.2550 52.2390];
D = 0;
[b,a] = ss2tf(A,B,C,D)
Phat = tf(b,a);
sys = ss(A,B,C,D);
G = eye(3);
n = size(A,1); %number of components in x
m = size(B,2); %number of components in u
p = size(G,1); %number of components in z (= Gx)
xinit = zeros(n,m);

%% reference signal(desired set point)
ref_type = 'step_zero21_at5_back20_at15';
t = [0:0.1:30]';
r = F_ref_at_t(t,ref_type);
ref = [t,r];

%% --- get equilibirum values of x and u
M = [A, B; C, D];
zm = [zeros(3,length(r));r'];
xstar = M\zm;
ustar =xstar(n+1,1:end);
ustar_t = [t,ustar'];

%% System parameters that simulink model uses   
%------- tuning knobs---- (in fact, gamma is one too)
q1 = 1;
r1 = .9;
R = r1*eye(m);
Q1 = q1*eye(p);
Q = G'*Q1*G;

%% Pole Calculation
[K,P,ctrl_poles] = lqr(A,B,Q,R,zeros(n,m));
obsv_poles = 8*ctrl_poles; %choose eigenvalues of A-LC
%compute L by using 'place' command here
L = place(A',C',obsv_poles);
L = L';

%% define Observer parameters
A_obs = A - L*C;
B_obs = [L B - (L*D)];
C_obs = eye(size(A_obs));
D_obs = zeros(size(B_obs));
xhat_init = zeros(n,1);

%% call simulink file
simOut = sim(modelName);
uout =simOut.uout;
yout = simOut.yout;
xhatout = simOut.xhatout;
tout = simOut.tout;

%% plot data
figure
plot(yout,'r--');
hold on
plot(uout,'b');
legend('Output','Input');
ylabel('Output / Input'), xlabel('Time(t)');
