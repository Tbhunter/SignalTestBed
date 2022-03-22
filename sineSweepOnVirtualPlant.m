clear
close all
clc

modelName = 'openLoopTestBed';
conductSineSweepExp = 1;
A_u = 5;
j = sqrt(-1)
numFreqRes = 2^11;

%% user inputs (for sine sweep only)
if conductSineSweepExp==1
    %-------------------------
    omega_array = [ 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10]; %it needs to be a row or column vecor of frequencies in 
    % "radians/time unit"
    %--------------------------
    
    %--------------------------
    decay_time = [1]; %in "time units"
    %--------------------------
    
    %--------------------------
    num_cycles2average_nom = [20]; %needs to be a positive integer
    %--------------------------
    
    Fs = [10]; %sampling period in "samples/time unit"
    % choose it large enough to give you enough samples in one period of the sinusoid
    % (depends on the highest frequency you use, of course), but not so
    % high that the simulation takes foreever.
    Ts = 1/Fs;

end

%% ----   SINE SWEEP ---
if conductSineSweepExp==1
    g_hat_array = nan*ones(length(omega_array),1);
    theta_hat_array = nan*ones(length(omega_array),1);
    
    for omega_index = 1:length(omega_array)
        omega = omega_array(omega_index);
        
        num_cycles2average = num_cycles2average_nom+5*ceil(omega);%this increases 
        % the value of N for higher frequencies. 
        
        t_final = decay_time+(2*pi/omega)*num_cycles2average;        
        time = [0:Ts:t_final]';
        uValues = A_u*sin(omega*time);
       
        u = [];
        u.time = time;
        u.signals.values = uValues;
        u.signals.dimensions = 1;
        
        %% do experiment to collect data
        Out = sim(modelName,'StopTime',num2str(t_final));
        timeOut = Out.y.time;
        y = Out.y.signals.values;
        
        %cut the transient period out, so that the data is consistent with the
        %theory of sine in sine out
        inds2average = [ceil(decay_time/Ts):1:length(time)]';
        N = length(inds2average); %number of samples to average over
        
        cosine_vector = cos(omega*time);
        sine_vector = sin(omega*time);
        ZcN = y(inds2average)'*cosine_vector(inds2average);
        ZsN = y(inds2average)'*sine_vector(inds2average);
        g_hat_omega = 2/A_u/N*sqrt(ZcN^2+ZsN^2); %gain est
        theta_hat_omega = atan2(ZcN,ZsN); %phase est, in rad
        
        %save estimates
        g_hat_array(omega_index) = g_hat_omega;
        theta_hat_array(omega_index) = theta_hat_omega;
        
        disp(['done with freq = ',num2str(omega),' rad/sec']);
    end
    
    %% compute the true magnitude and phase at a large number of frequencies: 
    
     w         = logspace(-2,3,1000);
     Mag       = g_hat_array;
     phase     = theta_hat_array;
     H         = Mag.*(exp(j*phase));  
    [B,A]      = invfreqs(H,omega_array,2,3);    
     G         = tf(B,A);
    [a,b,c,d] = tf2ss(B,A);
     e         = eig(a);
    [G_fit]    = freqresp(G,w);
     G_fit     = G_fit(:);
        
    %% superimpose the estimate and the true frequency response on the same Bode plot
    bodefig = figure
    ax1 = axes('position',[0.1300 0.55 0.7750 0.4])
    semilogx(w,20*log10(abs(G_fit)),'b-');
    hold on;
    semilogx(omega_array,20*log10(g_hat_array),'ro');
    ylabel('gain, dB');
    legend('true','est');
    
    ax2 = axes('position',[0.1300 0.1 0.7750 0.4])
    semilogx(w,angle(G_fit)*180/pi,'b-');
    hold on;
    semilogx(omega_array,theta_hat_array*180/pi,'ro');
    xlabel('\omega (rad/sec)');
    ylabel('Phase, degree');
    legend('true','est');
    
    
end
