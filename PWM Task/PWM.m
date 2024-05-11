clc;
clear all;
inputStr = input('Enter the value of f0: ', 's');                                %input f0 from user
f0 = str2double(inputStr);
inputStr = input('Enter the two values of fs in the format fs=[...]: ', 's');      %input fs from user
Switching_frequency = str2num(inputStr);
%%inputStr = input('Enter the Mode:(0-->Symmetric, 1-->Asymmetric): ', 's');       %input the mode from user
%%Mode = str2double(inputStr);
matrix_str = input('Enter the modulation index in the format M=[...]: ', 's');   %input modulation index
Mi = str2num(matrix_str); 
model_name = 'Task1_version5';
simulation_time = 0.05;                                                          %Specify the simulation time in seconds
open_system(model_name);
%%Mi =[0.1 0.3 0.5 0.6 0.8 0.9 1 1.1 1.3];                                       %Modulation index
for s=1:2
    if s==1
    Mode=0;
for j=1:2
    fs=Switching_frequency(j);
    
    if j==1
for i = 1:length(Mi)                                                             %for loop to run every mudulation index
    M = Mi(i);
   simOut = sim(model_name, 'StopTime', num2str(simulation_time));               %simulate the model and save the data to workspace

   %SPWM signal
   calcFFT(simOut.V_SPWM_Symm.time, simOut.V_SPWM_Symm.signals.values, 0, 0.05,0); %call the calcfft function to get amplitude and freuqency
   [S_thdSPWM_1000(i),S_wthdSPWM_1000(i),S_shdSPWM_1000(i),S_SPWM_voltage_1000(i)] = calculate_all_out(frequency, AMP, f0);     %calculate the THD & WTHD
  
   
   %SVM signal
   calcFFT(simOut.V_SVM_Symm.time, simOut.V_SVM_Symm.signals.values, 0, 0.05,0); 
   [S_thdSVM_1000(i),S_wthdSVM_1000(i),S_shdSVM_1000(i),S_SVM_voltage_1000(i)] = calculate_all_out(frequency, AMP, f0);
 
    
   %DPWM signal
   calcFFT(simOut.V_DWPM1_Symm.time , simOut.V_DWPM1_Symm.signals.values , 0, 0.05,0);
   [S_thdDPWM_1000(i),S_wthdDPWM_1000(i),S_shdDPWM_1000(i),S_DPWM_voltage_1000(i)] = calculate_all_out(frequency, AMP, f0);

 
end


    elseif j==2
     
for i = 1:length(Mi)                                                 %for loop to run every mudulation index
    M = Mi(i);
   simOut = sim(model_name, 'StopTime', num2str(simulation_time));  %simulate the model and save the data to workspace

   %SPWM signal
   calcFFT(simOut.V_SPWM_Symm.time, simOut.V_SPWM_Symm.signals.values, 0, 0.05,0);  %call the calcfft function to get amplitude and freuqency
   [S_thdSPWM_5000(i),S_wthdSPWM_5000(i),S_shdSPWM_5000(i),S_SPWM_voltage_5000(i)] = calculate_all_out(frequency, AMP, f0);     %calculate the THD & WTHD
  
   
   %SVM signal
   calcFFT(simOut.V_SVM_Symm.time, simOut.V_SVM_Symm.signals.values, 0, 0.05,0); 
   [S_thdSVM_5000(i),S_wthdSVM_5000(i),S_shdSVM_5000(i),S_SVM_voltage_5000(i)] = calculate_all_out(frequency, AMP, f0);
 
    
   %DPWM signal
   calcFFT(simOut.V_DWPM1_Symm.time , simOut.V_DWPM1_Symm.signals.values , 0, 0.05,0);
   [S_thdDPWM_5000(i),S_wthdDPWM_5000(i),S_shdDPWM_5000(i),S_DPWM_voltage_5000(i)] = calculate_all_out(frequency, AMP, f0);

 
end

end
        
        
    end 


%%
%%Asymmetric
    elseif s==2
    
    Mode=1;
    for j=1:2
    fs=Switching_frequency(j);
    
    if j==1
for i = 1:length(Mi)                                                             %for loop to run every mudulation index
    M = Mi(i);
   simOut = sim(model_name, 'StopTime', num2str(simulation_time));               %simulate the model and save the data to workspace

   %SPWM signal
   calcFFT(simOut.V_SPWM_Symm.time, simOut.V_SPWM_Symm.signals.values, 0, 0.05,0); %call the calcfft function to get amplitude and freuqency
   [A_thdSPWM_1000(i),A_wthdSPWM_1000(i),A_shdSPWM_1000(i),A_SPWM_voltage_1000(i)] = calculate_all_out(frequency, AMP, f0);     %calculate the THD & WTHD
  
   
   %SVM signal
   calcFFT(simOut.V_SVM_Symm.time, simOut.V_SVM_Symm.signals.values, 0, 0.05,0); 
   [A_thdSVM_1000(i),A_wthdSVM_1000(i),A_shdSVM_1000(i),A_SVM_voltage_1000(i)] = calculate_all_out(frequency, AMP, f0);
 
    
   %DPWM signal
   calcFFT(simOut.V_DWPM1_Symm.time , simOut.V_DWPM1_Symm.signals.values , 0, 0.05,0);
   [A_thdDPWM_1000(i),A_wthdDPWM_1000(i),A_shdDPWM_1000(i),A_DPWM_voltage_1000(i)] = calculate_all_out(frequency, AMP, f0);

 
end


    elseif j==2
     
for i = 1:length(Mi)                                                 %for loop to run every mudulation index
    M = Mi(i);
   simOut = sim(model_name, 'StopTime', num2str(simulation_time));  %simulate the model and save the data to workspace

   %SPWM signal
   calcFFT(simOut.V_SPWM_Symm.time, simOut.V_SPWM_Symm.signals.values, 0, 0.05,0);  %call the calcfft function to get amplitude and freuqency
   [A_thdSPWM_5000(i),A_wthdSPWM_5000(i),A_shdSPWM_5000(i),A_SPWM_voltage_5000(i)] = calculate_all_out(frequency, AMP, f0);     %calculate the THD & WTHD
  
   
   %SVM signal
   calcFFT(simOut.V_SVM_Symm.time, simOut.V_SVM_Symm.signals.values, 0, 0.05,0); 
   [A_thdSVM_5000(i),A_wthdSVM_5000(i),A_shdSVM_5000(i),A_SVM_voltage_5000(i)] = calculate_all_out(frequency, AMP, f0);
 
    
   %DPWM signal
   calcFFT(simOut.V_DWPM1_Symm.time , simOut.V_DWPM1_Symm.signals.values , 0, 0.05,0);
   [A_thdDPWM_5000(i),A_wthdDPWM_5000(i),A_shdDPWM_5000(i),A_DPWM_voltage_5000(i)] = calculate_all_out(frequency, AMP, f0);
end
end
    end 
    end
end

%% plotting section:
Mideal=Mi;
Mode=0;
%Plot FFT for DPWM1 and SVM at modulation=1:
M=1;
simOut = sim(model_name, 'StopTime', num2str(simulation_time));
figure('Name', 'FFT DPWM1');
calcFFT(simOut.V_DWPM1_Symm.time , simOut.V_DWPM1_Symm.signals.values , 0, 0.05,1);
figure('Name', 'FFT SVM');
calcFFT(simOut.V_SVM_Symm.time, simOut.V_SVM_Symm.signals.values, 0, 0.05,1); 

% Plot Symmetric output:

Draw(Mi, S_SPWM_voltage_5000, S_SVM_voltage_5000, S_DPWM_voltage_5000,'Fundamental Voltage vs Min 5000' ,Mode, 'Modulation index', 'Fundamental Voltage', {'SPWM','SVM','DPWM'}, 'Voltage5000 vs Min Sym');plot(Mi,Mideal,'--');
Draw(Mi, S_SPWM_voltage_1000, S_SVM_voltage_1000, S_DPWM_voltage_1000,'Fundamental Voltage vs Min 1000' ,Mode, 'Modulation index', 'Fundamental Voltage', {'SPWM','SVM','DPWM'}, 'Voltage1000 vs Min Sym');plot(Mi,Mideal,'--');
Draw2(Mi, S_SPWM_voltage_5000, S_SPWM_voltage_1000,'SPWM vs Min 5000vs1000' ,Mode, 'Modulation index', 'Fundamental Voltage', {'SPWM 5000','SPWM100'}, 'SPWM vs Min 5000vs1000 Sym');plot(Mi,Mideal,'--');
Draw2(Mi, S_SVM_voltage_5000, S_SVM_voltage_1000,'SVM vs Min 5000vs1000' ,Mode, 'Modulation index', 'Fundamental Voltage', {'SVM 5000','SVM 1000'}, 'SVM vs Min 5000vs1000 Sym');plot(Mi,Mideal,'--');
Draw2(Mi, S_DPWM_voltage_5000, S_DPWM_voltage_1000,'DPWM vs Min 5000vs1000' ,Mode, 'Modulation index', 'Fundamental Voltage', {'DPWM 5000','DPWM1000'}, 'DPWM vs Min 5000vs1000 Sym');plot(Mi,Mideal,'--');

Draw(Mi, S_thdSPWM_5000, S_thdSVM_5000, S_thdDPWM_5000,'THD vs Min 5000' ,Mode, 'Modulation index', 'THD 5000', {'SPWM','SVM','DPWM'}, 'THD 5000 vs Min Sym');
Draw(Mi, S_thdSPWM_1000, S_thdSVM_1000, S_thdDPWM_1000,'THD vs Min 1000' ,Mode, 'Modulation index', 'THD 1000', {'SPWM','SVM','DPWM'}, 'THD 1000 vs Min Sym');
Draw2(Mi, S_thdSPWM_5000, S_thdSPWM_1000,'THD vs Min 5000vs1000' ,Mode, 'Modulation index', 'THD', {'SPWM 5000','SPWM 1000'}, 'SPWM THD vs Min 5000vs1000 Sym');
Draw2(Mi, S_thdSVM_5000, S_thdSVM_1000,'THD vs Min 5000vs1000' ,Mode, 'Modulation index', 'THD', {'SVM 5000','SVM 1000'}, 'SVM_THD vs Min 5000vs1000 Sym');
Draw2(Mi, S_thdDPWM_5000, S_thdDPWM_1000,'THD vs Min 5000vs1000' ,Mode, 'Modulation index', 'THD', {'DPWM 5000','DPWM 1000'}, 'DPWM_THD vs Min 5000vs1000 Sym');

Draw(Mi, S_wthdSPWM_5000, S_wthdSVM_5000, S_wthdDPWM_5000,'WTHD vs Min 5000' ,Mode, 'Modulation index', 'WTHD 5000', {'SPWM','SVM','DPWM'}, 'WTHD 5000 vs Min Sym');
Draw(Mi, S_wthdSPWM_1000, S_wthdSVM_1000, S_wthdDPWM_1000,'WTHD vs Min 1000' ,Mode, 'Modulation index', 'WTHD 1000', {'SPWM','SVM','DPWM'}, 'WTHD 1000 vs Min Sym');
Draw2(Mi, S_wthdSPWM_5000, S_wthdSPWM_1000,'WTHD vs Min 5000vs1000' ,Mode, 'Modulation index', 'WTHD', {'SPWM 5000','SPWM 1000'}, 'SPWM_WTHD vs Min 5000vs1000 Sym');
Draw2(Mi, S_wthdSVM_5000, S_wthdSVM_1000,'WTHD vs Min 5000vs1000' ,Mode, 'Modulation index', 'WTHD', {'SVM 5000','SVM 1000'}, 'SVM_WTHD vs Min 5000vs1000 Sym');
Draw2(Mi, S_wthdDPWM_5000, S_wthdDPWM_1000,'WTHD vs Min 5000vs1000' ,Mode, 'Modulation index', 'WTHD', {'DPWM 5000','DPWM 1000'}, 'DPWM_WTHD vs Min 5000vs1000 Sym');

Draw(Mi, S_shdSPWM_5000, S_shdSVM_5000, S_shdDPWM_5000,'SHD vs Min 5000' ,Mode, 'Modulation index', 'SHD 5000', {'SPWM','SVM','DPWM'}, 'SHD 5000 vs Min Sym');
Draw(Mi, S_shdSPWM_1000, S_shdSVM_1000, S_shdDPWM_1000,'SHD vs Min 1000' ,Mode, 'Modulation index', 'SHD 1000', {'SPWM','SVM','DPWM'}, 'SHD 1000 vs Min Sym');
Draw2(Mi, S_shdSPWM_5000, S_shdSPWM_1000,'SHD vs Min 5000vs1000' ,Mode, 'Modulation index', 'SHD', {'SPWM 5000','SPWM 1000'}, 'SPWM_SHD vs Min 5000vs1000 Sym');
Draw2(Mi, S_shdSVM_5000, S_shdSVM_1000,'SHD vs Min 5000vs1000' ,Mode, 'Modulation index', 'SHD', {'SVM 5000','SVM 1000'}, 'SVM_SHD vs Min 5000vs1000 Sym');
Draw2(Mi, S_shdDPWM_5000, S_shdDPWM_1000,'SHD vs Min 5000vs1000' ,Mode, 'Modulation index', 'SHD', {'DPWM 5000','DPWM 1000'}, 'DPWM_SHD vs Min 5000vs1000 Sym');


% Plot Asymmetric output
Mode=1;
Draw(Mi, A_SPWM_voltage_5000, A_SVM_voltage_5000, A_DPWM_voltage_5000,'Fundamental Voltage vs Min 5000' ,Mode, 'Modulation index', 'Fundamental Voltage', {'SPWM','SVM','DPWM'}, 'Voltage5000 vs Min Asym');plot(Mi,Mideal,'--');
Draw(Mi, A_SPWM_voltage_1000, A_SVM_voltage_1000, A_DPWM_voltage_1000,'Fundamental Voltage vs Min 1000' ,Mode, 'Modulation index', 'Fundamental Voltage', {'SPWM','SVM','DPWM'}, 'Voltage1000 vs Min Asym');plot(Mi,Mideal,'--');
Draw2(Mi, A_SPWM_voltage_5000, A_SPWM_voltage_1000,'SPWM vs Min 5000vs1000' ,Mode, 'Modulation index', 'Fundamental Voltage', {'SPWM 5000','SPWM100'}, 'SPWM vs Min 5000vs1000 Asym');plot(Mi,Mideal,'--');
Draw2(Mi, A_SVM_voltage_5000, A_SVM_voltage_1000,'SVM vs Min 5000vs1000' ,Mode, 'Modulation index', 'Fundamental Voltage', {'SVM 5000','SVM 1000'}, 'SVM vs Min 5000vs1000 Asym');plot(Mi,Mideal,'--');
Draw2(Mi, A_DPWM_voltage_5000, A_DPWM_voltage_1000,'DPWM vs Min 5000vs1000' ,Mode, 'Modulation index', 'Fundamental Voltage', {'DPWM 5000','DPWM1000'}, 'DPWM vs Min 5000vs1000 Asym');plot(Mi,Mideal,'--');

Draw(Mi, A_thdSPWM_5000, A_thdSVM_5000, A_thdDPWM_5000,'THD vs Min 5000' ,Mode, 'Modulation index', 'THD 5000', {'SPWM','SVM','DPWM'}, 'THD 5000 vs Min Asym');
Draw(Mi, A_thdSPWM_1000, A_thdSVM_1000, A_thdDPWM_1000,'THD vs Min 1000' ,Mode, 'Modulation index', 'THD 1000', {'SPWM','SVM','DPWM'}, 'THD 1000 vs Min Asym');
Draw2(Mi, A_thdSPWM_5000, A_thdSPWM_1000,'THD vs Min 5000vs1000' ,Mode, 'Modulation index', 'THD', {'SPWM 5000','SPWM 1000'}, 'SPWM THD vs Min 5000vs1000 Asym');
Draw2(Mi, A_thdSVM_5000, A_thdSVM_1000,'THD vs Min 5000vs1000' ,Mode, 'Modulation index', 'THD', {'SVM 5000','SVM 1000'}, 'SVM_THD vs Min 5000vs1000 Asym');
Draw2(Mi, A_thdDPWM_5000, A_thdDPWM_1000,'THD vs Min 5000vs1000' ,Mode, 'Modulation index', 'THD', {'DPWM 5000','DPWM 1000'}, 'DPWM_THD vs Min 5000vs1000 Asym');

Draw(Mi, A_wthdSPWM_5000, A_wthdSVM_5000, A_wthdDPWM_5000,'WTHD vs Min 5000' ,Mode, 'Modulation index', 'WTHD 5000', {'SPWM','SVM','DPWM'}, 'WTHD 5000 vs Min Asym');
Draw(Mi, A_wthdSPWM_1000, A_wthdSVM_1000, A_wthdDPWM_1000,'WTHD vs Min 1000' ,Mode, 'Modulation index', 'WTHD 1000', {'SPWM','SVM','DPWM'}, 'WTHD 1000 vs Min Asym');
Draw2(Mi, A_wthdSPWM_5000, A_wthdSPWM_1000,'WTHD vs Min 5000vs1000' ,Mode, 'Modulation index', 'WTHD', {'SPWM 5000','SPWM 1000'}, 'SPWM_WTHD vs Min 5000vs1000 Asym');
Draw2(Mi, A_wthdSVM_5000, A_wthdSVM_1000,'WTHD vs Min 5000vs1000' ,Mode, 'Modulation index', 'WTHD', {'SVM 5000','SVM 1000'}, 'SVM_WTHD vs Min 5000vs1000 Asym');
Draw2(Mi, A_wthdDPWM_5000, A_wthdDPWM_1000,'WTHD vs Min 5000vs1000' ,Mode, 'Modulation index', 'WTHD', {'DPWM 5000','DPWM 1000'}, 'DPWM_WTHD vs Min 5000vs1000 Asym');

Draw(Mi, A_shdSPWM_5000, A_shdSVM_5000, A_shdDPWM_5000,'SHD vs Min 5000' ,Mode, 'Modulation index', 'SHD 5000', {'SPWM','SVM','DPWM'}, 'SHD 5000 vs Min Asym');
Draw(Mi, A_shdSPWM_1000, A_shdSVM_1000, A_shdDPWM_1000,'SHD vs Min 1000' ,Mode, 'Modulation index', 'SHD 1000', {'SPWM','SVM','DPWM'}, 'SHD 1000 vs Min Asym');
Draw2(Mi, A_shdSPWM_5000, A_shdSPWM_1000,'SHD vs Min 5000vs1000' ,Mode, 'Modulation index', 'SHD', {'SPWM 5000','SPWM 1000'}, 'SPWM_SHD vs Min 5000vs1000 Asym');
Draw2(Mi, A_shdSVM_5000, A_shdSVM_1000,'SHD vs Min 5000vs1000' ,Mode, 'Modulation index', 'SHD', {'SVM 5000','SVM 1000'}, 'SVM_SHD vs Min 5000vs1000 Asym');
Draw2(Mi, A_shdDPWM_5000, A_shdDPWM_1000,'SHD vs Min 5000vs1000' ,Mode, 'Modulation index', 'SHD', {'DPWM 5000','DPWM 1000'}, 'DPWM_SHD vs Min 5000vs1000 Asym');
        
%%
function[f, Amp] = calcFFT(tiVect, sig, tiStart, tiEnd,draw_status)
if (~ isempty(tiStart))&&(~ isempty(tiEnd))
    % Trim signal if start and end time is given
    idxVld       = find ((tiVect>=tiStart) & (tiVect<=tiEnd));
    
    tiVect = tiVect(idxVld)- tiVect((idxVld(1)));
    sig    = sig(idxVld);
end
dtMin = mean(diff(tiVect(2:end)));         % Calculate the minimum time interval between data points
tiVectNew = min(tiVect):dtMin:max(tiVect); % Create a new time vector with equidistant time points
sigNew    = interp1(tiVect, sig, tiVectNew);

T = mean(diff(tiVectNew));                 % Sampling period
Fs = 1/T;                                  % Sampling frequency
L = length(tiVectNew);                     % Calculate the length of the new signal
Y = fft(sigNew);                           
P2 = abs(Y/L);                             %Calculate the one-sided amplitude spectrum
P1 = P2(1:(L/2+1));
P1(2:end-1) = 2*P1(2:end-1);
Amp = P1;
assignin('base', 'AMP', Amp);              % Assign the amplitude values at workspace
f = Fs*(0:(L/2))/L;                        %Calculate the corresponding frequencies for the amplitude values
assignin('base', 'frequency', f);          % Assign the frequency values at workspace
     if draw_status==1
      sigName = 'Udc [V]';
    subplot(2,1,1);
    plot(tiVectNew,  sigNew);
    grid on; grid minor;
    ylabel(sigName);
    xlabel('Time[s]')
    
    ax(1) = subplot(2,1,2);
    plot(f,P1) 
    set(gca, 'XScale', 'log')
    grid on; grid minor;
    ylabel('Amplitude');
    xlabel('frequency [Hz]')
    title([sigName, ', DC value = ', num2str(P1(1))]);
    ylim([0 1.1*max(P1(2:end))]);
    linkaxes(ax,'x');
        
    end
end
function [thd, wthd,shd,fundamental_value] = calculate_all_out(frequencies, amplitudes, fundamental_frequency)

    [maxValue, ~] = max(amplitudes);             %fundamental amplitude equal to the maximum number of AMP
    fundamental_value=(maxValue);                %save fundamental amplitude of signal
    
    fundamentalIndex = find(round(frequencies) == fundamental_frequency);     % Find the index of the fundamental and harmonics frequency
    second_harmonic_Index = find(round(frequencies) == 2*fundamental_frequency); 
    l=second_harmonic_Index-fundamentalIndex;
    assignin('base', 'l', l); 
    harmonic_index=round(frequencies((fundamentalIndex+l):l:end)/fundamental_frequency);
    assignin('base', 'harmonic_indexx', harmonic_index); 
    % Extract fundamental and harmonic amplitudes
    fundamental_amplitude = amplitudes(fundamentalIndex);
    harmonic_amplitudes = amplitudes((fundamentalIndex+l):l:end);
    assignin('base', 'harmonic_amddd', harmonic_amplitudes); 
    % Calculate RMS amplitudes
    rms_fundamental = fundamental_amplitude / sqrt(2);
    rms_harmonics = harmonic_amplitudes / sqrt(2);
    % Calculate THD (Percentage)
    thd = sqrt(sum(rms_harmonics.^2)) / rms_fundamental * 100;  
    % Calculate WTHD (Percentage)
    wthd = sqrt(sum((rms_harmonics./harmonic_index).^2)) / rms_fundamental * 100;
    % Calulate SHD (Percentage)
    SHD_frequency(1) = find(round(frequencies) == 5*fundamental_frequency);
    SHD_frequency(2) = find(round(frequencies) == 7*fundamental_frequency);
    SHD_frequency(3) = find(round(frequencies) == 11*fundamental_frequency);
    SHD_frequency(4) = find(round(frequencies) == 13*fundamental_frequency);
    harmonic_STD(1)=amplitudes(SHD_frequency(1))/sqrt(2);
    harmonic_STD(2)=amplitudes(SHD_frequency(2))/sqrt(2);
    harmonic_STD(3)=amplitudes(SHD_frequency(3))/sqrt(2);
    harmonic_STD(4)=amplitudes(SHD_frequency(4))/sqrt(2);
    assignin('base', 'harmonic_STD', harmonic_STD); 
    shd = (sqrt(sum(harmonic_STD.^2)) / rms_fundamental) * 100;

end
function Draw2(X_signal, Y_signal_1, Y_signal_2,title_name ,title_mode, x_label, y_label, Legend, figure_name)
    figure('Name', figure_name);
    
    % Plot the signals
    hold on;
    plot(X_signal, Y_signal_1, 'r');
    plot(X_signal, Y_signal_2, 'g');
    grid on;
    % Set plot labels and title based on mode
    switch title_mode
        case 0
            title_text =[title_name, '(Symmetric)'];
        case 1
            title_text =[title_name, '(Asymmetric)'];

    end
    
    xlabel(x_label);
    ylabel(y_label);
    title(title_text);
    
    % Set legend
    if ~isempty(Legend)
        legend(Legend);
    end
end
function Draw(X_signal, Y_signal_1, Y_signal_2, Y_signal_3,title_name ,title_mode, x_label, y_label, Legend, figure_name)
    figure('Name', figure_name);
    
    % Plot the signals
    hold on;
    plot(X_signal, Y_signal_1, 'r');
    plot(X_signal, Y_signal_2, 'g');
    plot(X_signal, Y_signal_3, 'b');

    grid on;
    % Set plot labels and title based on mode
    switch title_mode
        case 0
            title_text =[title_name, '(Symmetric)'];
        case 1
            title_text =[title_name, '(Asymmetric)'];

    end
    
    xlabel(x_label);
    ylabel(y_label);
    title(title_text);
    
    % Set legend
    if ~isempty(Legend)
        legend(Legend);
    end
end

