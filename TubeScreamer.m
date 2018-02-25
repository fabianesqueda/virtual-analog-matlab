% This script emulates a guitar distortion stage reminiscent to that of the
% Ibanez Tube Screamer.
% 
% Input values
% 
% Vin - Input Signal 
%  Fs - Sampling Rate (approx 8x44.1 kHz)
%   k - "Distortion" level [0,1]
%
% Code written by F. Esqueda 1/10/2017
% based on the schematic found at: https://www.electrosmash.com/tube-screamer-analysis
%
%% Load input signal

[Vin, Fs] = audioread('TestGuitarPhrase.wav');
Vin = Vin(:,1);             % Convert to mono
OS = 8;                     % Oversampling factor
Vin = resample(Vin,OS,1);   % Upsample
Vin = Vin/max(abs(Vin));    % Normalize
Vout = zeros(size(Vin));    % Output vector

%% User Parameters

G = 0.1;      % Input gain (in V)  [0,1]
k = 1;          % "Distortion" level [0,1]
if (k<0) || (k>1)
    error('Value of k should be between 0–1');
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%$NO MORE USER-ADJUSTABLE PARAMETERS PAST THIS POINT%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Circuit Simulation Parameters

Is = 2.52e-9;   % Saturation Current    (N914)
N = 1.752;      % Diode ideality factor (N914)
VT = 25.864e-3; % Thermal Voltage

% Component Values
R1 = 4.7e3;
R2 = 51e3 + k*500e3;
C1 = 47e-9;
C2 = 51e-12;

Ts = 1/(OS*Fs); % Sampling period

%% Processing Loop

% State variables
Vin_n1 = 0; V_HP_n1 = 0; V_n1 = 0;

% Processing Loop
for n=1:length(Vin)
   
    Vin(n) = G*Vin(n);  % Input Gain
    
    % High Pass Input Stage
    V_HP = ((1 - Ts/(2*R1*C1))*V_HP_n1 + Vin(n) - Vin_n1)/(1 + Ts/(2*R1*C1));
    I = V_HP/R1;
    
    % Newton Solver
    V = V_n1;
    for m=1:100
       
        F = V - V_n1 - (Ts/C2)*I + (Ts/(2*R2*C2))*(V+V_n1) + (Ts*Is/C2)*sinh(V/(N*VT)) + (Ts*Is/C2)*sinh(V_n1/(N*VT));
        dF = 1 + (Ts/(2*R2*C2)) + ((Ts*Is)/(N*VT*C2))*cosh(V/(N*VT));
        err = F/dF;
        
        if abs(err)<10e-10
            break;
        else
            V = V - err;
        end
                
    end
    
    % Read Output
    Vout(n) = Vin(n) - V;
    
    % Update states
    Vin_n1 = Vin(n);    V_HP_n1 = V_HP;    V_n1 = V;
    
end

Vout = resample(Vout,1,8);
soundsc(Vout,Fs);
