function [ Vout ] = SergeWavefolder( Vin )
% This function implements a single stage of the Serge wavefolder (aka the 
% middle section of the Serge VCM). A detailed description of the model and 
% its recommended implementation can be found in the accompanying artcile:
%
% Virtual Analog Models of the Lockhart and Serge Wavefolders by
% F. Esqueda, H. Pöntynen, J. D. Parker and S. Bilbao
% Available at: http://www.mdpi.com/2076-3417/7/12/1328
%
% This function requires the function "Lambert_W_Fritsch.m" to run.
% Alternatively, it can be replaced with MATLAB's (much slower) native
% "lambertw" function.
%
% Input values:
%   Vin (V) - Input signal
% Output:
%   Vout (V) - Folded signal
% 
% Code written by F. Esqueda 12/10/2017

%% Component Values
R = 33e3;
Is = 2.52e-9;   % Saturation current (N914)
N = 1.752;      % Diode ideality factor (N914)
VT = 0.026; % Thermal Voltage

Vout = zeros(size(Vin));

for n=1:length(Vin)

    lambda = sign(Vin(n));
        
    Vout(n) = Vin(n) + 2*(lambda*R*Is - lambda*N*VT*LambertWFritsch(((R*Is)/(N*VT))*exp( (lambda*(Vin(n) + lambda*R*Is))/(N*VT))));

end


end

