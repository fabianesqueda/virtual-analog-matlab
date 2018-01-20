function [ Vout ] = LockhartWavefolder( Vin, RL )
% This function implements a single stage of the Lockhart wavefolder. 
% A detailed description of the model and its recommended implementation
% can be found in the accompanying artcile:
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
%   RL (Ohm) - Load resistance (defaults to 7k5)
% Output:
%   Vout (V) - Folded signal
% 
% Code written by F. Esqueda 12/10/2017
%% Validate input parameters

if (size(Vin,2)>1)
    error('Please input a mono signal.')
end

if (nargin<2)
    RL = 7.5e3;
end

%% Component Values
R = 15e3;       % Supply resistor
Is = 1e-14;     % Reverse bias current (transistor p–n junction)
N = 1;          % Ideality factor      (transistor p–n junction)
VT = 25.864e-3; % Thermal Voltage

Vout = zeros(size(Vin));

for n=1:length(Vin)

    lambda = sign(Vin(n));
    Lu = LambertWFritsch((RL*Is/VT)*exp(lambda*((R+2*RL)/(VT*R))*Vin(n)));
    Vout(n) = (2*RL/R)*Vin(n) - lambda*VT*Lu;
   
end


end

