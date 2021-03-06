%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% InfluenzaC - Code to simulating the model of Hancioglu et. al. (2007)
% form "A dynamical model of human immune response to influenza A virus 
% infection" developed upon previosuly developed code. Code also models
% the spread of the Influenza virus over a user defined grid.
%
% Copyright (C) 2013  Jackie Ang, Jonny Brook-Bartlett, Alexander Erlich,
% James Mbewu and Robert Ross.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Note: this code was developed atop of a previously code developed
% by the following authors:
% Mike Boemo, Lukas Hutter, Noemi Picco, Alex Saunders, Huw
% Colin-York, Elizabeth McMillan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ dYdt ] = ODE_function(t,Y, state_matrix, param)
%System of ODEs taken from the paper below

% Variables, see B. Hancioglu et al. J. Theor. Biol., 
% 246, 2007, p70 - 86, Table1
V = Y(1);
H = Y(2);
I = Y(3);
M = Y(4);
F = Y(5);
R = Y(6);
E = Y(7);
P = Y(8);
A = Y(9);
S = Y(10);
D = 1 - H - R - I;




        %SCHEME FOR NORMAL AND VACCINE%
        if state_matrix == 0 || state_matrix == 1

        dVdt = param(1) * I - param(2) * S * A * V - param(3) * H * V - param(4) * V - (param(5) * V) / (1 + (param(6) * V));
        dHdt = param(7) * D * (H + R) + param(8) * R - param(9) * V * H - param(10) * F * H;
        dIdt = param(9) * V * H - param(11) * E * I - param(12) * I;
        dMdt = (param(13) * D + param(14) * V) * (1 - M) - param(15) * M;
        dFdt = param(16) * M + param(17) * I - param(18) * H * F - param(19) * F;
        dRdt = param(10) * F * H - param(8) * R;
        dEdt = param(20) * M * E - param(21) * E * I + param(22) * (1 - E);
        dPdt = param(23) * M * P + param(24) * (1 - P);
        dAdt = param(25) * P - param(26) * S * A * V - param(27) * A;
        dSdt = param(28) * P * (1 - S);

        dYdt = [dVdt ; dHdt ; dIdt ; dMdt ; dFdt ; dRdt ; dEdt ; dPdt ; dAdt ; dSdt];

        end
        
        
        %SCHEME FOR NORMAL/TREATMENT AND VACCINE/TREATMENT%
        if state_matrix == 2 || state_matrix == 3
   
        dVdt = param(1) * I - param(2) * S * A * V - param(3) * H * V - param(4) * V - (param(5) * V) / (1 + (param(6) * V)) - param(29)*V;
        dHdt = param(7) * D * (H + R) + param(8) * R - param(9) * V * H - param(10) * F * H;
        dIdt = param(9) * V * H - param(11) * E * I - param(12) * I;
        dMdt = (param(13) * D + param(14) * V) * (1 - M) - param(15) * M;
        dFdt = param(16) * M + param(17) * I - param(18) * H * F - param(19) * F;
        dRdt = param(10) * F * H - param(8) * R;
        dEdt = param(20) * M * E - param(21) * E * I + param(22) * (1 - E);
        dPdt = param(23) * M * P + param(24) * (1 - P);
        dAdt = param(25) * P - param(26) * S * A * V - param(27) * A;
        dSdt = param(28) * P * (1 - S);

        dYdt = [dVdt ; dHdt ; dIdt ; dMdt ; dFdt ; dRdt ; dEdt ; dPdt ; dAdt ; dSdt];
    
        end


end

