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


function [Vc_output] = SDEpredsolver(initial_value,state_matrix,times,param,random)

tBegin=times(1);
tEnd=times(2);

dt=param(30);
t=tBegin:dt:tEnd;
N=length(t);

%% Preallocate the predictor and corrector solution matrices 

Vc = zeros(10,N);
Vp = zeros(10,N);

D = zeros(1,N);

Vc(:,1) = initial_value;
D(1) = 1 - Vc(2,1) - Vc(3,1) - Vc(6,1);

%% Solve the ODES

%SCHEME FOR NORMAL AND VACCINE%
    if state_matrix == 0 || state_matrix == 1
        for i = 2:N
    
            Vp(1,i) = Vc(1,i-1)+dt*(param(1) * Vc(3,i-1) - param(2) * Vc(10,i-1) * Vc(9,i-1) * Vc(1,i-1) - param(3) * Vc(2,i-1) * Vc(1,i-1) - param(4) * Vc(1,i-1) - (param(5) * Vc(1,i-1)) / (1 + (param(6) * Vc(1,i-1))));
            Vp(2,i) = Vc(2,i-1)+dt*(param(7) * D(i-1) * (Vc(2,i-1) + Vc(6,i-1)) + param(8) * Vc(6,i-1) - param(9) * Vc(1,i-1) * Vc(2,i-1) - param(10) * Vc(5,i-1) * Vc(2,i-1));
            Vp(3,i) = Vc(3,i-1)+dt*(param(9) * Vc(1,i-1) * Vc(2,i-1) - param(11) * Vc(7,i-1) * Vc(3,i-1) - param(12)* Vc(3,i-1));
            Vp(4,i) = Vc(4,i-1)+dt*((param(13) * D(i-1) + param(14) * Vc(1,i-1)) * (1 - Vc(4,i-1)) - param(15) * Vc(4,i-1));
            Vp(5,i) = Vc(5,i-1)+dt*(param(16) * Vc(4,i-1) + param(17) * Vc(3,i-1) - param(18) * Vc(2,i-1) * Vc(5,i-1) - param(19) * Vc(5,i-1));
            Vp(6,i) = Vc(6,i-1)+dt*(param(10) * Vc(5,i-1) * Vc(2,i-1) - param(8) * Vc(6,i-1));
            Vp(7,i) = Vc(7,i-1)+dt*(param(20) * Vc(4,i-1) * Vc(7,i-1) - param(21) * Vc(7,i-1) * Vc(3,i-1) + param(22) * (1 - Vc(7,i-1)));
            Vp(8,i) = Vc(8,i-1)+dt*(param(23) * Vc(4,i-1) * Vc(8,i-1) + param(24) * (1 - Vc(8,i-1)));
            Vp(9,i) = Vc(9,i-1)+dt*(param(25) * Vc(8,i-1) - param(26) * Vc(10,i-1) * Vc(9,i-1) * Vc(1,i-1) - param(27) * Vc(9,i-1));
            Vp(10,i) = Vc(10,i-1)+dt*(param(28) * Vc(8,i-1) * (1 - Vc(10,i-1)));
            
            
            Vc(1,i) = Vc(1,i-1)+dt*(param(1) * Vp(3,i) - param(2) * Vp(10,i) * Vp(9,i) * Vp(1,i) - param(3) * Vp(2,i) * Vp(1,i) - param(4) * Vp(1,i) - (param(5) * Vp(1,i)) / (1 + (param(6) * Vp(1,i))));
            Vc(2,i) = Vc(2,i-1)+dt*(param(7) * D(i-1) * (Vp(2,i) + Vp(6,i)) + param(8) * Vp(6,i) - param(9) * Vp(1,i) * Vp(2,i) - param(10) * Vp(5,i) * Vp(2,i));
            Vc(3,i) = Vc(3,i-1)+dt*(param(9) * Vp(1,i) * Vp(2,i) - param(11) * Vp(7,i) * Vp(3,i) - param(12) * Vp(3,i));
            Vc(4,i) = Vc(4,i-1)+dt*((param(13) * D(i-1) + param(14) * Vp(1,i)) * (1 - Vp(4,i)) - param(15) * Vp(4,i));
            Vc(5,i) = Vc(5,i-1)+dt*(param(16) * Vp(4,i) + param(17) * Vp(3,i) - param(18) * Vp(2,i) * Vp(5,i) - param(19) * Vp(5,i));
            Vc(6,i) = Vc(6,i-1)+dt*(param(10)* Vp(5,i) * Vp(2,i) - param(8) * Vp(6,i));
            Vc(7,i) = Vc(7,i-1)+dt*(param(20) * Vp(4,i) * Vp(7,i) - param(21) * Vp(7,i) * Vp(3,i) + param(22) * (1 - Vp(7,i)));
            Vc(8,i) = Vc(8,i-1)+dt*(param(23) * Vp(4,i) * Vp(8,i) + param(24) * (1 - Vp(8,i)));
            Vc(9,i) = Vc(9,i-1)+dt*(param(25) * Vp(8,i) - param(26) * Vp(10,i) * Vp(9,i) * Vp(1,i) - param(27) * Vp(9,i));
            Vc(10,i) = Vc(10,i-1)+dt*(param(28) * Vp(8,i) * (1 - Vp(10,i))) + random*sqrt(dt)*randn(1)*Vc(10,i-1);
            
            D(i) = 1 - Vc(2,i) - Vc(6,i) - Vc(3,i);
            
        end
        
    end

    %SCHEME FOR NORMAL/TREATMENT AND VACCINE/TREATMENT%
    if state_matrix == 2  || state_matrix == 3
    
        for i = 2:N
   
            Vp(1,i) = Vc(1,i-1)+dt*(param(1) * Vc(3,i-1) - param(2) * Vc(10,i-1) * Vc(9,i-1) * Vc(1,i-1) - param(3) * Vc(2,i-1) * Vc(1,i-1) - param(4) * Vc(1,i-1) - (param(5) * Vc(1,i-1)) / (1 + (param(6) * Vc(1,i-1))) - param(29)*Vc(1,i-1));
            Vp(2,i) = Vc(2,i-1)+dt*(param(7) * D(i-1) * (Vc(2,i-1) + Vc(6,i-1)) + param(8) * Vc(6,i-1) - param(9) * Vc(1,i-1) * Vc(2,i-1) - param(10) * Vc(5,i-1) * Vc(2,i-1));
            Vp(3,i) = Vc(3,i-1)+dt*(param(9) * Vc(1,i-1) * Vc(2,i-1) - param(11) * Vc(7,i-1) * Vc(3,i-1) - param(12)* Vc(3,i-1));
            Vp(4,i) = Vc(4,i-1)+dt*((param(13) * D(i-1) + param(14) * Vc(1,i-1)) * (1 - Vc(4,i-1)) - param(15) * Vc(4,i-1));
            Vp(5,i) = Vc(5,i-1)+dt*(param(16) * Vc(4,i-1) + param(17) * Vc(3,i-1) - param(18) * Vc(2,i-1) * Vc(5,i-1) - param(19) * Vc(5,i-1));
            Vp(6,i) = Vc(6,i-1)+dt*(param(10) * Vc(5,i-1) * Vc(2,i-1) - param(8) * Vc(6,i-1));
            Vp(7,i) = Vc(7,i-1)+dt*(param(20) * Vc(4,i-1) * Vc(7,i-1) - param(21) * Vc(7,i-1) * Vc(3,i-1) + param(22) * (1 - Vc(7,i-1)));
            Vp(8,i) = Vc(8,i-1)+dt*(param(23) * Vc(4,i-1) * Vc(8,i-1) + param(24) * (1 - Vc(8,i-1)));
            Vp(9,i) = Vc(9,i-1)+dt*(param(25) * Vc(8,i-1) - param(26) * Vc(10,i-1) * Vc(9,i-1) * Vc(1,i-1) - param(27) * Vc(9,i-1));
            Vp(10,i) = Vc(10,i-1)+dt*(param(28) * Vc(8,i-1) * (1 - Vc(10,i-1)));
            
            Vc(1,i) = Vc(1,i-1)+dt*(param(1) * Vp(3,i) - param(2) * Vp(10,i) * Vp(9,i) * Vp(1,i) - param(3) * Vp(2,i) * Vp(1,i) - param(4) * Vp(1,i) - (param(5) * Vp(1,i)) / (1 + (param(6) * Vp(1,i)))- param(29)*Vc(1,i-1));
            Vc(2,i) = Vc(2,i-1)+dt*(param(7) * D(i-1) * (Vp(2,i) + Vp(6,i)) + param(8) * Vp(6,i) - param(9) * Vp(1,i) * Vp(2,i) - param(10) * Vp(5,i) * Vp(2,i));
            Vc(3,i) = Vc(3,i-1)+dt*(param(9) * Vp(1,i) * Vp(2,i) - param(11) * Vp(7,i) * Vp(3,i) - param(12) * Vp(3,i));
            Vc(4,i) = Vc(4,i-1)+dt*((param(13) * D(i-1) + param(14) * Vp(1,i)) * (1 - Vp(4,i)) - param(15) * Vp(4,i));
            Vc(5,i) = Vc(5,i-1)+dt*(param(16) * Vp(4,i) + param(17) * Vp(3,i) - param(18) * Vp(2,i) * Vp(5,i) - param(19) * Vp(5,i));
            Vc(6,i) = Vc(6,i-1)+dt*(param(10)* Vp(5,i) * Vp(2,i) - param(8) * Vp(6,i));
            Vc(7,i) = Vc(7,i-1)+dt*(param(20) * Vp(4,i) * Vp(7,i) - param(21) * Vp(7,i) * Vp(3,i) + param(22) * (1 - Vp(7,i)));
            Vc(8,i) = Vc(8,i-1)+dt*(param(23) * Vp(4,i) * Vp(8,i) + param(24) * (1 - Vp(8,i)));
            Vc(9,i) = Vc(9,i-1)+dt*(param(25) * Vp(8,i) - param(26) * Vp(10,i) * Vp(9,i) * Vp(1,i) - param(27) * Vp(9,i));
            Vc(10,i) = Vc(10,i-1)+dt*(param(28) * Vp(8,i) * (1 - Vp(10,i))) + random*sqrt(dt)*randn(1)*Vc(10,i-1);
    
            D(i) = 1 - Vc(2,i) - Vc(6,i) - Vc(3,i);
            
        end
        
    end

    Vc_output = Vc(:,N);
    
end

