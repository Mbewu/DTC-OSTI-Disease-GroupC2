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


function [ final_value ] = SolveODESystem( initial_value, state_matrix, param, doStochastic )
% Input an initial vector for passing to the ODE solver
%   Returns a vector final_value of same length as initial_value
final_value = 0;
%parameter to control the amount of stochasticity
randomParam = 0.1;

if(doStochastic)
    % Run the stochastic ODE solver             
    [ Y ] = SDEpredsolver(initial_value, state_matrix, [0,1], param,randomParam);
    %Set the final_value struct to the final results from the ODE solver
    final_value = Y';    
else
    % Run the ODE solver   
    [ x, Y ] = ode15s(@(x,Y)ODE_function(x,Y,state_matrix,param),[0,1],initial_value);
    %Set the final_value struct to the final results from the ODE solver
    final_value = Y(end,:);
end

end

