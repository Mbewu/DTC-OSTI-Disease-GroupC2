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


function SetParams(treatment)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global param;

%NORMAL PARAMETERS%

param = zeros(1,29);


param(1) = 510;
param(2) = 619.2;
param(3) = 1.02;
param(4) = 1.7;
param(5) = 100;
param(6) = 23000;
param(7) = 4;
param(8) = 1;
param(9) = 0.34;
param(10) = 0.01;
param(11) = 0.066;
param(12) = 1.5;
param(13) = 1;
param(14) = 0.0037;
param(15) = 1;
param(16) = 250000;
param(17) = 2000;
param(18) = 17;
param(19) = 8;
param(20) = 8.3;
param(21) = 2.72;
param(22) = 0.4;
param(23) = 11.5;
param(24) = 0.4;
param(25) = 0.043;
param(26) = 146.2;
param(27) = 0.043;
param(28) = (3 * 10 ^ (-5));
param(29) = treatment;
param(30) = 0.001;  %deltaT



end



