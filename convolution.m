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



%do convolution with randomised average kernel, scaled by scalar
function [ outputMatrix ] = convolution( matrix, susceptibility)
%convolution - does convolution of 

   N = length(matrix);
   outputMatrix = zeros(N);
   matrix_calculate = zeros(N+2);
   matrix_calculate(2:N+1,2:N+1) = matrix;

   for i = 2:N+1
       
       for j = 2:N+1
           
           infectionKernel = zeros(3);
           y = rand(1,8);
           S = sum(y);
           y = y/S;
           infectionKernel(1,1) = y(1);
           infectionKernel(1,2) = y(2);
           infectionKernel(1,3) = y(3);
           infectionKernel(2,1) = y(4);
           infectionKernel(2,3) = y(5);
           infectionKernel(3,1) = y(6);
           infectionKernel(3,2) = y(7);
           infectionKernel(3,3) = y(8);
           
           infectionKernel = infectionKernel*susceptibility;
       
           outputMatrix(i-1,j-1) = sum(sum(matrix_calculate(i-1:i+1,j-1:j+1).*infectionKernel));
   
       end
   end


end


