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



function SolveODESystem()

global param;
treatmentParam = 5;
SetParams(treatmentParam);

%0 - normal
%1 - normal and vaccine
%2 - normal and treatment
%3 - vaccination and treatment
state = 0;

%y0 = [V,     H,  I,  M,  F,  R,  E,  P,  A,  S]
y0 = [ 0.1 , 1 , 0 , 0 , 0 , 0 , 1 , 1 , 1 , 0.1 ];
%y0 = [ 0.1 , 1 , 0 , 0 , 0 , 0 , 1 , 1 , 1 , 0.35 ];
figure
endTime = 15;

tic;
[x,Y]=ode15s(@(x,Y)ODE_function(x,Y,state,param),[0,endTime],y0);
toc;

% Run the stochastic ODE solver             
%[ Y ] = SDEpredsolver(y0, 0,[0,endTime], param);
%Set the final_value struct to the final results from the ODE solver
%final_value = Y(:,end)';    
D= 1-Y(:,2)-Y(:,6)-Y(:,3);
Y2= [Y'; D']';

% 
% 
% 
h1 = subplot(4,3,1), plot(x,Y2(:,1));
title('V: Viral Load');
xlabel('Days');
ylabel('V');
h2 = subplot(4,3,2), plot(x,Y2(:,2));
title('H: Healthy Cells');
xlabel('Days');
ylabel('H');
h3 = subplot(4,3,3), plot(x,Y2(:,3));
title('I: Infected Cells');
xlabel('Days');
ylabel('I');
h4 = subplot(4,3,4), plot(x,Y2(:,4));
title('M: Activated Antigen-presenting cells');
xlabel('Days');
ylabel('M');
h5 = subplot(4,3,5), semilogy(x,Y2(:,5));
title('F: Interferons');
xlabel('Days');
ylabel('F');
h6 = subplot(4,3,6), plot(x,Y2(:,6));
title('R: Resistant Cells');
xlabel('Days');
ylabel('R');
h7 = subplot(4,3,7), semilogy(x,Y2(:,7));
title('E: Effector Cells');
xlabel('Days');
ylabel('E');
h8 = subplot(4,3,8), semilogy(x,Y2(:,8));
title('P: Plasma Cells');
xlabel('Days');
ylabel('P');
h9 = subplot(4,3,9), semilogy(x,Y2(:,9));
title('A: Antibodies');
xlabel('Days');
ylabel('A');
h10 = subplot(4,3,11), plot(x,Y2(:,10));
title('S: Antigenic Distance');
xlabel('Days');
ylabel('S');
h11 = subplot(4,3,10), plot(x,Y2(:,11));
title('D: Dead Cells');
xlabel('Days');
ylabel('D');

axis([h1 h5 h7 h8 h9],'tight');
axis([h1],[0 endTime 0.00001 150]);
axis([h5],[0 endTime 0.1 15000]);
axis([h7],[0 endTime 0.01 150]);
axis([h8],[0 endTime 0.1 15000]);
axis([h9],[0 endTime 0.00001 1000]);
axis([h2 h3 h4 h6 h10 h11],[0 endTime 0 1]);

%legend('V','H','I','M','F','R','E','P','A','S','D')
end

