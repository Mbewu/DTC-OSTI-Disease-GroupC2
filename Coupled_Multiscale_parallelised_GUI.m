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


%ONLYGUI
function [variable healthy state_matrix] = Coupled_Multiscale_parallelised_GUI(parameters,coordinates)

%##############Section.1##################
%###INITIALISATION OF THE DATASTRUCTURE###
%#########################################

global h;

global param;

%treatmentAmount
treatment = 5;
SetParams(treatment);

gridDim = parameters(1);
endTime = parameters(2);
viralStart = parameters(3);
percentageQuarantined = parameters(4);
quarantineStartThreshold = parameters(5);
quarantineEndThreshold = parameters(6);
percentageTreated = parameters(7);
percentageVaccinated = parameters(8);
infectivity = parameters(9);        % proportion viral load that gets from surrounding people
randomInfectiousness = parameters(10);
stochastic = parameters(11);
parallel = parameters(12);
%GUI tunable parameters
quarantine = 0;
if(percentageQuarantined>0)
    quarantine = 1;
end

%Initial Conditions for the immune response model ( B. Hancioglu J. Theo.Biol. 2006)
%y0= [  V,    H,  I,  M,  F,  R,  E,  P,  A,   S  ]
y0 = [ 0.1 , 1 , 0 , 0 , 0 , 0 , 1 , 1 , 1 , 0.1 ];
y0vaccine = [ 0.1 , 1 , 0 , 0 , 0 , 0 , 1 , 1 , 1 , 0.35 ];





%state matrix has 0 - normal, 1 - vaccine and normal, 2 - treatment and
%normal - 3 - vaccine and treatment
%set treated values of statematrix
state_matrix = zeros(gridDim,gridDim);
Nel = numel(state_matrix);
Rindices = randperm(Nel);
if(percentageVaccinated > 0)
    vaccineIndices = floor(Nel*percentageVaccinated/100);
    state_matrix(Rindices(1:vaccineIndices)) = 1;    
end

%Establishing a routine for accessing elements of the struct 'variables'

varNames = {'V',...         %Viral load per epithelial cell
            'H',...         %Proportion of healthy cells
            'I',...         %Proportion of infected cells
            'M',...         %Activated antigen-presenting Cells (APC) per homeostatic level
            'F',...         %Interferons per homeostatic level of macrophages
            'R',...         %Proportion of resistant cells         
            'E',...         %Effector cells per homeostatic level
            'P',...         %Plasma cells per homeostatic level
            'A',...         %Antibodies per homeostatic level
            'S'};           %Antigenic distance
                            %...Proportion of dead cells (D) can be calculated
                            %as D = 1 - H - R - I

%Initialises a structure containing the  initial values of each variable except V for each position in the grid.
%+ Memory preallocation

for i = 2:numel(varNames)
    variable.(varNames{i}) = repmat(y0(i), [gridDim gridDim endTime]);
end

%set initial condition of vaccinated
for index_x=1:gridDim
    for index_y= 1:gridDim
        if state_matrix(index_x,index_y)
            for time=1:endTime
                variable.H(index_x,index_y,time) = y0vaccine(2);
                variable.I(index_x,index_y,time) = y0vaccine(3);
                variable.M(index_x,index_y,time) = y0vaccine(4);
                variable.F(index_x,index_y,time) = y0vaccine(5);
                variable.R(index_x,index_y,time) = y0vaccine(6);
                variable.E(index_x,index_y,time) = y0vaccine(7);
                variable.P(index_x,index_y,time) = y0vaccine(8);
                variable.A(index_x,index_y,time) = y0vaccine(9);
                variable.S(index_x,index_y,time) = y0vaccine(10);
            end
        end
    end
end

%Initialises variable V as zero for all elements in the grid
variable.V = zeros([gridDim gridDim endTime]);

% creates and initialises 3D matrix for a person's state: HEALTHY =0,
% INFECTED = 1... not used...
healthy = zeros([gridDim gridDim endTime]);


%Assigning the value for V taken from the vector holding the initial conditions
%to individuals specified in input by user.
%ONLY GUI
for i=1:length(coordinates(1,:))
    variable.V(coordinates(1,i),coordinates(2,i)) = y0(1);
end


%set treated values of statematrix
Rindices = randperm(Nel);
if(percentageTreated > 0)
    treatedIndices = floor(Nel*percentageTreated/100)
    state_matrix(Rindices(1:treatedIndices)) = state_matrix(Rindices(1:treatedIndices)) + 2;
end


%Create quarantine false
%quarantineTrue is true is candidate for quarantine

quarantineTrue = zeros(gridDim,gridDim);
Nel = numel(quarantineTrue);
Rindices = randperm(Nel);
quarantinedIndices = 0;


if(quarantine)
    quarantinedIndices = floor(Nel*percentageQuarantined/100);
    quarantineTrue(Rindices(1:quarantinedIndices)) = 1;
end    


%quarantineFalse is false if actively quarantined
quarantineFalse = ones(gridDim,gridDim);

transmitV = zeros(gridDim,gridDim);
stochasticTimer = 6*ones(gridDim,gridDim);
doStochastic = 0;

%Assigning the value for V taken from the vector holding the initial conditions
%to individuals specified in input by user.

%##############Section.2##################
%#########RUNNING THE SIMULATION##########
%#########################################

for index_time = 2: endTime
    t = tic;
    numInfectedCells = 0;
    %viralStart just zero
    [i,j]=find(variable.V(:,:,index_time-1) ~= viralStart);
    
    
    for n = 1:numel(i)
        temp1(n) = variable.V(i(n),j(n),index_time - 1);
        temp2(n) = variable.H(i(n),j(n),index_time - 1);
        temp3(n) = variable.I(i(n),j(n),index_time - 1);
        temp4(n) = variable.M(i(n),j(n),index_time - 1);
        temp5(n) = variable.F(i(n),j(n),index_time - 1);
        temp6(n) = variable.R(i(n),j(n),index_time - 1);
        temp7(n) = variable.E(i(n),j(n),index_time - 1);
        temp8(n) = variable.P(i(n),j(n),index_time - 1);
        temp9(n) = variable.A(i(n),j(n),index_time - 1);
        temp10(n) = variable.S(i(n),j(n),index_time - 1);
        temp11(n) = state_matrix(i(n),j(n));
        stochasticTimer(i(n),j(n)) = stochasticTimer(i(n),j(n)) - 1;
    end
        %parallel or not
        if(parallel)
            parfor n = 1:numel(i)                              
                
                if(stochastic == 1  && stochasticTimer(i(n),j(n)) > 0)
                    doStochastic = 1;
                else
                    doStochastic = 0;
                end
                 
                     
                 %run Ode and update variables for this person     
                 [tempVec] = SolveODESystem ([temp1(n), ...
                                              temp2(n), ...
                                              temp3(n), ...
                                              temp4(n), ...
                                              temp5(n), ...
                                              temp6(n), ...
                                              temp7(n), ...
                                              temp8(n), ...
                                              temp9(n), ...
                                              temp10(n)], ...
                                              temp11(n), ...
                                              param, ...
                                              doStochastic);
                                          
                     temp1(n) = tempVec(1);                                                      
                     temp2(n) = tempVec(2);
                     temp3(n) = tempVec(3);
                     temp4(n) = tempVec(4);
                     temp5(n) = tempVec(5);
                     temp6(n) = tempVec(6);
                     temp7(n) = tempVec(7);
                     temp8(n) = tempVec(8);
                     temp9(n) = tempVec(9);
                     temp10(n) = tempVec(10);

                     
            end%if...
        else
            for n = 1:numel(i)  
                stochasticTimer(i(n),j(n)) = stochasticTimer(i(n),j(n)) - 1;
                if(stochastic == 1 && stochasticTimer(i(n),j(n)) > 0)
                    doStochastic = 1;
                else
                    doStochastic = 0;
                end
                 %run Ode and update variables for this person     
                 [tempVec] = SolveODESystem ([temp1(n), ...
                                              temp2(n), ...
                                              temp3(n), ...
                                              temp4(n), ...
                                              temp5(n), ...
                                              temp6(n), ...
                                              temp7(n), ...
                                              temp8(n), ...
                                              temp9(n), ...
                                              temp10(n)], ...
                                              temp11(n), ...
                                              param, ...
                                              doStochastic);
                                          
                     temp1(n) = tempVec(1);                                                      
                     temp2(n) = tempVec(2);
                     temp3(n) = tempVec(3);
                     temp4(n) = tempVec(4);
                     temp5(n) = tempVec(5);
                     temp6(n) = tempVec(6);
                     temp7(n) = tempVec(7);
                     temp8(n) = tempVec(8);
                     temp9(n) = tempVec(9);
                     temp10(n) = tempVec(10);
            end
        end

                  
            
            
        for n = 1:numel(i)    
            variable.V(i(n),j(n),index_time) = temp1(n);
            variable.H(i(n),j(n),index_time) = temp2(n);
            variable.I(i(n),j(n),index_time) = temp3(n);
            variable.M(i(n),j(n),index_time) = temp4(n);
            variable.F(i(n),j(n),index_time) = temp5(n);
            variable.R(i(n),j(n),index_time) = temp6(n);
            variable.E(i(n),j(n),index_time) = temp7(n);
            variable.P(i(n),j(n),index_time) = temp8(n);
            variable.A(i(n),j(n),index_time) = temp9(n);
            variable.S(i(n),j(n),index_time) = temp10(n);
            %counter for num infected cells, although could 
            %just count from healthy matrix in clever way
            numInfectedCells = numInfectedCells + 1;
            
        end
                
        clear i j;
                   
        %ONLY GUI
        waitbar(index_time/endTime,h,sprintf('timestep %d / %d...',index_time,endTime)) %sends progress to the gui
        %Infection of other people based on convolution of kernel and array of
        %variable V
        %generate random numbers and scale
        
        
        infectionKernel = 1/8*ones(3);        
        infectionKernel(2,2) = 0;
        
        infectionKernel = infectionKernel*infectivity; %needed to infect from one reaching max
                
        if(quarantine)
            if randomInfectiousness
                susceptibility = infectivity * abs((1 + 3*(randn(1))));
                transmitV = convolution(variable.V(:,:,index_time).*quarantineFalse,susceptibility);
            else                
                transmitV = convn(variable.V(:,:,index_time).*quarantineFalse,infectionKernel,'same');
            end
            transmitV = transmitV.*quarantineFalse;
        else
            if randomInfectiousness
                susceptibility = abs(infectivity *(1 + 3*(randn(1))));
                transmitV = convolution(variable.V(:,:,index_time).*quarantineFalse,susceptibility); 
            else
                transmitV = convn(variable.V(:,:,index_time),infectionKernel,'same');
            end
        end
        
        variable.V(:,:,index_time) = variable.V(:,:,index_time) + transmitV;            
                
        % update quarantine
        for index_x=1:gridDim
                for index_y= 1:gridDim
                    
                    %if above quarantine threshold and not quarantined then
                    %quarantine
                    %if below quarantin threshold and quarantined then
                    %remove from quarantine
                    if (quarantineTrue(index_x,index_y) && quarantineFalse(index_x,index_y) == 1 ...
                            && variable.V(index_x,index_y,index_time)>quarantineStartThreshold)
                        quarantineFalse(index_x,index_y)=0;
                    elseif (quarantineTrue(index_x,index_y) && quarantineFalse(index_x,index_y) == 0 ...
                            && variable.V(index_x,index_y,index_time)<quarantineEndThreshold)
                        quarantineFalse(index_x,index_y)=1;
                    end
                    
                end
        end
        
end %t

end %func