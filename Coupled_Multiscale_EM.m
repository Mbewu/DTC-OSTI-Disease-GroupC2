%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%    Authors: Mike Boemo, Lukas Hutter, Noemi Picco, Alex Saunders, Huw
%    Colin-York, Elizabeth McMillan

% to run it type:
% [variable healthy] = Coupled_Multiscale_EM(10, [7, 4; 4 6],510,10)
function [variable healthy] = Coupled_Multiscale_EM(gridDim_input,coordinates,gammaV_input,endTime_input)
% Implementation of standard forward Euler. As vectorised as it can be.
% Very poor for our coupled stiff system.
% 
%##############Section.1##################
%###INITIALISATION OF THE DATASTRUCTURE###
%#########################################

% (gridDim_input,coordinates,gammaV_input) Must be set in the following GUI
% script
global h;

%Model parameters:

endTime= endTime_input;           %Number of timesteps (in days) to be performed
gridDim = gridDim_input;           %Size of the grid/lattice for the large-scale infection modelling

global param;
SetParams(gammaV_input);

%Definition of the location of patient zero within the grid
%         xPatientZero = 2;
%         yPatientZero = 5;


%Initial Conditions for the immune response model ( B. Hancioglu J. Theo.Biol. 2006)

%y0= [  V,    H,  I,  M,  F,  R,  E,  P,  A,   S  ]

% ##########################################
% change the last value: 0.1 to YOUR VALUE
% ##########################################
y0 = [ 0.01 , 1 , 0 , 0 , 0 , 0 , 1 , 1 , 1 , 0.1 ];



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

%Initialises variable V as zero for all elements in the grid
variable.V = zeros([gridDim gridDim endTime]);

% creates and initialises 3D matrix for a person's state: HEALTHY =0, INFECTED = 1
healthy = zeros([gridDim gridDim endTime]);

%Assigning the value for V taken from the vector holding the initial conditions
%to individuals specified in input by user.
for i=1:length(coordinates(1,:))
    variable.V(coordinates(1,i),coordinates(2,i)) = y0(1);
    healthy(coordinates(1,i),coordinates(2,i),1)=1;
end

Vtemp = zeros(gridDim);
Htemp = zeros(gridDim);
Itemp = zeros(gridDim);
Mtemp = zeros(gridDim);
Ftemp = zeros(gridDim);
Rtemp = zeros(gridDim);
Etemp = zeros(gridDim);
Ptemp = zeros(gridDim);
Atemp = zeros(gridDim);
Stemp = zeros(gridDim);

dt = 1e-2;


for index_time = 2: endTime
    t = tic;
    

           
                Vtemp = variable.V(:,:,index_time-1);
                Htemp = variable.H(:,:,index_time-1);
                Itemp = variable.I(:,:,index_time-1);
                Mtemp = variable.M(:,:,index_time-1);
                Ftemp = variable.F(:,:,index_time-1);
                Rtemp = variable.R(:,:,index_time-1);
                Etemp = variable.E(:,:,index_time-1);
                Ptemp = variable.P(:,:,index_time-1);
                Atemp = variable.A(:,:,index_time-1);
                Stemp = variable.S(:,:,index_time-1);
                
                figure(5);
                pcolor(Vtemp);
                colorbar;
                %pause(0.1);
                Vtemp
                
                % sickos is a vector (not matrix) of unknown length
                % V, H, I etc will have the same dimension as sickos
                sickos = find(Vtemp>0); 
                
                V = Vtemp(sickos);
                H = Htemp(sickos);
                I = Itemp(sickos);
                M = Mtemp(sickos);
                F = Ftemp(sickos);
                R = Rtemp(sickos);
                E = Etemp(sickos);
                P = Ptemp(sickos);
                A = Atemp(sickos);
                S = Stemp(sickos);
                D = 1 - H - R - I;
                
                for i=0:dt:1
                   
                    
                    Va = V + dt*(param.gammaV * I - param.gammaVA * S .*A .* V - param.gammaVH * H .* V - param.alphaV * V - (param.aV1 * V) ./ (1 + (param.aV2 * V))); % ... + sqrt(4) * ones(size(sickos)) * 
                    Ha = H + dt*(param.bHD * D .* (H + R) + param.aR * R - param.gammaHV * V .* H - param.bHF * F .* H);
                    Ia = I + dt*(param.gammaHV * V .* H - param.bIE * E .* I - param.aI * I);
                    Ma = M + dt*(param.bMD * D + param.bMV * V .* (1 - M) - param.aM * M);
                    Fa = F + dt*(param.bF * M + param.cF * I - param.bFH * H .* F - param.aF * F);
                    Ra = R + dt*(param.bHF * F .* H - param.aR * R);
                    Ea = E + dt*(param.bEM * M .* E - param.bEI * E .* I + param.aE * (1 - E));
                    Pa = P + dt*(param.bPM * M .* P + param.aP * (1 - P));
                    Aa = A + dt*(param.bA * P - param.gammaAV * S .* A .* V - param.aA * A);
                    Sa = S + dt*(param.r * P .* (1 - S));
                    Da = 1 - H - R - I;
                    
                V = Va;
                H = Ha;
                I = Ia;
                M = Ma;
                F = Fa;
                R = Ra;
                E = Ea;
                P = Pa;
                A = Aa;
                S = Sa;
                D = 1 - Ha - Ra - Ia;
                    
                end
                
                Vtemp = zeros(gridDim);
                Htemp = zeros(gridDim);
                Itemp = zeros(gridDim);
                Mtemp = zeros(gridDim);
                Ftemp = zeros(gridDim);
                Rtemp = zeros(gridDim);
                Etemp = zeros(gridDim);
                Ptemp = zeros(gridDim);
                Atemp = zeros(gridDim);
                Stemp = zeros(gridDim);
                
                Vtemp(sickos) = V;
                Htemp(sickos) = H;
                Itemp(sickos) = I;
                Mtemp(sickos) = M;
                Ftemp(sickos) = F;
                Rtemp(sickos) = R;
                Etemp(sickos) = E;
                Ptemp(sickos) = P;
                Atemp(sickos) = A;
                Stemp(sickos) = S;
                
                
                variable.V(:,:,index_time) = Vtemp;
                variable.H(:,:,index_time) = Htemp;
                variable.I(:,:,index_time) = Itemp;
                variable.M(:,:,index_time) = Mtemp;
                variable.F(:,:,index_time) = Ftemp;
                variable.R(:,:,index_time) = Rtemp;
                variable.E(:,:,index_time) = Etemp;
                variable.P(:,:,index_time) = Ptemp;
                variable.A(:,:,index_time) = Atemp;
                variable.S(:,:,index_time) = Stemp;
                
        
        % waitbar(index_time/endTime,h,sprintf('timestep %d / %d...',index_time,endTime)) %sends progress to the gui
        %Infection of other people based on convolution of kernel and array of
        %variable V

        infectionKernel = ones(3)/8;
        infectionKernel(2,2) = 0;

        transmitV = convn(variable.V(:,:,index_time),infectionKernel,'same');
        variable.V(:,:,index_time) = variable.V(:,:,index_time) + transmitV;


        timePerIteration = zeros(endTime);        
        timePerIteration(index_time-1) = toc(t);

        % update healthy
        for index_x=1:gridDim
                for index_y= 1:gridDim
                    if (variable.V(index_x,index_y,index_time)>0)
                        healthy(index_x,index_y,index_time)=1;
                    end
                end
        end
        
        current_slice = variable.V(:,:,index_time);
        clims = [ 0 0.000001 ];
        imagesc(squeeze(current_slice),clims);
        set(gca,'YDir','normal')
        colorbar
        pause(0.1)

end %t

end% function 
