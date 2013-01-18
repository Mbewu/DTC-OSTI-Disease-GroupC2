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



function varargout = gui2(varargin)
% GUI2 MATLAB code for gui2.fig
%      GUI2, by itself, creates a new GUI2 or raises the existing
%      singleton*.
%
%      H = GUI2 returns the handle to a new GUI2 or the handle to
%      the existing singleton*.
%
%      GUI2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI2.M with the given input arguments.
%
%      GUI2('Property','Value',...) creates a new GUI2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui2

% Last Modified by GUIDE v2.5 16-Jan-2013 21:21:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui2_OpeningFcn, ...
                   'gui_OutputFcn',  @gui2_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before gui2 is made visible.
function gui2_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui2 (see VARARGIN)
handles.stop = false;   %added
% Choose default command line output for gui2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.stop = false;
% Get default command line output from handles structure
varargout{1} = handles.output;

guidata(hObject, handles);



function pop_size_input_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function pop_size_input_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double



% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double



% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double




% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1




% --- Executes on button press in pop_size_set.
function pop_size_set_Callback(hObject, eventdata, handles)
% hObject    handle to pop_size_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pop_size; % grid size
global frames_per_sec; % grid size

pop_size = str2double(get(handles.pop_size_input,'String')); %get input from 'enter grid size'
frames_per_sec = str2double(get(handles.edit4,'String')); %get input from 'enter grid size'

fprintf('pop_size = %d',pop_size)
axis([handles.matrix_plot],[0 pop_size 0 pop_size]); %set axis to  pop_size

global counter; %define counter to count number of points selected on axis 
global coordinates_input; %define vector to store x,y coordiantes of points
coordinates_input = 0;
counter = 1;



% --- Executes on button press in run_simulation. slow 
function run_simulation_Callback(hObject, eventdata, handles)
% hObject    handle to run_simulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global pop_size;
global quarantine_end_threshold; % grid size
global frames_per_sec; % grid size
global quarantine_start_threshold; % grid size
global percentage_quarantine; % grid size
global percentage_treated; % grid size
global viral_start; % grid size
global parallel;
global susceptibility; % grid size
global percentage_vaccine; % grid size
global randomInfection; % grid size
global stochastic; % grid size
global coordinates_input;
global runtime; %number of days to run simulation 
global matrix_values; %takes output from simulation (3D matrix of structs) 
global healthy_output;  %takes output from simulation (1D matrix) 
global h;
global parameters;
global state_matrix;
global wait_time;

quarantine_end_threshold = str2double(get(handles.edit3,'String')); %get input from 'enter grid size'
frames_per_sec = str2double(get(handles.edit4,'String')); %get input from 'enter grid size'
wait_time = 1/frames_per_sec;
quarantine_start_threshold = str2double(get(handles.edit6,'String')); %get input from 'enter grid size'
percentage_quarantine = str2double(get(handles.edit7,'String')); %get input from 'enter grid size'
percentage_treated = str2double(get(handles.edit8,'String')); %get input from 'enter grid size'
viral_start = str2double(get(handles.edit10,'String')); %get input from 'enter grid size'

if get(handles.checkbox1,'Value') == get(handles.checkbox1,'Max')
    parallel = 1; %get input from 'enter grid size' 
else
    parallel = 0;
end
susceptibility = str2double(get(handles.edit11,'String')); %get input from 'enter grid size'
percentage_vaccine = str2double(get(handles.edit13,'String')); %get input from 'enter grid size' 
if get(handles.checkbox2,'Value') == get(handles.checkbox2,'Max')
    randomInfection= 1; %get input from 'enter grid size'
else
    randomInfection = 0;
end
if get(handles.checkbox3,'Value') == get(handles.checkbox3,'Max')
    stochastic= 1; %get input from 'enter grid size' 
else
    stochastic = 0;
end
runtime = str2double(get(handles.runtime_input,'String')); %get input from 'enter run time'

parameters = zeros(13,1);
parameters(1) = pop_size;
parameters(2) = runtime;
parameters(3) = viral_start;
parameters(4) = percentage_quarantine;
parameters(5) = quarantine_start_threshold;
parameters(6) = quarantine_end_threshold;
parameters(7) = percentage_treated;
parameters(8) = percentage_vaccine;
maxViralLoad = 80;
parameters(9) = 0.0023*8/maxViralLoad/susceptibility;
parameters(10) = randomInfection;
parameters(11) = stochastic;
parameters(12) = parallel;

h = waitbar(0,'Initializing waitbar...'); %progress bar for simulation
waitbar(0,h,'0%...')
[matrix_values, healthy_output, state_matrix] = Coupled_Multiscale_parallelised_GUI(parameters,coordinates_input); %run full simulation
waitbar(1,h,'Complete...')
pause(1)
close(h)

%global counter;
%counter = 1; %reset counter to 1


% --- Executes on mouse press over axes background.
function matrix_plot_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to matrix_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global counter;
global coordinates_input; 
global x1;
global y1;

handles.xy1 = round(get(handles.matrix_plot,'Currentpoint')); % gets x,y values of point on mousebuttondown 
x1 = handles.xy1(1,1);
y1 = handles.xy1(1,2);

hold on
plot(x1,y1,'ro'); % plots red * at selected point
hold off

coordinates_input (1,counter) = y1; %fills coordiantes vector with x y values
coordinates_input (2,counter) = x1;

counter = counter + 1;




function runtime_input_Callback(hObject, eventdata, handles)
% hObject    handle to runtime_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of runtime_input as text
%        str2double(get(hObject,'String')) returns contents of runtime_input as a double


% --- Executes during object creation, after setting all properties.
function runtime_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to runtime_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% Hints: contents = cellstr(get(hObject,'String')) returns varibale_input contents as cell array
%        contents{get(hObject,'Value')} returns selected item from varibale_input

% --- Executes on button press in play_sim.
function play_sim_Callback(hObject, eventdata, handles)
% hObject    handle to play_sim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global matrix_values;
global runtime;
global variable_index;
global wait_time;
%global healthy_output;
fps = round(1/wait_time);

if variable_index==1 %if user has selected V
maximum = max(max(max(matrix_values.V)));
for index_time=2:runtime
slice = matrix_values.V(:,:,index_time);
imagesc(squeeze(slice))
caxis([min(min(slice)) maximum])
cb = colorbar;
TL = get(cb, 'YTick');
%set(cb, 'YTickLabel', sprintf('10^%g|', TL));
set(gca,'YDir','normal')
ylabel(cb,'Viral load','fontsize',12)
%Mov(index_time-1) = getframe(handles.figure1);
pause(wait_time);
end
elseif variable_index==2 %if user has selected H (does not work with Run Simulation 2)
for index_time=2:runtime
slice = matrix_values.H(:,:,index_time);
imagesc(squeeze(slice))
caxis([0 1])
cb = colorbar;
TL = get(cb, 'YTick');
%set(cb, 'YTickLabel', sprintf('10^%g|', TL));
set(gca,'YDir','normal')
ylabel(cb,'Healthy cells','fontsize',12)
pause(wait_time);
%Mov(index_time-1) = getframe(handles.figure1);
end
elseif variable_index==3 %if user has selected H (does not work with Run Simulation 2)
for index_time=2:runtime
slice = matrix_values.I(:,:,index_time);
imagesc(squeeze(slice))
caxis([0 1])
cb = colorbar;
TL = get(cb, 'YTick');
%set(cb, 'YTickLabel', sprintf('10^%g|', TL));
set(gca,'YDir','normal')
ylabel(cb,'Infected cells','fontsize',12)
pause(wait_time);
%Mov(index_time-1) = getframe(handles.figure1);
end
elseif variable_index==4 %if user has selected H (does not work with Run Simulation 2)
for index_time=2:runtime
slice = matrix_values.M(:,:,index_time);
imagesc(squeeze(slice))
caxis([0 1])
cb = colorbar;
TL = get(cb, 'YTick');
%set(cb, 'YTickLabel', sprintf('10^%g|', TL));
set(gca,'YDir','normal')
ylabel(cb,'Activated antigen presenting cells','fontsize',12)
pause(wait_time);
%Mov(index_time-1) = getframe(handles.figure1);
end
elseif variable_index==5 %if user has selected H (does not work with Run Simulation 2)
maximum = max(max(max(matrix_values.F)));
for index_time=2:runtime
slice = matrix_values.F(:,:,index_time);
imagesc(log(squeeze(slice)))
caxis([min(min(slice)) maximum])
cb = colorbar;
TL = get(cb, 'YTick');
set(cb, 'YTickLabel', sprintf('10^%g|', TL));
set(gca,'YDir','normal')
ylabel(cb,'Inteferons','fontsize',12)
pause(wait_time);
%Mov(index_time-1) = getframe(handles.figure1);
end
elseif variable_index==6 %if user has selected H (does not work with Run Simulation 2)
for index_time=2:runtime
slice = matrix_values.R(:,:,index_time);
imagesc(squeeze(slice))
caxis([0 1])
cb = colorbar;
TL = get(cb, 'YTick');
%set(cb, 'YTickLabel', sprintf('10^%g|', TL));
set(gca,'YDir','normal')
ylabel(cb,'Resistant cells','fontsize',12)
pause(wait_time);
%Mov(index_time-1) = getframe(handles.figure1);
end
elseif variable_index==7 %if user has selected H (does not work with Run Simulation 2)
maximum = max(max(max(matrix_values.E)));
for index_time=2:runtime
slice = matrix_values.E(:,:,index_time);
imagesc(log(squeeze(slice)))
caxis([min(min(slice)) maximum])
cb = colorbar;
TL = get(cb, 'YTick');
set(cb, 'YTickLabel', sprintf('10^%g|', TL));
set(gca,'YDir','normal')
ylabel(cb,'Effector cells','fontsize',12)
pause(wait_time);
%Mov(index_time-1) = getframe(handles.figure1);
end
elseif variable_index==8 %if user has selected H (does not work with Run Simulation 2)
maximum = max(max(max(matrix_values.P)));
for index_time=2:runtime
slice = matrix_values.P(:,:,index_time);
imagesc(log(squeeze(slice)))
caxis([min(min(slice)) maximum])
cb = colorbar;
TL = get(cb, 'YTick');
set(cb, 'YTickLabel', sprintf('10^%g|', TL));
set(gca,'YDir','normal')
ylabel(cb,'Plasma cells','fontsize',12)
pause(wait_time);
%Mov(index_time-1) = getframe(handles.figure1);
end
elseif variable_index==9 %if user has selected H (does not work with Run Simulation 2)
maximum = max(max(max(matrix_values.A)));
for index_time=2:runtime
slice = matrix_values.A(:,:,index_time);
imagesc(log(squeeze(slice)))
caxis([min(min(slice)) maximum])
cb = colorbar;
TL = get(cb, 'YTick');
set(cb, 'YTickLabel', sprintf('10^%g|', TL));
set(gca,'YDir','normal')
ylabel(cb,'Antibodies','fontsize',12)
pause(wait_time);
%Mov(index_time-1) = getframe(handles.figure1);
end
elseif variable_index==10 %if user has selected H (does not work with Run Simulation 2)
for index_time=2:runtime
slice = matrix_values.S(:,:,index_time);
imagesc(squeeze(slice))
caxis([0 1])
cb = colorbar;
TL = get(cb, 'YTick');
%set(cb, 'YTickLabel', sprintf('10^%g|', TL));
set(gca,'YDir','normal')
ylabel(cb,'Antigenic distance','fontsize',12)
pause(wait_time);
%Mov(index_time-1) = getframe(handles.figure1);
end
end

% 
% if variable_index==1
%     movie2avi(Mov,'fluMovieViralLoad','compression','None','fps',fps);
% elseif variable_index==2
%     movie2avi(Mov,'fluMovieHealthy','compression','None','fps',fps);
% elseif variable_index==3
%     movie2avi(Mov,'fluMovieInfected','compression','None','fps',fps);
% elseif variable_index==4
%     movie2avi(Mov,'fluMovieActivatedAntigen','compression','None','fps',fps);
% elseif variable_index==5
%     movie2avi(Mov,'fluMovieInteferonsDistance','compression','None','fps',fps);
% elseif variable_index==6
%     movie2avi(Mov,'fluMovieResistantDistance','compression','None','fps',fps);
% elseif variable_index==7
%     movie2avi(Mov,'fluMovieEffectorDistance','compression','None','fps',fps);
% elseif variable_index==8
%     movie2avi(Mov,'fluMoviePlasma','compression','None','fps',fps);
% elseif variable_index==9
%     movie2avi(Mov,'fluMovieAntibodies','compression','None','fps',fps);
% elseif variable_index==10
%     movie2avi(Mov,'fluMovieAntigenDistance','compression','None','fps',fps);
% end




% --- Executes on button press in play_3D_sim.
function play_3D_sim_Callback(hObject, eventdata, handles)
% hObject    handle to play_3D_sim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%makes movie of surfaces
global matrix_values;
global runtime;
global healthy_output;
global variable_index;
global wait_time;

figure; 
if variable_index==1 %if user has selected V
for index_time=2:runtime
    surf(matrix_values.V(:,:,index_time),healthy_output(:,:,index_time));
    zlabel('V')
    title('Viral load per epithelial cell for each individual')
    pause(wait_time);
end
elseif variable_index==2%if user has selected H (does not work with Run Simulation 2)
for index_time=2:runtime
    surf(matrix_values.H(:,:,index_time),healthy_output(:,:,index_time));
    zlabel('H')
    title('Proportion healthy cells')
    pause(wait_time);
end
elseif variable_index==3%if user has selected H (does not work with Run Simulation 2)
for index_time=2:runtime
    surf(matrix_values.I(:,:,index_time),healthy_output(:,:,index_time));
    zlabel('V')
    title('Proportion infected cells')
    pause(wait_time);
end
elseif variable_index==4%if user has selected H (does not work with Run Simulation 2)
for index_time=2:runtime
    surf(matrix_values.M(:,:,index_time),healthy_output(:,:,index_time));
    zlabel('M')
    title('Antigen presenting cells')
    pause(wait_time);
end
elseif variable_index==5%if user has selected H (does not work with Run Simulation 2)
for index_time=2:runtime
    surf(matrix_values.F(:,:,index_time),healthy_output(:,:,index_time));
    zlabel('F')
    title('Inteferons')
    pause(wait_time);
end
elseif variable_index==6%if user has selected H (does not work with Run Simulation 2)
for index_time=2:runtime
    surf(matrix_values.R(:,:,index_time),healthy_output(:,:,index_time));
    zlabel('R')
    title('Proportion resistant cells')
    pause(wait_time);
end
elseif variable_index==7%if user has selected H (does not work with Run Simulation 2)
for index_time=2:runtime
    surf(matrix_values.E(:,:,index_time),healthy_output(:,:,index_time));
    zlabel('E')
    title('Effector cells')
    pause(wait_time);
end
elseif variable_index==8%if user has selected H (does not work with Run Simulation 2)
for index_time=2:runtime
    surf(matrix_values.P(:,:,index_time),healthy_output(:,:,index_time));
    zlabel('P')
    title('Plasma cells')
    pause(wait_time);
end
elseif variable_index==9%if user has selected H (does not work with Run Simulation 2)
for index_time=2:runtime
    surf(matrix_values.A(:,:,index_time),healthy_output(:,:,index_time));
    zlabel('A')
    title('Antibodies')
    pause(wait_time);
end
elseif variable_index==10%if user has selected H (does not work with Run Simulation 2)
for index_time=2:runtime
    surf(matrix_values.S(:,:,index_time),healthy_output(:,:,index_time));
    zlabel('S')
    title('Antigenic distance')
    pause(wait_time);
end
end





% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global pop_size; 
global coordinates_input;
global counter;
global runtime; %number of days to run simulation 
global matrix_values; %takes output from simulation (3D matrix of structs) 
global healthy_output;  %takes output from simulation (1D matrix) 
global variable_index;

coordinates_input
counter
runtime

for i=1:counter-1
    x = 1:runtime;
    index_x = coordinates_input(1,i);
    index_y = coordinates_input(2,i);
    figure; 
    squeeze(matrix_values.V(index_x,index_y,:))
    
    h1 = subplot(4,3,1), plot(x,squeeze(matrix_values.V(index_x,index_y,:)));
    title('V: Viral Load');
    xlabel('Days');
    ylabel('V');
    h2 = subplot(4,3,2), plot(x,squeeze(matrix_values.H(index_x,index_y,:)));
    title('H: Healthy Cells');
    xlabel('Days');
    ylabel('H');
    h3 = subplot(4,3,3), plot(x,squeeze(matrix_values.I(index_x,index_y,:)));
    title('I: Infected Cells');
    xlabel('Days');
    ylabel('I');
    h4 = subplot(4,3,4), plot(x,squeeze(matrix_values.M(index_x,index_y,:)));
    title('M: Activated Antigen-presenting cells');
    xlabel('Days');
    ylabel('M');
    h5 = subplot(4,3,5), semilogy(x,squeeze(matrix_values.F(index_x,index_y,:)));
    title('F: Interferons');
    xlabel('Days');
    ylabel('F');
    h6 = subplot(4,3,6), plot(x,squeeze(matrix_values.R(index_x,index_y,:)));
    title('R: Resistant Cells');
    xlabel('Days');
    ylabel('R');
    h7 = subplot(4,3,7), semilogy(x,squeeze(matrix_values.E(index_x,index_y,:)));
    title('E: Effector Cells');
    xlabel('Days');
    ylabel('E');
    h8 = subplot(4,3,8), semilogy(x,squeeze(matrix_values.P(index_x,index_y,:)));
    title('P: Plasma Cells');
    xlabel('Days');
    ylabel('P');
    h9 = subplot(4,3,9), semilogy(x,squeeze(matrix_values.A(index_x,index_y,:)));
    title('A: Antibodies');
    xlabel('Days');
    ylabel('A');
    h10 = subplot(4,3,11), plot(x,squeeze(matrix_values.S(index_x,index_y,:)));
    title('S: Antigenic Distance');
    xlabel('Days');
    ylabel('S');
    h11 = subplot(4,3,10), plot(x,1-squeeze(matrix_values.H(index_x,index_y,:))-squeeze(matrix_values.R(index_x,index_y,:))-squeeze(matrix_values.I(index_x,index_y,:)));
    title('D: Dead Cells');
    xlabel('Days');
    ylabel('D');

    axis([h1 h5 h7 h8 h9],'tight');
    axis(h1,[0 runtime 0.00001 150]);
    axis(h5,[0 runtime 0.1 15000]);
    axis(h7,[0 runtime 0.01 150]);
    axis(h8,[0 runtime 0.1 15000]);
    axis(h9,[0 runtime 0.00001 1000]);
    axis([h2 h3 h4 h6 h10 h11],[0 runtime 0 1]);
end


counter = 1; %reset counter to 1






% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global counter; %define counter to count number of points selected on axis 
global coordinates_input; %define vector to store x,y coordiantes of points
global pop_size;

axis([handles.matrix_plot],[0 pop_size 0 pop_size]); %set axis to  pop_size

coordinates_input = 0;
counter = 1;


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5
global variable_index;

variable_index = 1;
variable_index = get(hObject,'Value')

% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
