
function varargout = BlackHoleSimulator(varargin)
% BLACKHOLESIMULATOR MATLAB code for BlackHoleSimulator.fig
%      BLACKHOLESIMULATOR, by itself, creates a new BLACKHOLESIMULATOR or raises the existing
%      singleton*.
%
%      H = BLACKHOLESIMULATOR returns the handle to a new BLACKHOLESIMULATOR or the handle to
%      the existing singleton*.
%
%      BLACKHOLESIMULATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BLACKHOLESIMULATOR.M with the given input arguments.
%
%      BLACKHOLESIMULATOR('Property','Value',...) creates a new BLACKHOLESIMULATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BlackHoleSimulator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BlackHoleSimulator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BlackHoleSimulator

% Last Modified by GUIDE v2.5 06-Apr-2017 09:28:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BlackHoleSimulator_OpeningFcn, ...
                   'gui_OutputFcn',  @BlackHoleSimulator_OutputFcn, ...
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


% --- Executes just before BlackHoleSimulator is made visible.
function BlackHoleSimulator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BlackHoleSimulator (see VARARGIN)

% Choose default command line output for BlackHoleSimulator
handles.output = hObject;
handles.isMassive = 0;
handles.isAnimated = 0;
handles.particleAngularMomentum = 0;
handles.BHAngularMomentum = 0;
handles.CIRCULAR_ORBIT_POTENTIAL = 9;
handles.CIRCULAR_ORBIT_RADIAL_DISTANCE = 7;
handles.STABLE_ORBIT_POTENTIAL = 4;
handles.STABLE_ORBIT_RADIAL_DISTANCE = 5;
handles.blackHoleMass = 0;
handles.blackHoleType = 1; % Enumeration
                           %    Classical Newtonian = 0
                           %    Schwarzchild = 1
                           %    Kerr = 2
handles.particleIndex = 6; % keeps track of which particle is selected
handles.plotRadialDistance = 0;
handles.plotEnergy = 0;
%store values in particleMatrix for each particle:
%in order: mass, black hole type, black hole mass, angular momentum, plot
%flag, potential, radial distance
handles.particleMatrix = zeros(6, 7);
handles.particleMatrix(1,2) = 1;
handles.particleMatrix(1,5) = 1;
handles.data = cell(6,7);
handles.data(:) = {''};
handles.data(:,1) = {'no'};
handles.data(2:6,5) = {'no'};
handles.data(1,5) = {'yes'};
handles.data(:, 3:4) = {'0'};
handles.data(:, 6:7) = {'0'};
handles.data(1, 2) = {'Schwarzschild'};
handles.data(2:6, 2) = {'Newtonian'};
set(handles.listOfParticles,'Data',handles.data);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BlackHoleSimulator wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BlackHoleSimulator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function txtParticleAngularMomentum_Callback(hObject, eventdata, handles)
% hObject    handle to txtParticleAngularMomentum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtParticleAngularMomentum as text
%        str2double(get(hObject,'String')) returns contents of txtParticleAngularMomentum as a double
handles.particleAngularMomentum = str2double(get(hObject,'String'));
handles.particleMatrix(handles.particleIndex, 4) = handles.particleAngularMomentum ;
handles.data(handles.particleIndex, 4) = {handles.particleAngularMomentum};
set(handles.listOfParticles,'Data',handles.data)
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txtParticleAngularMomentum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtParticleAngularMomentum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderRadialDistance_Callback(hObject, eventdata, handles)
% hObject    handle to sliderRadialDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.plotRadialDistance = get(hObject, 'Value');
plot(handles.plotRadialDistance, handles.plotEnergy, 'go');
handles.particleMatrix(handles.particleIndex, 7) = handles.plotRadialDistance ;
axis([0 1 0 1]);

handles.data(handles.particleIndex, 7) = {handles.plotRadialDistance};
set(handles.listOfParticles,'Data',handles.data)
guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sliderRadialDistance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderRadialDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on slider movement.
function sliderPotentialEnergy_Callback(hObject, eventdata, handles)
% hObject    handle to sliderPotentialEnergy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.plotEnergy = get(hObject, 'Value');
plot(handles.plotRadialDistance , handles.plotEnergy, 'ro');
axis([0 1 0 1]); % [xmin xmax ymin ymax]

handles.particleMatrix(handles.particleIndex, 6) = handles.plotEnergy ;
handles.data(handles.particleIndex, 6) = {handles.plotEnergy};
set(handles.listOfParticles,'Data',handles.data)
guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sliderPotentialEnergy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderPotentialEnergy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in isAnimated.
function isAnimated_Callback(hObject, eventdata, handles)
% hObject    handle to isAnimated (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.isAnimated = get(hObject, 'Value');

guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of isAnimated


% --- Executes on button press in btnRunSimulation.
function btnRunSimulation_Callback(hObject, eventdata, handles)
% hObject    handle to btnRunSimulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.isAnimated == 1 
    plot()
else
    plot()%% can't use plot function
end    


% --- Executes on button press in btnAddParticle.
function btnAddParticle_Callback(hObject, eventdata, handles)
% hObject    handle to btnAddParticle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.particleMatrix(handles.particleIndex, 5) = 1;
handles.data(handles.particleIndex, 5) = {'yes'};
set(handles.listOfParticles,'Data',handles.data)
guidata(hObject, handles);

% --- Executes on button press in btnRemoveParticle.
function btnRemoveParticle_Callback(hObject, eventdata, handles)
% hObject    handle to btnRemoveParticle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.particleMatrix(handles.particleIndex, 5) = 0;
handles.data(handles.particleIndex, 5) = {'no'};
set(handles.listOfParticles,'Data',handles.data)
guidata(hObject, handles);

% --- Executes on button press in btnCreateStableOrbit.
function btnCreateStableOrbit_Callback(hObject, eventdata, handles)
% hObject    handle to btnCreateStableOrbit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.blackHoleMass = handles.particleMatrix(handles.particleIndex, 3);
handles.particleAngularMomentum = handles.particleMatrix(handles.particleIndex, 4);
handles.STABLE_ORBIT_RADIAL_DISTANCE = handles.particleAngularMomentum + handles.blackHoleMass; %function of angular momentum, mass etc.
handles.STABLE_ORBIT_POTENTIAL= handles.particleAngularMomentum + handles.blackHoleMass;%function of angular momentum, mass etc.

handles.particleMatrix(handles.particleIndex, 6) = handles.STABLE_ORBIT_POTENTIAL ;
handles.data(handles.particleIndex, 6) = {handles.STABLE_ORBIT_POTENTIAL};
handles.particleMatrix(handles.particleIndex, 7) = handles.STABLE_ORBIT_RADIAL_DISTANCE ;
handles.data(handles.particleIndex, 7) = {handles.STABLE_ORBIT_RADIAL_DISTANCE};
set(handles.listOfParticles,'Data',handles.data)
guidata(hObject, handles);


function txtBlackHoleMass_Callback(hObject, eventdata, handles)
% hObject    handle to txtBlackHoleMass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtBlackHoleMass as text
%        str2double(get(hObject,'String')) returns contents of txtBlackHoleMass as a double
handles.blackHoleMass = str2double(get(hObject,'String'));
handles.particleMatrix(handles.particleIndex, 3) = handles.blackHoleMass ;
handles.data(handles.particleIndex, 3) = {handles.blackHoleMass};
set(handles.listOfParticles,'Data',handles.data)
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txtBlackHoleMass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtBlackHoleMass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
%if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%    set(hObject,'BackgroundColor','white');
%end

function txtBHAngularMomentum_Callback(hObject, eventdata, handles)
% hObject    handle to txtBHAngularMomentum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtBHAngularMomentum as text
%        str2double(get(hObject,'String')) returns contents of txtBHAngularMomentum as a double


% --- Executes during object creation, after setting all properties.
function txtBHAngularMomentum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtBHAngularMomentum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CreateCircularOrbit.
function CreateCircularOrbit_Callback(hObject, eventdata, handles)
% hObject    handle to CreateCircularOrbit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.blackHoleMass = particleMatrix(handles.particleIndex, 3);
handles.particleAngularMomentum = particleMatrix(handles.particleIndex, 4);
handles.CIRCULAR_ORBIT_RADIAL_DISTANCE = handles.particleAngularMomentum + handles.blackHoleMass; %function of angular momentum, mass etc.
handles.CIRCULAR_ORBIT_POTENTIAL= handles.particleAngularMomentum + handles.blackHoleMass;%function of angular momentum, mass etc.

handles.particleMatrix(handles.particleIndex, 6) = handles.CIRCULAR_ORBIT_POTENTIAL ;
handles.data(handles.particleIndex, 6) = {handles.CIRCULAR_ORBIT_POTENTIAL};
handles.particleMatrix(handles.particleIndex, 7) = handles.CIRCULAR_ORBIT_RADIAL_DISTANCE ;
handles.data(handles.particleIndex, 7) = {handles.CIRCULAR_ORBIT_RADIAL_DISTANCE};
set(handles.listOfParticles,'Data',handles.data)
guidata(hObject, handles);


% --- Executes on button press in optClassicalNewtonian.
function optClassicalNewtonian_Callback(hObject, eventdata, handles)
% hObject    handle to optClassicalNewtonian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.blackHoleType = 0;
handles.particleMatrix(handles.particleIndex, 2) = 0;
handles.data(handles.particleIndex, 2) = {'Newtonian'};
set(handles.listOfParticles,'Data',handles.data)
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of optClassicalNewtonian


% --- Executes on button press in optSchwarzchild.
function optSchwarzchild_Callback(hObject, eventdata, handles)
% hObject    handle to optSchwarzchild (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.blackHoleType = 1;
handles.particleMatrix(handles.particleIndex, 2) = 1;
handles.data(handles.particleIndex, 2) = {'Schwarzchild'};
set(handles.listOfParticles,'Data',handles.data)
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of optSchwarzchild


% --- Executes on button press in optKerr.
function optKerr_Callback(hObject, eventdata, handles)
% hObject    handle to optKerr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.blackHoleType = 2;
handles.particleMatrix(handles.particleIndex, 2) = 2;
handles.data(handles.particleIndex, 2) = {'Kerr'};
set(handles.listOfParticles,'Data',handles.data)
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of optKerr


% --- Executes on button press in btnDebug.
function btnDebug_Callback(hObject, eventdata, handles)
% hObject    handle to btnDebug (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles


% --- Executes on button press in particle1.
function particle1_Callback(hObject, eventdata, handles)
% hObject    handle to particle1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.particleIndex = 1  ;
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of particle1


% --- Executes on button press in particle2.
function particle2_Callback(hObject, eventdata, handles)
% hObject    handle to particle2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.particleIndex = 2;
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of particle2


% --- Executes on button press in particle3.
function particle3_Callback(hObject, eventdata, handles)
% hObject    handle to particle3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.particleIndex = 3;
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of particle3


% --- Executes on button press in particle4.
function particle4_Callback(hObject, eventdata, handles)
% hObject    handle to particle4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.particleIndex = 4 ;
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of particle4


% --- Executes on button press in particle5.
function particle5_Callback(hObject, eventdata, handles)
% hObject    handle to particle5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.particleIndex = 5;
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of particle5


% --- Executes on button press in particle6.
function particle6_Callback(hObject, eventdata, handles)
% hObject    handle to particle6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.particleIndex = 6;
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of particle6


% --- Executes on button press in isMassive.
function isMassive_Callback(hObject, eventdata, handles)
% hObject    handle to isMassive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.massData = {'no'};
handles.isMassive = get(hObject, 'Value');
handles.particleMatrix(handles.particleIndex, 1) = handles.isMassive;
if handles.isMassive == 1
    handles.massData = {'yes'};
else
    handles.massData = {'no'};
end
handles.data(handles.particleIndex, 1) = handles.massData
set(handles.listOfParticles,'Data',handles.data)
guidata(hObject, handles);
% Hint: get(yhObject,'Value') returns toggle state of isMassive

% run whenever fields are changed, bh type is changed
% Kerr case: run when slider moved, bhAngMom changed, isMassive changed
function updateEnergyPlot(handles)
    z = findMaxMins(handles)
    if min(z) == 0 
        radDist = linspace(0.1, max(z)*1.5, 1000);
    else
        radDist = linspace(min(z)*.5, max(z)*1.5, 1000);
    end
    
    Ueff = getPotential(handles, radDist);
    
    hold off;
    plot(radDist, Ueff, '-b');

% finds potential energy as a function of parameter r, the radial distance
function pot = getPotential(handles, r)
    type = handles.blackHoleType;
    m = handles.blackHoleMass;
    J = handles.BHangularMomentum;
    L = handles.particleAngularMomentum;
    E = handles.plotEnergy;
    K = handles.isMassive;
    
    if type == 0
        pot = L^2./(2*m*r.^2) - m./r;
    elseif type == 1
        pot = 1 - 2*m./r + L^2./r.^2 - 2*m*L^2./r.^3;
    else
        pot = -K*m./r + L^2./(2*r.^2) + (K - E^2)*(1 + J^2./(2*r.^2)) - m*(L - J*E/m)^2./r.^3;
    end
    
% Finds the maximum and/or the minimum of the effective potential
% depending on the type of black hole and field paramaters
function z = findMaxMins(handles)
    type = handles.blackHoleType;
    m = handles.blackHoleMass;
    J = handles.BHangularMomentum;
    L = handles.particleAngularMomentum;
    E = handles.plotEnergy;
    K = handles.isMassive;
    
    if type == 0
        z = L^2/m;
    elseif type == 1
        potRoots = roots([m, -L^2, 3*m*L^2]);
        z = [];
        
        for n = 1:length(potRoots)
            if imag(potRoots(n))
                z = [z, potRoots(n)];
            end
        end
    else type == 2
        potRoots = roots([K*m, -L^2 - (K - E^2)*(J/m)^2, 3*m*(L-J*E/m)^2]);
        z = [];
        
        for n = 1:length(potRoots)
            if imag(potRoots(n)) == 0
                z = [z, potRoots(n)];
            end
        end
    end
    
%TODO: set slider limits, xlim instead of axis,
