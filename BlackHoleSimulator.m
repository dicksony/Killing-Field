
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

% Last Modified by GUIDE v2.5 07-Apr-2017 10:48:51

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
handles.particleAngularMomentum = 4;
handles.BHangularMomentum = 0.5;
handles.CIRCULAR_ORBIT_POTENTIAL = 9;
handles.CIRCULAR_ORBIT_RADIAL_DISTANCE = 7;
handles.STABLE_ORBIT_POTENTIAL = 4;
handles.STABLE_ORBIT_RADIAL_DISTANCE = 5;
handles.blackHoleMass = 1;
handles.blackHoleType = 1; % Enumeration
                           %    Classical Newtonian = 0
                           %    Schwarzchild = 1
                           %    Kerr = 2
handles.particleIndex = 6; % keeps track of which particle is selected
handles.plotRadialDistance = 0;
handles.plotEnergy = 1;
handles.energy = 0;
%store values in particleMatrix for each particle:
%in order: 1 mass, 2 black hole type, 3 black hole mass, 4 particle angular momentum, 5 plot
%flag, 6 potential, 7 radial distance, 8 black hole angular momentum 
%9 SIZE(for plotting, not user input)
handles.particleMatrix = zeros(6, 10);
handles.particleMatrix(1,2) = 1;
handles.particleMatrix(1,5) = 1;
handles.particleMatrix(:,9) = 30000;
handles.particleMatrix(:,10) = 0.03;
handles.data = cell(6, 8);
handles.data(:) = {''};
handles.data(:,1) = {'no'};
handles.data(2:6,5) = {'no'};
handles.data(1,5) = {'yes'};
handles.data(:, 3:4) = {'0'};
handles.data(:, 6:8) = {'0'};
handles.data(1, 2) = {'Schwarzschild'};
handles.data(2:6, 2) = {'Newtonian'};
set(handles.listOfParticles,'Data',handles.data);
handles.plot_cell_labels = cellstr('labels');
handles.all_labels = ['Particle 1'; 'Particle 2'; 'Particle 3'; 'Particle 4'; 'Particle 5'; 'Particle 6'];
handles.cell_labels = cellstr(handles.all_labels);
 
% Update handles structure
guidata(hObject, handles);
updateEnergyPlot(handles);

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

updateEnergyPlot(handles);

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
minPot = getPotential(handles, [handles.plotRadialDistance]);
if handles.plotEnergy < minPot
    set(handles.sliderEnergy, 'Value', minPot);
    handles.plotEnergy = minPot;
end

plot(handles.plotRadialDistance, handles.plotEnergy, 'ro');
hold off
handles.particleMatrix(handles.particleIndex, 7) = handles.plotRadialDistance ;

handles.data(handles.particleIndex, 7) = {handles.plotRadialDistance};
set(handles.listOfParticles,'Data',handles.data)
guidata(hObject, handles);

updateEnergyPlot(handles);
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
function sliderEnergy_Callback(hObject, eventdata, handles)
% hObject    handle to sliderEnergy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
minPot = getPotential(handles, [handles.plotRadialDistance]);
if get(hObject, 'Value') < minPot
    set(hObject, 'Value', minPot);
end
handles.plotEnergy = get(hObject, 'Value');

plot(handles.plotRadialDistance , handles.plotEnergy, 'ro');
hold off;

handles.particleMatrix(handles.particleIndex, 6) = handles.plotEnergy ;
handles.data(handles.particleIndex, 6) = {handles.plotEnergy};
set(handles.listOfParticles,'Data',handles.data)
guidata(hObject,handles);

updateEnergyPlot(handles);
    

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sliderEnergy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderEnergy (see GCBO)
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
%in order: 1 mass, 2 black hole type, 3 black hole mass, 4 particle angular momentum, 5 plot
%flag, 6 potential, 7 radial distance, 8 black hole angular momentum
figure(1)
    handles.plot_cell_labels = cellstr('labels');
    for i = 1:6
        if handles.particleMatrix(i, 5) == 1
            handles.plot_cell_labels = [handles.plot_cell_labels; handles.cell_labels(i)];
            if handles.particleMatrix(i, 2) == 1 || handles.particleMatrix(i,2) == 2%if Swarzchild or Kerr
                if handles.particleMatrix(i, 1) == 1 %if massive
                    handles.energy = sqrt(handles.particleMatrix(i, 6)*2 +1);
                else
                    handles.energy = sqrt(handles.particleMatrix(i, 6)*2);
                end
            else
                handles.energy = handles.particleMatrix(i, 6);
            end    
            [radius, angle, time] = BH(handles.particleMatrix(i, 2), ...
                                       handles.particleMatrix(i, 8), ...
                                       handles.particleMatrix(i, 3), ...
                                       handles.particleMatrix(i, 7), ...
                                       handles.energy, ...
                                       handles.particleMatrix(i, 4), ...
                                       handles.particleMatrix(i, 1), ...
                                       handles.particleMatrix(i, 9), ...
                                       handles.particleMatrix(i, 10), 1);
            if ~isempty(time)
                polarplot(angle,radius);%Display data
                hold on
            end
        end
    end 
    legend('show')
    handles.plot_cell_labels = handles.plot_cell_labels(2:length(handles.plot_cell_labels));
    legend(handles.plot_cell_labels)


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

updateEnergyPlot(handles);

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

handles.BHangularMomentum = str2double(get(hObject,'String'));
handles.particleMatrix(handles.particleIndex, 8) = handles.BHangularMomentum ;
handles.data(handles.particleIndex, 8) = {handles.BHangularMomentum};
set(handles.listOfParticles,'Data',handles.data)
guidata(hObject, handles);

    if handles.blackHoleType == 2
        updateEnergyPlot(handles);
    end 
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


% --- Executes on button press in DemoOrbit2.
function DemoOrbit2_Callback(hObject, eventdata, handles)
% hObject    handle to DemoOrbit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%in order: 1 mass, 2 black hole type, 3 black hole mass, 4 particle angular momentum, 5 plot
%flag, 6 potential, 7 radial distance, 8 black hole angular momentum
handles.particleMatrix(handles.particleIndex, 1) = 1;
handles.particleMatrix(handles.particleIndex, 2) = 1;
handles.particleMatrix(handles.particleIndex, 3) = 1;
handles.particleMatrix(handles.particleIndex, 4) = 4.2;
handles.particleMatrix(handles.particleIndex, 5) = 1;
handles.particleMatrix(handles.particleIndex, 6) = (((1)^2)-1)/2;
handles.particleMatrix(handles.particleIndex, 7) = 40;
handles.particleMatrix(handles.particleIndex, 8) = 0;
handles.particleMatrix(handles.particleIndex, 9) = 20000;
handles.particleMatrix(handles.particleIndex, 10) = 0.03;

handles.data(handles.particleIndex, 1) = {'yes'};
handles.data(handles.particleIndex, 2) = {'Schwarzschild'};
handles.data(handles.particleIndex, 3) = {handles.particleMatrix(handles.particleIndex, 3)};
handles.data(handles.particleIndex, 4) = {handles.particleMatrix(handles.particleIndex, 4)};
handles.data(handles.particleIndex, 5) = {'yes'};
handles.data(handles.particleIndex, 6) = {handles.particleMatrix(handles.particleIndex, 6)};
handles.data(handles.particleIndex, 7) = {handles.particleMatrix(handles.particleIndex, 7)};
handles.data(handles.particleIndex, 8) = {handles.particleMatrix(handles.particleIndex, 8)};
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

updateEnergyPlot(handles);
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

updateEnergyPlot(handles);
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

updateEnergyPlot(handles);
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
handles.data(handles.particleIndex, 1) = handles.massData;
set(handles.listOfParticles,'Data',handles.data)
guidata(hObject, handles);

updateEnergyPlot(handles);
% Hint: get(yhObject,'Value') returns toggle state of isMassive

function updateEnergyPlot(handles)
    numPoints = 20000;
    
    z = findMaxMins(handles);
    
    if length(z) ~= 0
        radDist = linspace(min(z)/2, max(z)*5, numPoints);
        rmin = radDist(1);
        rmax = radDist(numPoints);
    else
       radDist = linspace(0.1, 10, numPoints);
       rmin = radDist(1);
       rmax = radDist(numPoints);
    end
    sliderRad = get(handles.sliderRadialDistance, 'Value');
    set(handles.sliderRadialDistance, 'min', rmin);
    set(handles.sliderRadialDistance, 'max', rmax);
    if sliderRad < rmin || sliderRad > rmax
        set(handles.sliderRadialDistance, 'Value', rmax)
        handles.plotRadialDistance = get(handles.sliderRadialDistance, 'Value');
    end
    
    Ueff = getPotential(handles, radDist);
    Emin = min(Ueff)-.1;
    Emax = max(Ueff)+.1;
    if Emax < .05
        Emax = .05;
    end
    sliderE = get(handles.sliderEnergy, 'Value');
    set(handles.sliderEnergy, 'min', Emin);
    set(handles.sliderEnergy, 'max', Emax);
    minPot = getPotential(handles, [handles.plotRadialDistance]);
    if sliderE < minPot || sliderE > Emax
        set(handles.sliderEnergy, 'Value', minPot);
    end
    handles.plotEnergy = get(handles.sliderEnergy,'Value');
    
    plot(handles.plotRadialDistance, handles.plotEnergy, 'ro');
    hold on;
    plot(radDist, Ueff, '-b');
    xlim([rmin, rmax]);
    ylim([Emin, Emax]);
    hold off;

% finds potential energy as a function of parameter r, the radial distance
function pot = getPotential(handles, r)
    type = handles.blackHoleType;
    m = handles.blackHoleMass;
    J = handles.BHangularMomentum;
    L = handles.particleAngularMomentum;
    E = handles.plotEnergy;
    K = handles.isMassive;
    
    if type == 0
        pot = K*(L^2./(2*m*r.^2) - m./r);
    elseif type == 1
        pot = -K*m./r + L^2./(2*r.^2) - m*L^2./r.^3;
    else
        pot = -K*m./r + L^2./(2*r.^2) + (K - E^2)*(1 + J^2./(m^2*r.^2))/2 - m*(L - J*E/m)^2./r.^3;
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
        potRoots = roots([K*m, -L^2, 3*m*L^2]);
        z = [];
        for n = 1:length(potRoots)
            if imag(potRoots(n)) == 0
                z = [z, potRoots(n)];
            end
        end
    else
        potRoots = roots([K*m, -L^2 - (K - E^2)*(J/m)^2, 3*m*(L-J*E/m)^2]);
        z = [];
        for n = 1:length(potRoots)
            if imag(potRoots(n)) == 0
                z = [z, potRoots(n)];
            end
        end
    end


% --- Executes on button press in DemoOrbit3.
function DemoOrbit3_Callback(hObject, eventdata, handles)
% hObject    handle to DemoOrbit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%in order: 1 mass, 2 black hole type, 3 black hole mass, 4 particle angular momentum, 5 plot
%flag, 6 potential, 7 radial distance, 8 black hole angular momentum
handles.particleMatrix(handles.particleIndex, 1) = 1;
handles.particleMatrix(handles.particleIndex, 2) = 0;
handles.particleMatrix(handles.particleIndex, 3) = 1;
handles.particleMatrix(handles.particleIndex, 4) = 1;
handles.particleMatrix(handles.particleIndex, 5) = 1;
handles.particleMatrix(handles.particleIndex, 6) = 1;
handles.particleMatrix(handles.particleIndex, 7) = 1;
handles.particleMatrix(handles.particleIndex, 8) = 0;
handles.particleMatrix(handles.particleIndex, 9) = 30000;
handles.particleMatrix(handles.particleIndex, 10) = 0.01;

handles.data(handles.particleIndex, 1) = {'yes'};
handles.data(handles.particleIndex, 2) = {'Newtonian'};
handles.data(handles.particleIndex, 3) = {handles.particleMatrix(handles.particleIndex, 3)};
handles.data(handles.particleIndex, 4) = {handles.particleMatrix(handles.particleIndex, 4)};
handles.data(handles.particleIndex, 5) = {'yes'};
handles.data(handles.particleIndex, 6) = {handles.particleMatrix(handles.particleIndex, 6)};
handles.data(handles.particleIndex, 7) = {handles.particleMatrix(handles.particleIndex, 7)};
handles.data(handles.particleIndex, 8) = {handles.particleMatrix(handles.particleIndex, 8)};
set(handles.listOfParticles,'Data',handles.data)
guidata(hObject, handles);


% --- Executes on button press in DemoOrbit5.
function DemoOrbit5_Callback(hObject, eventdata, handles)
% hObject    handle to DemoOrbit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%in order: 1 mass, 2 black hole type, 3 black hole mass, 4 particle angular momentum, 5 plot
%flag, 6 potential, 7 radial distance, 8 black hole angular momentum
handles.particleMatrix(1, 1) = 1;
handles.particleMatrix(1, 2) = 2;
handles.particleMatrix(1, 3) = 1;
handles.particleMatrix(1, 4) = 4;
handles.particleMatrix(1, 5) = 1;
handles.particleMatrix(1, 6) = (((1.001)^2)-1)/2;
handles.particleMatrix(1, 7) = 9.606;
handles.particleMatrix(1, 8) = 0.05;
handles.particleMatrix(1, 9) = 7000;
handles.particleMatrix(1, 10) = 0.08;

handles.data(1, 1) = {'yes'};
handles.data(1, 2) = {'Kerr'};
handles.data(1, 3) = {handles.particleMatrix(1, 3)};
handles.data(1, 4) = {handles.particleMatrix(1, 4)};
handles.data(1, 5) = {'yes'};
handles.data(1, 6) = {handles.particleMatrix(1, 6)};
handles.data(1, 7) = {handles.particleMatrix(1, 7)};
handles.data(1, 8) = {handles.particleMatrix(1, 8)};
set(handles.listOfParticles,'Data',handles.data)

handles.particleMatrix(2, 1) = 1;
handles.particleMatrix(2, 2) = 2;
handles.particleMatrix(2, 3) = 1;
handles.particleMatrix(2, 4) = 4;
handles.particleMatrix(2, 5) = 1;
handles.particleMatrix(2, 6) = (((1.001)^2)-1)/2;
handles.particleMatrix(2, 7) = 9.606;
handles.particleMatrix(2, 8) = -0.05;
handles.particleMatrix(2, 9) = 7000;
handles.particleMatrix(2, 10) = 0.08;

handles.data(2, 1) = {'yes'};
handles.data(2, 2) = {'Kerr'};
handles.data(2, 3) = {handles.particleMatrix(2, 3)};
handles.data(2, 4) = {handles.particleMatrix(2, 4)};
handles.data(2, 5) = {'yes'};
handles.data(2, 6) = {handles.particleMatrix(2, 6)};
handles.data(2, 7) = {handles.particleMatrix(2, 7)};
handles.data(2, 8) = {handles.particleMatrix(2, 8)};

handles.particleMatrix(3, 5) = 0;
handles.particleMatrix(4, 5) = 0;
handles.particleMatrix(5, 5) = 0;
handles.particleMatrix(6, 5) = 0;
set(handles.listOfParticles,'Data',handles.data)
guidata(hObject, handles);


% --- Executes on button press in DemoOrbit4.
function DemoOrbit4_Callback(hObject, eventdata, handles)
% hObject    handle to DemoOrbit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%in order: 1 mass, 2 black hole type, 3 black hole mass, 4 particle angular momentum, 5 plot
%flag, 6 potential, 7 radial distance, 8 black hole angular momentum
handles.particleMatrix(handles.particleIndex, 1) = 1;
handles.particleMatrix(handles.particleIndex, 2) = 0;
handles.particleMatrix(handles.particleIndex, 3) = 1;
handles.particleMatrix(handles.particleIndex, 4) = sqrt(120);
handles.particleMatrix(handles.particleIndex, 5) = 1;
handles.particleMatrix(handles.particleIndex, 6) = 0.55;
handles.particleMatrix(handles.particleIndex, 7) = 10;
handles.particleMatrix(handles.particleIndex, 8) = 0;
handles.particleMatrix(handles.particleIndex, 9) = 30000;
handles.particleMatrix(handles.particleIndex, 10) = 0.01;

handles.data(handles.particleIndex, 1) = {'yes'};
handles.data(handles.particleIndex, 2) = {'Newtonian'};
handles.data(handles.particleIndex, 3) = {handles.particleMatrix(handles.particleIndex, 3)};
handles.data(handles.particleIndex, 4) = {handles.particleMatrix(handles.particleIndex, 4)};
handles.data(handles.particleIndex, 5) = {'yes'};
handles.data(handles.particleIndex, 6) = {handles.particleMatrix(handles.particleIndex, 6)};
handles.data(handles.particleIndex, 7) = {handles.particleMatrix(handles.particleIndex, 7)};
handles.data(handles.particleIndex, 8) = {handles.particleMatrix(handles.particleIndex, 8)};
set(handles.listOfParticles,'Data',handles.data)
guidata(hObject, handles);


% --- Executes on button press in DemoOrbit6.
function DemoOrbit6_Callback(hObject, eventdata, handles)
% hObject    handle to DemoOrbit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%in order: 1 mass, 2 black hole type, 3 black hole mass, 4 particle angular momentum, 5 plot
%flag, 6 potential, 7 radial distance, 8 black hole angular momentum
handles.particleMatrix(handles.particleIndex, 1) = 1;
handles.particleMatrix(handles.particleIndex, 2) = 1;
handles.particleMatrix(handles.particleIndex, 3) = 1;
handles.particleMatrix(handles.particleIndex, 4) = 4;
handles.particleMatrix(handles.particleIndex, 5) = 1;
handles.particleMatrix(handles.particleIndex, 6) = -0.0333;
handles.particleMatrix(handles.particleIndex, 7) = 8.6001;
handles.particleMatrix(handles.particleIndex, 8) = 0.5;

handles.data(handles.particleIndex, 1) = {'yes'};
handles.data(handles.particleIndex, 2) = {'Schwarzschild'};
handles.data(handles.particleIndex, 3) = {handles.particleMatrix(handles.particleIndex, 3)};
handles.data(handles.particleIndex, 4) = {handles.particleMatrix(handles.particleIndex, 4)};
handles.data(handles.particleIndex, 5) = {'yes'};
handles.data(handles.particleIndex, 6) = {handles.particleMatrix(handles.particleIndex, 6)};
handles.data(handles.particleIndex, 7) = {handles.particleMatrix(handles.particleIndex, 7)};
handles.data(handles.particleIndex, 8) = {handles.particleMatrix(handles.particleIndex, 8)};
set(handles.listOfParticles,'Data',handles.data)
guidata(hObject, handles);


% --- Executes on button press in DemoOrbit1.
function DemoOrbit1_Callback(hObject, eventdata, handles)
% hObject    handle to DemoOrbit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%in order: 1 mass, 2 black hole type, 3 black hole mass, 4 particle angular momentum, 5 plot
%flag, 6 potential, 7 radial distance, 8 black hole angular momentum
%9 size 10 timestep
handles.particleMatrix(handles.particleIndex, 1) = 1;
handles.particleMatrix(handles.particleIndex, 2) = 1;
handles.particleMatrix(handles.particleIndex, 3) = 1;
handles.particleMatrix(handles.particleIndex, 4) = 80;
handles.particleMatrix(handles.particleIndex, 5) = 1;
handles.particleMatrix(handles.particleIndex, 6) = (((7)^2)-1)/2;
handles.particleMatrix(handles.particleIndex, 7) = 50;
handles.particleMatrix(handles.particleIndex, 8) = 0;
handles.particleMatrix(handles.particleIndex, 9) = 700;
handles.particleMatrix(handles.particleIndex, 10) = 0.03;

handles.data(handles.particleIndex, 1) = {'yes'};
handles.data(handles.particleIndex, 2) = {'Schwarzschild'};
handles.data(handles.particleIndex, 3) = {handles.particleMatrix(handles.particleIndex, 3)};
handles.data(handles.particleIndex, 4) = {handles.particleMatrix(handles.particleIndex, 4)};
handles.data(handles.particleIndex, 5) = {'yes'};
handles.data(handles.particleIndex, 6) = {handles.particleMatrix(handles.particleIndex, 6)};
handles.data(handles.particleIndex, 7) = {handles.particleMatrix(handles.particleIndex, 7)};
handles.data(handles.particleIndex, 8) = {handles.particleMatrix(handles.particleIndex, 8)};
set(handles.listOfParticles,'Data',handles.data)
guidata(hObject, handles);
