
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

% Last Modified by GUIDE v2.5 30-Mar-2017 11:46:40

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
handles.angularMomentum = 0;
handles.blackHoleMass = 0;
handles.blackHoleType = 1; % Enumeration
                           %    Classical Newtonian = 0
                           %    Schwarzchild = 1
                           %    Kerr = 2
handles.particleIndex = 0; % keeps track of the total number of particles in existence
handles.plotRadialDistance = 0;
handles.plotPotentialEnergy = 0;

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


% --- Executes on selection change in listParticles.
function listParticles_Callback(hObject, eventdata, handles)
% hObject    handle to listParticles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listParticles contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listParticles


% --- Executes during object creation, after setting all properties.
function listParticles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listParticles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in isMassive.
function isMassive_Callback(hObject, eventdata, handles)
% hObject    handle to isMassive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.isMassive = get(hObject, 'Value');

guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of isMassive



function txtAngularMomentum_Callback(hObject, eventdata, handles)
% hObject    handle to txtAngularMomentum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtAngularMomentum as text
%        str2double(get(hObject,'String')) returns contents of txtAngularMomentum as a double


% --- Executes during object creation, after setting all properties.
function txtAngularMomentum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtAngularMomentum (see GCBO)
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
plot(handles.plotRadialDistance, handles.plotPotentialEnergy, 'go');

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
handles.plotPotentialEnergy = get(hObject, 'Value');
plot(handles.plotRadialDistance , handles.plotPotentialEnergy, 'ro');

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


% --- Executes on button press in btnAddParticle.
function btnAddParticle_Callback(hObject, eventdata, handles)
% hObject    handle to btnAddParticle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.particleIndex = handles.particleIndex + 1

guidata(hObject, handles);

% --- Executes on button press in btnClearParticles.
function btnClearParticles_Callback(hObject, eventdata, handles)
% hObject    handle to btnClearParticles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.particleIndex = 0;

guidata(hObject, handles);

% --- Executes on button press in btnAddStableOrbit.
function btnAddStableOrbit_Callback(hObject, eventdata, handles)
% hObject    handle to btnAddStableOrbit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function txtBlackHoleMass_Callback(hObject, eventdata, handles)
% hObject    handle to txtBlackHoleMass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtBlackHoleMass as text
%        str2double(get(hObject,'String')) returns contents of txtBlackHoleMass as a double


% --- Executes during object creation, after setting all properties.
function txtBlackHoleMass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtBlackHoleMass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AddCircularOrbit.
function AddCircularOrbit_Callback(hObject, eventdata, handles)
% hObject    handle to AddCircularOrbit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in optClassicalNewtonian.
function optClassicalNewtonian_Callback(hObject, eventdata, handles)
% hObject    handle to optClassicalNewtonian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.blackHoleType = 0;

guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of optClassicalNewtonian


% --- Executes on button press in optSchwarzchild.
function optSchwarzchild_Callback(hObject, eventdata, handles)
% hObject    handle to optSchwarzchild (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.blackHoleType = 1;

guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of optSchwarzchild


% --- Executes on button press in optKerr.
function optKerr_Callback(hObject, eventdata, handles)
% hObject    handle to optKerr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.blackHoleType = 2;

guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of optKerr


% --- Executes on button press in btnDebug.
function btnDebug_Callback(hObject, eventdata, handles)
% hObject    handle to btnDebug (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles
