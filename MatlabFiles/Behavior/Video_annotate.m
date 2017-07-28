function varargout = Video_annotate(varargin)
% VIDEO_ANNOTATE M-file for Video_annotate.fig
%      VIDEO_ANNOTATE, by itself, creates a new VIDEO_ANNOTATE or raises the existing
%      singleton*.
%
%      H = VIDEO_ANNOTATE returns the handle to a new VIDEO_ANNOTATE or the handle to
%      the existing singleton*.
%
%      VIDEO_ANNOTATE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIDEO_ANNOTATE.M with the given input arguments.
%
%      VIDEO_ANNOTATE('Property','Value',...) creates a new VIDEO_ANNOTATE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before confocal_idsG_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Video_annotate_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Video_annotate

% Last Modified by GUIDE v2.5 26-Oct-2014 18:54:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Video_annotate_OpeningFcn, ...
                   'gui_OutputFcn',  @Video_annotate_OutputFcn, ...
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


% --- Executes just before Video_annotate is made visible.
function Video_annotate_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Video_annotate (see VARARGIN)
clearvars DATA SSWSG
% Choose default command line output for Video_annotate
handles.output = hObject;

handles.path = 'C:\Users\Greenspan Lab\AcquisitionData';

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes Video_annotate wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% --- Outputs from this function are returned to the command line.
function varargout = Video_annotate_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in open_file.
function open_file_Callback(hObject, eventdata, handles)
% hObject    handle to open_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% free some space

clear DATA handles.N handles.Tserie handles.filename handles.nom handles.path

clearvars DATA SSWSG


global DATA 

%light on/off a, Walk w, Idle e, Groom g, goes right r, goes left d, goes
%straight s, Other t

global Behavior

[nom,chemin] = uigetfile('*.avi','pick a avi file',handles.path, 'MultiSelect', 'on');
handles.path = chemin;
handles.filename = strcat(chemin,nom);
handles.nom = nom;
fichier = strcat(chemin,nom)

mov = VideoReader(fichier);

DATA=[];

ii = 1;

while hasFrame(mov)
   A=readFrame(mov); 
   DATA(:,:,ii)=A(:,:,1);
   ii = ii+1;
end

num_images=ii-1;

% for i=1:num_images
% A=imread(fichier,i);
% DATA(:,:,i)=A(:,:,1); 
% end

Behavior=zeros(7,num_images);

S2=size(DATA);
Ximage=S2(2);
Yimage=S2(1);


handles.N = 1;
handles.Yimage = Yimage;
handles.Ximage = Ximage;

handles.Tserie = num_images;

handles.timeserie = 1;

imagesc(DATA(2:Yimage,2:Ximage,1))
colormap(gray)
set(handles.slider2,'Max',num_images,'Min',1,'Value',1,'SliderStep',[1/(num_images-1) 1/(num_images-1)]);

    
%% Update handles structure

guidata(hObject, handles);



%% --- Executes on button press in save_pushbutton.
function save_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to save_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global DATA Behavior

[name,path] = uiputfile('*.dat','nom du film',char(handles.filename(1)));%
nom = strcat(path,name);

save (nom,'Behavior','-ascii','-tabs')




%% --- Executes on button press in green_checkbox.
function green_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to green_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA

 handles.green_checkbox_channel = get(hObject,'Value') ;
 

% Update handles structure
guidata(hObject, handles);




% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
%light on/off a, Walk w, Idle e, Groom g, goes right r, goes left d, goes
%straight s, t other

global Behavior

switch eventdata.Key
  case 's'
    Behavior(1,handles.N)=1;
    case 'd'
    Behavior(2,handles.N)=1; 
    case 'w'
    Behavior(3,handles.N)=1;     
    case 'e'
    Behavior(4,handles.N)=1;    
    case 'g'
    Behavior(5,handles.N)=1;    
    case 'r'
    Behavior(6,handles.N)=1;    
    case 't'
    Behavior(7,handles.N)=1;
    case 'c'
    Behavior(:,handles.N)=0;
        end


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

%light on/off a, Walk w, Idle e, Groom g, goes right r, goes left d, goes
%straight s, t other
global Behavior

switch eventdata.Key
  case 's'
    Behavior(1,handles.N)=1;
    case 'd'
    Behavior(2,handles.N)=1; 
    case 'w'
    Behavior(3,handles.N)=1;     
    case 'e'
    Behavior(4,handles.N)=1;    
    case 'g'
    Behavior(5,handles.N)=1;    
    case 'r'
    Behavior(6,handles.N)=1;    
    case 't'
    Behavior(7,handles.N)=1;
    case 'c'
    Behavior(:,handles.N)=0;
        end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global DATA

handles.N = round(get(hObject,'Value'));
S2=size(DATA);
Ximage=S2(2);
Yimage=S2(1);


 imagesc(DATA(2:Yimage,2:Ximage,handles.N ))
     colormap(gray)


 set(handles.text5,'String',strcat(num2str(handles.N),' / ',num2str(handles.Tserie)));
 % Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
