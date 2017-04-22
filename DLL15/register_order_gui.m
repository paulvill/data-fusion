function varargout = register_order_gui(varargin)
% REGISTER_ORDER_GUI MATLAB code for register_order_gui.fig
%      REGISTER_ORDER_GUI, by itself, creates a new REGISTER_ORDER_GUI or raises the existing
%      singleton*.
%
%      H = REGISTER_ORDER_GUI returns the handle to a new REGISTER_ORDER_GUI or the handle to
%      the existing singleton*.
%
%      REGISTER_ORDER_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REGISTER_ORDER_GUI.M with the given input arguments.
%
%      REGISTER_ORDER_GUI('Property','Value',...) creates a new REGISTER_ORDER_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
% %      applied to the GUI before register_order_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to register_order_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help register_order_gui

% Last Modified by GUIDE v2.5 02-Feb-2015 10:19:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @register_order_gui_OpeningFcn, ...
    'gui_OutputFcn',  @register_order_gui_OutputFcn, ...
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


% --- Executes just before register_order_gui is made visible.
function register_order_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to register_order_gui (see VARARGIN)

% Choose default command line output for register_order_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes register_order_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = register_order_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function output_image_dir_Callback(hObject, eventdata, handles)
% hObject    handle to output_image_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output_image_dir as text
%        str2double(get(hObject,'String')) returns contents of output_image_dir as a double

handles.output_image_dir = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function output_image_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_image_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_images.
function save_images_Callback(hObject, eventdata, handles)
% hObject    handle to save_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'images_analyzed')
    msgbox('Images have not been analyzed.')
    return
end

if handles.reread_images
    choice = questdlg('Initial image settings have changed, but images have not been reread. Are you sure you want to continue?', ...
        'Reread images', ...
        'Yes','No','No');
    % Handle response
    switch choice
        case 'Yes'
            
        case 'No'
            return
    end
end
if handles.reapply_image_functions
    choice = questdlg('Image function parameters have changed, but have not been reapplied. Are you sure you want to continue?', ...
        'Repply image functions', ...
        'Yes','No','No');
    % Handle response
    switch choice
        case 'Yes'
            
        case 'No'
            return
    end
end
if handles.recalc_pairwise_rotations
    choice = questdlg('Number of pairwise rotations or image set/parameters have changed, but images have not been reanalyzed. Are you sure you want to continue?', ...
        'Reanalyze images', ...
        'Yes','No','No');
    % Handle response
    switch choice
        case 'Yes'
            
        case 'No'
            return
    end
end
if handles.rerun_dmaps
    choice = questdlg('Kernel scale or image set/parameters have changed, but images have not been reanalyzed. Are you sure you want to continue?', ...
        'Reanalyze images', ...
        'Yes','No','No');
    % Handle response
    switch choice
        case 'Yes'
            
        case 'No'
            return
    end
end
if handles.reregister_image_set
    choice = questdlg('The selected images to analyze (raw vs. preprocessed) has changed, but images have not been reanalyzed. Are you sure you want to continue?', ...
        'Reread images', ...
        'Yes','No','No');
    % Handle response
    switch choice
        case 'Yes'
            
        case 'No'
            return
    end
end

try
    save_images(handles.images_analyzed, handles.dim, handles.output_image_dir, handles.image_name, handles.image_ext, handles.stack_name)
catch ME
    if strcmp(ME.identifier, 'MATLAB:nonExistentField')
        msgbox('Some required fields are not populated. Please check inputs.')
        return
        %     elseif strcmp(ME.identifier, 'MATLAB:imagesci:imread:fileDoesNotExist')
    else
        msgbox(sprintf('ERROR: %s', ME.identifier))
        return
    end
end

function number_images_Callback(hObject, eventdata, handles)
% hObject    handle to number_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_images as text
%        str2double(get(hObject,'String')) returns contents of number_images as a double

handles.reread_images = true;

nimages = str2double(get(hObject,'String'));
handles.nimages = nimages;
if isnan(handles.nimages) || handles.nimages < 1 || mod(handles.nimages, 1) ~= 0
    msgbox('Invalid number of images.')
    handles = rmfield(handles, 'nimages');
    set(hObject, 'String', '');
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function number_images_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function image_extension_Callback(hObject, eventdata, handles)
% hObject    handle to image_extension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of image_extension as text
%        str2double(get(hObject,'String')) returns contents of image_extension as a double

handles.reread_images = true;

handles.image_ext = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function image_extension_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_extension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in read_images.
function read_images_Callback(hObject, eventdata, handles)
% hObject    handle to read_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%h = msgbox('Please wait...');

try
    [handles.images_raw, handles.nchannels] = read_images(handles.image_dir, handles.image_name, handles.image_ext, handles.stack_name, handles.nimages, handles.nstack, handles.dim);
catch ME
    if strcmp(ME.identifier, 'MATLAB:nonExistentField')
        msgbox('Some required fields are not populated. Please check inputs.')
        return
        %     elseif strcmp(ME.identifier, 'MATLAB:imagesci:imread:fileDoesNotExist')
    else
        msgbox(sprintf('ERROR: %s', ME.identifier))
        return
    end
end

handles.reread_images = false;
handles.reapply_image_functions = true;
handles.recalc_pairwise_rotations = true;
handles.rerun_dmaps = true;

handles.channel_weight = zeros(handles.nchannels, 1);
handles.channel_blur = zeros(handles.nchannels, 1);
handles.channel_mean_center = false(handles.nchannels, 1);
handles.channel_normalize = false(handles.nchannels, 1);

for i=1:handles.nchannels
    handles.channel_weight(i) = get(handles.weight_slider(i), 'Value');
    handles.channel_blur(i) = get(handles.blur_slider(i), 'Value');
    handles.channel_mean_center(i) = get(handles.mean_center_checkbox(i), 'Value');
    handles.channel_normalize(i) = get(handles.normalize_checkbox(i), 'Value');
    
    set(handles.weight_slider(i), 'enable','on')
    set(handles.blur_slider(i), 'enable','on')
    set(handles.weight_numbers(i), 'enable','on')
    set(handles.blur_numbers(i), 'enable','on')
    set(handles.normalize_checkbox(i), 'enable','on')
    set(handles.mean_center_checkbox(i), 'enable','on')
    set(handles.color_label(i), 'enable','on')
end
for i=handles.nchannels+1:3
    set(handles.weight_slider(i), 'enable','off')
    set(handles.blur_slider(i), 'enable','off')
    set(handles.weight_numbers(i), 'enable','off')
    set(handles.blur_numbers(i), 'enable','off')
    set(handles.normalize_checkbox(i), 'enable','off')
    set(handles.mean_center_checkbox(i), 'enable','off')
    set(handles.normalize_checkbox(i), 'value',0)
    set(handles.mean_center_checkbox(i), 'value',0)
    set(handles.color_label(i), 'enable','off')
end

guidata(hObject,handles);


% --- Executes on selection change in image_dim.
function image_dim_Callback(hObject, eventdata, handles)
% hObject    handle to image_dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns image_dim contents as cell array
%        contents{get(hObject,'Value')} returns selected item from image_dim

str = get(hObject, 'String');
val = get(hObject,'Value');

handles.reread_images = true;

% Set current data to the selected data set.
switch str{val};
    case '2D'
        handles.dim = 2;
        set(handles.stack_prefix, 'enable','off')
        set(handles.number_stack, 'enable','off')
        set(handles.stack_prefix, 'String','')
        set(handles.number_stack, 'String','')
        set(handles.register_button, 'enable', 'on')
        set(handles.order_button, 'enable', 'on')
        set(handles.register_order_button, 'enable', 'on')
        set(handles.register_order_zstack_button, 'enable', 'off')
    case '3D'
        handles.dim = 3;
        set(handles.stack_prefix, 'enable','on')
        set(handles.number_stack, 'enable','on')
        set(handles.register_button, 'enable', 'off')
        set(handles.order_button, 'enable', 'on')
        set(handles.register_order_button, 'enable', 'off')
        set(handles.register_order_zstack_button, 'enable', 'on')
end
% Save the handles structure.
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function image_dim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_dim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.dim = 2;
guidata(hObject,handles)


function image_dir_Callback(hObject, eventdata, handles)
% hObject    handle to image_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of image_dir as text
%        str2double(get(hObject,'String')) returns contents of image_dir as a double

handles.reread_images = true;

handles.image_dir = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function image_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function image_prefix_Callback(hObject, eventdata, handles)
% hObject    handle to image_prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of image_prefix as text
%        str2double(get(hObject,'String')) returns contents of image_prefix as a double

handles.reread_images = true;

handles.image_name = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function image_prefix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function number_stack_Callback(hObject, eventdata, handles)
% hObject    handle to number_stack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_stack as text
%        str2double(get(hObject,'String')) returns contents of number_stack as a double

handles.reread_images = true;

handles.nstack = str2double(get(hObject,'String'));
if isnan(handles.nstack) || handles.nstack < 1 || mod(handles.nstack, 1) ~= 0
    msgbox('Invalid number of images in stack.')
    handles = rmfield(handles, 'nstack');
    set(hObject, 'String', '');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function number_stack_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_stack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'enable','off');
handles.number_stack = hObject;
handles.nstack = 0;
guidata(hObject,handles);


function stack_prefix_Callback(hObject, eventdata, handles)
% hObject    handle to stack_prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stack_prefix as text
%        str2double(get(hObject,'String')) returns contents of stack_prefix as a double

handles.reread_images = true;

handles.stack_name = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function stack_prefix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stack_prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'enable','off');
handles.stack_prefix = hObject;
handles.stack_name = '';
guidata(hObject,handles);



function number_pixels_Callback(hObject, eventdata, handles)
% hObject    handle to number_pixels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_pixels as text
%        str2double(get(hObject,'String')) returns contents of number_pixels as a double

handles.reapply_image_functions = true;

handles.npixels = str2double(get(hObject,'String'));
if isnan(handles.npixels) || handles.npixels < 1 || mod(handles.npixels, 1) ~= 0
    msgbox('Invalid number of pixels.')
    handles = rmfield(handles, 'npixels');
    set(hObject, 'String', '');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function number_pixels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_pixels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in show_images1.
function show_images1_Callback(hObject, eventdata, handles)
% hObject    handle to show_images1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'images_raw')
    msgbox('Images have not been read.')
    return
end

plot_images(handles.images_raw, handles.dim)

% --- Executes on button press in apply_image_functions.
function apply_image_functions_Callback(hObject, eventdata, handles)
% hObject    handle to apply_image_functions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'images_raw')
    msgbox('No images have been read.')
    return
end

if handles.reread_images
    choice = questdlg('Initial image settings have changed, but images have not been reread. Are you sure you want to continue?', ...
        'Reread images', ...
        'Yes','No','No');
    % Handle response
    switch choice
        case 'Yes'
            
        case 'No'
            return
    end
end

try
    handles.images = apply_image_functions(handles.images_raw, handles.npixels, handles.dim, handles.channel_weight, handles.channel_blur, handles.channel_normalize, handles.channel_mean_center, handles.resize_image);
    
    handles.reapply_image_functions = false;
    handles.recalc_pairwise_rotations = true;
    handles.rerun_dmaps = true;
    
    guidata(hObject, handles);
catch ME
    if strcmp(ME.identifier, 'MATLAB:nonExistentField')
        msgbox('Some required fields are not populated. Please check inputs.')
        return
    else
        msgbox(sprintf('ERROR: %s', ME.identifier))
        return
    end
end

% --- Executes on button press in show_images2.
function show_images2_Callback(hObject, eventdata, handles)
% hObject    handle to show_images2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'images')
    msgbox('Images have not been preprocessed.')
    return
end

plot_images(handles.images, handles.dim)

function kernel_scale_Callback(hObject, eventdata, handles)
% hObject    handle to kernel_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kernel_scale as text
%        str2double(get(hObject,'String')) returns contents of kernel_scale as a double

handles.rerun_dmaps = true;

handles.eps_scale = str2double(get(hObject,'String')) ;

if isnan(handles.eps_scale)
    msgbox('Invalid kernel scale.')
    handles = rmfield(handles, 'eps_scale');
    set(hObject, 'String', '');
elseif handles.eps_scale < 0
    msgbox('Invalid kernel scale: kernel scale must be greater than 0.')
    handles = rmfield(handles, 'eps_scale');
    set(hObject, 'String', '');
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function kernel_scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kernel_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', '0.5');
handles.eps_scale = 0.5;
handles.rerun_dmaps = true;

guidata(hObject, handles);


% --- Executes on button press in register_order_vdm.
function register_order_vdm_Callback(hObject, eventdata, handles)
% hObject    handle to register_order_vdm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'images')
    msgbox('Images have not been preprocessed.')
    return
end

if handles.reread_images
    choice = questdlg('Initial image settings have changed, but images have not been reread. Are you sure you want to continue?', ...
        'Reread images', ...
        'Yes','No','No');
    % Handle response
    switch choice
        case 'Yes'
            
        case 'No'
            return
    end
end
if handles.reapply_image_functions
    choice = questdlg('Image function parameters have changed, but have not been reapplied. Are you sure you want to continue?', ...
        'Repply image functions', ...
        'Yes','No','No');
    % Handle response
    switch choice
        case 'Yes'
            
        case 'No'
            return
    end
end


if handles.recalc_pairwise_rotations
    [handles.R, handles.W] = compute_pairwise_alignments(handles.images, handles.ang_dis);
    handles.recalc_pairwise_rotations = false;
end

try
    h = multi_waitbar(0, 'Registering and ordering images....') ;
    
    % compute optimal rotations + embedding coordinates using vector diffusion maps
    ncomps = 1;
    [R_opt, embed_coord, D2] = vdm(handles.R, handles.W, handles.eps_scale, ncomps);
    
    multi_waitbar(0.3, h);
    
    % register images using optimal rotations
    if handles.use_raw
        images_registered = register_all_images(handles.images_raw, R_opt);
    else
        images_registered = register_all_images(handles.images, R_opt);
    end
    
    multi_waitbar(0.6, h);
    handles.images_analyzed = order_all_images(images_registered, embed_coord);
    
    handles.reregister_image_set = false;
    handles.rerun_dmaps = false;
    
    multi_waitbar(1, h);
    multi_waitbar(Inf, h);
    
    guidata(hObject, handles)
    
catch ME
    if strcmp(ME.identifier, 'MATLAB:nonExistentField')
        msgbox('Some required fields are not populated. Please check inputs.')
        return
        %     elseif strcmp(ME.identifier, 'MATLAB:imagesci:imread:fileDoesNotExist')
    else
        msgbox(sprintf('ERROR: %s', ME.identifier))
        return
    end
end


function number_rotations_Callback(hObject, eventdata, handles)
% hObject    handle to number_rotations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of number_rotations as text
%        str2double(get(hObject,'String')) returns contents of number_rotations as a double
handles.ang_dis = str2double(get(hObject,'String'));
handles.recalc_pairwise_rotations = true;
if isnan(handles.ang_dis) || handles.ang_dis < 0
    msgbox('Invalid angular discretization.')
    handles = rmfield(handles, 'ang_dis');
    set(hObject, 'String', '');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function number_rotations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_rotations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String','10');
handles.ang_dis = 10;
handles.recalc_pairwise_rotations = true;
guidata(hObject, handles)


% --- Executes on button press in order.
function order_Callback(hObject, eventdata, handles)
% hObject    handle to order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'images')
    msgbox('Images have not been preprocessed.')
    return
end

if handles.reread_images
    choice = questdlg('Initial image settings have changed, but images have not been reread. Are you sure you want to continue?', ...
        'Reread images', ...
        'Yes','No','No');
    % Handle response
    switch choice
        case 'Yes'
            
        case 'No'
            return
    end
end
if handles.reapply_image_functions
    choice = questdlg('Image function parameters have changed, but have not been reapplied. Are you sure you want to continue?', ...
        'Repply image functions', ...
        'Yes','No','No');
    % Handle response
    switch choice
        case 'Yes'
            
        case 'No'
            return
    end
end

try
    h = multi_waitbar(0, 'Ordering images....') ;
    
    W = squareform(pdist(reshape(double(handles.images), [], handles.nimages)')).^2;
    ncomps = 1;
    [embed_coord, D2] = dm(W, handles.eps_scale, ncomps);
    
    multi_waitbar(0.5, h);
    
    if handles.use_raw
        handles.images_analyzed = order_all_images(handles.images_raw, embed_coord);
    else
        handles.images_analyzed = order_all_images(handles.images, embed_coord);
    end
    
    handles.rerun_dmaps = false;
    handles.reregister_image_set = false;
    
    multi_waitbar(1, h);
    multi_waitbar(Inf, h);
    
    guidata(hObject, handles)
    
catch ME
    if strcmp(ME.identifier, 'MATLAB:nonExistentField')
        msgbox('Some required fields are not populated. Please check inputs.')
        return
        %     elseif strcmp(ME.identifier, 'MATLAB:imagesci:imread:fileDoesNotExist')
    else
        msgbox(sprintf('ERROR: %s', ME.identifier))
        return
    end
end

% --- Executes on button press in register.
function register_Callback(hObject, eventdata, handles)
% hObject    handle to register (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'images')
    msgbox('Images have not been preprocessed.')
    return
end

if handles.reread_images
    choice = questdlg('Initial image settings have changed, but images have not been reread. Are you sure you want to continue?', ...
        'Reread images', ...
        'Yes','No','No');
    % Handle response
    switch choice
        case 'Yes'
            
        case 'No'
            return
    end
end
if handles.reapply_image_functions
    choice = questdlg('Image function parameters have changed, but have not been reapplied. Are you sure you want to continue?', ...
        'Repply image functions', ...
        'Yes','No','No');
    % Handle response
    switch choice
        case 'Yes'
            
        case 'No'
            return
    end
end

try
    if handles.recalc_pairwise_rotations
        [handles.R, handles.W] = compute_pairwise_alignments(handles.images, handles.ang_dis);
        handles.recalc_pairwise_rotations = false;
    end
    
    h = multi_waitbar(0, 'Registering images....') ;
    % compute optimal rotations + embedding coordinates using vector diffusion maps
    ncomps = 1;
    [R_opt, embed_coord, D2] = vdm(handles.R, handles.W, handles.eps_scale, ncomps);
    
    multi_waitbar(0.5, h);
    
    % register images using optimal rotations
    if handles.use_raw
        handles.images_analyzed = register_all_images(handles.images_raw, R_opt);
    else
        handles.images_analyzed = register_all_images(handles.images, R_opt);
    end
    
    handles.rerun_dmaps = false;
    handles.reregister_image_set = false;
    
    multi_waitbar(1, h);
    multi_waitbar(Inf, h);
    
    guidata(hObject, handles)
    
catch ME
    if strcmp(ME.identifier, 'MATLAB:nonExistentField')
        msgbox('Some required fields are not populated. Please check inputs.')
        return
        %     elseif strcmp(ME.identifier, 'MATLAB:imagesci:imread:fileDoesNotExist')
    else
        msgbox(sprintf('ERROR: %s', ME.identifier))
        return
    end
end


% --- Executes on button press in register_order_zstacks.
function register_order_zstacks_Callback(hObject, eventdata, handles)
% hObject    handle to register_order_zstacks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'images')
    msgbox('Images have not been preprocessed.')
    return
end

if handles.reread_images
    choice = questdlg('Initial image settings have changed, but images have not been reread. Are you sure you want to continue?', ...
        'Reread images', ...
        'Yes','No','No');
    % Handle response
    switch choice
        case 'Yes'
            
        case 'No'
            return
    end
end
if handles.reapply_image_functions
    choice = questdlg('Image function parameters have changed, but have not been reapplied. Are you sure you want to continue?', ...
        'Repply image functions', ...
        'Yes','No','No');
    % Handle response
    switch choice
        case 'Yes'
            
        case 'No'
            return
    end
end

try
    max_proj_images = squeeze(max(handles.images, [], ndims(handles.images)-1));
    
    if handles.recalc_pairwise_rotations
        [handles.R, handles.W] = compute_pairwise_alignments(max_proj_images, handles.ang_dis);
        handles.recalc_pairwise_rotations = false;
    end
    
    % compute optimal rotations + embedding coordinates using vector diffusion maps
    ncomps = 1;
    [R_opt, embed_coord, D2] = vdm(handles.R, handles.W, handles.eps_scale, ncomps);
    
    if handles.use_raw
        images_registered = register_all_images(handles.images_raw, R_opt);
    else
        images_registered = register_all_images(handles.images_raw, R_opt);
    end
    
    h = multi_waitbar(0, 'Registering and ordering images....') ;
    
    W = squareform(pdist(reshape(double(images_registered), [], handles.nimages)')).^2;
    [embed_coord, D2] = dm(W, handles.eps_scale, ncomps);
    
    multi_waitbar(0.5, h);
    
    handles.images_analyzed = order_all_images(images_registered, embed_coord);
    
    handles.reregister_image_set = false;
    handles.rerun_dmaps = false;
    
    multi_waitbar(1, h);
    multi_waitbar(Inf, h);
    
    guidata(hObject, handles)
    
catch ME
    if strcmp(ME.identifier, 'MATLAB:nonExistentField')
        msgbox('Some required fields are not populated. Please check inputs.')
        return
        %     elseif strcmp(ME.identifier, 'MATLAB:imagesci:imread:fileDoesNotExist')
    else
        msgbox(sprintf('ERROR: %s', ME.identifier))
        return
    end
end


% --- Executes on slider movement.
function red_weight_Callback(hObject, eventdata, handles)
% hObject    handle to red_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.reapply_image_functions = true;

handles.channel_weight(1) = get(hObject,'Value');
set(handles.weight_numbers(1), 'String', sprintf('%2.2f',get(hObject,'Value')));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function red_weight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject, 'Value', 1)
handles.channel_weight(1) = 1;
handles.weight_slider(1) = hObject;
guidata(hObject, handles);


% --- Executes on slider movement.
function green_weight_Callback(hObject, eventdata, handles)
% hObject    handle to green_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.reapply_image_functions = true;


handles.channel_weight(2) = get(hObject,'Value');
set(handles.weight_numbers(2), 'String', sprintf('%2.2f',get(hObject,'Value')));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function green_weight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject, 'Value', 1)
handles.channel_weight(2) = 1;
handles.weight_slider(2) = hObject;
guidata(hObject, handles);

% --- Executes on slider movement.
function blue_weight_Callback(hObject, eventdata, handles)
% hObject    handle to blue_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.reapply_image_functions = true;

handles.channel_weight(3) = get(hObject,'Value');
set(handles.weight_numbers(3), 'String', sprintf('%2.2f',get(hObject,'Value')));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function blue_weight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue_weight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject, 'Value', 1)
handles.channel_weight(3) = 1;

handles.weight_slider(3) = hObject;
guidata(hObject, handles);

% --- Executes on slider movement.
function red_blur_Callback(hObject, eventdata, handles)
% hObject    handle to red_blur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.reapply_image_functions = true;

handles.channel_blur(1) = get(hObject,'Value');
set(handles.blur_numbers(1), 'String', sprintf('%2.0f%%',100*get(hObject,'Value')));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function red_blur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red_blur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject, 'Value', 0)
handles.channel_blur(1) = 0;

handles.blur_slider(1) = hObject;
guidata(hObject, handles);

% --- Executes on slider movement.
function green_blur_Callback(hObject, eventdata, handles)
% hObject    handle to green_blur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.reapply_image_functions = true;

handles.channel_blur(2) = get(hObject,'Value');
set(handles.blur_numbers(2), 'String', sprintf('%2.0f%%',100*get(hObject,'Value')));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function green_blur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green_blur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject, 'Value', 0)
handles.channel_blur(2) = 0;

handles.blur_slider(2) = hObject;
guidata(hObject, handles);

% --- Executes on slider movement.
function blue_blur_Callback(hObject, eventdata, handles)
% hObject    handle to blue_blur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.reapply_image_functions = true;

handles.channel_blur(3) = get(hObject,'Value');
set(handles.blur_numbers(3), 'String', sprintf('%2.0f%%',100*get(hObject,'Value')));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function blue_blur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue_blur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject, 'Value', 0)
handles.channel_blur(3) = 0;

handles.blur_slider(3) = hObject;
guidata(hObject, handles);

% --- Executes on button press in red_normalize.
function red_normalize_Callback(hObject, eventdata, handles)
% hObject    handle to red_normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of red_normalize

handles.reapply_image_functions = true;

handles.channel_normalize(1) = get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in green_normalize.
function green_normalize_Callback(hObject, eventdata, handles)
% hObject    handle to green_normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of green_normalize

handles.reapply_image_functions = true;

handles.channel_normalize(2) = get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in blue_normalize.
function blue_normalize_Callback(hObject, eventdata, handles)
% hObject    handle to blue_normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of blue_normalize

handles.reapply_image_functions = true;

handles.channel_normalize(3) = get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in red_mean_center.
function red_mean_center_Callback(hObject, eventdata, handles)
% hObject    handle to red_mean_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of red_mean_center

handles.reapply_image_functions = true;

handles.channel_mean_center(1) = get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in green_mean_center.
function green_mean_center_Callback(hObject, eventdata, handles)
% hObject    handle to green_mean_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of green_mean_center

handles.reapply_image_functions = true;

handles.channel_mean_center(2) = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in blue_mean_center.
function blue_mean_center_Callback(hObject, eventdata, handles)
% hObject    handle to blue_mean_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of blue_mean_center

handles.reapply_image_functions = true;

handles.channel_mean_center(3) = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function red_weight_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red_weight_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.weight_numbers(1) = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function green_weight_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green_weight_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.weight_numbers(2) = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function blue_weight_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue_weight_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.weight_numbers(3) = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function red_blur_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red_blur_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.blur_numbers(1) = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function green_blur_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green_blur_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.blur_numbers(2) = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function blue_blur_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue_blur_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.blur_numbers(3) = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function red_normalize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red_normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.normalize_checkbox(1) = hObject;
set(hObject, 'Value',0)
handles.channel_normalize(1) = false;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function green_normalize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green_normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.normalize_checkbox(2) = hObject;
set(hObject, 'Value',0)
handles.channel_normalize(2) = false;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function blue_normalize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue_normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.normalize_checkbox(3) = hObject;
set(hObject, 'Value',0)
handles.channel_normalize(3) = false;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function red_mean_center_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red_mean_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.mean_center_checkbox(1) = hObject;
set(hObject, 'Value',1)
handles.channel_mean_center(1) = true;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function green_mean_center_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green_mean_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.mean_center_checkbox(2) = hObject;
set(hObject, 'Value',1)
handles.channel_mean_center(2) = true;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function blue_mean_center_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue_mean_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.mean_center_checkbox(3) = hObject;
set(hObject, 'Value',1)
handles.channel_mean_center(3) = true;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function green_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.color_label(2) = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function blue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


handles.color_label(3) = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function register_CreateFcn(hObject, eventdata, handles)
% hObject    handle to register (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.register_button = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - hand

handles.order_button = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function register_order_vdm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to register_order_vdm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.register_order_button = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function register_order_zstacks_CreateFcn(hObject, eventdata, handles)
% hObject    handle to register_order_zstacks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.register_order_zstack_button = hObject;
set(hObject,'enable','off');
guidata(hObject, handles);


% --- Executes on button press in show_images3.
function show_images3_Callback(hObject, eventdata, handles)
% hObject    handle to show_images3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'images_analyzed')
    msgbox('Images have not been analyzed.')
    return
end

plot_images(handles.images_analyzed, handles.dim)


% --- Executes on button press in resize_image.
function resize_image_Callback(hObject, eventdata, handles)
% hObject    handle to resize_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of resize_image

handles.reapply_image_functions = true;

handles.resize_image = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function resize_image_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resize_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject, 'Value', 0)
handles.resize_image = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function red_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.color_label(1) = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function read_images_CreateFcn(hObject, eventdata, handles)
% hObject    handle to read_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.reread_images = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function apply_image_functions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to apply_image_functions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.reapply_image_functions = true;
guidata(hObject, handles);


% --- Executes when selected object is changed in select_raw_or_preprocessed.
function select_raw_or_preprocessed_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in select_raw_or_preprocessed
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

handles.reregister_image_set = true;

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'raw_button'
        handles.use_raw = true;
    case 'preprocessed_button'
        handles.use_raw = false;
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function raw_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to raw_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'Value',get(hObject,'Min'));
handles.use_raw = false;
handles.reregister_image_set = true;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function preprocessed_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to preprocessed_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'Value',get(hObject,'Max'));
handles.use_raw = false;
handles.reregister_image_set = true;
guidata(hObject, handles);



function nimages_avg_Callback(hObject, eventdata, handles)
% hObject    handle to nimages_avg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nimages_avg as text
%        str2double(get(hObject,'String')) returns contents of nimages_avg as a double

nimages = str2double(get(hObject,'String'));
handles.nimages_avg = nimages;
if isnan(handles.nimages_avg) || handles.nimages_avg < 1 || mod(handles.nimages_avg, 1) ~= 0
    msgbox('Invalid number of images.')
    handles = rmfield(handles, 'nimages_avg');
    set(hObject, 'String', '');
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function nimages_avg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nimages_avg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function avg_width_Callback(hObject, eventdata, handles)
% hObject    handle to avg_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of avg_width as text
%        str2double(get(hObject,'String')) returns contents of avg_width as a double

handles.avg_width = str2double(get(hObject,'String')) ;

if isnan(handles.avg_width)
    msgbox('Invalid width of averaging window.')
    handles = rmfield(handles, 'avg_width');
    set(hObject, 'String', '');
elseif handles.avg_width < 0
    msgbox('Invalid width of averaging window: width must be greater than 0.')
    handles = rmfield(handles, 'avg_width');
    set(hObject, 'String', '');
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function avg_width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to avg_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in show_avg_traj.
function show_avg_traj_Callback(hObject, eventdata, handles)
% hObject    handle to show_avg_traj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'images_analyzed')
    msgbox('Images have not been analyzed.')
    return
end

try
    avg_images = compute_average_trajectory(handles.images_analyzed, handles.nimages_avg, handles.avg_width);
    plot_images(avg_images, handles.dim)
catch ME
    if strcmp(ME.identifier, 'MATLAB:nonExistentField')
        msgbox('Some required fields are not populated. Please check inputs.')
        return
        %     elseif strcmp(ME.identifier, 'MATLAB:imagesci:imread:fileDoesNotExist')
    else
        msgbox(sprintf('ERROR: %s', ME.identifier))
        return
    end
end

% --- Executes on button press in save_avg_traj.
function save_avg_traj_Callback(hObject, eventdata, handles)
% hObject    handle to save_avg_traj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'images_analyzed')
    msgbox('Images have not been analyzed.')
    return
end

try
    avg_images = compute_average_trajectory(handles.images_analyzed, handles.nimages_avg, handles.avg_width);
    
    save_images(avg_images, handles.dim, handles.output_image_dir, strcat('avg_',handles.image_name), handles.image_ext, handles.stack_name);
catch ME
    if strcmp(ME.identifier, 'MATLAB:nonExistentField')
        msgbox('Some required fields are not populated. Please check inputs.')
        return
        %     elseif strcmp(ME.identifier, 'MATLAB:imagesci:imread:fileDoesNotExist')
    else
        msgbox(sprintf('ERROR: %s', ME.identifier))
        return
    end
end
