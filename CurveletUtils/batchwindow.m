function varargout = batchwindow(varargin)
% BATCHWINDOW M-file for batchwindow.fig
%      BATCHWINDOW, by itself, creates a new BATCHWINDOW or raises the existing
%      singleton*.
%
%      H = BATCHWINDOW returns the handle to a new BATCHWINDOW or the handle to
%      the existing singleton*.
%
%      BATCHWINDOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BATCHWINDOW.M with the given input arguments.
%
%      BATCHWINDOW('Property','Value',...) creates a new BATCHWINDOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before batchwindow_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to batchwindow_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help batchwindow

% Last Modified by GUIDE v2.5 01-Jul-2008 20:32:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @batchwindow_OpeningFcn, ...
                   'gui_OutputFcn',  @batchwindow_OutputFcn, ...
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


% --- Executes just before batchwindow is made visible.
function batchwindow_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to batchwindow (see VARARGIN)

% Choose default command line output for batchwindow
handles.output = hObject;

handles.thresh = varargin{1};
handles.thrshtype = varargin{2};
handles.levels = varargin{3};
handles.dirs = varargin{4};
handles.totdirs = varargin{5};

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes batchwindow wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = batchwindow_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if ishandle(hObject),
    varargout{1} = handles.output;
end



function output = analyzeimage(im, handles)

output.im = im;
C0 = fdct_wrapping(im, 0);

% remove unselected levels and dirs
for ll=1:length(C0)
    levok = ismember(ll,handles.levels);
    ndir = length(C0{ll});
    if ndir ==1,
        if ~levok,
            C0{ll}{1} = zeros(size(C0{ll}{1}));
        end
    else
        step = round(ndir/handles.totdirs);
        for dd=1:ndir,
            dirok = ismember(floor((dd-1)/step) + 1, handles.dirs);
            if ~(levok && dirok)
                C0{ll}{dd} = zeros(size(C0{ll}{dd}));
            end
        end
    end
end

% threshold
[Cmod, output.nnzcf, output.smallestcf] = threshold_coeffs(C0, handles);

output.rim = ifdct_wrapping(Cmod, 0);
output.errim = output.rim - output.im;

output.levels = handles.levels;
output.dirs = handles.dirs;

if get(handles.check_writecurvelets, 'Value')
    output.Cmod = Cmod;
end

if get(handles.check_writefft, 'Value'),
    output.fim = fftshift(fft2(im));
    output.frim = fftshift(fft2(output.rim));
    output.ferrim = output.frim - output.fim;
end


% Threshold coefficients according to settings
function [Cth, nnzcf, smcoef] = threshold_coeffs(C, handles)
% C         curvelet coeffs to be thresholded
% handles   figure handles

val = handles.thresh;
if isnan(val) || (val < 0),
    errordlg('Bad threshold value','Curvelet viewer error!');
    return
end
switch handles.thrshtype,
    case 1, % do nothing
        Cth = C;
        smcoef = 0.0;
        nnzcf = crvlt_countnnz(Cth);
    case 2,  % threshold value
        Cth = crvlt_thresh(C, val);
        nnzcf =  crvlt_countnnz(Cth);
        smcoef = val;
    case 3,  % number of coeffs
        if val <= 0,
            Cth = zerodct(C);
            smcoef = 0.0;
            nnzcf = crvlt_countnnz(Cth);
        else
            [Cth, smcoef] = crvlt_keeplargest(C, round(val));
            nnzcf = round(val);
        end
    case 4,  % percentage (assumes val is in [0,1])
        if (val > 1 || val <= 0)
            errordlg('Bad percentage value. Must be between 0 and 1.','Curvelet viewer error!')
        else
            [Cth, smcoef] = crvlt_keeplargest(C, val);
            nnzcf = crvlt_countnnz(Cth);
        end
end



function edit_dir_Callback(hObject, eventdata, handles)

handles.dirname = get(handles.edit_dir, 'String');
guidata(hObject, handles);


function edit_dir_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function button_browse_Callback(hObject, eventdata, handles)


handles.dirname = uigetdir('', 'Choose image directory');
if ~isequal(handles.dirname, 0),
    handles.dirname = [handles.dirname filesep];
    set(handles.edit_dir, 'String', handles.dirname);
    guidata(hObject, handles);
end



function button_go_Callback(hObject, eventdata, handles)

if isempty(handles.dirname),
    errordlg('No directory specified', 'CViewer error');
    return
end
files = dir(handles.dirname);

if ~exist(fullfile(handles.dirname, 'Output'), 'dir'),
    success = mkdir(fullfile(handles.dirname, 'Output'));
    if ~success,
        errordlg('Error creating Output directory', 'CViewer error');
        return
    end
end

cnt = 0;
for k=1:length(files),
    if files(k).isdir,
        continue
    end
    try
        im = double(imread(fullfile(handles.dirname, files(k).name)));
        
    catch
        continue
    end
    
    % choose non-zero color channel
    if (ndims(im) > 2),
        nzv = zeros(1,min(size(im,3),3));
        for m=1:min(size(im,3),3),
            nzv(m) = nnz(im(:,:,m));
        end
        [M,I] = max(nzv);
        im = im(:,:,I);
    end
    set(handles.text_currentfile, 'String', ['Analyzing: ' files(k).name])
    drawnow('update');
    cviewdata = analyzeimage(im, handles);

    [pt, nm, ext] = fileparts(files(k).name);
    matfilename = fullfile(handles.dirname, 'Output', [nm '_analyzed.mat']);
    save(matfilename, 'cviewdata');
    if get(handles.check_writeimages, 'Value'),
        imfilename = fullfile(handles.dirname, 'Output', [nm '_analyzed' ext]);
        tmpim = real(cviewdata.rim) - min(real(cviewdata.rim(:)));
        imwrite(tmpim/max(tmpim(:)) * 256, repmat(linspace(0,1,256)',[1 3]), imfilename);
    end

    cnt = cnt + 1;

end

if cnt == 0,
    errordlg('No images found', 'CViewer error');
else
    msgbox(sprintf('Analyzed %d images.', cnt));
end

uiresume(handles.figure1);
close(handles.figure1);



function check_writeimages_Callback(hObject, eventdata, handles)


function check_writecurvelets_Callback(hObject, eventdata, handles)


function check_writefft_Callback(hObject, eventdata, handles)


