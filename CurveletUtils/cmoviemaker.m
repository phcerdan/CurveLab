function varargout = cmoviemaker(varargin)
% CMOVIEMAKER M-file for cmoviemaker.fig
%      CMOVIEMAKER, by itself, creates a new CMOVIEMAKER or raises the existing
%      singleton*.
%
%      H = CMOVIEMAKER returns the handle to a new CMOVIEMAKER or the handle to
%      the existing singleton*.
%
%      CMOVIEMAKER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CMOVIEMAKER.M with the given input arguments.
%
%      CMOVIEMAKER('Property','Value',...) creates a new CMOVIEMAKER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cmoviemaker_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cmoviemaker_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cmoviemaker

% Last Modified by GUIDE v2.5 28-Nov-2007 11:42:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cmoviemaker_OpeningFcn, ...
                   'gui_OutputFcn',  @cmoviemaker_OutputFcn, ...
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


% --- Executes just before cmoviemaker is made visible.
function cmoviemaker_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cmoviemaker (see VARARGIN)

% Choose default command line output for cmoviemaker
handles.output = hObject;

% set default values for handle structure
handles.filenames = {};
handles.pathname = '';
handles.lastsavefilename = '';
handles.lastsavepath = '';
handles.hfig = [];
handles.computeinfo = struct('xlim',[],'ylim',[],'invert',0,'thrshtype',1,'thrshval',0,'levels',[],'dirs',{{}});
handles.layout = [1 1];
handles.plotinfo = struct('type', {1}, 'pos', {1}, 'transform', {1}, 'value', {1}, 'colorbar', {1});

% Update handles structure
guidata(hObject, handles);

setpath;

% UIWAIT makes cmoviemaker wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cmoviemaker_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_loadimages.
function button_loadimages_Callback(hObject, eventdata, handles)
% hObject    handle to button_loadimages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile('*.*','Select images',handles.pathname,'MultiSelect','on');
if ~(isequal(filename,0) || isequal(pathname,0)),
    if ischar(filename),
        handles.filenames = {filename};
    else
        handles.filenames = filename;
    end
    handles.pathname = pathname;
    
    % init filename listbox
    set(handles.listbox_images, 'String', handles.filenames);
    
    % show first image
    im1 = displayfile(1, handles);

    % initialize levels and dirs listboxes
    handles = initlistboxes(size(im1), handles);
    
    % save handles
    guidata(gcbo, handles);
end



% --- Executes on selection change in listbox_images.
function listbox_images_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nr = get(hObject, 'Value');
if nr > 0,
    im = displayfile(nr, handles);
    nlev = floor(log2(min(size(im)))) - 3;
    handles = reinitlistboxes(nlev, handles);
    guidata(hObject, handles);
end



% --- Executes during object creation, after setting all properties.
function listbox_images_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_xlim_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xlim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = get(hObject, 'String');
xlim = eval(str);
if (length(xlim) == 2 && xlim(2) > xlim(1) && xlim(1) >= 1) ...
        || isempty(xlim),
    handles.computeinfo.xlim = xlim;
    im = displaycurrentfile(handles);
    nlev = floor(log2(min(size(im)))) - 3;
    handles = reinitlistboxes(nlev, handles);
    guidata(hObject, handles)
else
    errordlg('Bad values for x-limits','CMovieMaker error')
end


% --- Executes during object creation, after setting all properties.
function edit_xlim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_xlim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ylim_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ylim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = get(hObject, 'String');
ylim = eval(str);
if (length(ylim) == 2 && ylim(2) > ylim(1) && ylim(1) >= 1) ...
        || isempty(ylim),
    handles.computeinfo.ylim = ylim;
    im = displaycurrentfile(handles);
    nlev = floor(log2(min(size(im)))) - 3;
    handles = reinitlistboxes(nlev, handles);
    guidata(hObject, handles)
else
    errordlg('Bad values for y-limits','CMovieMaker error')
end


% --- Executes during object creation, after setting all properties.
function edit_ylim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ylim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_levels.
function listbox_levels_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_levels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.computeinfo.levels = get(handles.listbox_levels, 'Value');
nlev = length(get(handles.listbox_levels, 'String'));
handles = reinitlistboxes(nlev, handles);
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function listbox_levels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_levels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_dirs.
function listbox_dirs_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_dirs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nrdirs = length(get(hObject,'String'));
if get(handles.checkbox_symmdirs,'Value') == get(handles.checkbox_symmdirs,'Max') && ...
         nrdirs > 1,
    topitem = get(hObject, 'ListBoxTop');
    dirsel = get(hObject,'Value');
    dirsel = unique([dirsel mod(dirsel + nrdirs/2 - 1, nrdirs)+1]);
    set(hObject,'Value',dirsel)    
    set(hObject, 'ListBoxTop', topitem)
end

% update computeinfo
for k=handles.computeinfo.levels,
    if ismember(k, [1 length(get(handles.listbox_levels, 'String'))]),
        thisnrdirs = 1;
    else
        thisnrdirs = 16 * (2.^floor((k-1)/2));
    end
    step = round(thisnrdirs / nrdirs);
    dirs = 1:thisnrdirs;
    dirs = dirs(ismember(floor((dirs-1)/step) + 1, dirsel));
    handles.computeinfo.dirs{k} = dirs;
end
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function listbox_dirs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_dirs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_symmdirs.
function checkbox_symmdirs_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_symmdirs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_symmdirs


% --- Executes on button press in button_makemovie.
function button_makemovie_Callback(hObject, eventdata, handles)
% hObject    handle to button_makemovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

curdir = pwd;
if ~isempty(handles.lastsavepath)
    cd(handles.lastsavepath);
end

[filename, pathname, index] = uiputfile({'*.avi', 'AVI files (*.avi)'}, 'Save movie as...', handles.lastsavefilename);
if index ==0,
   return; % user pressed cancel
end

handles.lastsavefilename = filename;
handles.lastsavepath = pathname;
cd(curdir);

if isempty(handles.hfig) || ~ishandle(handles.hfig),
    handles.hfig = figure;
    guidata(hObject, handles);
end

if get(handles.popup_thresh, 'Value') == 4,   % get threshold from first image
   handles = getthrshfromfirstimage(handles);
end

clear M;
for k=1:length(handles.filenames),
    hf = drawmovieframe(handles.hfig, fullfile(handles.pathname, handles.filenames{k}), ...
            handles.layout, handles.plotinfo, handles.computeinfo);
    pause(0.1);
    M(k) = getframe(hf);
end
movie2avi(M, fullfile(pathname, filename));




function edit_thrshvalue_Callback(hObject, eventdata, handles)
% hObject    handle to edit_thrshvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = str2double(get(hObject, 'String'));
if isnan(val) || val <= 0,
    errordlg('Bad value for threshold/number of coefficients', 'CMovieMaker error');
    return;
end
handles.computeinfo.thrshval = val;
guidata(hObject, handles)




% --- Executes during object creation, after setting all properties.
function edit_thrshvalue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_thrshvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_thresh.
function popup_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to popup_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.computeinfo.thrshtype = get(hObject, 'Value');
guidata(hObject, handles)
edit_thrshvalue_Callback(handles.edit_thrshvalue, eventdata, handles);



% --- Executes during object creation, after setting all properties.
function popup_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_thrshdisp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_thrshdisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes during object creation, after setting all properties.
function edit_thrshdisp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_thrshdisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% ------------------------------
% custom helper functions
%

% load, display and return an image from the list
function im = displayfile(nr, handles)

im = imread(fullfile(handles.pathname, handles.filenames{nr}));
if get(handles.checkbox_invert, 'Value'),
    im = max(im(:)) - im;
end
if isempty(handles.computeinfo.xlim),
    xvals = 1:size(im,2);
else
    xvals = floor(handles.computeinfo.xlim(1):min(handles.computeinfo.xlim(2), size(im,2)));
end
if isempty(handles.computeinfo.ylim),
    yvals = 1:size(im,1);
else
    yvals = floor(handles.computeinfo.ylim(1):min(handles.computeinfo.ylim(2),size(im,1)));
end
im = im(yvals, xvals);

axes(handles.axes_top);
imagesc([xvals(1) xvals(end)], [yvals(1) yvals(end)], im);
axis equal tight; colormap gray; colorbar


% displays the currently selected file in listbox_images
function im = displaycurrentfile(handles)

nr = get(handles.listbox_images, 'Value');
im = displayfile(nr, handles);




% Get selected values from a listbox containing values as entries
function vals = getlistboxvals(lbhandle)
% lbhandle      handle to a listbox
sel = get(lbhandle,'Value');
strs = get(lbhandle,'String');
vals = str2double({strs{sel}});



% initialize level and dir listboxes
function handles = initlistboxes(dims, handles)
% dims      dimensions of the current image
% handles   figure handles

im0 = zeros(dims);
C0 = fdct_wrapping(im0,0);
str = cell(1,length(C0));
for ll=1:length(C0)
    str{ll} = sprintf('%d',ll);
end
set(handles.listbox_levels,'String',str)
set(handles.listbox_levels,'Value',1:length(C0))

str=cell(1,length(C0{2}));
for dd=1:length(C0{2})
    str{dd} = sprintf('%d',dd);
end
set(handles.listbox_dirs,'String',str)
set(handles.listbox_dirs,'Value',1:length(C0{2}))

handles.computeinfo.levels = 1:length(C0);
for k=1:length(C0),
    handles.computeinfo.dirs{k} = 1:length(C0{k});
end
guidata(handles.listbox_levels, handles);



% call to reinitialize listboxes for new image or new level selection
function handles = reinitlistboxes(nlev, handles)

oldlevs = getlistboxvals(handles.listbox_levels);
olddirs = getlistboxvals(handles.listbox_dirs);
oldnrdirs = length(get(handles.listbox_dirs,'String'));
firstclev = min(oldlevs(oldlevs > 1));

% update levels
if nlev ~= length(get(handles.listbox_levels,'String')),
    str = cell(1, nlev);
    for ll=1:nlev,
        str{ll} = sprintf('%d',ll);
    end
    set(handles.listbox_levels,'String',str);
    if firstclev > nlev,
        % there were more levels in previous image, and only finer levels
        % selected => select finest level
        set(handles.listbox_levels,'Value',nlev);
    else
        set(handles.listbox_levels,'Value',oldlevs(oldlevs <= nlev));
    end
end

% update dirs
set(handles.listbox_dirs,'ListboxTop',1);
if isempty(firstclev) || firstclev >= nlev,
    % only wavelet levels selected, only one dir
    str = {'1'};
    nrdirs = 1;
    dirsel = 1;
else
    nrdirs = 16 * 2.^(floor((firstclev-1)/2));
    step = nrdirs/oldnrdirs;
    str = cell(1,nrdirs);
    dirsel = [];
    for dd=1:nrdirs,
        str{dd} = sprintf('%d',dd);
        if ismember(floor((dd-1)/step) + 1, olddirs),
            dirsel = [dirsel dd];
        end
    end
end
set(handles.listbox_dirs,'String',str);
set(handles.listbox_dirs,'Value',dirsel);
 




% resizes the plotinfo structure and shows correct data
function resizeplotinfo(nrows, ncols, handles)

if ~isempty(nrows) && ~isempty(ncols),
    % update listbox_subplotidx
    str = cell(nrows, ncols);
    for kr = 1:nrows,
        for kc = 1:ncols,
            str{kr, kc} = sprintf('[%d,%d]', kr, kc);
        end
    end
    set(handles.listbox_subplotidx, 'String', {str{:}});
    
    % resize plotinfo structure (saving old data), and layout member
    nel = nrows*ncols;
    sz = min(length(handles.plotinfo), nel);
    [ptypes{1:nel}] = deal(1);
    [ptypes{1:sz}] = deal(handles.plotinfo(1:sz).type);
    [ppos{1:nel}] = deal(1);
    [ppos{1:sz}] = deal(handles.plotinfo(1:sz).pos);
    [ptrans{1:nel}] = deal(1);
    [ptrans{1:sz}] = deal(handles.plotinfo(1:sz).transform);
    [pvals{1:nel}] = deal(1);
    [pvals{1:sz}] = deal(handles.plotinfo(1:sz).value);
    [pcbar{1:nel}] = deal(1);
    [pcbar{1:sz}] = deal(handles.plotinfo(1:sz).colorbar);
    handles.layout = [nrows ncols];
    handles.plotinfo = struct('type',ptypes,'pos',ppos,'transform',ptrans,'value',pvals,'colorbar',pcbar);
    guidata(gcbo, handles)
end


% sets the threshold using a specified number of coefficients from the
% first image
function handles = getthrshfromfirstimage(handles)

im = imread(fullfile(handles.pathname, handles.filenames{1}));
if get(handles.checkbox_invert, 'Value'),
    im = max(im(:)) - im;
end
if isempty(handles.computeinfo.xlim),
    xvals = 1:size(im,2);
else
    xvals = floor(handles.computeinfo.xlim(1):min(handles.computeinfo.xlim(2), size(im,2)));
end
if isempty(handles.computeinfo.ylim),
    yvals = 1:size(im,1);
else
    yvals = floor(handles.computeinfo.ylim(1):min(handles.computeinfo.ylim(2),size(im,1)));
end
im = im(yvals, xvals);

ncoeff = str2double(get(handles.edit_thrshvalue, 'String'));
if isnan(ncoeff) || isempty(ncoeff) || ncoeff <= 0,
    errordlg('Bad number of coefficients','Error in CMovieMaker!')
    return
end
wtypes = {'db3','db5','dmey'};
for k=unique(cat(2,handles.plotinfo.transform)),
    if k==1,
        C = fdct_wrapping(im, 0);
        [Cmod, thval{k}] = crvlt_keeplargest(C, round(ncoeff));
    else
        nlev = max([1 (wmaxlev(size(im), wtypes{k-1}) - 1)]);   % estimate number of useful levels
        [D, S] = wavedec2(im, nlev, wtypes{k-1});
        Ds = sort(abs(D),2,'descend');
        thval{k} = Ds(round(ncoeff)+1);
    end
end
handles.computeinfo.thrshval = thval;
handles.computeinfo.thrshtype = 2;   % set type to thresholding
guidata(gcbo, handles);


% --- Executes on button press in checkbox_invert.
function checkbox_invert_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_invert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.computeinfo.invert = get(hObject, 'Value');
guidata(hObject, handles);
displaycurrentfile(handles);


function edit_nrows_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nrows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nrows = str2num(get(handles.edit_nrows, 'String'));
ncols = str2num(get(handles.edit_ncols, 'String'));
resizeplotinfo(nrows, ncols, handles);


% --- Executes during object creation, after setting all properties.
function edit_nrows_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nrows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ncols_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ncols (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nrows = str2num(get(handles.edit_nrows, 'String'));
ncols = str2num(get(handles.edit_ncols, 'String'));
resizeplotinfo(nrows, ncols, handles);



% --- Executes during object creation, after setting all properties.
function edit_ncols_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ncols (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_subplotidx.
function listbox_subplotidx_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_subplotidx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

curplot = get(hObject, 'Value');
set(handles.popup_plottype, 'Value', handles.plotinfo(curplot).type);
set(handles.popup_posplot, 'Value', handles.plotinfo(curplot).pos);
set(handles.popup_transform, 'Value', handles.plotinfo(curplot).transform);
set(handles.popup_values, 'Value', handles.plotinfo(curplot).value);
set(handles.checkbox_colorbar, 'Value', handles.plotinfo(curplot).colorbar);






% --- Executes during object creation, after setting all properties.
function listbox_subplotidx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_subplotidx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popup_posplot.
function popup_posplot_Callback(hObject, eventdata, handles)
% hObject    handle to popup_posplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

curplot = get(handles.listbox_subplotidx, 'Value');
handles.plotinfo(curplot).pos = get(hObject, 'Value');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popup_posplot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_posplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_preview.
function button_preview_Callback(hObject, eventdata, handles)
% hObject    handle to button_preview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.hfig) || ~ishandle(handles.hfig),
    handles.hfig = figure;
    guidata(hObject, handles);
end

if get(handles.popup_thresh, 'Value') == 4,   % get threshold from first image
   handles = getthrshfromfirstimage(handles);
end

nr = get(handles.listbox_images, 'Value');
fullfilename = fullfile(handles.pathname, handles.filenames{nr});
drawmovieframe(handles.hfig, fullfilename, handles.layout, handles.plotinfo, handles.computeinfo);


% --- Executes on selection change in popup_transform.
function popup_transform_Callback(hObject, eventdata, handles)
% hObject    handle to popup_transform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

curplot = get(handles.listbox_subplotidx, 'Value');
handles.plotinfo(curplot).transform = get(hObject, 'Value');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popup_transform_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_transform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_plottype.
function popup_plottype_Callback(hObject, eventdata, handles)
% hObject    handle to popup_plottype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

curplot = get(handles.listbox_subplotidx, 'Value');
handles.plotinfo(curplot).type = get(hObject, 'Value');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popup_plottype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_plottype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_values.
function popup_values_Callback(hObject, eventdata, handles)
% hObject    handle to popup_values (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

curplot = get(handles.listbox_subplotidx, 'Value');
handles.plotinfo(curplot).value = get(hObject, 'Value');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popup_values_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_values (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_colorbar.
function checkbox_colorbar_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_colorbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

curplot = get(handles.listbox_subplotidx, 'Value');
handles.plotinfo(curplot).colorbar = get(hObject, 'Value');
guidata(hObject, handles);

