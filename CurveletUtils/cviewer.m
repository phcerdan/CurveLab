function varargout = cviewer(varargin)
% CVIEWER M-file for cviewer.fig
%      CVIEWER, by itself, creates a new CVIEWER or raises the existing
%      singleton*.
%
%      H = CVIEWER returns the handle to a new CVIEWER or the handle to
%      the existing singleton*.
%
%      CVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CVIEWER.M with the given input arguments.
%
%      CVIEWER('Property','Value',...) creates a new CVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cviewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cviewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cviewer

% Last Modified by GUIDE v2.5 01-Jul-2008 20:36:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cviewer_OpeningFcn, ...
                   'gui_OutputFcn',  @cviewer_OutputFcn, ...
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


% --- Executes just before cviewer is made visible.
function cviewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cviewer (see VARARGIN)

% Choose default command line output for cviewer
handles.output = hObject;

% Set paths to curvelet functions
setpath;

% set some default values
handles.imdims = [256 256];
handles.Cmod = {};
handles.loadedim = zeros(handles.imdims);
handles.im = zeros(handles.imdims);
handles.fim = fftshift(fft2(handles.im));
handles.rim = zeros(handles.imdims);
handles.frim = fftshift(fft2(handles.rim));
handles.errim = handles.rim - handles.im;
handles.ferrim = handles.frim - handles.fim;
handles.C = fdct_wrapping(handles.im,0);
handles.stickyaxes = 0;
handles.lastdisptype = [1 1];
handles.curdirdata = struct('levels',[],'dirs',{},'values',{},'dothandle',[],'dotcoords',[]);
handles.lastpath = '';


% initialize the levels and dirs listboxes
handles = initlistboxes(handles.imdims, handles);

% Update handles structure
guidata(hObject, handles);



% UIWAIT makes cviewer wait for user response (see UIRESUME)
% uiwait(handles.cviewerfig);


% --- Outputs from this function are returned to the command line.
function varargout = cviewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox_levels.
function listbox_levels_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_levels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = reinitlistboxes(handles);
guidata(hObject, handles);

% if we view coefficients - plot directly
disptype_top = get(handles.displaymenu_top, 'Value');
disptype_bottom = get(handles.displaymenu_bottom, 'Value');
if any(disptype_top == [7 8]) || any(disptype_bottom == [7 8]),  % view coefficients in top or bottom
    displaydata(handles)
end


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
if get(handles.chk_symmetric,'Value') == get(handles.chk_symmetric,'Max') && ...
         nrdirs > 1,
    topitem = get(hObject, 'ListBoxTop');
    dirsel = get(hObject,'Value');
    dirsel = unique([dirsel mod(dirsel + nrdirs/2 - 1, nrdirs)+1]);
    set(hObject,'Value',dirsel)    
    set(hObject, 'ListBoxTop', topitem)
end

% if we view coefficients - plot directly
disptype_top = get(handles.displaymenu_top, 'Value');
disptype_bottom = get(handles.displaymenu_bottom, 'Value');
if any(disptype_top == [7 8]) || any(disptype_bottom == [7 8]),  % view coefficients in top or bottom
    displaydata(handles)
end


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


% --- Executes on selection change in displaymenu_bottom.
function displaymenu_bottom_Callback(hObject, eventdata, handles)
% hObject    handle to displaymenu_bottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dtype = get(handles.displaymenu_bottom,'Value');
handles.stickyaxes = (isphysspace(dtype) == isphysspace(handles.lastdisptype(2)));

showinaxis(handles.axes_bottom, dtype, 1, ...
    get(handles.plottype_menu_bottom,'Value'), handles.stickyaxes, handles);
handles.lastdisptype(2) = dtype;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function displaymenu_bottom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to displaymenu_bottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in displaymenu_top.
function displaymenu_top_Callback(hObject, eventdata, handles)
% hObject    handle to displaymenu_top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dtype = get(handles.displaymenu_top,'Value');
handles.stickyaxes = (isphysspace(dtype) == isphysspace(handles.lastdisptype(1)));

showinaxis(handles.axes_top, dtype, get(handles.menu_plot_curveletpos_top,'Value'), ...
    get(handles.plottype_menu_top,'Value'), handles.stickyaxes, handles);
handles.lastdisptype(1) = dtype;
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function displaymenu_top_CreateFcn(hObject, eventdata, handles)
% hObject    handle to displaymenu_top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function cmax = getlargestcoeff(C)

cmax = 0;
for j=1:length(C),
    for l = 1:length(C{j})
        cmax = max([cmax; abs(C{j}{l}(:))]);
    end
end


% --- Load image from file 
function loadimage(hObject, handles)
% hObject    handle to button_loadimage (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

curdir = pwd;
if ~isempty(handles.lastpath)
    cd(handles.lastpath)
end
[filename, pathname]=uigetfile({'*' 'All files (*.*)'}, 'Open image');
if ~(isequal(filename,0) || isequal(pathname,0))
    [pstr,nm,ext,ver] = fileparts(fullfile(pathname,filename));
    if strcmp(ext,'.dcm')
        im = dicomread(fullfile(pathname,filename));
    elseif strcmp(ext,'.dat')
        im = load(fullfile(pathname,filename));
    else
        im = imread(fullfile(pathname,filename));
    end
    if (ndims(im) > 2),
        nzv = zeros(1,min(size(im,3),3));
        for k=1:min(size(im,3),3),
            nzv(k) = nnz(im(:,:,k));
        end
        [M,I] = max(nzv);
        im = im(:,:,I);
    end
    handles.loadedim = double(im);
    
    set(handles.cviewerfig, 'Name', ['CViewer - ' filename]);
    
    if get(handles.radio_showimage,'Value') == get(handles.radio_showimage,'Max'),
        handles.im = handles.loadedim;
        handles.C = fdct_wrapping(handles.im, 0);
        handles.fim = fftshift(fft2(handles.im));
     
        handles.lastpath = pathname;
        handles = reinitlistboxes(handles);
        handles.stickyaxes = 0;
        button_go_Callback(hObject, [], handles);  % do calculations
        handles = guidata(hObject);  % get new handles, since button_go_Callback() changes them
    end
end
guidata(hObject, handles)
cd(curdir);




function edit_thrsh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_thrsh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_thrsh as text
%        str2double(get(hObject,'String')) returns contents of edit_thrsh as a double


% --- Executes during object creation, after setting all properties.
function edit_thrsh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_thrsh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in menu_thrshtype.
function menu_thrshtype_Callback(hObject, eventdata, handles)
% hObject    handle to menu_thrshtype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns menu_thrshtype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menu_thrshtype


% --- Executes during object creation, after setting all properties.
function menu_thrshtype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menu_thrshtype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in menu_plot_curveletpos_top.
function menu_plot_curveletpos_top_Callback(hObject, eventdata, handles)
% hObject    handle to menu_plot_curveletpos_top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

showinaxis(handles.axes_top, get(handles.displaymenu_top,'Value'), get(handles.menu_plot_curveletpos_top,'Value'), ...
    get(handles.plottype_menu_top,'Value'), 1, handles);


% --- Executes during object creation, after setting all properties.
function menu_plot_curveletpos_top_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menu_plot_curveletpos_top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_go.
function button_go_Callback(hObject, eventdata, handles)
% hObject    handle to button_go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mbh = msgbox('Working...');

if get(handles.radio_showimage,'Value') == get(handles.radio_showimage,'Max'),
    
    [nnzcf, nel] = crvlt_countnnz(handles.C);
    set(handles.text_nrcoeffs, 'String', sprintf('Contains %d coeffs.', nel));
    cmax = getlargestcoeff(handles.C);
    set(handles.text_maxvalue, 'String', sprintf('Max coeff value: %.2g', cmax));
    
    C0 = handles.C;
    levels = getlistboxvals(handles.listbox_levels);
    dirs = getlistboxvals(handles.listbox_dirs);
    totdirs = length(get(handles.listbox_dirs,'String'));
    
    for ll=1:length(C0)
        levok = ismember(ll,levels);
        ndir = length(C0{ll});
        if ndir ==1,
            if ~levok,
                C0{ll}{1} = zeros(size(C0{ll}{1}));
            end
        else
            step = round(ndir/totdirs);
            for dd=1:ndir,
                dirok = ismember(floor((dd-1)/step) + 1, dirs); 
                if ~(levok && dirok)
                    C0{ll}{dd} = zeros(size(C0{ll}{dd}));
                end
            end
        end
    end
    
    handles.Cmod = threshold_coeffs(C0, handles);
    handles.rim = ifdct_wrapping(handles.Cmod, 0);
    handles.frim = fftshift(fft2(handles.rim));
    handles.errim = handles.rim - handles.im;
    handles.ferrim = handles.frim - handles.fim;
    
    if strcmp(get(handles.showdirdata_tool, 'State'), 'on') && ~isempty(handles.curdirdata),
        % update directional data
        levels = getlistboxvals(handles.listbox_levels);
        levels = levels(levels > 1 & levels < length(get(handles.listbox_levels,'String')));
        [dirs, vals] = crvlt_getdirdata(handles.C, handles.curdirdata.dotcoords, levels);
        handles.curdirdata.levels = levels;
        handles.curdirdata.dirs = dirs;
        handles.curdirdata.values = vals;
    end    
    
elseif get(handles.radio_showcurvelet,'Value') == get(handles.radio_showcurvelet,'Max'),
    im0 = zeros(handles.imdims);
    C0 = fdct_wrapping(im0,0);
    levels = getlistboxvals(handles.listbox_levels);
    dirs = getlistboxvals(handles.listbox_dirs);
    totdirs = length(get(handles.listbox_dirs,'String'));
    for ll=levels
        ndir = length(C0{ll});
        if ndir ==1,
            cx = floor((size(C0{ll}{1},1)-1)/2) + 1;
            cy = floor((size(C0{ll}{1},2)-1)/2) + 1;
            C0{ll}{1}(cx,cy) = 1;
        else
            step = round(ndir/totdirs);
            for dd=(step*(dirs-1)+1),
                for dofs=0:step-1,
                    didx = dd + dofs;
                    cx = floor((size(C0{ll}{didx},1)-1)/2) + 1;
                    cy = floor((size(C0{ll}{didx},2)-1)/2) + 1;
                    C0{ll}{didx}(cx,cy) = 1;
                end
            end
        end
    end
    
    handles.C = C0;
    handles.Cmod = C0;
    handles.im = ifdct_wrapping(C0, 0);
    handles.fim = fftshift(fft2(handles.im));
    handles.rim = handles.im;
    handles.frim = handles.fim;
    handles.errim = handles.rim - handles.im;
    handles.ferrim = handles.frim - handles.fim;
    
end

close(mbh);

guidata(hObject, handles);
displaydata(handles)



% --- Executes when selected object is changed in displaychoice_panel.
function displaychoice_panel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in displaychoice_panel 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.stickyaxes = 0;
if get(handles.radio_showcurvelet,'Value') == get(handles.radio_showcurvelet,'Max'),
    handles.im = zeros(handles.imdims);
    handles.C = fdct_wrapping(handles.im, 0);
    handles.fim = fftshift(fft2(handles.im));
    handles = reinitlistboxes(handles);
    button_go_Callback(get(hObject,'Parent'), eventdata, handles);
else
    handles.im = handles.loadedim;
    handles.C = fdct_wrapping(handles.im, 0);
    handles.fim = fftshift(fft2(handles.im));
    handles = reinitlistboxes(handles);
    button_go_Callback(get(hObject,'Parent'), eventdata, handles);
end





%
% --- Own helper functions  --- %
%

% True if disptype is the index of a physical space display type
function isphys = isphysspace(disptype)
isphys = ismember(disptype, [1 3 5]);


% Get selected values from a listbox containing values as entries
function vals = getlistboxvals(lbhandle)
% lbhandle      handle to a listbox
sel = get(lbhandle,'Value');
strs = get(lbhandle,'String');
vals = str2double({strs{sel}});


% Display the data stored in handles
function displaydata(handles)


keeplims = handles.stickyaxes;
showinaxis(handles.axes_top, get(handles.displaymenu_top,'Value'), ...
    get(handles.menu_plot_curveletpos_top,'Value'), get(handles.plottype_menu_top,'Value'), keeplims, handles);
if strcmp(get(handles.showdirdata_tool, 'State'), 'on'),
    makedirplot(handles.axes_bottom, handles)
elseif strcmp(get(handles.anisomeas_tool, 'State'), 'on'),
    makeanisoplot(handles.axes_bottom, handles)
else
    showinaxis(handles.axes_bottom, get(handles.displaymenu_bottom,'Value'), 1, ...
        get(handles.plottype_menu_bottom,'Value'), keeplims, handles);
end

l2err = norm(handles.errim(:)/sqrt(numel(handles.errim)));
set(handles.edit_l2error, 'String', sprintf('%.5g', l2err))
psnr = 20 * log10((max(handles.im(:)) - min(handles.im(:))) / l2err);
set(handles.edit_PSNR, 'String', sprintf('%.5f', psnr));




% Display a type of data in an axis
function showinaxis(ax, type, postype, valuetype, keeplims, handles)
% ax        handle to axis to plot in
% type      the plot types in displaymenu_top/bottom
% postype   the type of position plot, as in menu_plot_curveletpos_top/bottom
% handles   figure handles

physspace = isphysspace(type);

if ~physspace,
    % compute frequency axes
    fsize = size(handles.fim);
    fdim = 2*pi*[-floor(fsize(1)/2) floor((fsize(1)-1)/2) ...
        -floor(fsize(2)/2) floor((fsize(2)-1)/2)];
end

curtag = get(ax, 'Tag');  % save to preserve tag
axes(ax)
curlims = axis(ax);
cla
hold off
switch type,
    case 1,
        imagesc(funcimag(handles.im, valuetype));
        axis equal tight;  colorbar
        title('Original image, physical space')
    case 2,
        imagesc(fdim(1:2),fdim(3:4),funcimag(handles.fim, valuetype));
        axis equal tight;  colorbar
        title('Original image, frequency space')
    case 3,
        imagesc(funcimag(handles.rim, valuetype));
        axis equal tight;  colorbar
        title('Reconstructed image, physical space')
    case 4,
        imagesc(fdim(1:2),fdim(3:4),funcimag(handles.frim, valuetype));
        axis equal tight;  colorbar
        title('Reconstructed image, frequency space')
    case 5,
        imagesc(funcimag(handles.errim, valuetype));
        axis equal tight;  colorbar
        title('Error, physical space')
    case 6,
        imagesc(fdim(1:2),fdim(3:4),funcimag(handles.ferrim, valuetype));
        axis equal tight;  colorbar
        title('Error, frequency space')
    case 7, 
        levs = getlistboxvals(handles.listbox_levels);
        dirs = getlistboxvals(handles.listbox_dirs);
        if levs(1) == 1 || levs(1) == length(handles.C),
            dirs = 1;
        end
        imagesc(funcimag(handles.C{levs(1)}{dirs(1)}, valuetype));
        axis square tight;  colorbar
        title(sprintf('Coefficients (level %d, dir %d) - original', levs(1), dirs(1)))
    case 8,
        levs = getlistboxvals(handles.listbox_levels);
        dirs = getlistboxvals(handles.listbox_dirs);
        if levs(1) == 1 || levs(1) == length(handles.Cmod),
            dirs = 1;
        end
        imagesc(funcimag(handles.Cmod{levs(1)}{dirs(1)}, valuetype));
        axis square tight;  colorbar
        title(sprintf('Coefficients (level %d, dir %d) - reconstruction', levs(1), dirs(1)))
end

if physspace,
    switch postype,
        case 1,
            % do nothing
        case 2,  % large dots
            plotcurveletpos(handles.Cmod, ax, 0, 6);
        case 3,  % small dots
            plotcurveletpos(handles.Cmod, ax, 0, 1);
        case 4,  % scaled arrows
            plotcurveletpos(handles.Cmod, ax, 1, 1);
        case 5,  % non-scaled arrows
            plotcurveletpos(handles.Cmod, ax, 1, 0);
    end
end

% set colormap
cmapstr = get(handles.popup_colormap,'String');
map = colormap(cmapstr{get(handles.popup_colormap,'Value')});
if get(handles.checkbox_invertcolormap,'Value') == 1,
    map = flipud(map);
    colormap(map);
end

zoom('reset')   % set maximum zoom out
if keeplims,
    % keep the previous limits
    axis(ax, curlims)
end

% set the callback function for clicking in image (set for all children to
% axes)
if physspace && isequal(ax, handles.axes_top),
    if ~isempty(handles.curdirdata) && strcmp(get(handles.showdirdata_tool,'State'), 'on'),
        hold(ax, 'on')
        crd = handles.curdirdata.dotcoords;
        handles.curdirdata.dothandle = plot(ax, crd(1), crd(2), 'r*');
    end
    hc = get(ax, 'Children');
    set(hc, 'ButtonDownFcn', 'cviewer(''axis_clickcallback'',gcbo,[],guidata(gcbo))');
end
set(ax, 'Tag', curtag)   % reset the axes tag (destroyed by plotting)

handles.stickyaxes = 1;  % make the axes remain the same next time we plot
guidata(gcbo, handles);


% apply function, do averages and sum for directional data
function [levs, dirs, vals] = preparedirdata(dirdata, valuetype, avgsize, sum)

levs = dirdata.levels; dirs = dirdata.dirs; vals = dirdata.values;
nlevs = length(levs);

% apply function
for k=1:nlevs,
    vals{k} = funcimag(vals{k},valuetype);
end

% sum over levels
if sum,
    maxdirs = max(cat(2,dirdata.dirs{:}));
    data = zeros(1,maxdirs);
    for k=1:nlevs,
        step = round(maxdirs/length(dirs{k}));
        data = data + reshape(repmat(vals{k},[step 1]),[1 maxdirs]);
    end
    levs = 1;
    nlevs = 1;
    vals = {data};
    dirs = {1:maxdirs};
end

if avgsize > 1
    % compute sliding average
    for k=1:nlevs,
        lshift = -floor(avgsize/2);
        rshift = avgsize + lshift - 1;
        avgdata = zeros(size(vals{k}));
        for sh = lshift:rshift,
            avgdata = avgdata + circshift(vals{k}, [0 sh]);
        end
        vals{k} = avgdata / avgsize;
    end
end


% plots the coefficient values for all dirs on specified levels
function plotdirections(ax, dirdata, valuetype, plottype, avgsize, sum)
% ax        axes to plot in
% dirdata   levels, dirs, values for directional data - see
%           handles.curdirdata
% valuetype type of value to plot, see funcimag() and plottype_menu_top
% plottype  0 for rectangular, 1 for polar

colors = {'b*-', 'r*-', 'c*-', 'g*-', 'm*-'}; 
curtag = get(ax,'Tag');
axes(ax)
cla
hold off

% apply transformations to data (function, circular average, sum)
[levs, dirs, vals] = preparedirdata(dirdata, valuetype, avgsize, sum);
maxdirs = max(cat(2,dirs{:}));
nlev = length(levs);

if plottype
    % sort levels according to max value to plot largest values first
    % (important for polar plot)
    maxval = zeros(1,nlev);
    for k=1:nlev,
        maxval(k) = max(funcimag(vals{k}(:), valuetype));
    end
    % sort both columns in descending order after first column
    levorder = sortrows([maxval' (1:nlev)'], -1);
else
    levorder = repmat((1:nlev)',[1 2]);
end

legstrs = cell(1,nlev);
cnt = 1;
for k=levorder(:,2)',
    ndirs = length(dirs{k});

    % plot (rectangular or polar)
    if plottype == 0,
        step = maxdirs / ndirs;
        X = step*(0.5:ndirs-0.5) + 0.5;
        plot(ax, X, vals{k}, colors{min(levs(k), length(colors))})
    else
        angles = crvlt_getangle(1:ndirs, ndirs);
        polar(ax, [angles angles(1)], [vals{k} vals{k}(1)], colors{min(levs(k), length(colors))})
    end
    hold on
    if sum,
        legstrs{cnt} = ['Sum of levels ' sprintf('%d,',dirdata.levels)];
    else
        legstrs{cnt} = sprintf('Level %d', levs(k));
    end
    cnt = cnt + 1;
end
legend(ax, legstrs{:},'Location','Best')
zoom('reset')
set(ax,'Tag',curtag)



% plot directions in axes ax using current settings
function makedirplot(ax, handles)

functype = get(handles.plottype_menu_bottom, 'Value');
plottype = (get(handles.polarplot_checkbox, 'Value') == get(handles.polarplot_checkbox, 'Max'));
if get(handles.diraverage_checkbox, 'Value') == get(handles.diraverage_checkbox, 'Max'),
    avgsize = str2double(get(handles.diraverage_edit,'String'));
    if isnan(avgsize)
        errordlg('Not a number in direction plot average size edit box','Curvelet viewer error');
        avgsize = 1;
    end
else
    avgsize = 1;
end
sum = (get(handles.checkbox_sumlevels, 'Value') == get(handles.checkbox_sumlevels, 'Max'));
plotdirections(ax, handles.curdirdata, functype, plottype, avgsize, sum);



% compute and plot anisotropy measure
function makeanisoplot(ax, handles)

A = crvlt_anisomeas(handles.C);
aimg = crvlt_vizaniso(A);
cla(ax)
axes(ax)
imagesc(aimg)
axis equal tight, axis off
colorbar
caxis([0 max(aimg(:))])




% evaluate image, according to choice
function imval = funcimag(im, valuetype)

switch valuetype,
    case 1, % real part
        imval = real(im);
    case 2, % imag part
        imval = imag(im);
    case 3, % absolute value
        imval = abs(im);
    case 4, % log of absolute value
        imval = log(abs(im));
end


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


% call to reinitialize listboxes for new image or new level selection
% assumes that handles.im is set to the new/current image
% and handles.C to its curvelet transform
function handles = reinitlistboxes(handles)

oldlevs = getlistboxvals(handles.listbox_levels);
olddirs = getlistboxvals(handles.listbox_dirs);
oldnrdirs = length(get(handles.listbox_dirs,'String'));
firstclev = min(oldlevs(oldlevs > 1));

% update levels
if length(handles.C) ~= length(get(handles.listbox_levels,'String')),
    str = cell(1,length(handles.C));
    for ll=1:length(handles.C)
        str{ll} = sprintf('%d',ll);
    end
    set(handles.listbox_levels,'String',str);
    if firstclev > length(handles.C),
        % there were more levels in previous image, and only finer levels
        % selected => select finest level
        set(handles.listbox_levels,'Value',length(handles.C));
    else
        set(handles.listbox_levels,'Value',oldlevs(oldlevs <= length(handles.C)));
    end
end

% update dirs
set(handles.listbox_dirs,'ListboxTop',1);
if isempty(firstclev) || firstclev >= length(handles.C),
    % only wavelet levels selected, only one dir
    str = {'1'};
    set(handles.listbox_dirs,'String',str);
    set(handles.listbox_dirs,'Value',1);    
else
    nrdirs = length(handles.C{firstclev});
    step = nrdirs/oldnrdirs;
    str = cell(1,nrdirs);
    dirsel = [];
    for dd=1:nrdirs,
        str{dd} = sprintf('%d',dd);
        if ismember(floor((dd-1)/step) + 1, olddirs),
            dirsel = [dirsel dd];
        end
    end
    set(handles.listbox_dirs,'String',str);
    set(handles.listbox_dirs,'Value',dirsel);
end



% Threshold coefficients according to settings
function Cth = threshold_coeffs(C, handles)
% C         curvelet coeffs to be thresholded
% handles   figure handles

val = str2double(get(handles.edit_thrsh,'String'));
if isnan(val) || (val < 0),
    errordlg('Bad threshold value','Curvelet viewer error!');
    return
end
switch get(handles.menu_thrshtype,'Value'),
    case 1, % do nothing
        Cth = C;
    case 2,  % threshold value
        Cth = crvlt_thresh(C, val);
        set(handles.edit_curthrsh, 'String', sprintf('# nnz: %d', crvlt_countnnz(Cth)));
    case 3,  % number of coeffs
        if val <= 0,
            Cth = zerodct(C);
            smcoef = 0.0;
        else
            [Cth, smcoef] = crvlt_keeplargest(C, round(val));
        end
        set(handles.edit_curthrsh, 'String', sprintf('smallest coeff: %f', smcoef)); 
    case 4,  % percentage (assumes val is in [0,1])
        if (val > 1 || val <= 0)
            errordlg('Bad percentage value. Must be between 0 and 1.','Curvelet viewer error!')
        else
            [Cth, smcoef] = crvlt_keeplargest(C, val);
            set(handles.edit_curthrsh, 'String', sprintf('smallest coeff: %f', smcoef)); 
        end
end
    


% --- Executes on selection change in plottype_menu_top.
function plottype_menu_top_Callback(hObject, eventdata, handles)
% hObject    handle to plottype_menu_top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

showinaxis(handles.axes_top, get(handles.displaymenu_top,'Value'), get(handles.menu_plot_curveletpos_top,'Value'), ...
    get(handles.plottype_menu_top,'Value'), 1, handles);


% --- Executes during object creation, after setting all properties.
function plottype_menu_top_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plottype_menu_top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in plottype_menu_bottom.
function plottype_menu_bottom_Callback(hObject, eventdata, handles)
% hObject    handle to plottype_menu_bottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(get(handles.showdirdata_tool, 'State'), 'off')
    showinaxis(handles.axes_bottom, get(handles.displaymenu_bottom,'Value'), 1, ...
        get(handles.plottype_menu_bottom,'Value'), 1, handles);
else
    makedirplot(handles.axes_bottom, handles)
end


% --- Executes during object creation, after setting all properties.
function plottype_menu_bottom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plottype_menu_bottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in menu_imageselect.
function menu_imageselect_Callback(hObject, eventdata, handles)
% hObject    handle to menu_imageselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

N = 512;
x=0:1/N:1-10*eps;
[X,Y]=meshgrid(x,x);

switch get(hObject,'Value'),
    case 1, % do nothing
        
    case 2, % load image
        loadimage(hObject, handles);
        return;
    case 3  % disc 
        Rlg = 0.4;
        handles.im = double((X-0.5).^2 + (Y-0.5).^2 < Rlg^2);
    case 4 % circle
        Rsm = 0.3;
        im0 = double((X-0.5).^2 + (Y-0.5).^2 < Rsm^2);
        handles.im = ([diff(im0,1,1); zeros(1,N)] + [diff(im0,1,2)  zeros(N,1)]);
        handles.im = double(abs(handles.im) > 0);
    case 5  % drop
        handles.im = eggshape(2, N, [0.5; 0.35], [0.55; 0.7], 0.25, 0.01);
    case 6  % front
        sl = 2.134;
        handles.im = double(sl*X - Y < 0.5);
    case 7  % line
        sl = 2.134;
        width = 0.04;
        handles.im = double((sl*X - Y < 0.5+(width/2)) & (sl*X - Y > 0.5-(width/2)));
    case 8  % c2-discontinuity
        handles.im = exp(-20*((X-0.5).^2 + (Y-0.5).^2));
        handles.im(X.^2 + (Y-1).^2 < 0.66^2) = 0.1 * handles.im(X.^2 + (Y-1).^2 < 0.66^2);
end

if get(handles.checkbox_smooth, 'Value') == get(handles.checkbox_smooth,'Max') ...
        && get(hObject,'Value') > 2,
    sig2 = (1/N)^2;
    gauss=exp(-((X-0.5).^2 + (Y-0.5).^2)/sig2);
    handles.im = ifftshift(ifft2(fft2(gauss).*fft2(handles.im))/sum(gauss(:)));
end

handles.C = fdct_wrapping(handles.im,0);
handles.fim = fftshift(fft2(handles.im));
handles = reinitlistboxes(handles);
guidata(hObject, handles);
handles.stickyaxes = 0;
button_go_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function menu_imageselect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menu_imageselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chk_symmetric.
function chk_symmetric_Callback(hObject, eventdata, handles)
% hObject    handle to chk_symmetric (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_symmetric



function edit_curthrsh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_curthrsh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_curthrsh as text
%        str2double(get(hObject,'String')) returns contents of edit_curthrsh as a double


% --- Executes during object creation, after setting all properties.
function edit_curthrsh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_curthrsh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_l2error_Callback(hObject, eventdata, handles)
% hObject    handle to edit_l2error (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_l2error as text
%        str2double(get(hObject,'String')) returns contents of edit_l2error as a double


% --- Executes during object creation, after setting all properties.
function edit_l2error_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_l2error (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_PSNR_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PSNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PSNR as text
%        str2double(get(hObject,'String')) returns contents of edit_PSNR as a double


% --- Executes during object creation, after setting all properties.
function edit_PSNR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PSNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_smooth.
function checkbox_smooth_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.menu_imageselect, 'Value') > 2,
    % predefined image, rerun selection with new smoothing option
    menu_imageselect_Callback(handles.menu_imageselect, [], handles);
end

    

% executes on button click in top axes
function axis_clickcallback(hObject, eventdata, handles)

if strcmp(get(handles.showdirdata_tool, 'State'), 'on') && strcmp(get(get(hObject, 'Parent'),'Tag'), 'axes_top'),
    clickdata = get(handles.axes_top, 'CurrentPoint');
    pt = clickdata(1,[2 1]);  % swap x and y because y is first index in matrix

    % plot clicked position
    if ~isempty(handles.curdirdata) && ishandle(handles.curdirdata.dothandle)
        delete(handles.curdirdata.dothandle);
    end
    hold(handles.axes_top, 'on')
    hdot = plot(handles.axes_top, pt(2), pt(1), 'r*');
    set(hdot, 'ButtonDownFcn', 'cviewer(''axis_clickcallback'',gcbo,[],guidata(gcbo))');
    
    % get and save directional data
    levels = getlistboxvals(handles.listbox_levels);
    levels = levels(levels > 1 & levels < length(get(handles.listbox_levels,'String')));
    [dirs, vals] = crvlt_getdirdata(handles.C, pt, levels);
    handles.curdirdata = struct('levels',levels,'dirs',{dirs},'values',{vals},'dothandle',hdot,'dotcoords',pt([2 1]));
    
    % plot the directional data in bottom axes
    makedirplot(handles.axes_bottom, handles)

    guidata(handles.axes_top, handles);
end


% --------------------------------------------------------------------
function linkaxestool_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to linkaxestool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Don't allow linked axes when showing directions
if strcmp(get(handles.showdirdata_tool, 'State'), 'on')
    set(hObject, 'State', 'off')
end

% turn on or off linked axes
if strcmp(get(hObject, 'State'),'on')
    linkaxes([handles.axes_top handles.axes_bottom])
else
    linkaxes([handles.axes_top handles.axes_bottom], 'off')
end


% --------------------------------------------------------------------
function showdirdata_tool_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to showdirdata_tool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(get(hObject, 'State'), 'on'),
    % turn off zoom and pan functions, and axis link
    linkaxes([handles.axes_top handles.axes_bottom], 'off')
    set(handles.linkaxestool, 'State', 'off')
    zoom(handles.axes_top, 'off')
    zoom(handles.axes_bottom, 'off')
    pan(handles.axes_top, 'off')
    pan(handles.axes_bottom, 'off')
    set(handles.anisomeas_tool, 'State', 'off')
    
    % make the bottom plot a directional plot
    set(handles.displaymenu_bottom, 'Enable', 'off')
    set(handles.plottype_menu_bottom, 'Value', 3)  % plot absolute values as default
    %set(handles.polarplot_checkbox, 'Visible', 'on')  % show polar/rectangular checkbox
    set(handles.dirdata_panel, 'Visible', 'on')  % show dirdata plot tools
else
    % make the bottom plot a normal plot
    set(handles.displaymenu_bottom, 'Enable', 'on')    
    %set(handles.polarplot_checkbox, 'Visible', 'off')  % hide polar/rectangular checkbox    
    set(handles.dirdata_panel, 'Visible', 'off')  % hide dirdata plot tools
end


% --- Executes on button press in polarplot_checkbox.
function polarplot_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to polarplot_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

makedirplot(handles.axes_bottom, handles)


% --- Executes on button press in showseparate_button_top.
function showseparate_button_top_Callback(hObject, eventdata, handles)
% hObject    handle to showseparate_button_top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% draw the top image in a new figure
hf = figure;
showinaxis(gca, get(handles.displaymenu_top,'Value'), get(handles.menu_plot_curveletpos_top,'Value'), ...
    get(handles.plottype_menu_top,'Value'), 0, handles);
axis(axis(handles.axes_top));


% --- Executes on button press in showseparate_button_bottom.
function showseparate_button_bottom_Callback(hObject, eventdata, handles)
% hObject    handle to showseparate_button_bottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% draw the bottom image in a new figure
hf = figure;
if strcmp(get(handles.showdirdata_tool, 'State'), 'off'),
    showinaxis(gca, get(handles.displaymenu_bottom,'Value'), 1, ...
        get(handles.plottype_menu_bottom,'Value'), 0, handles);
else
    makedirplot(gca, handles)
end
axis(axis(handles.axes_bottom));


% --- Executes on button press in exportdata_button.
function exportdata_button_Callback(hObject, eventdata, handles)
% hObject    handle to exportdata_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cviewdata = struct('im',[],'fim',[],'rim',[],'frim',[],'errim',[],'ferrim',[],'C',[],'Cmod',[],'dirdata',[]);
cviewdata.im = handles.im;
cviewdata.fim = handles.fim;
cviewdata.rim = handles.rim;
cviewdata.frim = handles.frim;
cviewdata.errim = handles.errim;
cviewdata.ferrim = handles.ferrim;
cviewdata.C = handles.C;
cviewdata.Cmod = handles.Cmod;
cviewdata.dirdata = handles.curdirdata;
assignin('base', 'cviewdata', cviewdata)
msgbox('Data has been exported to the variable ''cviewdata'' in base workspace.','Export data')


% --- Executes on button press in diraverage_checkbox.
function diraverage_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to diraverage_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

makedirplot(handles.axes_bottom, handles)


function diraverage_edit_Callback(hObject, eventdata, handles)
% hObject    handle to diraverage_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of diraverage_edit as text
%        str2double(get(hObject,'String')) returns contents of diraverage_edit as a double


% --- Executes during object creation, after setting all properties.
function diraverage_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diraverage_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menu_tools_Callback(hObject, eventdata, handles)
% hObject    handle to menu_tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuitem_makemovie_Callback(hObject, eventdata, handles)
% hObject    handle to menuitem_makemovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cmoviemaker;


% --- Executes on button press in checkbox_sumlevels.
function checkbox_sumlevels_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_sumlevels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

makedirplot(handles.axes_bottom, handles)


% --- Executes on selection change in popup_colormap.
function popup_colormap_Callback(hObject, eventdata, handles)
% hObject    handle to popup_colormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

displaydata(handles);


% --- Executes during object creation, after setting all properties.
function popup_colormap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_colormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_invertcolormap.
function checkbox_invertcolormap_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_invertcolormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

displaydata(handles);


% --------------------------------------------------------------------
function anisomeas_tool_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to anisomeas_tool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(get(hObject, 'State'), 'on'),
    % turn off zoom and pan functions, and axis link
    linkaxes([handles.axes_top handles.axes_bottom], 'off')
    set(handles.linkaxestool, 'State', 'off')
    zoom(handles.axes_top, 'off')
    zoom(handles.axes_bottom, 'off')
    pan(handles.axes_top, 'off')
    pan(handles.axes_bottom, 'off')
    set(handles.showdirdata_tool, 'State', 'off')
    
    % make the bottom plot a directional plot
    set(handles.displaymenu_bottom, 'Enable', 'off')
    set(handles.plottype_menu_bottom, 'Enable', 'off')  % disable plot type
    
    handles.stickyaxes = 0;
    displaydata(handles);
else
    % make the bottom plot a normal plot
    set(handles.displaymenu_bottom, 'Enable', 'on')    
    set(handles.plottype_menu_bottom, 'Enable', 'on')  % enable plot type
end


function menu_batch_Callback(hObject, eventdata, handles)

thresh = str2double(get(handles.edit_thrsh,'String'));
if isnan(thresh) || (thresh < 0),
    errordlg('Bad threshold value','Curvelet viewer error!');
    return
end
thtype = get(handles.menu_thrshtype,'Value');
levels = getlistboxvals(handles.listbox_levels);
dirs = getlistboxvals(handles.listbox_dirs);
totdirs = length(get(handles.listbox_dirs,'String'));
batchwindow(thresh, thtype, levels, dirs, totdirs);

