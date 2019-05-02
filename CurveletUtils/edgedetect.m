function varargout = edgedetect(varargin)
% EDGEDETECT M-file for edgedetect.fig
%      EDGEDETECT, by itself, creates a new EDGEDETECT or raises the existing
%      singleton*.
%
%      H = EDGEDETECT returns the handle to a new EDGEDETECT or the handle to
%      the existing singleton*.
%
%      EDGEDETECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EDGEDETECT.M with the given input arguments.
%
%      EDGEDETECT('Property','Value',...) creates a new EDGEDETECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before edgedetect_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to edgedetect_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help edgedetect

% Last Modified by GUIDE v2.5 02-Jul-2008 09:58:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @edgedetect_OpeningFcn, ...
                   'gui_OutputFcn',  @edgedetect_OutputFcn, ...
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


% --- Executes just before edgedetect is made visible.
function edgedetect_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to edgedetect (see VARARGIN)

% Choose default command line output for edgedetect
handles.output = hObject;

handles.lastpath = [];
handles.im = [];
handles.C = {};
handles.fld = {};
handles.totmag = [];
handles.e = [];
handles.E = [];
handles.sheet = [];
handles.prevvmode = 0;
handles.prevlevs = [];

% Update handles structure
guidata(hObject, handles);

setpath


% UIWAIT makes edgedetect wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = edgedetect_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuitem_loadfile_Callback(hObject, eventdata, handles)
% hObject    handle to menuitem_loadfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

loadimage(hObject, handles);



% --- Load image from file 
function loadimage(hObject, handles)
% hObject    handle to button_loadimage (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

curdir = pwd;
if ~isempty(handles.lastpath) && isdir(handles.lastpath),
    cd(handles.lastpath)
end
[filename, pathname]=uigetfile({'*' 'All files (*.*)'}, 'Open image');
if ~(isequal(filename,0) || isequal(pathname,0))
    [pstr,nm,ext] = fileparts(fullfile(pathname,filename));
    if strcmp(ext,'.dcm')
        im = dicomread(fullfile(pathname,filename));
    else
        im = imread(fullfile(pathname,filename));
    end
    
    set(handles.figure1, 'Name', ['EdgeDetect - ' filename]);
    
    % get one color channel
    if (ndims(im) > 2),
        nzv = zeros(1,min(size(im,3),3));
        for k=1:min(size(im,3),3)
            nzv(k) = nnz(im(:,:,k));
        end
        [M,I] = max(nzv);
        im = im(:,:,I);
    end

    handles.im = double(im);
    
    hmb = msgbox('Working...');
    
    handles.C = applycurvelettransform(handles.im, get(handles.check_allcrvlt, 'Value'));
    
    curlevs = get(handles.listbox_levels, 'Value');
    strs = cell(1,length(handles.C));
    for k=1:length(handles.C), strs{k} = sprintf('%d',k); end
    set(handles.listbox_levels, 'String', strs);
    if isempty(handles.lastpath),
        set(handles.listbox_levels, 'Value', length(handles.C)-1);
    else
        set(handles.listbox_levels, 'Value', min(curlevs, length(handles.C)-1));
    end
    
    handles.lastpath = pathname;
    handles.prevvmode = 0;
    findedges(handles, 1);
    handles = guidata(hObject);  % get new handles, since findedges() changes them

    close(hmb)
end
guidata(hObject, handles)
cd(curdir);


function C = applycurvelettransform(im, allcurvelet)

[m,n] = size(im);
nbscales = floor(log2(min(m,n)))-3;
nbangles_coarse = 16;

%call mex function
C = fdct_wrapping_mex(m,n,nbscales, nbangles_coarse, allcurvelet, im);


% Get the current value from a popup menu with numerical values
function cval = getpopupvalue(hpopup)

strs = get(hpopup, 'String');
cval = str2double(strs{get(hpopup, 'Value')});


% --- Do the computations
function findedges(handles, step)
% handles    structure with handles and user data

if nargin < 2,
    step = 1;
end

hbox = msgbox('Working...','non-modal');

drad = str2double(get(handles.edit_tophatradius, 'String'));
if isempty(drad),
    close(hbox)
    errordlg('Bad disk radius value!', 'EdgeDetect error')
    return
end

if step == 1,
    cavg = getpopupvalue(handles.popup_nrcircavg);
    crvltsize = getpopupvalue(handles.popup_crvltsize);
    nrdirfld = getpopupvalue(handles.popup_nrfields);
    dirspace = getpopupvalue(handles.popup_dirspace);
    fltrtype = get(handles.popup_filtertype, 'Value') - 1;
    levs = [get(handles.listbox_levels, 'Value')];

    handles.fld = crvlt_extractdirs(handles.C, levs, nrdirfld, cavg, crvltsize, dirspace, fltrtype);
    
    handles.totmag = crvlt_getmagnitude(handles.C, levs, crvltsize);
    
end

thlow = str2double(get(handles.edit_thrshlow, 'String'));
thhigh = str2double(get(handles.edit_thrshhigh, 'String'));
if isempty(thlow) || isempty(thhigh),
    close(hbox)
    errordlg('Bad threshold value!', 'EdgeDetect error!')
    return
elseif thlow >= thhigh,
    close(hbox)
    errordlg('Low threshold must be smaller than high threshold!','EdgeDetect error!')
    return
end
extlen = get(handles.popup_extendlen, 'Value') - 1;
extglth = str2double(get(handles.edit_global_extth, 'String'));
extlcth = str2double(get(handles.edit_local_extth, 'String'));
if isempty(extglth) || isempty(extlcth),
    close(hbox)
    errordlg('Bad extension threshold value!', 'EdgeDetect error')
    return
end
thinedge = getpopupvalue(handles.popup_thinedges);
dilatethin = getpopupvalue(handles.popup_dilatethin);
despur = getpopupvalue(handles.popup_despur);

%do tophat transformation
if get(handles.check_tophat, 'Value') > 0,
    mag = sqrt(handles.fld{1}{1}.^2 + handles.fld{1}{2}.^2);
    newmag = imtophat(mag, strel('disk', drad)); 
    newfld{1} = {handles.fld{1}{1} .* newmag ./ (mag + 10*eps), handles.fld{1}{2} .* newmag ./ (mag + 10*eps)};
    handles.e = curvecanny_multi(newfld, [thlow thhigh], extlen, [extglth extlcth], thinedge, dilatethin, despur);
    %handles.fld = newfld;
    
    %figure(1), clf
    %subplot(1,3,1), imagesc(mag), axis equal tight; colormap gray; colorbar, hold on
    %subplot(1,3,2), imagesc(mag-newmag > 0.9*mean(mag(:))), axis equal tight; colormap gray; colorbar, hold on
    %title(sprintf('%f', nnz(mag-newmag > 0.9*mean(mag(:)))/numel(mag)))
    %subplot(1,3,3), imagesc(newmag), axis equal tight; colormap gray; colorbar, hold on
    %quiver(newfld{1}{1}, -newfld{1}{2}, 'r')
else
    handles.e = curvecanny_multi(handles.fld, [thlow thhigh], extlen, [extglth extlcth], thinedge, dilatethin, despur);
end

% compute sheets
shtthrsh = str2double(get(handles.edit_sheetthrsh, 'String'));
if isempty(shtthrsh),
    close(hbox)
    errordlg('Bad value for sheet threshold!', 'EdgeDetect error')
    return
end

opentotmag = imopen(handles.totmag, strel('disk', drad));
handles.sheet = (opentotmag > max(opentotmag(:)) * shtthrsh);
%handles.sheet = (opentotmag > mean(handles.totmag(:)) * shtthrsh);
%figure(1), clf
%imagesc(opentotmag), axis equal tight, colormap gray

% exclude sheets from edges
if get(handles.check_excludesheets, 'Value') > 0,
    handles.e = handles.e & ~handles.sheet;
end

% get edges in image resolution
eh = 1 ./ size(handles.e);
ih = 1 ./ size(handles.im);
handles.E = interp2((-0:eh(2):1-eh(2)/2)', -0:eh(1):1-eh(1)/2, double(handles.e), ...
    (0:ih(2):1-ih(2)/2)', 0:ih(1):1-ih(1)/2, 'linear', 0);
handles.E = handles.E > 0.4;

% compute and show statistics
eskel = bwmorph(handles.E, 'thin', Inf);
[L, ncmp] = bwlabel(eskel, 8);
totlen = 0;
for k=1:ncmp,
    totlen = totlen + edgelength(L == k);
end
totlen = totlen/size(handles.im,2);
sha = nnz(handles.sheet) / numel(handles.sheet);

set(handles.text_numcmp, 'String', sprintf('# components: %d', ncmp));
set(handles.text_totlen, 'String', sprintf('total length: %f * imwidth', totlen));
set(handles.text_sheetarea, 'String', sprintf('sheet area: %f %%', sha*100));


% save skeletonized edges
if get(handles.check_skeletonize, 'Value') > 0
    handles.E = eskel;
end

% save and display data
guidata(gcbo, handles);

displaydata(handles)
close(hbox)



% Get current view mode (image (=1) or magnitudes (=2))
function vm = getviewmode(viewtype)

vm = int32(ismember(viewtype, 2)) + 1;



% --- display results
function displaydata(handles, keepax)

type = get(handles.popup_displaytype, 'Value');
vmode = getviewmode(type);
clevs = get(handles.listbox_levels, 'Value');
fldnr = get(handles.popup_fieldsel, 'Value');
dispedge = get(handles.check_displayedges, 'Value');
dispsheet = get(handles.check_displaysheets, 'Value');
dispfield = get(handles.check_displayfield, 'Value');

if nargin < 2,
    keepax = (handles.prevvmode == vmode && isequal(clevs,handles.prevlevs));  % keep the axes if we haven't changed mode
end

axes(handles.axes_main)
if keepax,
    curax = axis;
end
cla

if getviewmode(type) == 1,
    imagesc(handles.im);
    axis equal tight
    colormap gray
    %colorbar off
else
    mag = sqrt(handles.fld{fldnr}{1}.^2 + handles.fld{fldnr}{2}.^2);
    imagesc(mag)
    axis equal tight
    colormap gray
    colorbar
end
hold on

switch type,
    case 1,
        scx = size(handles.im, 2) / size(handles.e, 2);
        scy = size(handles.im, 1) / size(handles.e, 1);
        if dispsheet,
            shr = imresize(handles.sheet, size(handles.im));
            hsh = image(255*shr);
            alpha(hsh, 0.3);
        end
        if dispedge
            [Ex,Ey] = find(handles.E > 0.5);
            plot(Ey, Ex, 'g.');
        end
        if dispfield,
            [X,Y] = meshgrid(1 + (0:size(handles.e,2)-1) * scx, 1 + (0:size(handles.e,1)-1) * scy);
            quiver(X, Y, handles.fld{fldnr}{1}, - handles.fld{fldnr}{2}, 'r')
        end
  case 2,
        if dispsheet,
            hsh = image(255*handles.sheet);
            alpha(hsh, 0.3);
        end
        if dispedge,
            [ex,ey] = find(handles.e);
            plot(ey,ex, 'g.');
        end
        if dispfield,
            [X,Y] = meshgrid(1:size(handles.e,2), 1:size(handles.e,1));
            quiver(X, Y, handles.fld{fldnr}{1}, - handles.fld{fldnr}{2}, 'r')
        end
end

if keepax,
    axis(curax)
else 
    zoom reset
end

handles.prevvmode = vmode;
handles.prevlevs = clevs;
guidata(gcbo, handles);


% --- Executes on selection change in popup_nrcircavg.
function popup_nrcircavg_Callback(hObject, eventdata, handles)
% hObject    handle to popup_nrcircavg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

findedges(handles, 1);


% --- Executes during object creation, after setting all properties.
function popup_nrcircavg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_nrcircavg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_crvltsize.
function popup_crvltsize_Callback(hObject, eventdata, handles)
% hObject    handle to popup_crvltsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

findedges(handles, 1);


% --- Executes during object creation, after setting all properties.
function popup_crvltsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_crvltsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_nrfields.
function popup_nrfields_Callback(hObject, eventdata, handles)
% hObject    handle to popup_nrfields (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% update the fieldsel popup menu
strs = {};
nflds = get(hObject, 'Value');
for k = 1:nflds, strs{k} = sprintf('%d',k); end
set(handles.popup_fieldsel, 'String', strs);
set(handles.popup_fieldsel, 'Value', 1);
guidata(hObject, handles)

findedges(handles, 1);


% --- Executes during object creation, after setting all properties.
function popup_nrfields_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_nrfields (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_levels.
function listbox_levels_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_levels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

findedges(handles, 1)


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



function edit_thrshlow_Callback(hObject, eventdata, handles)
% hObject    handle to edit_thrshlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

findedges(handles, 2);


% --- Executes during object creation, after setting all properties.
function edit_thrshlow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_thrshlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_thrshhigh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_thrshhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

findedges(handles, 2);


% --- Executes during object creation, after setting all properties.
function edit_thrshhigh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_thrshhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_extendlen.
function popup_extendlen_Callback(hObject, eventdata, handles)
% hObject    handle to popup_extendlen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

findedges(handles, 2);


% --- Executes during object creation, after setting all properties.
function popup_extendlen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_extendlen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in popup_thinedges.
function popup_thinedges_Callback(hObject, eventdata, handles)
% hObject    handle to popup_thinedges (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

findedges(handles, 2);



% --- Executes on button press in popup_dilatethin.
function popup_dilatethin_Callback(hObject, eventdata, handles)
% hObject    handle to popup_dilatethin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

findedges(handles, 2);


% --- Executes on selection change in popup_dirspace.
function popup_dirspace_Callback(hObject, eventdata, handles)
% hObject    handle to popup_dirspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

findedges(handles, 1);


% --- Executes during object creation, after setting all properties.
function popup_dirspace_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_dirspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_displaytype.
function popup_displaytype_Callback(hObject, eventdata, handles)
% hObject    handle to popup_displaytype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

displaydata(handles)


% --- Executes during object creation, after setting all properties.
function popup_displaytype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_displaytype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_fieldsel.
function popup_fieldsel_Callback(hObject, eventdata, handles)
% hObject    handle to popup_fieldsel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

displaydata(handles, 1)

% --- Executes during object creation, after setting all properties.
function popup_fieldsel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_fieldsel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in popup_despur.
function popup_despur_Callback(hObject, eventdata, handles)
% hObject    handle to popup_despur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

findedges(handles, 2);



function edit_global_extth_Callback(hObject, eventdata, handles)
% hObject    handle to edit_global_extth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

findedges(handles, 2);


% --- Executes during object creation, after setting all properties.
function edit_global_extth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_global_extth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_local_extth_Callback(hObject, eventdata, handles)
% hObject    handle to edit_local_extth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

findedges(handles, 2);


% --- Executes during object creation, after setting all properties.
function edit_local_extth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_local_extth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in check_tophat.
function check_tophat_Callback(hObject, eventdata, handles)
% hObject    handle to check_tophat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

findedges(handles, 2);

% --- Executes on button press in check_excludesheets.
function check_excludesheets_Callback(hObject, eventdata, handles)
% hObject    handle to check_excludesheets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

findedges(handles, 2);



function edit_tophatradius_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tophatradius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

findedges(handles, 2);


% --- Executes during object creation, after setting all properties.
function edit_tophatradius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tophatradius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in check_skeletonize.
function check_skeletonize_Callback(hObject, eventdata, handles)
% hObject    handle to check_skeletonize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

findedges(handles, 2);


% --- Executes on button press in check_displayedges.
function check_displayedges_Callback(hObject, eventdata, handles)
% hObject    handle to check_displayedges (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

displaydata(handles)


% --- Executes on button press in check_displaysheets.
function check_displaysheets_Callback(hObject, eventdata, handles)
% hObject    handle to check_displaysheets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

displaydata(handles)


% --- Executes on button press in check_displayfield.
function check_displayfield_Callback(hObject, eventdata, handles)
% hObject    handle to check_displayfield (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

displaydata(handles)


function edit_sheetthrsh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sheetthrsh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

findedges(handles, 2);


% --- Executes during object creation, after setting all properties.
function edit_sheetthrsh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sheetthrsh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_export.
function button_export_Callback(hObject, eventdata, handles)
% hObject    handle to button_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

est = struct('image', handles.im, 'cc', {handles.C}, 'field', {handles.fld}, 'edges_fld', handles.e, 'edges', handles.E, 'sheets', handles.sheet);
assignin('base', 'edgedata', est);
msgbox('Data was exported to base workspace in struct "edgedata".');


% --- Executes on button press in check_allcrvlt.
function check_allcrvlt_Callback(hObject, eventdata, handles)
% hObject    handle to check_allcrvlt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.C),
    allcurvelets = get(hObject, 'Value');
    hb = msgbox('Working...');
    handles.C = applycurvelettransform(handles.im, allcurvelets);
    guidata(hObject, handles);
    close(hb)
    
    % recompute if finest level selected
    if ismember(length(get(handles.listbox_levels, 'String')), [get(handles.listbox_levels, 'Value')]),
        findedges(handles, 1);
    end
end


% --- Executes on selection change in popup_filtertype.
function popup_filtertype_Callback(hObject, eventdata, handles)
% hObject    handle to popup_filtertype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

findedges(handles, 1);

% --- Executes during object creation, after setting all properties.
function popup_filtertype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_filtertype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function menu_batch_Callback(hObject, eventdata, handles)

batch_edge(handles);



