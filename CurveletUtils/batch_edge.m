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

handles.mainhandles = varargin{1};

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


function cval = getpopupvalue(hpopup)

strs = get(hpopup, 'String');
cval = str2double(strs{get(hpopup, 'Value')});



function C = applycurvelettransform(im, allcurvelet)

[m,n] = size(im);
nbscales = floor(log2(min(m,n)))-3;
nbangles_coarse = 16;

%call mex function
C = fdct_wrapping_mex(m,n,nbscales, nbangles_coarse, allcurvelet, im);



% see edgedetect.m:findedges()
function edgedata = findedges(im, C, handles)
% handles    structure with handles and user data

drad = str2double(get(handles.edit_tophatradius, 'String'));
if isempty(drad),
    close(hbox)
    errordlg('Bad disk radius value!', 'EdgeDetect error')
    return
end


cavg = getpopupvalue(handles.popup_nrcircavg);
crvltsize = getpopupvalue(handles.popup_crvltsize);
nrdirfld = getpopupvalue(handles.popup_nrfields);
dirspace = getpopupvalue(handles.popup_dirspace);
fltrtype = get(handles.popup_filtertype, 'Value') - 1;
levs = [get(handles.listbox_levels, 'Value')];

edgedata.fld = crvlt_extractdirs(C, levs, nrdirfld, cavg, crvltsize, dirspace, fltrtype);

edgedata.totmag = crvlt_getmagnitude(C, levs, crvltsize);
    

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
    mag = sqrt(edgedata.fld{1}{1}.^2 + edgedata.fld{1}{2}.^2);
    newmag = imtophat(mag, strel('disk', drad)); 
    newfld{1} = {edgedata.fld{1}{1} .* newmag ./ (mag + 10*eps), edgedata.fld{1}{2} .* newmag ./ (mag + 10*eps)};
    edgedata.e = curvecanny_multi(newfld, [thlow thhigh], extlen, [extglth extlcth], thinedge, dilatethin, despur);
else
    edgedata.e = curvecanny_multi(edgedata.fld, [thlow thhigh], extlen, [extglth extlcth], thinedge, dilatethin, despur);
end

% compute sheets
shtthrsh = str2double(get(handles.edit_sheetthrsh, 'String'));
if isempty(shtthrsh),
    close(hbox)
    errordlg('Bad value for sheet threshold!', 'EdgeDetect error')
    return
end

opentotmag = imopen(edgedata.totmag, strel('disk', drad));
edgedata.sheet = (opentotmag > max(opentotmag(:)) * shtthrsh);

% exclude sheets from edges
if get(handles.check_excludesheets, 'Value') > 0,
    edgedata.e = edgedata.e & ~edgedata.sheet;
end

% get edges in image resolution
eh = 1 ./ size(edgedata.e);
ih = 1 ./ size(im);
edgedata.E = interp2((-0:eh(2):1-eh(2)/2)', -0:eh(1):1-eh(1)/2, double(edgedata.e), ...
    (0:ih(2):1-ih(2)/2)', 0:ih(1):1-ih(1)/2, 'linear', 0);
edgedata.E = edgedata.E > 0.4;

% compute and show statistics
eskel = bwmorph(edgedata.E, 'thin', Inf);
[L, ncmp] = bwlabel(eskel, 8);
totlen = 0;
for k=1:ncmp,
    totlen = totlen + edgelength(L == k);
end
edgedata.totlen = totlen/size(im,2);
edgedata.sheetarea = nnz(edgedata.sheet) / numel(edgedata.sheet);

% save skeletonized edges
if get(handles.check_skeletonize, 'Value') > 0
    edgedata.E = eskel;
end



function outim = makeimage(im, edgedata, mainhandles)

type = get(mainhandles.popup_displaytype, 'Value');
fldnr = get(mainhandles.popup_fieldsel, 'Value');
dispedge = get(mainhandles.check_displayedges, 'Value');
dispsheet = get(mainhandles.check_displaysheets, 'Value');
%dispfield = get(mainhandles.check_displayfield, 'Value');

switch type,
    case 1,
        outim = im - min(im(:));
    case 2,
        mag = sqrt(edgedata.fld{fldnr}{1}.^2 + edgedata.fld{fldnr}{2}.^2);
        outim = mag - min(mag(:));
end
if dispsheet,
    outim = uint8(repmat(outim/max(outim(:)) * 191, [1 1 3]));
else
    outim = uint8(repmat(outim/max(outim(:)) * 255, [1 1 3]));
end

switch type,
    case 1,
%        scx = size(handles.im, 2) / size(handles.e, 2);
%        scy = size(handles.im, 1) / size(handles.e, 1);
        if dispsheet,
            shr = imresize(edgedata.sheet, size(im));
            outim = outim + repmat(64*uint8(shr), [1 1 3]);
        end
        if dispedge
            EI = find(edgedata.E);
            outim(EI + numel(im)) = 255;
            outim(EI) = 0;
            outim(EI + 2*numel(im)) = 0;
        end
%         if dispfield,
%             [X,Y] = meshgrid(1 + (0:size(handles.e,2)-1) * scx, 1 + (0:size(handles.e,1)-1) * scy);
%             quiver(X, Y, handles.fld{fldnr}{1}, - handles.fld{fldnr}{2}, 'r')
%         end
  case 2,
        if dispsheet,
            outim = outim + repmat(64*uint8(edgedata.sheet), [1 1 3]);
        end
        if dispedge,
            eI = find(edgedata.e);
            outim(eI + numel(edgedata.e)) = 255;
            outim(eI) = 0;
            outim(eI + 2*numel(edgedata.e)) = 0;
        end
%         if dispfield,
%             [X,Y] = meshgrid(1:size(handles.e,2), 1:size(handles.e,1));
%             quiver(X, Y, handles.fld{fldnr}{1}, - handles.fld{fldnr}{2}, 'r')
%         end
end




function button_go_Callback(hObject, eventdata, handles)

if isempty(handles.dirname),
    errordlg('No directory specified', 'EdgeDetect error');
    return
end
files = dir(handles.dirname);

if ~exist(fullfile(handles.dirname, 'Output'), 'dir'),
    success = mkdir(fullfile(handles.dirname, 'Output'));
    if ~success,
        errordlg('Error creating Output directory', 'EdgeDetect error');
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
        continue;   % not an image file
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
    
    allcurvelets = get(handles.mainhandles.check_allcrvlt, 'Value');
    C = applycurvelettransform(im, allcurvelets);
    edgedata = findedges(im, C, handles.mainhandles);

    if get(handles.check_writecurvelets, 'Value')
        edgedata.C = C;
    end

    % print results
    disp(sprintf('%s: total edge length: %.3f (of image width); sheet area fraction: %.3f', files(k).name, edgedata.totlen, edgedata.sheetarea))
    
    % write data file
    [pt, nm, ext] = fileparts(files(k).name);
    matfilename = fullfile(handles.dirname, 'Output', [nm '_edges.mat']);
    save(matfilename, 'edgedata');
    
    % write image file
    if get(handles.check_writeimages, 'Value'),
        imfilename = fullfile(handles.dirname, 'Output', [nm '_edges' ext]);
        
        outim = makeimage(im, edgedata, handles.mainhandles);
        imwrite(outim, imfilename);
    end

    cnt = cnt + 1;

end

if cnt == 0,
    errordlg('No images found', 'EdgeDetect error');
else
    msgbox(sprintf('Analyzed %d images.', cnt));
end

uiresume(handles.figure1);
close(handles.figure1);



function check_writeimages_Callback(hObject, eventdata, handles)


function check_writecurvelets_Callback(hObject, eventdata, handles)



