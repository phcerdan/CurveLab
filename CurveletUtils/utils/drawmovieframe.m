function hfig = drawmovieframe(hfig, fullfilename, layout, plotinfo, computeinfo)

% Draws a movieframe in figure hfig, using plotinfo for layout and
% computeinfo for computations (see cmoviemaker for structure members)

wtypes = {'db3', 'db5', 'dmey'};


% load image file
[pstr,nm,ext,ver] = fileparts(fullfilename);
if strcmp(ext,'.dcm')
    im = dicomread(fullfilename);
else
    im = imread(fullfilename);
end
if (ndims(im) > 2),
    nzv = zeros(1,size(im,3));
    for k=1:size(im,3)
        nzv(k) = nnz(im(:,:,k));
    end
    [M,I] = max(nzv);
    im = im(:,:,I);
end
im = double(im);

% invert and crop if options set
if computeinfo.invert,
    im = max(im(:)) - im;
end
if isempty(computeinfo.xlim),
    xvals = 1:size(im,2);
else
    xvals = floor(computeinfo.xlim(1):min(computeinfo.xlim(2), size(im,2)));
end
if isempty(computeinfo.ylim),
    yvals = 1:size(im,1);
else
    yvals = floor(computeinfo.ylim(1):min(computeinfo.ylim(2),size(im,1)));
end
im = im(yvals, xvals);

% compute curvelet/wavelet transforms
if any(cat(2,plotinfo.type) > 2),  % do we need to compute transforms?
    if ismember(1,cat(2, plotinfo.transform)),
        C = fdct_wrapping(im, 0);
        Cmod = select_crvltlevels(C, computeinfo);
        Cmod = thrsh_crvltcoeffs(Cmod, computeinfo);
    end
    D = cell(1,max(cat(2, plotinfo.transform)));
    S = cell(1,max(cat(2, plotinfo.transform)));
    Dmod = cell(1,max(cat(2, plotinfo.transform)));
    for k = unique(cat(2, plotinfo.transform)),
        if k>1,
            nlev = max([1 (wmaxlev(size(im), wtypes{k-1}) - 1)]);   % estimate number of useful levels
            [D{k}, S{k}] = wavedec2(im, nlev, wtypes{k-1});
            Dmod{k} = thrsh_waveletcoeffs(D{k}, computeinfo, k);
        end
    end
end

fsize = size(im);
fdim = 2*pi*[-floor(fsize(1)/2) floor((fsize(1)-1)/2) ...
    -floor(fsize(2)/2) floor((fsize(2)-1)/2)];

figure(hfig)
clf
nplots = prod(layout);
for k=1:nplots,
    ax = subplot(layout(1), layout(2), k);
    
    physspace = ismember(plotinfo(k).type, [1 3 5]);
    
    if plotinfo(k).type > 2,
        % make image reconstruction if necessary
        ttype = plotinfo(k).transform;
        if ttype == 1, % curvelets
            rim = ifdct_wrapping(Cmod, 0);
        else
            rim = waverec2(Dmod{ttype}, S{ttype}, wtypes{ttype-1});
        end
        errim = rim - im;
    end
    
    valuetype = plotinfo(k).value;
    switch plotinfo(k).type,
        case 1,
        imagesc(funcimag(im, valuetype));
        title('Original image, physical space')
    case 2,
        imagesc(fdim(1:2),fdim(3:4),funcimag(fftshift(fft2(im)), valuetype));
        title('Original image, frequency space')
    case 3,
        imagesc(funcimag(rim, valuetype));
        title('Reconstructed image, physical space')
    case 4,
        imagesc(fdim(1:2),fdim(3:4),funcimag(fftshift(fft2(rim)), valuetype));
        title('Reconstructed image, frequency space')
    case 5,
        imagesc(funcimag(errim, valuetype));
        %set(gca, 'YDir', 'normal')
        title('Error, physical space')
    case 6,
        imagesc(fdim(1:2),fdim(3:4),funcimag(fftshift(fft2(errim)), valuetype));
        title('Error, frequency space')
    end
    axis equal tight
    colormap gray
    if plotinfo(k).colorbar,
        colorbar
    end
    
    hold on
    if physspace,
        switch plotinfo(k).pos,
            case 1,
                % do nothing
            case 2,  % large dots
                if plotinfo(k).transform == 1,
                    plotcurveletpos(Cmod, ax, 0, 6);
                else
                    plotwaveletpos(Dmod{ttype}, S{ttype}, ax, 6);
                end
            case 3,  % small dots
                if plotinfo(k).transform == 1,
                    plotcurveletpos(Cmod, ax, 0, 1);
                else
                    plotwaveletpos(Dmod{ttype}, S{ttype}, ax, 1);
                end
            case 4,  % scaled arrows
                plotcurveletpos(Cmod, ax, 1, 1);
            case 5,  % non-scaled arrows
                plotcurveletpos(Cmod, ax, 1, 0);
        end
    end
end




% select levels and dirs
function Cmod = select_crvltlevels(C, computeinfo)

Cmod = C;
for ll=1:length(Cmod)
    levok = ismember(ll,computeinfo.levels);
    ndir = length(Cmod{ll});
    if ndir ==1,
        if ~levok,
            Cmod{ll}{1} = zeros(size(Cmod{ll}{1}));
        end
    else
        dirok = ismember(1:ndir, computeinfo.dirs{ll});
        for dd=1:ndir,
            if ~(levok && dirok(dd))
                Cmod{ll}{dd} = zeros(size(Cmod{ll}{dd}));
            end
        end
    end
end


% threshold coefficients
function Cth = thrsh_crvltcoeffs(C, computeinfo)

val = computeinfo.thrshval;
if iscell(val)
    val = val{1};
end

switch computeinfo.thrshtype,
    case 1, % do nothing
        Cth = C;
    case 2,  % threshold value
        Cth = crvlt_thresh(C, val);
        %set(handles.edit_curthrsh, 'String', sprintf('# nnz: %d', crvlt_countnnz(Cth)));
    case 3,  % number of coeffs
        if val <= 0,
            Cth = zerodct(C);
        else
            [Cth, smcoef] = crvlt_keeplargest(C, round(val));
        end
        %set(handles.edit_curthrsh, 'String', sprintf('smallest coeff: %f', smcoef)); 
    case 4,  % get threshold from first image
        error('Should not occur here!')
end

% threshold wavelet coefficients
function Dth = thrsh_waveletcoeffs(D, computeinfo, ttype)

val = computeinfo.thrshval;
if iscell(val),       % multiple thresholds (one per transform type)
    val = val{ttype};
end

switch computeinfo.thrshtype,
    case 1, % do nothing
        Dth = D;
    case 2,  % threshold value
        Dth = wthresh(D, 'h', val);
    case 3,  % number of coeffs
        Ds = sort(abs(D),2,'descend');
        thrsh = Ds(round(val)+1);
        Dth = wthresh(D, 'h', thrsh);
    case 4,  % get threshold from first image
        error('Should not occur here!')
end



% evaluate function on image, according to choice
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


