function results_phasor = phasor_code(FileName, PathName)

close all
% -- set constants
thresh_low  = 0.02;         % Lower threshold [0,1]
thresh_high = 0.44;         % High threshold [0,1]
harmonic    = 1;            % harmonic number (1,2,3,4,...)
channel     = 1;            % dectector channel
shift_bins  = 0;            % number of bins shifted to right from max for minimizing effect of IRF
scale_bar   = 10;           % in microns

% ------ single calculation or batch processing ------ %

if nargin<2
    [FileName, PathName] = uigetfile({'*.mat'}, 'Select FLIM file');
    datadir = 1 ;     % dummy for for-loop
else
    folder_FLIM = [PathName '\*.mat'];
    datadir = dir(folder_FLIM);
end

% ------ begin for-loop for batch data processing ------ %

h = waitbar(0, 'Phasoring...','WindowStyle', 'normal');

for jj = 1:numel(datadir)
    
    if nargin<2
        decayfile = FileName;
    else
        decayfile = datadir(jj).name;
    end
    
    % ------ load data ------ %
    
    fidd = [PathName '\' decayfile];
    res = load(fidd);
    
    % ------ extract metadata ------ %
    
    freq0 = 1/res.head.MeasDesc_GlobalResolution;     % laser frequency
    delta_t = res.head.MeasDesc_Resolution;           % width of one time channel
    time_bins = 1 / (freq0 * delta_t) / harmonic;
    time_bins = round(time_bins);
    % number of time bins for one period (not 64!). 80 for harmonic no 1
    
    % ------ initial calculations ------ %
    
    freq = harmonic * freq0;        % virtual frequency for higher harmonics
    w = 2 * pi * freq ;             % angular frequency
    
    % Find Max and remove data before max
    [~,I] = max(res.decay(:,channel));
%     decdata = res.img(:,:,I+shift_bins:end-2,channel);
    decdata = cat(3,res.img(:,:,I+shift_bins:end-2,channel), res.img(:,:,1:round(2*I/3),channel));
    
    timechannels_data = size(decdata,3);
    decdatatooshort = 0;
    
    if timechannels_data < time_bins
        time_bins = timechannels_data;
        %         decdata(1,1,time_bins) = 0;
        %         decdatatooshort = 1;
        %     elseif  timechannels_data > time_bins
        %         decdata = decdata(:,:,1:time_bins);
    end
    
    %decdata = int32(medfilt1(double(decdata)));        %smoothing
    
    % Threshold (not peak vs peak. Used whole intensity vs max whole intensity)
    maxmax = max(max(sum(decdata,3)));
    update = false;
    
    pixels_low  = any(sum(decdata,3) < (thresh_low * maxmax), 3);
    pixels_high = any(sum(decdata,3) > (thresh_high * maxmax), 3);
    thresh_img = sum(decdata.*~pixels_low.*~pixels_high,3);
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0.1, 0.4, 0.74]);
    imagesc(thresh_img)
    a = gray;
    a(1,:) = [ 1 0.3 0.3];
    colormap(a)
    colorbar
    hold on
    
    %
    while (~update)
        
        %   numberofpixels = numel(sum(res.img(:,:,:,channel),3)) - sum(sum(pixels_low+pixels_high));
        prompt = {'Threshold Low','Threshold High'};
        dlgtitle = 'Thresholding';
        dims = [1 50];
        definput = {num2str(thresh_low), num2str(thresh_high)};
        options.WindowStyle = 'normal';
        answer = inputdlg(prompt,dlgtitle,dims,definput,options);
        if str2double(answer{1}) ~= thresh_low || str2double(answer{2}) ~= thresh_high
            thresh_low = str2double(answer{1});
            thresh_high = str2double(answer{2});
            pixels_low  = any(sum(decdata,3) < (thresh_low * maxmax), 3);
            pixels_high = any(sum(decdata,3) > (thresh_high * maxmax), 3);
            thresh_img = sum(decdata,3).*~pixels_low.*~pixels_high;
            drawnow;
            imagesc(thresh_img, [0 max(max(thresh_img))])
            colorbar
            numberofpixels = numel(sum(res.img(:,:,:,channel),3)) - sum(sum(pixels_low+pixels_high));
            hold on
        elseif thresh_low == 0
            pixels_low  = any(sum(decdata,3) < 1, 3);
            pixels_high = any(sum(decdata,3) > (thresh_high * maxmax), 3);
            thresh_img = sum(decdata,3).*~pixels_low.*~pixels_high;
            update = true;
            decdata = decdata.*~pixels_low.*~pixels_high;
            numberofpixels = numel(sum(res.img(:,:,:,channel),3)) - sum(sum(pixels_low+pixels_high));
            close
        else
            update = true;
            decdata = decdata.*~pixels_low.*~pixels_high;
            numberofpixels = numel(sum(res.img(:,:,:,channel),3)) - sum(sum(pixels_low+pixels_high));
            close
        end
    end
    %}

%%

pixels_low  = any(sum(decdata,3) < (thresh_low * maxmax), 3);
pixels_high = any(sum(decdata,3) > (thresh_high * maxmax), 3);
decdata = decdata.*~pixels_low.*~pixels_high;
numberofpixels = numel(sum(res.img(:,:,:,channel),3)) - sum(sum(pixels_low+pixels_high));

%%
    tic
    
    % remove offset from data
    decdata = reshape(decdata, [], size(decdata,3));
%     data_off = round(mean(res.img(:,:,round(I/3):round(2*I/3),channel),3));
%     data_off = reshape(data_off, [],1);
%     data_off = round(mean(decdata(:,end-2:end),2));
    data_off = min(decdata,[],2);
    data_off = repmat(data_off,[1 size(decdata,2)]);
    decdata = decdata - data_off ;
    decdata(decdata<0) = 0 ;    
    
    % calculate G/S-sin/cos-matrix
    tb_vec = linspace(1,time_bins,time_bins)'; % time bin indexing vector
    
    Gn_ma = cos(w * delta_t * (tb_vec - 0.5));
    Sn_ma = sin(w * delta_t * (tb_vec - 0.5));
    
    % calculate data phasor
    Gn = double(decdata) * double(Gn_ma) ;  %take decdata from preprocessing
    Sn = double(decdata) * double(Sn_ma) ;
    area = sum(decdata, 2) ;
    
    % normalization
    G_f = Gn ./ area;
    S_f = Sn ./ area;
    
    G_f_img = reshape(G_f,size(res.img,1),size(res.img,2));
    S_f_img = reshape(S_f,size(res.img,1),size(res.img,2));
    G_f_img(isnan(G_f_img))= 0;
    S_f_img(isnan(S_f_img))= 0;
    
    G_f(isnan(G_f))=[];
    S_f(isnan(S_f))=[];
    
    % data
    Z = [G_f,S_f];
    
    % ------ calculate lifetimes ------ %
    
    tau_phi = mean(1 / w .* (Z(:,2) ./ Z(:,1)) .* 1e9);
    tau_m = mean(1 / w .* sqrt(1./(Z(:,2).^2+Z(:,1).^2)-1) .* 1e9);
    tau_phi_std = std(1 / w .* (Z(:,2) ./ Z(:,1)) .* 1e9);
    tau_m_std = std(1 / w .* sqrt(1./(Z(:,2).^2+Z(:,1).^2)-1) .* 1e9);
    tau_mean = mean([1 / w .* (Z(:,2) ./ Z(:,1)) .* 1e9 1 / w .* sqrt(1./(Z(:,2).^2+Z(:,1).^2)-1) .* 1e9],'all');
    tau_mean_std = std([1 / w .* (Z(:,2) ./ Z(:,1)) .* 1e9 1 / w .* sqrt(1./(Z(:,2).^2+Z(:,1).^2)-1) .* 1e9],0,'all');
    tau_mean_img = mean(cat(3,1 / w .* (S_f_img ./ G_f_img) .* 1e9, 1 / w .* sqrt(1./(S_f_img.^2+G_f_img.^2)-1) .* 1e9),3);
    
    % ------ Binning for 2D histogram ------ %
    
    % bin centers
    steps = 0.0025;
    G_f_bins = 0:steps:1;
    S_f_bins = 0:steps:1;
    G_f_NumBins = numel(G_f_bins);
    S_f_NumBins = numel(S_f_bins);
    
    % map X/Y values to bin indices
    G_f_i = round( interp1(G_f_bins, 1:G_f_NumBins, G_f, 'linear', 'extrap') );
    S_f_i = round( interp1(S_f_bins, 1:S_f_NumBins, S_f, 'linear', 'extrap') );
    
    % limit indices to the range [1,numBins]
    G_f_i = max( min(G_f_i,G_f_NumBins), 1);
    S_f_i = max( min(S_f_i,S_f_NumBins), 1);
    
    % count number of elements in each bin
    H = accumarray([S_f_i(:) G_f_i(:)], 1, [S_f_NumBins G_f_NumBins]);
    H = medfilt2(H);
  
    % plot 2D histogram (w	 contour)
    figure;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.45, 0.6, 0.5]);
    subplot(1,2,1)
    contourf(G_f_bins, S_f_bins, H);
    imagesc(G_f_bins, S_f_bins, H)
    j = jet;
    j(1,:) = [ 0 0 0 ];     %comment that out for blue background in plot
    colormap(j);
    set(gca,'Ydir','Normal')
    text(0.01,0.57,['\bf Phase lifetime: ' num2str(tau_phi,'%.2f') ' +/- ' num2str(tau_phi_std,'%.2f') 'ns'], 'Color', 'w', 'FontSize', 12)
    hold on
%     plot(G_f, S_f, 'b.', 'MarkerSize',2);  % scatter plot of all phasor points
    
    % ------ uni. circle ------ %
    
    % Plotuniversal circle
    x = 0:0.005:1;
    circle = sqrt(0.25 - (x - 0.5) .^ 2);
    plot(x,circle,'Color',[0.8 0.8 0.8],'LineWidth',0.5);
    axis equal tight
    axis([0 1 0 0.6]);
    xlabel('G');
    ylabel('S');
    title([decayfile(1:end-8) ' - Phasor Plot / Thr- :' num2str(thresh_low) ' Thr+ :' num2str(thresh_high)]);
    
    % ------ Plot FLIM image ------ %
    x1 = axes;
    subplot(1,2,2,x1)
    x2 = axes;
    subplot(1,2,2,x2)
    thresh_n = max(max(thresh_img));
    imagesc(x1,tau_mean_img)
    axis off
    thresh_img(floor(end-(size(res.img,1)/20)):floor(end-(size(res.img,1)/30))...
        ,floor(size(res.img,2)/30):floor(size(res.img,2)/20)-1+ceil(scale_bar/res.head.ReqHdr_SpatialResolution)) = thresh_n;
    imagesc(x2,thresh_img);
%    caxis(x1,[tau_phi-tau_phi_std tau_phi+tau_phi_std])
    caxis(x1,[3 3.4]);
    al = 1-thresh_img/thresh_n;
    al(floor(end-(size(res.img,1)/20)):floor(end-(size(res.img,1)/30))...
        ,floor(size(res.img,2)/30):floor(size(res.img,2)/20)-1+ceil(scale_bar/res.head.ReqHdr_SpatialResolution)) = 1;
    alpha(x2, al)
    axis off
    colormap(x1,'parula')
    colormap(x2,'gray')
    %set(x1,'Position',[.555 .11 .3 .815]);
    %set(x2,'Position',[.555 .11 .3 .815]);
    set([x1,x2], 'visible', 'off')
    cb1 = colorbar(x1,'Position',[.55 .11 .015 .815]);
    cb2 = colorbar(x2,'Position',[.91 .11 .015 .815]);
    title([decayfile(1:end-8) ' - FLIM image / scale bar ' num2str(scale_bar) '\mum / Thr- :' num2str(thresh_low) ' Thr+ :' num2str(thresh_high)]);
    
    % ------ save temporarily variables to array ------ %
    
    kk = jj;
    decayfile_exp{kk}        = decayfile;
    tau_phi_exp(kk)          = tau_phi;
    tau_m_exp(kk)            = tau_m;
    tau_phi_std_exp(kk)      = tau_phi_std;
    tau_m_std_exp(kk)        = tau_m_std;
    tau_mean_exp(kk)         = tau_mean;
    tau_mean_std_exp(kk)     = tau_mean_std;
    thresh_low_exp(kk)       = thresh_low;
    thresh_high_exp(kk)      = thresh_high;
    harmonickk(kk)           = harmonic;
    tc_data_exp(kk)          = timechannels_data;
    time_bins_needed(kk)     = time_bins;
    if decdatatooshort == 1
        toolesstc{kk} = {'Too less time channels. Try again with higher harmonic'};
    else
        toolesstc{kk} = {'fine'};
    end
    exec_time1_exp(kk)       = toc;
    pixelnumber(kk)          = numberofpixels;
    
    waitbar(jj/(numel(datadir)),h)
    
    %     counter = [num2str(jj),' of ',num2str(numel(datadir)),' done'];
    %     disp(counter)
    %
    %     switch lower(singleorbatch)
    %         case {'batch'}
    %             close all
    %     end
    
    disp('Press a key !')  % Press a key here.You can see the message 'Paused: Press any key' in
    % the lower left corner of MATLAB window.
    pause;
    
    print([PathName '\' decayfile(1:end-8) '_P+F.png'],'-dpng','-r300')
    
    close
end

close(h)

h = waitbar(1,'Exporting Results...','WindowStyle', 'normal');

% ------ collect and export results ------ %

results_phasor = table(decayfile_exp', tau_phi_exp', tau_m_exp', tau_mean_exp', ...
    tau_phi_std_exp', tau_m_std_exp',tau_mean_std_exp', thresh_low_exp', thresh_high_exp', harmonickk',...
    exec_time1_exp', pixelnumber', tc_data_exp', time_bins_needed', toolesstc',...
    'VariableNames', {'filename' 'tau_phi' 'tau_m' 'tau_mean' 'tau_phi_std' ...
    'tau_m_std' 'tau_mean_std' 'thres_low' 'thres_high' 'harmonic' 'exec_time1' 'pixelnumber' ...
    'time_channels_data' 'time_bins_needed' 'no_of_tc'})  %#ok<NOPRT>


%------ write output ------ %

if nargin<2
    Filenameout_sug = [PathName '\' FileName(1:end-8) '_phasorresult.xlsx'];
else
    Filenameout_sug = [PathName '\Phasor_Result.xlsx'];
end

[FileNameout, PathNameout] = uiputfile(Filenameout_sug, 'Save phasor results as...' );
outputfile = [PathNameout '\' FileNameout];
writetable(results_phasor,outputfile);

close(h)

end
