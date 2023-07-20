function [head, im_sync, im_tcspc, im_chan, im_line, im_col, im_frame, num] = PTU_LineScanRead_Phasor(name, cnts)
% reads the photons in a ptu file
% input:
% name:     name of the ptu file
% cnts:     photons to read (cnts(1): photon number to start from; cnts(2): number of photons to read)
%
% output:
% head:     header of the ptu file
% im_sync:  sync number of photon (global time)
% etc for the other outputs

if strcmp(name(end-2:end),'pt3') || strcmp(name(end-2:end),'ptu') %confirming that it is a .ptu file
    if strcmp(name(end-2:end),'pt3')
        head = PT3_read(name); % Read PT3 Files
        anzch      = 32;
        Resolution = max(head.Resolution);
        chDiv      = Resolution/head.Resolution;
        Ngate      = ceil(head.Sync./Resolution)+1;
        
        LineStart = 1;
        LineStop  = 2;
        
        fprintf('\n\n');
        
        [im_sync, im_tcspc, im_chan, markers, num] = PT3_read(name, cnts, head);
        
        ind = (markers>0)|((im_chan<anzch)&(im_tcspc<Ngate*chDiv));
    
    else 
        head = PTU_Read_Head(name);
        anzch      = 32;
        Resolution = max(head.MeasDesc_Resolution);
        chDiv      = Resolution/head.MeasDesc_Resolution;
        Ngate      = ceil(head.MeasDesc_GlobalResolution./Resolution)+1;
        
        LineStart = 1;
        LineStop  = 2;
        FrameMark = 4;
        
        fprintf('\n\n');
        
        [im_sync, im_tcspc, im_chan, markers, num] = PTU_Read(name, cnts, head);
        
        ind = (markers>0)|((im_chan<anzch)&(im_tcspc<Ngate*chDiv));
    
    end
    
    nx = head.ImgHdr_PixX;
    
    %        MultiFrame Scan or Line-scanning for Leica SP8
    %                 im_sync    = [im_sync ; y(ind)];
    %                 im_tcspc   = uint16([im_tcspc ; floor(tmpx(ind)./chDiv)+1]);
    %                 im_chan    = uint8([im_chan ; chan(ind)+1]);
    %                 markers    = [markers ; mark(ind)];
    
    im_sync    = im_sync(ind);
    im_tcspc   = uint16(floor(im_tcspc(ind)./chDiv)+1);
    im_chan    = uint8(im_chan(ind));
    markers    = markers(ind);

    
    Turns1d = im_sync.*(markers==LineStart);
    Turns2d = flipud(im_sync.*(markers==LineStop));
    Turns1d = interp1(1:nnz(Turns1d), Turns1d(Turns1d ~= 0), cumsum(Turns1d ~= 0));
    Turns2d = flipud(interp1(1:nnz(Turns2d), Turns2d(Turns2d ~= 0), cumsum(Turns2d ~= 0)));
    Norm = Turns2d-Turns1d;
    im_col = uint16(ceil(nx*(im_sync-Turns1d)./Norm));
    im_col(isnan(im_col)) = 0;
    indmark = interp1(1:nnz(markers), markers(markers ~= 0), cumsum(markers ~= 0));
    
    ind = markers~=0 | indmark==2 | im_col==0 | im_col == nx+1;
    
    im_line = uint16(ceil(cumsum(markers==1)));
    im_frame = uint16(ceil(cumsum(markers==FrameMark)));
    im_frame(ind)   = [];
    im_line(ind)    = [];
    im_col(ind)     = [];
    im_chan(ind)    = [];
    im_tcspc(ind)   = [];
    im_sync(ind)    = [];
    Norm(ind)       = [];
    
%     [Au,~,ic] = unique(im_sync,'stable');
%     fr = accumarray(ic, 1);
%     ind = logical(im_sync .* ~ismember(im_sync, Au(fr > 1)));
%     im_frame(~ind)   = [];
%     im_line(~ind)    = [];
%     im_col(~ind)     = [];
%     im_chan(~ind)    = [];
%     im_tcspc(~ind)   = [];
%     im_sync(~ind)    = [];
%     Norm(~ind)       = [];
    
    dt       = mean(Norm);
    %         im_frame = markers==uint8(Frame);
    
    im_line = mod(im_line,nx);
    im_line(im_line==0) = nx;
    %                 if sum(im_frame) == 0
    %                     im_frame = floor(im_line./(nx*2))+1;
    %                     im_line = mod(im_line,nx*2);
    %                     im_line(im_line==0) = nx*2;
    %                 else
    %                     im_frame      = single(ceil(cumsum(markers==Frame)+1));
    %                     im_frame(ind) = [];
    %                 end
    
    if strcmp(name(end-2:end),'pt3')
        head.ImgHdr_PixelTime = dt/nx/head.CntRate0;
        head.ImgHdr_DwellTime = head.ImgHdr_PixelTime;
    else
        head.ImgHdr_PixelTime = dt/nx/head.TTResult_SyncRate;
        head.ImgHdr_DwellTime = head.ImgHdr_PixelTime;
    end
else
    error('not a ptu file')
end
end

function [tag, tau, tcspc_pix, time] = Process_Frame(im_sync,im_col,im_line,im_chan,im_tcspc,head) %#ok<DEFNU>

Resolution = max([head.MeasDesc_Resolution*1e9 0.256]); % resolution of 1 ns to calculate average lifetimes
chDiv      = ceil(1e-9*Resolution./head.MeasDesc_Resolution);
SyncRate   = 1./head.MeasDesc_GlobalResolution;
nx = head.ImgHdr_PixX;
ny = head.ImgHdr_PixY;
dind    = double(unique(im_chan, 'legacy'));
Ngate   = round(head.MeasDesc_GlobalResolution./head.MeasDesc_Resolution*(head.MeasDesc_Resolution/Resolution)*1e9);
maxch_n = numel(dind);

tcspc_pix = zeros(nx,ny,Ngate,maxch_n);
% time = zeros(numel(im_sync),maxch_n);
time = {};
tag = zeros(nx,ny,maxch_n);
tau = tag;
bin = permute(repmat((1:Ngate)',[1 nx,ny]),[2,3,1])*Resolution; % 3D time axis
for ch = 1:maxch_n
    ind = im_chan==dind(ch);
    tcspc_pix(:,:,:,ch) =  mHist3(double(im_line(ind)),double(im_col(ind)),double(im_tcspc(ind)./chDiv),1:nx,1:ny,1:Ngate); % tcspc histograms for all the pixels at once!
    time{ch} = round(im_sync(ind)*1/SyncRate/Resolution/1e-9)+double(im_tcspc(ind)); %#ok<AGROW> % in tcspc bins
    tag(:,:,ch) = sum(tcspc_pix(:,:,:,ch),3);
    tau(:,:,ch) = real(sqrt((sum(bin.^2.*tcspc_pix(:,:,:,ch),3)./tag(:,:,ch))-(sum(bin.*tcspc_pix(:,:,:,ch),3)./tag(:,:,ch)).^2));
end
end
