
function res = PTU_to_Prehist(name, P, C)

photons = P; % Photons batch size for processing.
Cores = C; % Cores in you PC processor

file = name(end-2:end);

if strcmp(file,'ptu')
    head=PTU_Read_Head(name);
    Timeunit=1/(head.TTResult_SyncRate);
    Resolution=head.MeasDesc_Resolution;
    NCounts=head.TTResult_NumberOfRecords;
else
    disp('Not a valid file!')
    return
end

dnum = head.HW_InpChannels;
NChannels=ceil(Timeunit/Resolution);

bin = 0:NChannels-1;
time = single(bin).*head.MeasDesc_Resolution.*1e9;
Bunches = ceil(NCounts/photons);
photons = ceil(NCounts/Bunches);

if Bunches < Cores
    Cores = Bunches;
end

im_tcspc = [];
im_chan  = [];
im_line  = [];
im_pixel = [];
im_frame = [];

i = 0;

%Read the image

h = waitbar(0, 'Reading file ...','WindowStyle','normal');

while (i <  Bunches)
    
    spmd(Cores)
        
        switch labindex
            
            case 1
                cnt = (labindex+i-1)*photons;
                [~, ~, tcspc, chan, line, pixel, frame] = PTU_LineScanRead_Phasor(name, [1+cnt,photons]);
                
            case 2
                cnt = (labindex+i-1)*photons;
                [~, ~, tcspc, chan, line, pixel, frame] = PTU_LineScanRead_Phasor(name, [1+cnt, photons]);
                
            case 3
                cnt = (labindex+i-1)*photons;
                [~, ~, tcspc, chan, line, pixel, frame] = PTU_LineScanRead_Phasor(name, [1+cnt,photons]); 
                
            case 4
                cnt = (labindex+i-1)*photons;
                [~, ~, tcspc, chan, line, pixel, frame] = PTU_LineScanRead_Phasor(name, [1+cnt,photons]); 
                
            case 5
                cnt = (labindex+i-1)*photons;
                [~, ~, tcspc, chan, line, pixel, frame] = PTU_LineScanRead_Phasor(name, [1+cnt,photons]); 
                
            otherwise 
                cnt = (labindex+i-1)*photons;
                [~, ~, tcspc, chan, line, pixel, frame] = PTU_LineScanRead_Phasor(name, [1+cnt,photons]);
                
        end
    end
    
    for j = 1:Cores                                
            im_tcspc = [im_tcspc ; tcspc{j}];
            im_chan  = [im_chan ; chan{j}];
            im_line  = [im_line ; line{j}];
            im_pixel = [im_pixel ; pixel{j}];
            im_frame = [im_frame ; frame{j}];
    end  
    
    i = i + Cores;
    waitbar(i/Bunches, h);
    
    if i+Cores > Bunches
        Cores = Bunches - i;
    end 
    
    clear tcspc
    clear chan
    clear line
    clear pixel
    clear frame
    
    delete(gcp('nocreate'));
    
end 

im_frame = [1 ; cumsum(diff(im_frame))+1];
frames = max(im_frame);
img = zeros(head.ImgHdr_PixY,head.ImgHdr_PixX,NChannels,dnum, 'single');
decay = zeros(NChannels,dnum);

close(h)

p = 1;
D = parallel.pool.DataQueue;
h = waitbar(0, 'Reconstructing image ...','WindowStyle','normal');
afterEach(D, @nUpdateWaitbar);
lines_cor = [];


for i = 1:frames
    idx = (im_frame == i);
    line1 =unique(im_line(idx));
    lines = [line1(1) ; cumsum(diff(im_line(idx)))+line1(1)];
    lines_cor = [lines_cor ; lines]; %#ok<AGROW>
end

waitbar(0.2,h, 'Reconstructing image ...','WindowStyle','normal');

im_line = lines_cor ; clear lines_cor ; clear lines

    parfor t = 1:NChannels
        for c = 1:dnum
            ind = (im_tcspc == t & im_chan == c);
            decay(t,c) = sum(ind);
            
            y = double(im_line(ind));
            x = im_pixel(ind);
            
            cbin = histogram2(y,x,1:head.ImgHdr_PixY+1,1:head.ImgHdr_PixX+1);
            binC = cbin.BinCounts;
            
            img(:,:,t,c) = img(:,:,t,c) + binC;
            send (D,t);
            
        end
    end

    delete(gcp('nocreate'));
    close(h)
    im_tcspc = []; im_chan = []; im_line = []; im_pixel = [];

figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.45, 0.6, 0.5]);
subplot(1,2,1)
semilogy(time, sum(decay,2))
subplot(1,2,2)
imagesc(sum(sum(img,3),4))

res.img   = img;
res.decay = decay;
res.time  = time;
res.head  = head;

save([name(1:end-4) ' bin.mat'],'img','decay','time','head','-v7.3');


function nUpdateWaitbar(~)
waitbar(0.2+(p/NChannels*0.8), h);
p = p + 1;
end


end