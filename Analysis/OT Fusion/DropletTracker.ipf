#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <Image Common>
#include <ImageSlider>
#pragma ModuleName=DropletTracker

Menu "Functions"
	SubMenu "Droplets"
		"Droplet Video Tracker", DropletTrackerPanel()
	End
End



Function/Wave LoadVideo(string vidname)
	ControlInfo/W=DropletTrackerPanel popVideos
	String graph = S_value
	if(WinType(graph)==0)
	endif
	Wave v = ImageNameToWaveRef(graph, StringFromList(0,ImageNameList(graph, ";")))
	Killwindow/Z $graph
	KillWaves/Z v
	
	String newVideoName = UniqueName("video", 1, 0)

	DFREF tagDF = :Tag0
	if(datafolderRefStatus(tagDF)!=0)
		abort "Tag folder :Tag0 still exists. Cannot read FPS. Exiting."
	endif

	ImageLoad/T=tiff/O/LTMD/S=0/C=-1/LR3D/Q/N=$newVideoName vidname
	Wave video = $newVideoName
	if(!WaveExists(video))
		return $""
	endif
	String newGraphName = "Graph"+newVideoName
	NewImage/N=$newGraphName video
	ModifyGraph/W=$newGraphName height={Plan, 1, left, top}
	WMAppend3DImageSlider()
	ShowInfo

	//--Find FPS
	Variable CALIB_fps = 1
	Wave/T T_Tags = :Tag0:T_Tags
	if(WaveExists(T_Tags))
		Variable i
		for(i=0;i<dimsize(T_Tags,0);i+=1)
			if(Stringmatch(T_Tags[i], "ImageDescription:*"))
				String tagImageDescription = T_tags[i]
				Variable sstart = strsearch(tagImageDescription, "\nfps=", 0) + strlen("\nfps=")
				Variable send = strsearch(tagImageDescription, "\n", sstart) - 1
				CALIB_fps = str2num(tagImageDescription[sstart,send])
				print "FPS:"+num2str(CALIB_fps)
			endif
		endfor
	else
		print "Could not read FPS!"
	endif
	KillDataFolder/Z :Tag0

	//Variable CALIB_nmPerpixel = 89.324
	Variable CALIB_nmPerpixel =   89.47
	Note video, "CALIB_nmPerPixel:"+num2str(CALIB_nmPerpixel)
	Note video, "CALIB_fps:"+num2str(CALIB_fps)
	Note video, "Path:"+S_path
	Note video, "Filename:"+S_fileName
	
	
	Setscale/P x, 0, CALIB_nmPerpixel/1000, video
	Setscale/P y, 0, CALIB_nmPerpixel/1000, video
	Setscale/P z, 0, 1/CALIB_fps, video
	ControlUpdate/W=DropletTrackerPanel popVideos
	return video
End


//************************************
// Generate radial template from video at a given pixel (pixel value can be non-integer)
//************************************
Function/Wave GenerateRadialTemplate(Wave video, Variable fr, Variable px, Variable py, Variable rMax [,Variable fromAngle, Variable toAngle, Variable rPntsPerPixel, Variable removeAverage])
	Imagetransform/P=(fr) getplane video
	Wave m_imagePlane
	Redimension/D M_imageplane
	Setscale/P x, 0, 1, M_imageplane
	Setscale/P y, 0, 1, M_imageplane

	fromAngle = paramisdefault(fromAngle) ? 0 : fromAngle
	toAngle = paramisdefault(toAngle) ? 360 : toAngle
	rPntsPerPixel = paramisdefault(rPntsPerPixel ) ? 1 : rPntsPerPixel
	removeAverage = paramisdefault(removeAverage) ? 1 : removeAverage

	Variable method = 3
	
	Wave profile = CreateRadialProfile(m_imageplane, px, py, rMax, method, fromAngle=fromAngle, toAngle=toAngle, rPntsPerPixel=rPntsPerPixel)
	
	//Remove average
	if(removeAverage)
		Wavestats/Q profile
		profile -= v_avg
	endif
	
	Variable size = numpnts(profile)*2
	if(mod(size,2)==1) //ENSURE EVEN NUMBER OF ROWS (FOR FFT), BUT NO WELL-DEFINED CENTER POINT!
		size -= 1
	endif
	Variable mid = size/2
	
	Make/O/N=(size, size) Template
		
	CopyScales/P video, Template
	Setscale/P x, -indexToScale(Template,size/2/rPntsPerPixel,0), dimdelta(template,0)/rPntsPerPixel, template
	Setscale/P y, -indexToScale(Template,size/2/rPntsPerPixel,1), dimdelta(template,1)/rPntsPerPixel, template
	
	Template[][] = sqrt(x^2+y^2) < rightx(profile) ? profile(sqrt(x^2+y^2)) : 0


	if(fromAngle>toAngle)
		Template[][] = sqrt(x^2+y^2) >= 0 && (atan360(y,x) >= fromAngle || atan360(y,x) <= toAngle) ? Template : 0
	else
		Template[][] = sqrt(x^2+y^2) >= 0 && atan360(y,x) >= fromAngle && atan360(y,x) <= toAngle ? Template : 0
	endif
	
	return Template
End

Static Function atan360(Variable y, Variable x)
	if(y>=0)
		return atan2(y,x)*180/pi
	else
		return atan2(y,x)*180/pi+360
	endif
End

//************************************
// Calculate cross-correlation between an image and a template
//************************************
#define WINDOW
Function/Wave CalcCCF3(Wave imageFiltered, Wave template, Variable px, Variable py, Variable HalfBox)
	Duplicate/O/D/R=[px-Halfbox,px+HalfBox+1][py-Halfbox,py+HalfBox+1] imageFiltered, imageFilteredBox
#ifdef WINDOW	
	ImageWindow/O Hamming imageFilteredBox
#endif	
	FFT/DEST=ImageFilteredBoxFFT imageFilteredBox
	
	Duplicate/O imageFilteredBox, templateSameSize

	templateSameSize=0

	
	Setscale/P x, -(dimdelta(templateSameSize,0)*(dimsize(templateSameSize,0)/2)), dimdelta(templateSameSize,0), templateSameSize
	Setscale/P y, -(dimdelta(templateSameSize,0)*(dimsize(templateSameSize,0)/2)), dimdelta(templateSameSize,1), templateSameSize

	templatesamesize = x>=dimOffset(template,0)&&y>=dimOffset(template,1)&&x<=dimOffset(template,0)+dimDelta(template,0)*(DimSize(template,0)-1)&&y<=dimOffset(template,1)+dimDelta(template,1)*(DimSize(template,1)-1) ? template(x)(y) : 0
	
	//Flip template (we're doing a convolution instead of a CCF)
	Imagetransform/O flipcols templatesamesize
	Imagetransform/O fliprows templatesamesize

#ifdef WINDOW
	ImageWindow/O Hamming templatesamesize
#endif	
	imagetransform swap templateSameSize
	FFT/Dest=templateFFT templateSameSize
	
	Duplicate/O/C ImageFilteredBoxFFT, CCF
	MatrixOp/O CCF = ImageFilteredBoxFFT * templateFFT
	IFFT CCF
	CopyScales/P imagefilteredbox, CCF
	return CCF
End



//************************************
// Fourier filter image (blur and remove offset)
//************************************
Function/Wave FilterImage(Wave image)
	Variable filterMin = 0.0005//*89.324
	Variable filterMax = 0.5//0.2//*89.324

	Variable old_offsetX = dimoffset(image,0)
	Variable old_offsetY = dimoffset(image,1)
	Variable old_deltaX = dimdelta(image,0)
	Variable old_deltaY = dimdelta(image,1)
	Setscale/P x, 0, 1, image
	Setscale/P y, 0, 1, image
	
	//ImageWindow /O/p=1.03 kaiser image
	fft/PAD={ceil(dimsize(image,0)/2)*2, ceil(dimsize(image,1)/2)*2}/DEST=image_fft image
	image_fft[][] = (x^2+y^2) > filterMax ? 0 : image_fft
	image_fft[][] = (x^2+y^2) < filterMin ? 0 : image_fft
	ifft/DEST=image image_fft

	Setscale/P x, old_offsetX, old_deltaX, image
	Setscale/P y, old_offsetY, old_deltaY, image
	
	return image
End


//************************************
// Find template in image by cross correlation
//************************************
Function/Wave FindByTemplate2(Wave srcImage, Wave template, Variable px, Variable py, Variable maxBoxHalf)

	Wave image = Filterimage(srcImage)
	Wave CCF = calcccf3(image, template, px, py, maxboxhalf)
	
	//--Find Maximum pixel of CCF---
	Variable dummy = 10
	Wavestats/Q/P/RMD=[maxBoxHalf-dummy,maxBoxHalf+dummy][maxboxHalf-dummy,maxBoxHalf+dummy] CCF
	Variable MaxOK = 0
	Variable locx = NaN, locy = NaN
	if(abs(maxBoxHalf-v_maxrowloc) < maxBoxHalf && abs(maxBoxHalf-v_maxcolloc) < maxBoxHalf)
		px = v_maxrowloc
		py = v_maxcolloc
		locX = indextoscale(CCF, px, 0)
		locY = indextoscale(CCF, py, 1)
		MaxOK = 1
	endif
	
	//---Then do 2D Gauss Fit around Maximum for sub-pixel resultion---
	Variable FitBoxHalf = 3 //pixels

	Variable ccfAmp = NaN
	if(MaxOK)
		Variable v_fiterror = 0
		Variable width = 0.3
		Make/O w_coef={0,v_max,IndexToScale(CCF,px,0),width,IndexToScale(CCF,py,1),width,0}
		Make/O/T/FREE T_Constr={"K1>0"}
		CurveFit/Q/H="1000001" Gauss2D,kwcwave=w_coef CCF[px-FitBoxHalf,px+FitBoxHalf][py-FitBoxHalf,py+FitBoxHalf]/C=T_constr
		Wave w_coef, w_sigma
		if(v_fiterror == 0 && w_sigma[2] < 1 && w_sigma[4] < 1 && w_coef[2] >= indextoScale(CCF, px-FitBoxHalf, 0) && w_coef[2] <= indextoScale(CCF, px+FitBoxHalf, 0) && w_coef[4] >= indextoScale(CCF, py-FitBoxHalf, 1) && w_coef[4] <= indextoScale(CCF, py+FitBoxHalf, 1))
			locX = w_coef[2]
			locY = w_coef[4]
			ccfAmp = w_coef[1]
		endif
	endif
	
	Make/O/FREE LOC = {locX, locY, ccfAmp}
	return LOC
End



Function UpdateLoc()
	ControlInfo/W=DropletTrackerPanel popVideos
	String graph = S_value
	Wave video = ImageNameToWaveRef(graph, StringFromList(0,ImageNameList(graph, ";")))

	Wave template=template

	ControlInfo/W=$graph WM3DVal
	Variable fr = v_value

	Imagetransform/P=(fr) getplane video
	Wave m_imagePlane
	Redimension/D M_imageplane
	Duplicate/O m_imagePlane, image_filtered
	FilterImage(image_filtered)

		
	Wave LOC = FindByTemplate2(image_filtered,  template,  pcsr(A,"Graphvideo9"), qcsr(A,"Graphvideo9"),  40)
	print LOC
	
	if(!numtype(LOC[2]))
		Cursor/W=$graph/I A $nameofwave(video) LOC[0], LOC[1]
		Make/O/N=1 LOCX, LOCY
		LOCY = LOC[1]
		LOCX = LOC[0]
	else
		LOCY = NaN
		LOCX = NaN
	endif
		
End


Window DropletTrackerPanel() : Panel
	PauseUpdate; Silent 1		// building window...
	DoWindow/F DropletTrackerPanel
	if(v_flag)
		return NaN
	endif
	
	NewPanel/K=1 /N=DropletTrackerPanel/W=(255,236,555,400)
	SetDrawLayer UserBack
	GroupBox grpTemplateRadialProfile,pos={7.00,10.00},size={288.00,121.00},title="Template from radial profile"
	Button btnGenerateTemplate,pos={206.00,100.00},size={86.00,28.00},title="Generate",proc=DropletTracker#procGenerateRadialTemplate
	SetVariable setVarTemplateRadius,pos={17.00,30.00},size={100.00,14.00},title="Radius (pxl)"
	SetVariable setVarTemplateRadius,format="%d",limits={1,inf,1},value= _NUM:30
	SetVariable setVarPntsPerPixelr,pos={17.00,50},size={100.00,14.00},title="radial pnts/pxl"
	SetVariable setVarPntsPerPixelr,format="%d",limits={1,inf,1},value= _NUM:1,disable=2
	SetVariable setVarFromAngle,pos={136.00,30},size={120.00,14.00},disable=0,title=">=Angle (°)"
	SetVariable setVarFromAngle,format="%d",limits={0,360,1},value= _NUM:0
	SetVariable setVarToAngle,pos={136.00,50},size={120.00,14.00},disable=0,title="<=Angle (°)"
	SetVariable setVarToAngle,format="%d",limits={0,360,1},value= _NUM:360
	SetVariable setVarIterateTemplate,pos={17.00,70.00},size={100.00,14.00},title="#Iterations"
	SetVariable setVarIterateTemplate,help={"Re-localize with template and re-generate # times (0: don't re-iterate)"}
	SetVariable setVarIterateTemplate,format="%d",limits={0,inf,1},value= _NUM:10
	PopupMenu popVideos,pos={17.00,100.00},size={150.00,23.00},bodyWidth=90,title="Selection"
	PopupMenu popVideos,mode=1,value=WinList("Graphvideo*", ";", "WIN:1"),proc=DropletTracker#procSelectVideo
EndMacro




Static Function procLoadVideo(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up		
			ControlUpdate/W=DropletTrackerPanel popVideos
			Wave video = LoadVideo("")
			TitleBox tbPath win=DropletTrackerPanel, title=StringByKey("Path", note(video), ":", "\r")+"\r"+StringByKey("Filename", note(video), ":", "\r")
			ControlUpdate/W=DropletTrackerPanel popVideos
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Static Function procDeleteVideo(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			ControlInfo/W=DropletTrackerPanel popVideos
			String graph = S_value
			if(WinType(graph)==0)
				abort "Graph "+graph+" not found."
			endif
						
			Wave video = ImageNameToWaveRef(graph, StringFromList(0,ImageNameList(graph, ";")))
			Killwindow/Z $graph
			KillWaves/Z video
			
			ControlUpdate/W=DropletTrackerPanel popVideos
			ControlUpdate/W=DropletTrackerPanel popTracks
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Static Function procSelectVideo(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr

			if(WinType(popStr)==0)
				break
			endif
			
			DoWindow/F $popStr			
			Wave video = ImageNameToWaveRef(popStr, StringFromList(0,ImageNameList(popStr, ";")))
			TitleBox tbPath win=DropletTrackerPanel, title=StringByKey("Path", note(video), ":", "\r")+"\r"+StringByKey("Filename", note(video), ":", "\r")
			ControlUpdate/W=DropletTrackerPanel popTracks
			
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End



Static Function procGenerateRadialTemplate(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			ControlInfo/W=DropletTrackerPanel popVideos
			String graph = S_value
			if(WinType(graph)==0)
				abort "Graph "+graph+" not found."
			endif
			
			if(strlen(CsrInfo(A, graph))==0)
				abort "Cursor A not in graph "+graph
			endif
			Wave video = ImageNameToWaveRef(graph, StringFromList(0,ImageNameList(graph, ";")))
			
			ControlInfo/W=$graph WM3DVal
			Variable fr = v_value
			
			ControlInfo/W=DropletTrackerPanel setVarTemplateRadius
			Variable rMax = v_value
			
			ControlInfo/W=DropletTrackerPanel setVarPntsPerPixelr
			Variable rPntsPerPixel = v_value
			
			ControlInfo/W=DropletTrackerPanel setVarFromAngle
			Variable fromAngle = v_value
			
			ControlInfo/W=DropletTrackerPanel setVarToAngle
			Variable toAngle = v_value
			
			ControlInfo/W=DropletTrackerPanel setVarIterateTemplate
			Variable numIterateTemplate = v_value

			get_radius(video, fr, rmax, rpntsperpixel, fromangle, toangle, numiteratetemplate, graph, "A")		
		   get_radius(video, fr, rmax, rpntsperpixel, fromangle, toangle, numiteratetemplate, graph, "B")				
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End




// create a radial profile in an image
// simg - 2D image wave
// cx, cy - center for circle in image
// rmax - maximum radius to take profile
// creates profile waves: M_RadialProfile, M_RadialProfileWeights, M_NormRadialProfile
// creates variables: V_TotalArea, V_TotalIntensity
//
//* Method = 1 - Pixel by Pixel with Linear Interpolation
//* Method = 2 - Pixel by Pixel with Quadratic Interpolation
//* Method = 3 - Radial Wave Profile with ImageLineProfile
//https://www.wavemetrics.com/code-snippet/radial-profiler
//Modified JS 2020
//
//Everything is calculated in pixel coordinates of the original image
Function/Wave CreateRadialProfile(simg, cx, cy, rmax, method [,fromAngle, toAngle, rPntsPerPixel])
    wave simg
    variable cx, cy, rmax, method
    Variable fromAngle, toAngle
    Variable rPntsPerPixel
    
    rPntsPerPixel = paramisdefault(rPntsPerPixel ) ? 1 : rPntsPerPixel // to achieve sub-pixel profiles
    fromAngle = paramisdefault(fromAngle) ? 0 : fromAngle
    toAngle = paramisdefault(toAngle) ? 360 : toAngle
    Variable angleRange = abs(toAngle-fromAngle)
    
    if(method != 3)
    	print "Only Method 3 is robust"
    	return $""
    endif
    
    variable ix, iy, ixmax, iymax, rpos, rposf, rposc, pmax, wfc
    variable ic, dnT, npts, r360
   
    Make/D/O/N=(rmax*rPntsPerPixel+1) M_RadialProfile=0, M_RadialProfileWeights=0, M_RadialFluxProfile=0
    Setscale/P x, 0, 1/rPntsPerPixel, M_RadialFluxProfile, M_RadialProfile, M_RadialProfileWeights
    variable/G V_TotalArea, V_TotalIntensity
    V_TotalArea = 0; V_TotalIntensity = 0
    
    if (method == 3)
        make/FREE/D/N=(1) xWave = 0, yWave = 0
        r360 = ceil(360/(2*Pi))
        for (ic=1;ic<rmax*rPntsPerPixel+1;ic+=1)
            npts = round(angleRange*ic/r360) + 1
            redimension/N=(npts) xWave, yWave
            dnT= angleRange/npts
            xWave = cx + ic*cos((fromAngle + p*dnT)*Pi/180) /rPntsPerPixel 
            yWave = cy + ic*sin((fromAngle + p*dnT)*Pi/180) /rPntsPerPixel 
            ImageLineProfile xWave=xWave, yWave=yWave, srcwave=simg
            wave W_ImageLineProfile
            M_RadialProfile[ic] = sum(W_ImageLineProfile)
            M_RadialProfileWeights[ic] = numpnts(W_ImageLineProfile)
        endfor
        M_RadialProfileWeights[0] = 1
        
        M_RadialProfile[0] = simg[cx][cy]
        //Bilinear center:
        Variable fx = cx-floor(cx)
        Variable fy = cy-floor(cy)
        M_RadialProfile[0] = simg[floor(cx)][floor(cy)] * (1-fx) * (1-fy) + simg[floor(cx)+1][floor(cy)] * fx*(1-fy) + simg[floor(cx)][floor(cy)+1] * (1-fx)*fy + simg[floor(cx)+1][floor(cy)+1] * fx*fy
        
        M_RadialFluxProfile = M_RadialProfile/M_RadialProfileWeights
        V_TotalArea = sum(M_RadialProfileWeights)
        V_TotalIntensity = sum(M_RadialProfile)
        KillWaves/Z W_ImageLineProfile, W_LineProfileDisplacement, W_LineProfileX, W_LineProfileY
        return M_RadialFluxProfile
    else    
       ixmax = rmax
       // outer sum over x
        for (ix=1-ixmax; ix<ixmax; ix+=1)
            // bound y by current x
            iymax = floor(sqrt(rmax^2 - ix^2))
            // inner sum over y
            for (iy=1-iymax; iy<iymax; iy+=1)
                rpos = sqrt(ix^2 + iy^2)
                // split value between the two closest points in the radial profile wave
                rposf = floor(rpos)
                rposc = ceil(rpos)
                switch(method)
                    case 1:     // linear interpolation
                        wfc = rposc - rpos
                        break
                    case 2:
                        wfc = (rposc^2 - rpos^2)/(rposc^2 - rposf^2)
                        wfc = numtype(wfc) == 2 ? 1 : wfc
                        break
                endswitch
                // fill the weighting array
                M_RadialProfileWeights[rposf]+= (1 - wfc)
                M_RadialProfileWeights[rposc]+= wfc
                // fill the profile array
             M_RadialProfile[rposf] += simg[cx+ix][cy+iy]*(1 - wfc)
             M_RadialProfile[rposc] += simg[cx+ix][cy+iy]*wfc                
             // add to the total area and total intensity
             V_TotalArea += 1
             V_TotalIntensity += simg[cx+ix][cy+iy]     
            endfor
        endfor
        // generate normalize profile
        M_RadialFluxProfile = M_RadialProfile/M_RadialProfileWeights
    endif
    return M_RadialFluxProfile
end

function get_radius(wave video, variable fr, variable rMax, variable rPntsPerPixel, variable fromAngle, variable toAngle, variable numIterateTemplate, string graph, string csr)
	nvar s1, s2
	wave template = GenerateRadialTemplate(video, fr, pcsr($csr, graph), qcsr($csr, graph), rMax, rPntsPerPixel=rPntsPerPixel, fromAngle=fromAngle, toAngle=toAngle)
	Variable it
			for(it=0;it<numIterateTemplate;it+=1)
				Imagetransform/P=(fr) getplane video
				Wave m_imagePlane
				Redimension/D M_imageplane
				Duplicate/O/FREE M_imageplane, image
				Wave LOC = FindByTemplate2(image,  template,  pcsr($csr, graph), qcsr($csr, graph),  40)
				Wave template = GenerateRadialTemplate(video, fr, (LOC[0] - DimOffset(video,0)) / DimDelta(video,0), (LOC[1] - DimOffset(video,1)) / DimDelta(video,1), rMax, rPntsPerPixel=rPntsPerPixel, fromAngle=fromAngle, toAngle=toAngle) //don't use scaletoindex (it rounds)
				Make/O/N=1 LOCX, LOCY
				LOCY = LOC[1]
				LOCX = LOC[0]
			endfor
			
			if(winType("WinTemplate")==0)
				NewImage/N=WinTemplate template
				ModifyGraph/W=WinTemplate height={Plan, 1, left, top}
			endif
				Wavestats/Q template
				Variable range = max(v_max,-v_min)
				ModifyImage/W=WinTemplate Template ctab={-range,range,BlueBlackRed,0}
	
	
			//Find bead radius (right now: right between white ring and zero transition)
			Wave profile = m_radialfluxprofile
			setscale/P x, 0, dimdelta(video,0)/rPntsperPixel, profile
		 	wavestats/Q/R=(0.5,) profile
		 	Variable maxProfile = v_maxloc
		 	CurveFit/Q poly 3, profile(v_maxloc-0.5, v_maxloc+0.5)
		 	Wave w_coef
		 	maxProfile = -w_coef[1]/(2*w_coef[2])
		 	FindLevel/Q/R=(maxProfile,0)/EDGE=1 profile, 0
		 	Variable zeroCrossProfile = v_levelx
		 	Variable radius= (maxProfile+zerocrossprofile)/2
		 	print "Cursor"+ csr + " Radius:",radius,"µm"
 			if(stringmatch(csr, "A") == 1)
 				s1 = radius
 			else
 				s2 = radius
 			endif
 			if(numIterateTemplate>0)
	 			DrawAction/W=$graph getgroup=$("circle"+csr), delete, begininsert
	 			SetDrawEnv/W=$graph gname=$("circle"+csr),gstart
	 			SetDrawEnv/W=$graph xcoord= top,ycoord= left,linefgc= (65535,32768,58981),fillpat= 0
	 			DrawOval/W=$graph LOC[0]-radius,LOC[1]-radius,LOC[0]+radius,LOC[1]+radius
	 			SetDrawEnv/W=$graph gstop
	 			DrawAction/W=$graph endinsert
			endif		
end


