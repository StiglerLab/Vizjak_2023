#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Menu "Functions"
	SubMenu "Droplets"
		"Droplet Fusion Fitter", droplet_window()
	End
End


function load_data(file_path)
	// Load selected file and get differential force + timestamps
	string file_path
	variable refNum
	KillWindow/Z FusFit#forceVtime
	HDF5OpenFile/R refNum as file_path
	if(v_flag == 0)
		print "Loaded " + file_path + " successfully."
	endif 
	HDF5ListGroup/R/Type=1 refNum, "/"
	string groups
	groups = S_HDF5ListGroup
	variable groupID
	HDF5OpenGroup refNum, "/", groupID
	// Load force data
	HDF5LoadData/Q/N=f1x/O refNum, "Force HF/Force 1x"
	wave f1x
	HDF5LoadData/Q/N=f2x/O refNum, "Force HF/Force 2x"
	wave f2x
	make/O/N=(numpnts(f1x)) force
	force = f1x-f2x  // Differential force
	// Get sample rate
	HDF5Loaddata/Q/A="Sample rate (Hz)"/N=sampleRate/O refNum, "Force HF/Force 1x"
	wave sampleRate
	variable sample_rate = sampleRate[0]
	Setscale/P x, 0, 1/sample_rate, "s", f1x, f2x, force
	
	Resample/DOWN=100 f1x, f2x, force
end


function normalise_force(force)
	// Normalise the loaded differential force
	wave force
	nvar finish, start
	variable posA, posB
	try
		force -= start
		force /= finish
	catch
		Variable err = GetRTError(1)
		print "Set Cursors"
	endtry
end


Function GeneralLogistic(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = finish+(start-finish)/(1+exp((x-x0)/tau))^gamma
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 5
	//CurveFitDialog/ w[0] = finish
	//CurveFitDialog/ w[1] = start
	//CurveFitDialog/ w[2] = x0
	//CurveFitDialog/ w[3] = tau
	//CurveFitDialog/ w[4] = gamma

	return w[0]+(w[1]-w[0])/(1+exp((x-w[2])/w[3]))^w[4]
End


function fit_logistic(force)
	Wave force
	
	FindLevel/EDGE=1/Q/R=(hcsr(A),hcsr(B)) force, (wavemax(force,xcsr(A),xcsr(B))+wavemin(force,xcsr(A),xcsr(B)))/2
	Variable x0start = v_levelx

	Make/O/D/N=5 w_coef
	w_coef[0] = {wavemax(force,xcsr(A),xcsr(B)),wavemin(force,xcsr(A),xcsr(B)),v_levelx,1,1}

	
	Make/O/FREE/T T_constr={"K4>0"}
	FuncFit/Q/TBOX=768 GeneralLogistic W_coef force[pcsr(A),pcsr(B)]/C=T_constr /D
	
	Variable/g finish = w_coef[0] //final level
	print(finish)
	Variable/g start = w_coef[1] //inital level
	Variable x0 = w_coef[2] //"mid" parameter
	Variable tau=w_coef[3] //"slope" parameter
	Variable gam = w_coef[4] //"shape" parameter
	
	Variable inflPoint = x0-tau*ln(gam) //Inflection point X
	print "Inflection point (s)", inflPoint
	
	Variable inflY = finish + (-finish + start)/(1 + 1/gam)^gam //Inflection point Y
	print "inflection level", inflY
	
	Variable vel = (gam*(finish - start))/((1 + 1/gam)^gam*((1 + gam)*tau)) //"Velocity" (dy/dx) at inflection point. Watch out for units!
	print "Velocity (dy/dx)", vel
	variable/G vel_norm = vel/(finish-start)
	print "Normalized velocity ((dy/dx)/(final-initial)", vel/(finish-start)
	
	Duplicate/O force, vel_extrapolation
	vel_extrapolation = vel*(x-inflPoint)+ inflY
	vel_extrapolation = vel_extrapolation>start && vel_extrapolation<finish ? vel_extrapolation : NaN

	//--Set zero-point to where velocity_extrapolation crosses the initial level
	Wave fit = $("fit_"+nameofwave(force))
	Wave w_coef, w_fitConstants
	Variable root = (start-inflY)/vel+inflPoint
	Setscale/P x, leftx(force)-root, deltax(force), "s", force
	Setscale/P x, leftx(fit)-root, deltax(fit), "s", fit
	Setscale/P x, leftx(vel_extrapolation)-root, deltax(vel_extrapolation), "s", vel_extrapolation
	Appendtograph /C=(0,0,0) vel_extrapolation
	return vel
end

function Select_file()
	// Select HDF5 file
	variable refnum
	string fileFilters = "HDF5 files (*.h5):.h5;"
	fileFilters += "All Files: .*;"
	open/D/R/F=fileFilters refnum
	string/G file_path
	file_path = S_fileName
	variable tmp_count = ItemsInList(file_path, ":")
	string/g vid_path = StringfromList(tmp_count-2, file_path, ":")
	vid_path = ReplaceString(vid_path, file_path, vid_path +":cropped")
end

function extract_name(file_path)
	string file_path
	string/g filename, folder_path
	filename = stringfromList(itemsInList(file_path,":")-1,file_path,":")
	folder_path = replaceString(filename, file_path, "")
	print folder_path
	filename = stringfromList(0, filename, ".")
	
	
end


function droplet_window()
   Execute "DropletTrackerPanel()"
	string/G category = ""
	
	wave force
	// Main function for starting the procedure
	// Creates the control window
	KillWindow/Z FusFit
	NewPanel/W= (100, 100, 1800, 800)
	DoWindow/C/T FusFit, "Fusion Fitting" 
	
	NewPanel/W=(0.745,0.01,0.999,0.3)/HOST=FusFit/N=controls
	Button bLoad pos={100, 20}, title="Load data", size={110, 20}, proc=LoadDataProc 	
	Button bNormalise pos={100, 45}, title="Normalise Force", size={110, 20}, proc=NormProc
	Button bLogFit pos={220, 20}, title="Fit logistic", size={110,20}, proc=LogFitProc
	Button bSaveFig pos={220, 45}, title="Save graph", size={110, 20}, proc=SaveGraphProc
	Button bSaveWaves pos={100, 70}, title="Save fit values", size={110, 20}, proc=SaveWavesProc
	Button bCleanLists pos={220, 140}, title="Clean-up fit waves", size={110, 20}, proc=CleanListsProc
	Setvariable strCat pos={100, 140}, size={110,20}, variable=category
	Button bLFC pos={220,70}, title="Show logfold change", size={110,20},proc=procLFC,fsize=10
	wave tau_logistic
	if (waveexists(tau_logistic) == 0)
		make/o/t/n=1 names
		make/o/n=1 tau_logistic
		make/o/n=1 gamma_logistic
		make/o/n=1 vel_logistic
		make/o/n=1 vel_n
		make/o/N=1 size1
		make/O/N=1 size2
		make/O/N=1 size
		size1 = nan
		size2 = nan
		size = Nan
		tau_logistic = NaN
		gamma_logistic = NaN
		vel_logistic = NaN
		vel_n = NaN
	endif
	edit/N=FusData names, tau_logistic, gamma_logistic, vel_logistic, size, size1, size2
end







function LogFitProc(ba) : ButtonControl
	// Button control to perform exponential fit to force data
	Struct WMButtonAction &ba
	switch(ba.eventCode)
		case 2:
			wave force
			variable vel,i
			for(i=0;i<3;i++)
				vel = fit_logistic(force)
				normalise_force(force)
			endfor			
			nvar vel_norm
			wave tau_logistic, tau, gamma_logistic, vel_logistic, w_coef, vel_n
			tau_logistic[numpnts(tau_logistic)-1] = w_coef[3] // Logistic tau value
			gamma_logistic[numpnts(gamma_logistic)-1] = w_coef[4] // Logistic gamma
			vel_logistic[numpnts(vel_logistic)-1] = vel // Logistic velocity
			vel_n[numpnts(vel_n)-1] = vel_norm
			break
		case -1:
		break
	endswitch
end




function LoadDataProc(ba) : ButtonControl
	// Button control to select and load HDF5 file
	Struct WMButtonAction &ba
	switch(ba.eventCode)
		case 2:
			select_file()
			svar file_path
			string/G vid_path			
			titleBox tbFileName pos={30,170}, variable=file_path, size={150,20}, win=FusFit#controls
			load_data(file_path)
			wave force, timestamps
			DoWindow/F FusFit
			Display/HOST=FusFit/N=forceVtime/W=(0.01,0.01,0.74,0.99) force
			ShowInfo/W=FusFit
			Label/W=FusFit#forceVtime bottom "Time"
			Cursor A force 0
			Cursor B force numpnts(force)
			get_vid_name()
			LoadVideo(vid_path)
			wave tau_logistic,  gamma_logistic, vel_logistic, vel_n, size1, size2, size
			wave/t names
			if (numtype(tau_logistic[numpnts(tau_logistic)-1]) != 2)
				redimension/N=(numpnts(names)+1) names
				redimension/N=(numpnts(tau_logistic)+1) tau_logistic
				redimension/N=(numpnts(gamma_logistic)+1) gamma_logistic
				redimension/N=(numpnts(vel_logistic)+1) vel_logistic
				redimension/N=(numpnts(vel_n)+1) vel_n
				redimension/N=(numpnts(size1)+1) size1
				redimension/N=(numpnts(size2)+1) size2
				redimension/N=(numpnts(size)+1) size
				tau_logistic[numpnts(tau_logistic)-1] = NaN
				gamma_logistic[numpnts(gamma_logistic)-1] = NaN
				vel_logistic[numpnts(vel_logistic)-1] = NaN
				vel_n[numpnts(vel_n)-1] = NaN
				size1[numpnts(size1)-1] = NaN
				size2[numpnts(size2)-1] = NaN
				size[numpnts(size)-1] = NaN
			endif	
			break
		case -1:
			break
	endswitch	
	
end


function NormProc(ba) : ButtonControl
	// Button procedure to normalise force
	Struct WMButtonAction &ba
	wave force
	wave timestamps
	switch(ba.eventCode)
		case 2: // mouse up
			normalise_force(force)
			break
		case -1: // control being killed
			break
	endswitch
	return 0
end

function SaveGraphProc(ba) : ButtonControl
	// Button procedure to save graph as .png (filename same as .h5 file)
	Struct WMButtonAction &ba
	switch(ba.eventCode)
		case 2:
			svar file_path
			nvar s1, s2
			wave size1, size2, size
			extract_name(file_path)
			svar filename, folder_path
			string dest = folder_path + filename
			print dest
			SavePICT/O/WIN=FusFit#forceVtime/E=-5 as dest+".png"
			wave/T names
			filename = Stringfromlist(ItemsinList(filename, " ")-1 ,filename, " ")
			print filename
			names[numpnts(names)-1] = filename
			size1[numpnts(size1)-1] = s1
			size2[numpnts(size2)-1] = s2
			size[numpnts(size)-1] = (s1+s2)/2
			break
		case -1:
			break
	endswitch
end





function CleanListsProc(ba) : ButtonControl
	Struct WMButtonAction &ba
	switch(ba.eventCode)
		case 2:
			cleanlists()
			break
		case -1:
			break
	endswitch
end



function CleanLists()  // Remove discarded events (those that don't have a name)
	wave tau_logistic, gamma_logistic, vel_logistic,vel_n, size, size1, size2
	wave/t names
	svar category
	variable i
	for(i=0;i<numpnts(tau_logistic);i+=1)	
		if(stringmatch(names[i], ""))
			DeletePoints i, 1, names, tau_logistic, gamma_logistic, vel_logistic, vel_n, size, size1, size2
			print "Removed point", i
			i--
		endif
	endfor	
	// Move data to category waves					
	duplicate/O names $(category+"_names")
	duplicate/O vel_n $(category)
	duplicate/FREE vel_n vel_nn
	duplicate/O size $(category+"_size")
	vel_nn *= size	// Size normalized velocity
	duplicate/O vel_nn $(category+"_n")  
	// Clear waves
	DeletePoints 0, numpnts(names), names, tau_logistic, gamma_logistic, vel_logistic, vel_n, size, size1, size2
	Redimension/N=1 names, tau_logistic,gamma_logistic, vel_logistic, vel_n, size, size1, size2
end


function save_waves()
	// Save fit value waves as tab-separated file
	svar folder_path, filename
	saveTableCopy/W=FusDATA/f=1/t=1 as folder_path+"fits.csv"
end

function SaveWavesProc(ba) : ButtonControl
	// Button procedure for saving waves of fitted values
	Struct WMButtonAction &ba
	switch(ba.eventCode)
		case 2:
			save_waves()
			break
		case -1:
			break
	endswitch
end




//e.g. SplitWaveByFactor(vel_logistic, names, "UTP;ATP;AMPPNP;")
//Returns a wave reference wave to factored data
Function/Wave SplitWaveByFactor(Wave w, Wave/T factors, String factorlist [, Variable verbose])
	if(numpnts(w) != numpnts(factors))
		abort "Waves need to have same number of points"
	endif
	if(!WaveExists(w) || !WaveExists(factors))
		abort "Wave not found"
	endif

	verbose = paramisdefault(verbose) ? 0 : verbose

	Make/O/N=(itemsinlist(factorlist))/FREE/WAVE resultwaves

	Variable i
	for(i=0;i<itemsinlist(factorlist);i+=1)
		String factor = stringfromlist(i, factorlist)
		String newWaveName = nameofwave(w)+"_"+factor
		
		Extract/O w, $newWaveName, stringmatch(factors, "*"+factor+"*")
		resultwaves[i] = $newWaveName
		
		if(verbose)
			printf "Created %s with %d points\r", newWaveName, numpnts($newWaveName)
		endif
	endfor
	
	return resultwaves
End

function normvel(string name)  // Perform size normalisation of velocity and rename result waves
	wave  r1, r2, vel_n
	duplicate vel_n $name
	duplicate r1 size
	size = (r1+r2)/2
	duplicate vel_n vel_nn
	vel_nn = vel_n * size
	duplicate size $(name+"_s")
	duplicate vel_nn $(name+"_n")
end


function normalise_to_chromatin(string chromatin)  // Normalise all normVel waves in CWD so that chromatin only sample == 1
	string waves = WaveList("*_n", ";","")
	variable normfactor = mean($(chromatin+"_n"))
	wave/T wavenames = listtoTextWave(waves, ";") 
	variable i
	for(i=0; i<numpnts(wavenames);i++)	
		duplicate/O $(wavenames[i]) $(wavenames[i]+"_norm")
		wave w = $(wavenames[i]+"_norm")
		w /= normfactor
	endfor	
end


function get_vid_name()  // Extracts name of cropped tiff file from selected .h5 file
	svar file_path
	svar vid_path
	string vid_file= ReplaceString(".h5", StringFromList(itemsInList(file_path," ")-1,file_path," "), ".tiff.tiff")
	vid_path = ReplaceString(StringFromList(itemsInList(file_path,":")-1,file_path,":"), file_path, "")+"cropped:"+vid_file
end


function logfoldChange()  // Get log10-fold change relative to chromatin condition
	string waves = WaveList("*_n",";","")
	wave/T wavenames = ListtoTextWave(waves,";")
	if(waveexists($("chrom_n")))
		wave chromatin = $("chrom_n")
	else
		wave chromatin = $("chromatin_n")
	endif
	variable meanChromatin = mean(chromatin)
	variable i
	for(i=0;i<numpnts(wavenames);i++)
		duplicate/O $(wavenames[i]) $(wavenames[i]+"_log")
		wave lfc = $(wavenames[i]+"_log")
		lfc /= meanChromatin
		lfc = log(lfc)
	endfor
end


function procLFC(ba) : ButtonControl
	Struct WMButtonAction &ba
	switch(ba.eventCode)
		case 2:
			wave chromatin_n, chrom_n
			if(!waveexists(chromatin_n))
				if(!waveexists(chrom_n))
					break
				endif
			endif
			logfoldChange()
			plotLFC()
			break
		case -1:
		break
	endswitch
end

function plotLFC()
	string waves = WaveList("*_n_log",";","")
	make/O/N=1 tickLocs
	make/O/N=1/T tickNames
	tickNames="Chromatin"
	wave chrom_n_log
	Display;AppendBoxPlot chrom_n_log
	variable i
	for(i=0;i<ItemsInList(waves,";");i++)
		string condition  = StringFromList(i,waves,";")
		if (!stringmatch(condition,"vel_n_log")&& !stringmatch(condition,"chrom_n_log"))
			AddwavestoBoxPlot $(condition)
			Redimension/N=(numpnts(ticklocs)+1) ticklocs, tickNames
			ticknames[numpnts(ticknames)-1] = StringfromList(0,condition,"_")
		endif	
	endfor
	ticklocs = x
	Modifygraph userticks(bottom)={ticklocs,ticknames}
	ModifyBoxPlot trace=chrom_n_log,markers={-1,8,8},markerSizes={3,3,3},lineThickness={1.5,1.5,1.5},markerThick={1.5,1.5,1.5},boxwidth=20
	ModifyGraph height=120
	ModifyGraph width={Aspect,0.25*itemsinList(waves,";")}
	SetAxis bottom -0.2,*
	Label left "log-fold change"
	SetAxis left -1.5,*
end