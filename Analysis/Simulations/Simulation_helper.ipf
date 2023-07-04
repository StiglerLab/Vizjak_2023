#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
function run()
	svar filename
	readparticlestrajectory("")
	makeparticlewin()
end






function filenames_to_condition()
	wave/T files
	variable i
	duplicate/T/O files, conditions
	for(i=0;i<numpnts(files);i++)
		conditions[i] = filename_to_condition(files[i]) 
	endfor
end

function/T filename_to_condition(string filename)
	variable n_parts = ItemsinList(filename,":")
	string condition = StringFromList(0, StringfromList(n_parts-1, filename, ":"),".")
	return condition
end


function normalize_func(string name)

duplicate/O $name, $(name+"_n")
wave n = $(name+"_n")
variable v_min = wavemin(n)
n -= v_min
variable v_max = wavemax(n)
n /= v_max
end



function xspread()
	wave particles_snapshot_0, particles_snapshot_1, particles_snapshot
	Concatenate/O/NP=0 {particles_snapshot_0, particles_snapshot_1}, histones
	Duplicate/O particles_snapshot histones
	Wavestats/PCST histones
	wave M_wavestats
	return M_wavestats[4][0]
end 

function xspread_trajectory()
	wave particles
	variable i
	make/O/N=(dimsize(particles,2)) spread_x, ratio, sphere
	setscale/P x, 0, dimdelta(particles,2), spread_x, ratio, sphere
	UpdateSnapshot(0)
	wave particles_snapshot_0, particles_snapshot_1
	wavestats/Q/RMD=[][0] particles_snapshot_0
	variable s0 = v_sdev
	wavestats/Q/RMD=[][0] particles_snapshot_1
	variable s1 = v_sdev
	variable size = (s0+s1)/2
	for(i=0;i<dimsize(particles,2);i++)
		UpdateSnapshot(i)
		spread_x[i] = xspread()
		wave M_wavestats
		ratio[i] = M_wavestats[4][0]/M_wavestats[4][1]
		sphere[i] = (M_wavestats[4][2]^2/(M_wavestats[4][0] * M_wavestats[4][1]))^(1/3)
	endfor
end


function fit_fusion()
	wave w_coef
	wave ratio
	Variable finish = w_coef[0] //final level
	Variable start = w_coef[1] //initial level
	Variable x0 = w_coef[2] //"mid" parameter
	Variable tau=w_coef[3] //"slope" parameter
	Variable gam = w_coef[4] //"shape" parameter
	Variable vel = (gam*(finish - start))/((1 + 1/gam)^gam*((1 + gam)*tau)) //"Velocity" (dy/dx) at inflection point. Watch out for units!
	print "Velocity (dy/dx)", vel
	variable/G vel_norm = vel/(finish-start)
	print "Normalized velocity ((dy/dx)/(final-initial)", vel/(finish-start)
	Variable inflPoint = x0-tau*ln(gam) //Inflection point X
	Variable inflY = finish + (-finish + start)/(1 + 1/gam)^gam
	Variable root = (start-inflY)/vel+inflPoint 
	Setscale/P x, leftx(ratio)-root, deltax(ratio), "s", ratio
end


function getREMinDroplet([variable n_frames])
	n_frames = paramIsDefault(n_frames) ? 100 : n_frames
	wave CM, particles_snapshot_SPECIAL, particles
	Variable fromT = 0
	nvar radius
	variable t
	make/O/N=(dimsize(particles_snapshot_special,0)) REMinDrop = 2
	for(t=0;t<dimsize(CM,0);t++)
		UpdateSnapshot(t)
		variable n
		for(n=0;n<dimsize(particles_snapshot_special,0);n++)
			if(REMinDrop[n] == 2)  // only compute distance if particle hasn't been eliminated yet
				variable distance_from_COM = sqrt((particles_snapshot_special[n][0]-CM[t][0])^2 + (particles_snapshot_special[n][1]-CM[t][1])^2 + (particles_snapshot_special[n][2]-CM[t][2])^2)
				if(distance_from_COM>radius+50)
					REMinDrop[n] = Nan
				endif
			endif
		endfor	
	endfor
end


function makeREMinDropletWave()  // Creates wave with remodeller coordinates in droplet
	wave REMinDrop, particles, particles_snapshot_special
	variable i
	for(i=0;i<numpnts(REMinDrop);i++)
		if(REMinDrop[i] != 2)
			particles_snapshot_special[i][*] = Nan
		endif
	endfor
	Splitwave/O/OREF=NaNRefs particles_snapshot_special
	Make/O/N=(dimsize(particles_snapshot_special,0))/FREE rowHasBlank = AnyColumnIsBlank(NaNRefs,p) // p is also table row number
	RemoveRowsWithBlanks(NaNRefs, rowHasBlank)
	wave wx = nanrefs[0]
	wave wy = nanrefs[1]
	wave wz = nanrefs[2]
	if(!exists("w_indices"))
		duplicate wx w_indices
	else
		wave w_indices
	endif
	concatenate/O {wx,wy,wz}, FRAP_remodellers
	sortcolumns keywaves=w_indices,sortWaves=FRAP_remodellers
	Duplicate/O/R=[0, dimsize(FRAP_remodellers,0)/2-1] FRAP_remodellers, FRAP_L
	Duplicate/O/R=[dimsize(FRAP_remodellers,0)/2, dimsize(FRAP_remodellers,0)-1] FRAP_remodellers, FRAP_R
	// Kill made waves
	Killwaves $(stringfromList(0, s_wavenames)),$(stringfromList(1, s_wavenames)),$(stringfromList(2, s_wavenames)) 
end



Function AnyColumnIsBlank(waveRefs,row)
    WAVE/WAVE waveRefs
    Variable row
   
    Variable i, nWaves= numpnts(waveRefs)
    for(i=0; i < nWaves; i+=1 )
        WAVE w = waveRefs[i]
        Variable type = numtype(w[row])
        if( type == 2 )     // NaN
            return 1            // at least one "column" is NaN
        endif
    endfor
   
    return 0 // none were blank if we got here
End


Function RemoveRowsWithBlanks(waveRefs, rowHasBlank)
    WAVE/WAVE waveRefs
    Wave rowHasBlank
   
    // work from the end to the start, because DeletePoint changes the effective row Number.
    Variable nRows = numpnts(rowHasBlank)
    Variable row= nRows-1
    do
        Variable hasBlank = rowHasBlank[row]
        if( hasBlank )
            Variable i, nWaves= numpnts(waveRefs)
            for(i=0; i < nWaves; i+=1 )
                WAVE w = waveRefs[i]
                DeletePoints/M=0 row, 1, w
            endfor
        endif
        row = row-1
    while( row >= 0 )
End


function ConcatenateCoordinates(waverefs)
	WAVE/WAVE waverefs
	wave FRAP_remodellers
	
	// waverefindices is used to sort the waves based on their initial x-position 
	// if it doesn't exist it is created and sorted according to the current remodeller wave
	// FRAP_remodellers can then be split in half to simulate FRAP
	wave wx = waverefs[0]
	wave wy = waverefs[1]
	wave wz = waverefs[2]
	Concatenate/O {wx, wy, wz}, FRAP_remodellers

	if (!Waveexists($"waverefindices"))
		duplicate wx waverefindices
		sort wx, waverefindices
	else
		wave waverefindices
	endif
	Sortcolumns keywaves=waverefindices,sortWaves=FRAP_remodellers
end



function RemodellerMixing()  // Get mixing index for REMODELLER in droplet in current snapshot 
	wave FRAP_remodellers
	Variable NHIST = 50
	Variable MAXDIST = 200
	Make/O/N=(NHIST) dist_all_hist_rem=0, dist_cross_hist_rem=0
	setscale/I x,0, MAXDIST, dist_all_hist_rem, dist_cross_hist_rem
	Make/O/N=(dimsize(FRAP_remodellers, 0))/FREE dist_rem = 0
	Make/O/N=(dimsize(FRAP_remodellers, 0)/2)/FREE dist_cross_rem = 0
	
	Variable i,j
	for(i=0;i<dimsize(FRAP_remodellers,0);i++)
		dist_rem = p==i ? NaN : sqrt( (FRAP_remodellers[i][0]-FRAP_remodellers[p][0])^2 + (FRAP_remodellers[i][1]-FRAP_remodellers[p][1])^2 + (FRAP_remodellers[i][2]-FRAP_remodellers[p][2])^2 ) 
		Histogram/A dist_rem, dist_all_hist_rem
	
		if(i<dimsize(FRAP_remodellers, 0)/2)
			dist_cross_rem = sqrt( (FRAP_remodellers[i][0]-FRAP_remodellers[p+dimsize(FRAP_remodellers,0)/2][0])^2 + (FRAP_remodellers[i][1]-FRAP_remodellers[p+dimsize(FRAP_remodellers,0)/2][1])^2 + (FRAP_remodellers[i][2]-FRAP_remodellers[p+dimsize(FRAP_remodellers,0)/2][2])^2 ) 
			Histogram/A dist_cross_rem, dist_cross_hist_rem
		endif
	endfor
	
	Variable ar = area(dist_cross_hist_rem)
	dist_cross_hist_rem /= ar
	ar = area(dist_all_hist_rem)
	dist_all_hist_rem /= ar
	Duplicate/O dist_all_hist_rem, tmp, weights
	tmp = (dist_cross_hist_rem - dist_all_hist_rem)^2
	weights = dist_all_hist_rem
	
	Duplicate/O/FREE tmp, tmp2
	tmp2 = tmp/weights
	
	wavestats/Q tmp2
	variable s = v_sum
	variable n = v_npnts
	
	return s/n
end


function RemodellerMixingTraj([variable stepsize, variable n_frames])
	stepsize = paramIsDefault(stepsize) ? 20 : stepsize
	n_frames = paramIsDefault(n_frames) ? 100 : n_frames
	ClearRemodellerMixing()
	wave particles 
	make/O/N=(dimsize(particles,2)) MixingIndexREMODELLER = NaN
	Setscale/P x, 0, dimdelta(particles, 2), MixingIndexREMODELLER
	RadialDensity(n_frames=n_frames)
	getREMinDroplet(n_frames=n_frames)
	variable t
	for(t=0;t<dimsize(particles,2);t+=stepsize)
		UpdateSnapshot(t)	
		MakeRemInDropletWave()
		DoUpdate
		MixingIndexREMODELLER[t] = RemodellerMixing()
	endfor
	wave FRAP_remodellers
	NVAR NSPPARTICLES
	print "Calculated mixing index for remodeller in droplet. " + num2str(dimsize(FRAP_remodellers,0)) + " out of " + num2str(NSPPARTICLES) + " (" + num2str(dimsize(FRAP_remodellers,0)/NSPPARTICLES*100) + "%) remodellers in droplet"
end


function ClearRemodellerMixing()
	wave CM, RemInDrop, w_indices
	KillWaves CM, RemInDrop, w_indices
end


function MakeMIWindow()
	wave MixingIndex, MixingIndexREMODELLER
	KillWindow/Z MIXINGINDEX0
	Display/N=MIXINGINDEX/W=(1000,50,1200,300) 
	if(waveexists(MIXINGINDEXREMODELLER))
		AppendtoGraph MixingIndexRemodeller
		Label left "Mixing Index Remodeller";DelayUpdate
	endif
	if(waveexists(MIXINGINDEX))
	AppendtoGraph/R MixingIndex
	endif
	Label bottom "time (s)";DelayUpdate
	Label right "Mixing Index Chromatin"
	Legend/C/N=text0/J/B=1/A=LC "\\s(MixingIndexREMODELLER) Remodeller\r\\s(MixingIndex) Chromatin"
	ModifyGraph gaps=0,rgb(MixingIndex)=(0,0,0),width=283.465,height={Aspect,1}
	Legend/C/N=text0/J/A=RT/X=9/Y=9
	
end


function MakeARWindow()
	wave ratio, sphere
	KillWindow/Z AR0
	Display/N=AR0/W=(1000,50,1200,300) ratio
	AppendtoGraph/W=AR0/R sphere
	ModifyGraph/W=AR0 rgb(sphere)=(0,0,0)
	Label bottom "time (s)";DelayUpdate
	Label left "Aspect Ratio";DelayUpdate
	ModifyGraph width=283.465,height={Aspect,0.75},gaps=0
	//SetAxis left 1,*
end

function make_chains()
	NEWDataFolder/O CHAINS
	NVAR NSPPARTICLES
	WAVE particlepID_rousebuddy1, particlepID_rousebuddy2
	duplicate/FREE/RMD=[][1][] particlepID_rousebuddy1, rouse1
	duplicate/FREE/RMD=[][1][] particlepID_rousebuddy2, rouse2
	redimension/N=(dimsize(rouse1,0)) rouse1, rouse2
	VARIABLE ROUSELENGTH = 13
	variable/G n_chains = (dimsize(rouse1, 0)-NSPPARTICLES)/ROUSELENGTH
	variable chain
	// Find beginning of all chains
	make/O/N=(n_chains) chain_starts
	findvalue/V=-1 rouse2
	chain_starts[0] = v_value
	for(chain=1;chain<n_chains;chain++)
		variable start = chain_starts[chain-1]+1
		findvalue/V=-1/S=(start) rouse2
		chain_starts[chain] = v_Value
	endfor
	for(chain=0;chain<n_chains;chain++)
		Make/N=(ROUSELENGTH)/O $("root:CHAINS:CHAIN"+num2str(chain))
		wave w = $("root:CHAINS:CHAIN"+num2str(chain))
		variable i
		w[0] = chain_starts[chain]
		for(i=1;i<ROUSELENGTH-1;i+=2)
			w[i] = rouse1[w[i-1]]
			w[i+1] = rouse2[w[i]]
		endfor
	endfor
end

function get_AVG_x()
	NVAR n_chains
	make/O/N=(n_chains) AVGX
	variable i
	for(i=0;i<n_chains; i++)
		AVGX[i] = AVG_X($("root:CHAINS:CHAIN"+num2str(i)))
	endfor
end


function get_AVG_3D()
	NVAR n_chains
	make/O/N=(n_chains) AVGX, AVGY, AVGZ
	variable i
	for(i=0;i<n_chains; i++)
		AVGX[i] = AVG_X($("root:CHAINS:CHAIN"+num2str(i)))
		AVGY[i] = AVG_Y($("root:CHAINS:CHAIN"+num2str(i)))
		AVGZ[i] = AVG_Z($("root:CHAINS:CHAIN"+num2str(i)))
	endfor
end

function AVG_X(wave chain)
	wave particles
	NVAR FRAPSTEP
	Wavestats/RMD=[chain[0], chain[numpnts(chain)-1]][0][FRAPSTEP]/Q particles
	return v_avg
end

function AVG_Y(wave chain)
	wave particles
	NVAR FRAPSTEP
	Wavestats/RMD=[chain[0], chain[numpnts(chain)-1]][1][FRAPSTEP]/Q particles
	return v_avg
end

function AVG_Z(wave chain)
	wave particles
	NVAR FRAPSTEP
	Wavestats/RMD=[chain[0], chain[numpnts(chain)-1]][2][FRAPSTEP]/Q particles
	return v_avg
end


function Make_Sort_Key()
	NVAR n_chains
	WAVE AVGX
	Duplicate/O AVGX X_SORT
	X_SORT = p
	Sort AVGX, X_SORT
	// Loop over X_SORT and concatenate waves
	variable i
	Make/N=0/O NEW_Indices 
	for(i=0;i<numpnts(X_SORT);i++)
		wave w = $("root:CHAINS:Chain"+num2str(X_sort[i]))
		Concatenate/NP=0 {w}, New_Indices
	endfor
	// Make Sort key for nucleosomes
	Duplicate/O New_Indices sort_key
	sort_key = p
	Sort/LOC New_Indices sort_key
	// Make Sort key for remodeller
	NVAR NSPPARTICLES
	Make/FREE/N=(NSPPARTICLES) sort_key_REM
	sort_key_REM = p + numpnts(sort_key) + 1
	Concatenate/NP=0 {sort_key_REM}, sort_key
end


function sort_particles()
	wave particles
	wave sort_key
	SortColumns keywaves=sort_key,sortwaves=particles

end

function FRAP_Chromatin()
	make_chains()
	get_avg_x()
	Make_Sort_key()
	sort_particles()
end






function plotFRAP()
	SetDataFolder root:FRAP
	string Chrom_waves = wavelist("*MIChrom",";","")
	string Rem_waves = wavelist("*MIRem",";","")
	variable i
	// Plot FRAP waves
	variable n = ItemsInList(Chrom_waves,";")
	Display
	for(i=0;i<n;i++)
		string chrom = StringfromList(i,Chrom_waves, ";")
		string rem = StringfromList(i, Rem_waves,";")
		string nucleotide = StringfromList(0,chrom, "_")
		wave chromatin = $chrom
		wave remodeler = $rem
		Appendtograph chromatin
		Appendtograph remodeler
		// Set colour
		if(stringmatch(nucleotide, "ATP"))
			ModifyGraph rgb($chrom)=(65535,0,0),rgb($rem)=(65535,0,0)						// Red
		elseif(stringmatch(nucleotide, "AMPPNP"))
			ModifyGraph rgb($chrom)=(0,0,65535),rgb($rem)=(0,0,65535)						// Blue
		else
			ModifyGraph rgb($chrom)=(34952,34952,34952),rgb($rem)=(34952,34952,34952)	// Grey
		endif
	endfor
	Modifygraph gaps=0
	ModifyGraph width=283.465,height={Aspect,1}
	SetDataFolder root:
end

function get_time(variable t)
	wave particles
	variable timer = t*dimdelta(particles,2)
	return timer
end



Function sphericity()
	wave particles
	variable i
	make/O/N=(dimsize(particles,2)) spread_x, ratio
	setscale/P x, 0, dimdelta(particles,2), spread_x, ratio
	UpdateSnapshot(0)
	wave particles_snapshot_0, particles_snapshot_1
	wavestats/Q/RMD=[][0] particles_snapshot_0
	variable s0 = v_sdev
	wavestats/Q/RMD=[][0] particles_snapshot_1
	variable s1 = v_sdev
	variable size = (s0+s1)/2
	for(i=0;i<dimsize(particles,2);i++)
		UpdateSnapshot(i)
		spread_x[i] = xspread()
		wave M_wavestats
		ratio[i] = (M_wavestats[4][2]^2/(M_wavestats[4][0] * M_wavestats[4][1]))^(1/3)
	endfor

End


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


function sortParticleByType(variable t)
	// Sort particles_snapshot into nucleosome and remodeler waves 
	// required if init files for simulation don't follow nuc->rem structure
	// e.g. if fused condensates are used as init for another fusion sim
	wave particles_snapshot, particlestates
	NVAR NSPParticles
	wave particles_snapshot_0, particles_snapshot_1, particles_snapshot_special
	variable n_nucleosomes = dimsize(particles_snapshot,0)-NSPParticles
	make/FREE/N=(0, 3) snap_sp, snap_0, snap_1
	variable i
	for(i=0;i<dimsize(particles_snapshot,0);i++)
		if(particlestates[i][0]!=0)
			Redimension/N=(dimsize(snap_sp,0)+1,3) snap_sp
			snap_sp[dimsize(snap_sp,0)-1][0] = particles_snapshot[i][0]
			snap_sp[dimsize(snap_sp,0)-1][1] = particles_snapshot[i][1]
			snap_sp[dimsize(snap_sp,0)-1][2] = particles_snapshot[i][2]
		else
			if(i<dimsize(particles_snapshot,0)/2)
				Redimension/N=(dimsize(snap_0,0)+1,3) snap_0
				snap_0[dimsize(snap_0,0)-1][0] = particles_snapshot[i][0]
				snap_0[dimsize(snap_0,0)-1][1] = particles_snapshot[i][1]
				snap_0[dimsize(snap_0,0)-1][2] = particles_snapshot[i][2]
			else
				Redimension/N=(dimsize(snap_1,0)+1,3) snap_1
				snap_1[dimsize(snap_1,0)-1][0] = particles_snapshot[i][0]
				snap_1[dimsize(snap_1,0)-1][1] = particles_snapshot[i][1]
				snap_1[dimsize(snap_1,0)-1][2] = particles_snapshot[i][2]
			endif
		endif
	endfor
	duplicate/O snap_sp particles_snapshot_special
	duplicate/O snap_0 particles_snapshot_0
	duplicate/O snap_1 particles_snapshot_1
end

function updateAR(variable i)
	wave ratio, AR2
	duplicate/FREE ratio AR_tmp
	duplicate/O/R=[0,i] AR_tmp AR2
end


function appendsim()
	// Combines 2 simulation runs into one trajectory
	// Adds runpt1 wave at the beginning of particles wave
	wave particles, runpt1, run_con
	Concatenate/NP/O {runpt1, particles}, run_con
	Duplicate/O particles, run2
	Duplicate/O run_con, particles
	Slider sldTimestep,limits={0,dimsize(particles,2),1}
end
