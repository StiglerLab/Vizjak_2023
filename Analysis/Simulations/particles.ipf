#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
//#include "Simulation_helper"


Structure float2
	float x
	float y
ENDStructure

Structure int2
	int32 x
	int32 y
ENDStructure

Structure uint2
	uint32 x
	uint32 y
ENDStructure

Structure float3
	float x
	float y
	float z
ENDStructure

Structure int3
	int32 x
	int32 y
	int32 z
ENDStructure

structure grid2D
  STRUCT float2 origin      //origin of grid
  STRUCT float2 boxsize     //nanometer size of full grid
  STRUCT int2 Ngrid         //MxN size of full grid

  //derived
  STRUCT float2 cellDelta   //nanometer size of one cell
EndStructure

structure grid3D
  STRUCT float3 origin      //origin of grid
  STRUCT float3 boxsize     //nanometer size of full grid
  STRUCT int3 Ngrid         //MxN size of full grid

  //derived
  STRUCT float3 cellDelta   //nanometer size of one cell
EndStructure


structure Pot_LJ_params
  double R0
  double epsilon
  double linCutoff //linear connecting line from this value to cutoff as in NAMD
  double cutoff //potential will be zero for this value
  double slopePoint    //linear interpolate from x=0 to x=slopePoint, matching potential slope at slopepoint
EndStructure

structure Pot_Hooke_params
  double R0
  double k
  double cutoff
EndStructure



structure SimParams2D
  uint32 version        //version of SimParams struct
  uint32 dimension      //2 or 3
  uint32 NPARTICLES     //number of particles in simulation
  uint32 DUMMY_
  STRUCT grid2D grid      //simulation grid

  uint32 outputInterval //number of delta-timesteps when outputs are written
  uint32 DUMMY2_
  
  double dt           //time step (seconds)
  double D            //diffusion coefficient (nm^2/second)

  uint32 potentialKind  //type of potential (enum)
  uint32 DUMMY3_
  STRUCT Pot_LJ_params LJ_params
  STRUCT Pot_Hooke_params hooke_params
  
  //derived
  double gam           //friction coefficient = kT/D
  double randFactor    //sqrt(2*kT*gam/dt)

  double dt_over_gamma
  
  uint32 PBC           //periodic boundary conditions
  uint32 DUMMY4_
  //uint32 DUMMY5_
EndStructure

structure SimParams3D
  uint32 version        //version of SimParams struct
  uint32 dimension      //2 or 3
  uint32 NPARTICLES     //total number of particles in simulation
  uint32 NSPPARTICLES   //number of special particles
  //uint32 DUMMY_
  //uint32 DUMMY2_
  STRUCT grid3D grid      //simulation grid

  uint32 outputInterval //number of delta-timesteps when outputs are written
  uint32 DUMMY2_
  
  double dt           //time step (seconds)
  double D            //diffusion coefficient (nm^2/second)

  uint32 potentialKind  //type of potential (enum)
  uint32 DUMMY3_
  STRUCT Pot_LJ_params LJ_params
  STRUCT Pot_LJ_params LJ_rem_params
  
  STRUCT Pot_Hooke_params hooke_params
  STRUCT Pot_Hooke_params hooke_rem_params
  
  
  //derived
  double gam           //friction coefficient = kT/D
  double randFactor    //sqrt(2*kT*gam/dt)

  double dt_over_gamma
  
  uint32 PBC           //periodic boundary conditions
  //uint32 DUMMY4_
  //uint32 DUMMY5_

	float atp
	//uint32 DUMMY6_
	float khyd
	float krelease
	
	char dummy0[256] //empty space for future variables
	char dummy1[256] //empty space for future variables
	char dummy2[256] //empty space for future variables
	char dummy3[256] //empty space for future variables
EndStructure




Function readParticlesTrajectory(string filename)
	Variable refNum = -1
	if(strlen(filename)==0)
		variable temp 
		open/D/R/F="*traj file: *.traj: LAMMPS file: *.lammpstrj" temp
		string/G filepath = s_filename
	else
		string/G filepath = filename
	endif
	print filepath
	SVAR file  
	file = StringFromList(0,StringFromList(ItemsinList(filepath,":")-1,filepath, ":"),".")
	try
		open/R/Z refNum as filepath
		
		if(v_flag == -1)
			return NaN
		endif
	
		FStatus refNum
		if(v_flag == 0)
			print "ERROR"
			return NaN
		endif
		
		Variable fileSizeBytes = V_logEOF
		printf "FileSize: %d bytes\r", fileSizeBytes
		
		//Read simulation parameters
		
		STRUCT simParams2D simParams
		
		Variable coordinate
		Variable type, state, status, pID_buddy1, pID_buddy2, pID_rousebuddy1, pID_rousebuddy2, pID

		Variable energy
		Variable dt, outputInterval
		
		FBinRead/B=0/F=0 refNum, simParams
		if(simParams.dimension==3)
			printf "This is a 3D file!"
			fsetpos refNum, 0
			STRUCT simParams3D simParams3D
			FBinRead/B=0/F=0 refNum, simParams3D
			print simParams3D
			dt = simParams3D.dt
			outputInterval = simParams3D.outputInterval
			
			Variable/G NSPPARTICLES = simParams3D.NSPPARTICLES
		else
			print simParams
			dt = simParams.dt
			outputInterval = simParams.outputInterval		
		endif		
		
	
		fgetPos refNum
		print "Header size", v_filePos
		
		
		//A String "ROUSECON" should be found at the start of rouse connectivity section. Make sure it's there.
		String rousecon = ""
		FreadLine/N=8/T="" refNum, rousecon
		if(stringmatch(rousecon, "ROUSECON")==0)
			abort "Wrong header. rousecon="+rousecon
		endif
		
		Variable NrouseConnections
		FbinRead/B=0/U/F=3 refNum, NrouseConnections
		Make/O/N=(NrouseConnections,2) s_RouseConnections = 0 //simplified rouse connections
		Variable j
		for(j=0;j<NrouseConnections;j+=1)
			STRUCT uint2 conn
			FbinRead/B=0/F=0 refNum, conn
			s_rouseConnections[j][0] = conn.x
			s_rouseConnections[j][1] = conn.y
		endfor
		
		//A String "DATADATA" should be found at the start of data section. Make sure it's there.
		String datadata = ""
		FreadLine/N=8/T="" refNum, datadata
		if(stringmatch(datadata, "DATADATA")==0)
			abort "Wrong header. datadata="+datadata
		endif
		
		
		Make/O/N=(simParams.NPARTICLES, simParams.dimension, 0) particles
		Setscale/P z, 0, dt * outputInterval, particles
		Make/O/N=(0) totalenergy
		Setscale/P x, 0, dt * outputInterval, totalenergy

		Make/O/N=(simParams.NPARTICLES, 0) particlestates
		Setscale/P y, 0, dt * outputInterval, particlestates
		
		Make/O/N=(simParams.NPARTICLES, 0) particlestatuses
		Setscale/P y, 0, dt * outputInterval, particlestatuses
		
		Make/O/N=(simParams.NPARTICLES, 0) particlepID_buddy1, particlepID_buddy2
		Setscale/P y, 0, dt * outputInterval, particlepID_buddy1, particlepID_buddy2
		
		
       Make/O/N=(simParams.NPARTICLES, 0) particlepID_rousebuddy1, particlepID_rousebuddy2
       Setscale/P y, 0, dt * outputInterval, particlepID_rousebuddy1, particlepID_rousebuddy2

       Make/O/N=(simParams.NPARTICLES, 0) particleID
       Setscale/P y, 0, dt * outputInterval, particleID

		
		Variable i, d
		do
			fgetPos refNum
			if(v_filePos >= fileSizeBytes)
				break
			endif
			
			Redimension/N=(-1,-1,dimsize(particles,2)+1) particles
			Redimension/N=(dimsize(totalenergy,0)+1) totalenergy
			Redimension/N=(-1,dimsize(particlestates,1)+1) particlestates
			Redimension/N=(-1,dimsize(particlestatuses,1)+1) particlestatuses
			Redimension/N=(-1,dimsize(particlepID_buddy1,1)+1) particlepID_buddy1
			Redimension/N=(-1,dimsize(particlepID_buddy2,1)+1) particlepID_buddy2
          Redimension/N=(-1,dimsize(particlepID_buddy1,1)+1) particlepID_rousebuddy1
          Redimension/N=(-1,dimsize(particlepID_buddy2,1)+1) particlepID_rousebuddy2
          Redimension/N=(-1,dimsize(particlepID_buddy1,1)+1) particleID
			
			for(i=0;i<simParams.Nparticles;i+=1)
				for(d=0;d<simParams.dimension;d+=1)
					FBinRead/B=0/F=4 refNum, coordinate //read single float [F=5: double, F=4: float]
					particles[i][d][dimsize(particles,2)-1] = coordinate
				endfor
				FBinRead/B=0/F=3 refNum, type //read single int (particle type) [F=3: int32]
				FBinRead/B=0/F=3 refNum, state //read single int (particle state) [F=3: int32]
				FBinRead/B=0/F=3 refNum, status //read single int (particle status: unpaired, paired) [F=3: int32]
				FBinRead/B=0/F=3 refNum, pID_buddy1 //read single int (buddy a remodeler is paired with) [F=3: int32]
				FBinRead/B=0/F=3 refNum, pID_buddy2 //read single int (buddy a remodeler is paired with) [F=3: int32]
				FBinRead/B=0/F=3 refNum, pID_rousebuddy1 //read single int (rouse buddy a nucleosome is paired with) [F=3: int32]
             FBinRead/B=0/F=3 refNum, pID_rousebuddy2 //read single int (rouse buddy a nucleosome is paired with) [F=3: int32]
             FBinRead/B=0/F=3 refNum, pID //read single int (particle ID) [F=3: int32]
				particlestatuses[i][dimsize(particlestatuses,1)-1] = status
				particlestates[i][dimsize(particlestates,1)-1] = state
				particlepID_buddy1[i][dimsize(particlepID_buddy1,1)-1] = pID_buddy1
				particlepID_buddy2[i][dimsize(particlepID_buddy2,1)-1] = pID_buddy2
				particlepID_rousebuddy1[i][dimsize(particlepID_rousebuddy1,1)-1] = pID_rousebuddy1
             particlepID_rousebuddy2[i][dimsize(particlepID_rousebuddy2,1)-1] = pID_rousebuddy2
             particleID[i][dimsize(particleID,1)-1] = pID

			endfor
			FBinRead/B=0/F=4 refNum, energy //read float energy
			totalenergy[dimsize(totalenergy,0)-1] = energy
		while(1)	
		
		duplicate/O totalenergy, logNegTotalEnergy
		logNegTotalEnergy = log(-totalEnergy)
		
		close refNum
		
		print "Number of time steps in file:", dimsize(particles,2)
	catch
		close /A
	endtry
	
	Slider sldTimestep,limits={0,dimsize(particles,2),1},value= 0,side= 0,vert= 0
	NVAR TWOPARTS
	if(TWOPARTS==0)
		Duplicate/O particles, particles_unsorted
		FRAP_Chromatin()
	endif
End





Function Play()
	Wave particles
	Variable i
	for(i=0;i<dimsize(particles,2);i+=1)
		UpdateSnapshot(i)
		doupdate
	endfor
End


Function UpdateSnapshot(Variable t)
	Wave particles
	
	Duplicate/O/R=[][][t] particles, particles_snapshot
	Redimension/N=(-1,-1) particles_snapshot
	
	
	NVAR NSPPARTICLES
	
	nvar TWOPARTS //Fusion
	
	if(twoparts)
		Duplicate/O/R=[,(dimsize(particles_snapshot,0)-NSPPARTICLES)/2-1][] particles_snapshot, particles_snapshot_0
		Duplicate/O/R=[dimsize(particles_snapshot,0)/2,dimsize(particles_snapshot,0)-NSPPARTICLES/2-1][] particles_snapshot, particles_snapshot_1
	
		Duplicate/O/R=[dimsize(particles_snapshot,0)/2-NSPPARTICLES/2,dimsize(particles_snapshot,0)/2-1][] particles_snapshot, particles_snapshot_SPECIAL_part1
		Duplicate/O/R=[dimsize(particles_snapshot,0)-NSPPARTICLES/2,][] particles_snapshot, particles_snapshot_SPECIAL_part2
		Concatenate/NP=0/O {particles_snapshot_SPECIAL_part1, particles_snapshot_SPECIAL_part2}, particles_snapshot_SPECIAL
	else
		Duplicate/O/R=[,(dimsize(particles_snapshot,0)-NSPPARTICLES)/2-1][] particles_snapshot, particles_snapshot_0
		Duplicate/O/R=[dimsize(particles_snapshot,0)/2-NSPPARTICLES-1,dimsize(particles_snapshot,0)-NSPPARTICLES-1][] particles_snapshot, particles_snapshot_1

		Duplicate/O/R=[dimsize(particles_snapshot,0)-NSPPARTICLES,][] particles_snapshot, particles_snapshot_SPECIAL
	endif
	
	updatechainsT(t)
	
	ValDisplay vdFrame,win=Panel0,value=_NUM:(t)
End



Function particlesSlider(sa) : SliderControl
	STRUCT WMSliderAction &sa
	NVAR TWOPARTS
	switch( sa.eventCode )
		case -3: // Control received keyboard focus
		case -2: // Control lost keyboard focus
		case -1: // Control being killed
			break
		default:
			if( sa.eventCode & 1 ) // value set
				Variable curval = sa.curval
				svar filetype
				if(stringmatch(filetype,"traj"))
					updateSnapshot(curval)
					if(TWOPARTS==0)
						MakeRemInDropletWave()
					endif
				endif
				ValDisplay vdFrame,win=Panel0, value=_NUM:curval	
			endif
			break
	endswitch

	return 0
End

Function MakeParticleWin() : Graph
	UpdateSnapshot(0)
	Wave particles_snapshot, particles

	PauseUpdate; Silent 1		// building window...
	Display /W=(899,338,1480,852)/N=ParticleWin particles_snapshot[*][1] vs particles_snapshot[*][0]
	ModifyGraph mode=2
	ModifyGraph lSize=3
	SetAxis left -23.1131031255738,17.8549273842367
	SetAxis bottom -33.3998553043557,21.8632860970295
	Cursor/P A particles_snapshot 0;Cursor/P B particles_snapshot 1
	ShowInfo
	ControlBar 40
	Slider sldTimestep,pos={1.00,1.00},size={568.00,16.00},proc=particlesSlider
	Slider sldTimestep,limits={0,dimsize(particles,2),1},value= 31,side= 0,vert= 0
	ValDisplay vdFrame,pos={1.00,22.00},size={50.00,13.00}
	ValDisplay vdFrame,limits={0,0,0},barmisc={0,1000},value= _NUM:400
	SetAxis/A
End



//PCA to measure spread of blob
Function PCASpread()
	Wave particles
	
	Make/O/N=(dimsize(particles,2), 3) spread = NaN
	Make/O/N=(dimsize(particles,2)) spreadTot = NaN
	duplicate/O spreadtot sphere
	setscale/P x, 0, dimdelta(particles,2), spread, spreadTot, sphere
	
	Variable i
	for(i=0;i<dimsize(particles,2);i+=1)
		UpdateSnapshot(i)
		Wave Particles_snapshot_0, Particles_snapshot_1
		Concatenate/NP=0/O/FREE {Particles_snapshot_0, Particles_snapshot_1}, particles_snapshot_bothParts
		
		//Duplicate/O/FREE/R=[][][i] particles, particles_snapshot

		pca/SEVC particles_snapshot_bothParts
		Wave w_eigen
		spread[i][] = w_eigen[q]

		wavestats/Q/RMD=[][0] particles_snapshot_bothParts //result: v_sdev
		spread[i][0] = v_sdev		
		wavestats/Q/RMD=[][1] particles_snapshot_bothParts //result: v_sdev
		spread[i][1] = v_sdev			
		wavestats/Q/RMD=[][2] particles_snapshot_bothParts //result: v_sdev
		spread[i][2] = v_sdev		
		
		spreadTot[i] = sqrt(spread[i][0]^2 + spread[i][1]^2 + spread[i][2]^2)
		sphere[i] = spread[i][0]
		doupdate
	endfor
End


//Radius of gyration over time for Rouse chains
Function RGyr()
	Wave particles
	
	Wave stretches = RouseConnectionStretches()
	
	
	Variable Nchains = dimsize(stretches,0)
	Make/O/N=(Nchains, dimsize(particles,2)) RadiusOfGyration //(chain, timestep)
	Make/O/N=(dimsize(particles,2)) RadiusOfGyration_avg //timestep
	Make/O/N=(dimsize(particles,2)) RadiusOfGyration_std //timestep
	
	
	Variable ch, t
	for(ch=0;ch<dimsize(stretches,0);ch+=1)
		for(t=0;t<dimsize(particles,2);t+=1)
			Duplicate/FREE/O/R=[stretches[ch][0],stretches[ch][1]][][t] particles, chain
			MatrixOp/O/FREE COM = averagecols(chain) // sumcols(chain) / numrows(chain)
			
			Variable NparticlesInChain = stretches[ch][1]-stretches[ch][0]+1
			MatrixOp/O/FREE RgSq = sum((chain - rowrepeat(COM,NparticlesInChain)) * (chain - rowrepeat(COM,NparticlesInChain)))
			RadiusOfGyration[ch][t] = sqrt(RgSq[0])
		endfor
	endfor
	
	for(t=0;t<dimsize(particles,2);t+=1)
		Wavestats/Q/RMD=[][t] RadiusOfGyration
		RadiusOfGyration_avg[t] = v_avg
		RadiusOfGyration_std[t] = v_sdev
	endfor
	
End


Function MakeChains()
	KillDataFolder/Z root:chains
	NewDataFolder/O root:Chains
	DFREF fld = root:Chains
	
	Wave rc = RouseConnectionStretches()
	Wave particles
	
	Variable i
	for(i=0;i<dimsize(rc,0);i+=1)
		Duplicate/O/R=[rc[i][0],rc[i][1]] particles, fld:$("chainT"+num2istr(i))/Wave=chainT
	endfor

End


Function UpdateChainsT(Variable t)
	DFREF fld = root:Chains
	
	DFREF cfld = getdatafolderdfr()
	try
		setdatafolder fld
		String wlist = wavelist("chainT*", ";", "")
		setdatafolder cfld
	catch
		setdatafolder cfld
	endtry
	
	Variable i
	for(i=0;i<itemsinlist(wlist);i+=1)
		Wave chainT = fld:$("chainT"+num2istr(i))
		Duplicate/O/R=[][][t] chainT, fld:$("chain"+num2istr(i))/Wave=chain
		redimension/N=(-1,-1) chain
	endfor
End

Function AppendChains()
	DFREF fld = root:Chains
	
	DFREF cfld = getdatafolderdfr()
	try
		setdatafolder fld
		String wlist = wavelist("chainT*", ";", "")
		setdatafolder cfld
	catch
		setdatafolder cfld
	endtry
	
	Variable i
	for(i=0;i<itemsinlist(wlist);i+=1)
		Wave chainT = fld:$("chainT"+num2istr(i))
		Duplicate/O/R=[][][0] chainT, fld:$("chain"+num2istr(i))/Wave=chain
		redimension/N=(-1,-1) chain
		AppendToGizmo/D/N=Gizmo0 path=chain,name=$("chain"+num2istr(i))
		ModifyGizmo modifyObject=$("chain"+num2istr(i)), objectType=path, property={drawTube, 1}, property={fixedRadius, 0.05},property={ calcNormals,1}
		if(i<itemsinlist(wlist)/2)
			//make green
			ModifyGizmo ModifyObject=$("chain"+num2istr(i)),objectType=path,property={ pathColor,3.0518e-05,0.6,1.5259e-05,1}
		else
			//make red
			ModifyGizmo ModifyObject=$("chain"+num2istr(i)),objectType=path,property={ pathColor,1,0,0,1}
		endif
	endfor
End	


Function RemoveChains()
	DFREF fld = root:Chains
	
	DFREF cfld = getdatafolderdfr()
	try
		setdatafolder fld
		String wlist = wavelist("chainT*", ";", "")
		setdatafolder cfld
	catch
		setdatafolder cfld
	endtry
	
	Variable i
	for(i=0;i<itemsinlist(wlist);i+=1)
		RemoveFromGizmo/Z/N=Gizmo0 object=$("chain"+num2istr(i))
	endfor
End	



Function/Wave RouseConnectionStretches()
	Wave w = s_RouseConnections
	
	if(mod(dimsize(w,0),2) != 0)
		abort "Rouse connections MUST come in pairs. Something is seriously wrong."
	endif
	
	Make/O/N=(0,2) stretches
	
	if(dimsize(w,0) == 0)
		return stretches
	endif
	
	sortColumns/KNDX={0,1} sortWaves=w
	
	//remove links that connect the same particles
	Duplicate/O/FREE w, wcpy
	
	Variable i
	for(i=dimsize(w,0)-1-1;i>=0;i-=1)
		if(w[i][0] == w[i+1][1] && w[i][1] == w[i+1][0])
			DeletePoints/M=0 i+1, 1, wcpy
		endif
	endfor
	sortColumns/KNDX={0,1} sortWaves=wcpy
	
	Variable sstretch = wcpy[0][0]
	for(i=0;i<dimsize(wcpy,0)-1;i+=1)
		if(wcpy[i][1] != wcpy[i+1][0])
			Redimension/N=(dimsize(stretches,0)+1,-1) stretches
			stretches[dimsize(stretches,0)-1][0] = sstretch
			stretches	[dimsize(stretches,0)-1][1] = wcpy[i][1] //end of stretch

			sstretch = wcpy[i][1]+1 //start next stretch
		endif
	endfor

	//deal with end of stretch list
	Redimension/N=(dimsize(stretches,0)+1,-1) stretches
	stretches[dimsize(stretches,0)-1][0] = sstretch
	stretches[dimsize(stretches,0)-1][1] = wcpy[dimsize(wcpy,0)-1][1] //end of stretch

	return stretches
End




//Find mixing index (redefined!) from trajectory
Function FindMixingTraj([Variable stepsize])
	Wave particles
	Variable Nparticles = dimsize(particles, 0)

	stepsize = paramisdefault(stepsize) ? 20 : stepsize

	make/O/N=(dimsize(particles,2)) MixingIndex = NaN
	Setscale/P x, 0, dimdelta(particles,2), mixingindex

	Variable t
	for(t=0;t<dimsize(particles,2);t+=stepsize)
		UpdateSnapshot(t)
		
		MixingIndex[t] = findMixing()
		DoUpdate
	endfor
End


Function FindMixing()
	Wave particles = particles_snapshot
	
	NVAR NSPPARTICLES
	Variable Nparticles = dimsize(particles, 0) - NSPPARTICLES
	
	Variable NHIST = 200
	Variable MAXDIST = 400
	Make/O/N=(NHIST) dist_ALL_hist = 0, dist_CROSS_hist = 0
	Setscale/I x, 0, MAXDIST, dist_ALL_hist, dist_CROSS_hist
	
	Make/O/N=(Nparticles)/FREE dist = 0
	Make/O/N=(Nparticles/2)/FREE distCross = 0

	
	Variable i, j
	for(i=0;i<Nparticles-NSPPARTICLES;i+=1)
		dist = p==i ? NaN : sqrt( (particles[i][0]-particles[p][0])^2 + (particles[i][1]-particles[p][1])^2 + (particles[i][2]-particles[p][2])^2 ) 
		Histogram/A dist, dist_ALL_hist
		
		if(i<(Nparticles-NSPPARTICLES)/2)
			distCross = sqrt( (particles[i][0]-particles[p+(Nparticles-NSPPARTICLES)/2][0])^2 + (particles[i][1]-particles[p+(Nparticles-NSPPARTICLES)/2][1])^2 + (particles[i][2]-particles[p+(Nparticles-NSPPARTICLES)/2][2])^2 ) 
			Histogram/A distCross, dist_CROSS_hist
		endif
	endfor
	Variable ar = area(dist_CROSS_hist)
	dist_CROSS_hist /= ar
	ar = area(dist_ALL_hist)
	dist_ALL_hist /= ar
	
	
	
	Duplicate/O dist_ALL_hist, tmp, weights
	tmp = (dist_CROSS_hist - dist_ALL_hist)^2
	
	weights = dist_ALL_hist
	
	Duplicate/O/FREE tmp, tmp2
	tmp2 = tmp/weights
	
	Wavestats/Q tmp2
	Variable s = v_sum
	Variable n = v_npnts
	
	return s/n	
End




Window Gizmo0() : GizmoPlot
	PauseUpdate; Silent 1		// building window...
	// Building Gizmo 8 window...
	NewGizmo/W=(100,300,400,600)
	ModifyGizmo aspectratio = 1
	ModifyGizmo startRecMacro=700
	ModifyGizmo scalingMode=8
	ModifyGizmo setOuterBox={-384,384,-384,384,-384,384}
	ModifyGizmo scalingOption=0
	AppendToGizmo Scatter=root:particles_snapshot_0,name=LeftDL
	ModifyGizmo ModifyObject=LeftDL,objectType=scatter,property={ scatterColorType,0}
	ModifyGizmo ModifyObject=LeftDL,objectType=scatter,property={ markerType,0}
	ModifyGizmo ModifyObject=LeftDL,objectType=scatter,property={ sizeType,0}
	ModifyGizmo ModifyObject=LeftDL,objectType=scatter,property={ rotationType,0}
	ModifyGizmo ModifyObject=LeftDL,objectType=scatter,property={ Shape,2}
	ModifyGizmo ModifyObject=LeftDL,objectType=scatter,property={ size,0.12}
	ModifyGizmo ModifyObject=LeftDL,objectType=scatter,property={ color,3.0518e-05,0.6,1.5259e-05,0.500008}
	AppendToGizmo Axes=boxAxes,name=axes0
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,lineWidth,3}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={-1,axisScalingMode,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={-1,axisColor,0,0,0,1}
	ModifyGizmo modifyObject=axes0,objectType=Axes,property={-1,Clipped,0}
	AppendToGizmo freeAxesCue={0,0,0,3},name=freeAxesCue0
	AppendToGizmo light=Directional,name=light0
	ModifyGizmo modifyObject=light0,objectType=light,property={ position,0.000000,0.000000,1.000000,0.000000}
	ModifyGizmo modifyObject=light0,objectType=light,property={ direction,0.000000,0.000000,1.000000}
	ModifyGizmo modifyObject=light0,objectType=light,property={ ambient,0.333333,0.333333,0.333333,1.000000}
	ModifyGizmo modifyObject=light0,objectType=light,property={ specular,1.000000,1.000000,1.000000,1.000000}
	AppendToGizmo Scatter=root:particles_snapshot_1,name=RightDL
	ModifyGizmo ModifyObject=RightDL,objectType=scatter,property={ scatterColorType,0}
	ModifyGizmo ModifyObject=RightDL,objectType=scatter,property={ markerType,0}
	ModifyGizmo ModifyObject=RightDL,objectType=scatter,property={ sizeType,0}
	ModifyGizmo ModifyObject=RightDL,objectType=scatter,property={ rotationType,0}
	ModifyGizmo ModifyObject=RightDL,objectType=scatter,property={ Shape,2}
	ModifyGizmo ModifyObject=RightDL,objectType=scatter,property={ size,0.12}
	ModifyGizmo ModifyObject=RightDL,objectType=scatter,property={ color,1,0,0,1}
	AppendToGizmo Path=root:Chains:chain0,name=chain0
	ModifyGizmo ModifyObject=chain0,objectType=path,property={ pathColorType,1}
	ModifyGizmo ModifyObject=chain0,objectType=path,property={ lineWidthType,0}
	ModifyGizmo ModifyObject=chain0,objectType=path,property={ pathColor,3.052e-05,0.6,1.526e-05,1}
	ModifyGizmo ModifyObject=chain0,objectType=path,property={ drawTube,1}
	ModifyGizmo ModifyObject=chain0,objectType=path,property={ fixedRadius,0.3500}
	ModifyGizmo ModifyObject=chain0,objectType=path,property={ calcNormals,1}
	ModifyGizmo modifyObject=chain0,objectType=Path,property={calcNormals,1}
	AppendToGizmo Path=root:Chains:chain1,name=chain1
	ModifyGizmo ModifyObject=chain1,objectType=path,property={ pathColorType,1}
	ModifyGizmo ModifyObject=chain1,objectType=path,property={ lineWidthType,0}
	ModifyGizmo ModifyObject=chain1,objectType=path,property={ pathColor,3.052e-05,0.6,1.526e-05,1}
	ModifyGizmo ModifyObject=chain1,objectType=path,property={ drawTube,1}
	ModifyGizmo ModifyObject=chain1,objectType=path,property={ fixedRadius,0.3500}
	ModifyGizmo ModifyObject=chain1,objectType=path,property={ calcNormals,1}
	ModifyGizmo modifyObject=chain1,objectType=Path,property={calcNormals,1}
	AppendToGizmo Path=root:Chains:chain2,name=chain2
	ModifyGizmo ModifyObject=chain2,objectType=path,property={ pathColorType,1}
	ModifyGizmo ModifyObject=chain2,objectType=path,property={ lineWidthType,0}
	ModifyGizmo ModifyObject=chain2,objectType=path,property={ pathColor,1,0,0,1}
	ModifyGizmo ModifyObject=chain2,objectType=path,property={ drawTube,1}
	ModifyGizmo ModifyObject=chain2,objectType=path,property={ fixedRadius,0.3500}
	ModifyGizmo ModifyObject=chain2,objectType=path,property={ calcNormals,1}
	ModifyGizmo modifyObject=chain2,objectType=Path,property={calcNormals,1}
	AppendToGizmo Path=root:Chains:chain3,name=chain3
	ModifyGizmo ModifyObject=chain3,objectType=path,property={ pathColorType,1}
	ModifyGizmo ModifyObject=chain3,objectType=path,property={ lineWidthType,0}
	ModifyGizmo ModifyObject=chain3,objectType=path,property={ pathColor,1,0,0,1}
	ModifyGizmo ModifyObject=chain3,objectType=path,property={ drawTube,1}
	ModifyGizmo ModifyObject=chain3,objectType=path,property={ fixedRadius,0.3500}
	ModifyGizmo ModifyObject=chain3,objectType=path,property={ calcNormals,1}
	ModifyGizmo modifyObject=chain3,objectType=Path,property={calcNormals,1}
	AppendToGizmo Scatter=root:particles_snapshot_SPECIAL,name=SpecialParticles
	ModifyGizmo ModifyObject=SpecialParticles,objectType=scatter,property={ scatterColorType,0}
	ModifyGizmo ModifyObject=SpecialParticles,objectType=scatter,property={ markerType,0}
	ModifyGizmo ModifyObject=SpecialParticles,objectType=scatter,property={ sizeType,0}
	ModifyGizmo ModifyObject=SpecialParticles,objectType=scatter,property={ rotationType,0}
	ModifyGizmo ModifyObject=SpecialParticles,objectType=scatter,property={ Shape,2}
	ModifyGizmo ModifyObject=SpecialParticles,objectType=scatter,property={ size,0.15}
	ModifyGizmo ModifyObject=SpecialParticles,objectType=scatter,property={ color,0,0,1,1}
	AppendToGizmo attribute diffuse={1,1,1,1,1032},name=diffuse0
	AppendToGizmo attribute shininess={42,42},name=shininess0
	AppendToGizmo attribute blendFunc={770,771},name=blendFunc0
	ModifyGizmo setDisplayList=0, object=light0
	ModifyGizmo setDisplayList=1, object=LeftDL
	ModifyGizmo setDisplayList=2, attribute=diffuse0
	ModifyGizmo setDisplayList=3, object=axes0
	ModifyGizmo setDisplayList=4, object=RightDL
	ModifyGizmo setDisplayList=5, object=chain0
	ModifyGizmo setDisplayList=6, object=chain1
	ModifyGizmo setDisplayList=7, object=SpecialParticles
	ModifyGizmo setDisplayList=8, object=chain2
	ModifyGizmo setDisplayList=9, object=chain3
	ModifyGizmo currentGroupObject=""
	ModifyGizmo zoomFactor=0.607000
	ModifyGizmo endRecMacro
	Execute/Q/Z "SetWindow kwTopWin sizeLimit={45.75,234,inf,inf}" // sizeLimit requires Igor 7 or later
EndMacro

Window Panel0() : Panel

	Variable/G STEPSIZE
	variable/G TWOPARTS
	variable/G MixRem
	Variable/G FRAPSTEP
	variable/G length
	
	if(exists("filepath"))
		string/G filepath = root:filepath
	else
		string/G filepath = ""
	endif
	string/G file  = StringFromList(0,StringFromList(ItemsinList(filepath,":")-1,filepath, ":"),".")
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(134,112,600,350)
	Button btnRead,pos={87.00,53.00},size={133.00,70.00},proc=ButtonProc,title="Read"
	Button btnRead,fColor=(52428,52425,1)
	Button btnPlay,pos={250,63},size={50,50},proc=ButtonProc,title="Play"
	Button btnAR, pos={87.00, 140}, size={133,25},proc=ButtonProc,title="Aspect Ratio"
	Button btnMix, pos={87.00, 170}, size={133,25},proc=ButtonProc,title="Mixing Index"
	CheckBox chkFusion,pos={250,25}, variable=Twoparts, title="Fusion"
	CheckBox chkRemodeller,pos={250,170}, variable=MixRem, value=0, title="Remodeller"
	Button btnResetGizmo, pos={87.00, 200}, size={133,25},proc=ButtonProc,title="Reset Windows"
	Titlebox ttlFilename, pos={250,200},size={113,25},variable=file
	Setvariable setStepsize, pos={250,140},size={113,25}, variable=STEPSIZE, title="Step size MI"
	Setvariable setFRAPstep, pos={85,25},size={113,25},variable=FRAPSTEP, title="FRAP at step"
	Slider sldTimestep,pos={1.00,1.00},size={450.00,16.00},proc=particlesSlider
	Slider sldTimestep,limits={0,100,1},value= 0,side= 0,vert= 0
	if(waveexists(root:particles))
		Slider sldTimestep,limits={0,dimsize(particles,2),1},value= 31,side= 0,vert= 0	
	endif
	ValDisplay vdFrame,pos={1.00,22.00},size={50.00,13.00}
	ValDisplay vdFrame,limits={0,0,0},barmisc={0,1000},value= _NUM:0
	Button btnShowAR, pos={30,140}, size={50,20}, proc=ButtonProc,title="Show"
	Button btnShowMI,pos={30,170},size={50,25},proc=ButtonProc,title="Show"
	ValDisplay vdLength,pos={320,22},size={90,13},value=length, title="length (s):"
EndMacro

Function ButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	NVAR MixRem, STEPSIZE, length
	wave particles
	variable n_frames = dimsize(particles, 2)
	switch( ba.eventCode )
		case 2: // mouse up
			strswitch(ba.ctrlName)
				case "btnRead":
					readParticlesTrajectory("")
						MakeChains()
						clean_windows()
						MakeParticleWin()
						SSQDistFromPos()
						ClearRemodellerMixing()
						RadialDensity(n_frames=n_frames-1)
						UpdateSnapshot(1)
						getREMinDroplet(n_frames=n_frames-1)
						length = dimsize(particles, 2) * dimdelta(particles, 2)
					break
				case "btnPlay":
					DoWindow/F ParticleWin 
					DoWindow/F Gizmo0
					Play()
					break
				case "btnMix":
					if(MixRem)
						RemodellerMixingTraj(stepsize=STEPSIZE)
					else
						FindMixingTraj(stepsize=STEPSIZE)
					endif
					break
				case "btnAR":
					xspread_trajectory()
					break
				case "btnResetGizmo":
					KillWindow/Z Gizmo0
					Execute "Gizmo0()"
					KillWindow/Z Panel0
					Execute "Panel0()"
					KillWindow/Z giz
					break
					break
				case "btnShowMI":
					makeMIWindow()
					break
				case "btnShowAR":
					makeARWindow()
				endswitch
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End




//Additional marker for fusion: Sum squared distances from final position
Function SSQDistFromPos()
	Wave particles
	
	Make/O/N=3 reference = {0,0,0} //Calculate reference from this point
	
	Make/O/N=(dimsize(particles,2)) SSQDist=NaN
	Setscale/P x, dimoffset(particles,2), dimdelta(particles,2), SSQdist
	
	NVAR NSPPARTICLES
	//Variable Nparticles = dimsize(particles,0)
	Variable Nparticles = dimsize(particles,0) -NSPPARTICLES
	
	Variable t
	for(t=0;t<dimsize(particles,2);t+=1)
		Duplicate/O/FREE/R=[0,Nparticles-1][][t] particles, snapshot
		Matrixop/O/FREE ssq = sumsqr(snapshot-rowrepeat(reference,Nparticles))
		SSQDist[t] = sqrt(ssq[0])
		//doupdate
	endfor
End







//find holes in chromatin (accessible volume)
Function FindHoles([Variable VERBOSE])
	VERBOSE = paramisdefault(VERBOSE) ? 0 : VERBOSE

	Variable nside = 100
	Variable boxhalfWidth = 200

	Variable radiusSq = 14^2 //square radius of a nucleosome
	if(VERBOSE)
		print "VolNucleosome (nm^3):", 4/3*pi*(sqrt(radiusSq))^3
	endif
	
	Wave coords = particles_snapshot
	NVAR NSPPARTICLES
		
	Variable npnts = dimsize(coords,0)-NSPPARTICLES
	
	if(VERBOSE)
		print "MINSIDELENGTH (nm):", 2*boxhalfwidth/nside
	endif
	Make/O/B/U/N=(nside, nside, nside) grid = 0
	
	setscale/I x, -boxhalfWidth, boxhalfWidth, grid
	setscale/I y, -boxhalfWidth, boxhalfWidth, grid
	setscale/I z, -boxhalfWidth, boxhalfWidth, grid
	
#if 0
	Variable i, j, k, n
	for(i=0;i<nside;i+=1)
		for(j=0;j<nside;j+=1)
			for(k=0;k<nside;k+=1)
				Variable cx = indexToScale(grid, i, 0)
				Variable cy = indexToScale(grid, j, 1)
				Variable cz = indexToScale(grid, k, 2)
				
				for(n=0;n<npnts;n+=1)
					if((cx-coords[n][0])^2+(cy-coords[n][1])^2+(cz-coords[n][2])^2 < radiusSq)
						grid[i][j][k] = 255 //chromatin
						break
					endif
				endfor
			endfor
		endfor
	endfor
#else
	Variable n
	for(n=0;n<npnts;n+=1)
		variable boxX_L = scaleToIndex(grid, coords[n][0]-sqrt(radiusSq)-1, 0)
		variable boxX_R = scaleToIndex(grid, coords[n][0]+sqrt(radiusSq)+1, 0)
		variable boxY_L = scaleToIndex(grid, coords[n][1]-sqrt(radiusSq)-1, 1)
		variable boxY_R = scaleToIndex(grid, coords[n][1]+sqrt(radiusSq)+1, 1)
		variable boxZ_L = scaleToIndex(grid, coords[n][2]-sqrt(radiusSq)-1, 2)
		variable boxZ_R = scaleToIndex(grid, coords[n][2]+sqrt(radiusSq)+1, 2)
		
		boxX_L = boxX_L < 0 ? 0 : boxX_L
		boxX_R = boxX_R < 0 ? 0 : boxX_R
		boxY_L = boxY_L < 0 ? 0 : boxY_L
		boxY_R = boxY_R < 0 ? 0 : boxY_R
		boxZ_L = boxZ_L < 0 ? 0 : boxZ_L
		boxZ_R = boxZ_R < 0 ? 0 : boxZ_R

		boxX_L = boxX_L > dimsize(grid,0)-1 ? dimsize(grid,0)-1 : boxX_L
		boxX_R = boxX_R > dimsize(grid,0)-1 ? dimsize(grid,0)-1 : boxX_R
		boxY_L = boxY_L > dimsize(grid,1)-1 ? dimsize(grid,1)-1 : boxY_L
		boxY_R = boxY_R > dimsize(grid,1)-1 ? dimsize(grid,1)-1 : boxY_R
		boxZ_L = boxZ_L > dimsize(grid,2)-1 ? dimsize(grid,2)-1 : boxZ_L
		boxZ_R = boxZ_R > dimsize(grid,2)-1 ? dimsize(grid,2)-1 : boxZ_R
		
		grid[boxX_L,boxX_R][boxY_L,boxY_R][boxZ_L,boxZ_R] = (x-coords[n][0])^2+(y-coords[n][1])^2+(z-coords[n][2])^2 < radiusSq ? 255 : grid
	endfor
#endif
	
	Wavestats/Q grid
	Variable volChromatin = v_sum/255
	if(VERBOSE)
		print "VolChromatin (nm^3):", volChromatin, "which is", (volChromatin)/(nside*nside*nside)*100, "%"
		//Meaningless: print "Nucleosome density (M): ", npnts/volChromatin * (1e-1/1e-9)^3 / 6.022e23
	endif
	
	//find contiguous stretches in chromatin and holes in between
	imageanalyzeparticles stats grid
	
	Wave grid, m_3dParticleInfo
	sortColumns/R/KNDX=(finddimlabel(M_3DParticleInfo,1,"area")) sortwaves=M_3DParticleInfo
	
	Duplicate/O grid, internalFill
	internalFill = 255
	
	Variable totVolumeHoles=0
	
	Variable i
	for(i=0;i<dimsize(m_3dparticleInfo,0);i+=1)
		if(i==0)
			continue //usually(!) biggest by far is the outside
		endif
	
		ImageSeedFill/B=255 min=0,max=1,seedP=m_3dParticleInfo[i][%seedX],seedQ=m_3dParticleInfo[i][%seedY],seedR=m_3dParticleInfo[i][%seedZ],target=0,srcWave=grid
		Wave m_seedfill
		Duplicate/O internalFill, tmp
		MatrixOp/O internalFill = bitand(M_seedFill, tmp)
		//Multithread internalFill = m_seedFill == 0 ? 0 : internalFill
		
		totVolumeHoles += m_3dParticleInfo[i][%volume]
	endfor
	Copyscales/P grid, internalFill
	
	if(VERBOSE)
		print "VolHoles (nm^3)", totVolumeHoles
		print "VolHoles/VolChromatin", totVolumeHoles/volChromatin*100, "%"
	endif
	
	return totVolumeHoles/volChromatin*100
End





Function radialDensity([variable n_frames]) //of all chained particles
	
	n_frames = paramIsDefault(n_frames) ? 100 : n_frames
	Wave particles
	
	NVAR NSPPARTICLES
	Variable N_noremodeler = dimsize(particles,0)-NSPPARTICLES
	print "Radial density for particles 0 to", N_noremodeler
	
	//Variable fromT = dimsize(particles,2)-n_frames //from this frame on
	variable fromT=0
	Make/O/N=100 radialHist = 0, radialHistRemodeler = 0
	Variable eps = 1e-9
	setscale/I x, 0, 300, radialHist, radialHistRemodeler

	Duplicate/O/FREE/R=[,N_noremodeler-1][][] particles, particles_noremodeler
	Duplicate/O/FREE/R=[N_noremodeler,dimsize(particles,0)-1][][] particles, particles_remodeler

	Make/N=(0,3)/O CM
	Variable N = 0
	Variable t
	for(t=fromT;t<dimsize(particles,2);t+=1)
		Duplicate/FREE/O/R=[][][t] particles_noremodeler, particles_noremodeler_t
		MatrixOp/O/FREE COMnoremodeler = averagecols(particles_noremodeler_t)
		Redimension/N=(dimsize(CM,0)+1,3) CM
		CM[dimsize(CM,0)-1][0] = COMnoremodeler[0]
		CM[dimsize(CM,0)-1][1] = COMnoremodeler[1]
		CM[dimsize(CM,0)-1][2] = COMnoremodeler[2]
		
		MatrixOp/O/FREE tmpnoremodeler = sqrt(sumrows((particles_noremodeler_t - rowrepeat(COMnoremodeler,N_noremodeler)) *  (particles_noremodeler_t - rowrepeat(COMnoremodeler,N_noremodeler))))
		
		Histogram/N/A/B=2 tmpnoremodeler, radialHist
		N+= 1
	endfor

	//radialHist /= (x+deltax(radialHist))^2
	radialHist /= (4*pi*(x+eps)^2*deltax(radialHist))
	Wave W_sqrtN
	Duplicate/O w_sqrtN, radialHist_err
	//radialHist_err /= (x+deltax(radialHist))^2
	radialHist_err /= (4*pi*(x+eps)^2*deltax(radialHist))

	radialHist /= N
	radialHist_err /= N

	//same for remodeler
	N=0
	for(t=fromT;t<dimsize(particles,2);t+=1)
		Duplicate/FREE/O/R=[][][t] particles_remodeler,   particles_remodeler_t
		MatrixOp/O/FREE COMremodeler = averagecols(particles_remodeler_t)
		
		MatrixOp/O/FREE tmpremodeler = sqrt(sumrows((particles_remodeler_t - rowrepeat(COMremodeler,NSPPARTICLES)) *  (particles_remodeler_t - rowrepeat(COMremodeler,NSPPARTICLES))))
		
		Histogram/N/A/B=2 tmpremodeler, radialHistRemodeler
		N+=1
	endfor

	//radialHistRemodeler /= (x+deltax(radialHistRemodeler))^2
	radialHistRemodeler /= (4*pi*(x+eps)^2*deltax(radialHistRemodeler))
	Wave W_sqrtN
	Duplicate/O w_sqrtN, radialHistRemodeler_err
	//radialHistRemodeler_err /= (x+deltax(radialHistRemodeler))^2
	radialHistRemodeler_err /= (4*pi*(x+eps)^2*deltax(radialHistRemodeler_err))

	radialHistRemodeler /= N
	radialHistRemodeler_err /= N


	//Get concentration
	CurveFit/Q/M=2/W=0 Sigmoid, radialHist[9,99]/D
	Wave w_coef
	
	print "Nucleosome concentration in core (µM):", w_coef[0]/(1e-8)^3/6.022e23 * 1e6
	print "Droplet radius (nm):", w_coef[2]
	Variable/G radius = w_coef[2]
End


function clean_windows()
	string windows = winlist("ParticleWin*",";","")
	variable i
	for(i=0;i<itemsinlist(windows);i++)
		string windowname = stringfromList(i, windows, ";")
		killwindow $windowname
	endfor
end