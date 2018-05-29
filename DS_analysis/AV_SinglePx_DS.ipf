#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

function SinglePx_DS()
//find the DSI of single pixels. Based on "GluSnFr_DS" from KF and on Baden et al. 2016

/// FLAGS
variable save_stuff = 0

wave wParamsNum
variable ZoomFactor = wParamsNum[30]
wave Conditions, Conditions_Offset
variable nConditions = Dimsize(Conditions, 0)
variable nIterations = 10000
variable Alpha = 0.05

/// IMPORTING DATA
string traces_name = "Traces0_znorm" // traces
string Roi_name = "ROIs" // rois
string triggerchannel_name = "wDataCh2" // triggers
string Averagesname = "Averages0" //average trace containing average across loops
duplicate /o $traces_name traces
duplicate /o $Roi_name Roi
duplicate /o $triggerchannel_name triggerchannel
duplicate /o $Averagesname Responses_average

/// GET DATA DIMENSIONS
variable nF = Dimsize(traces,0)
variable nRois = Dimsize(traces,1)
variable nLines = Dimsize(triggerchannel,1) // nY
variable nX = Dimsize(triggerchannel,0)
variable frame_duration = (nLines*2)/1000 // in s 


variable ff,rr,ll,tt, cc, ii

/// SNIPPET VARIABLES
wave Triggertimes_Frame
wave Triggertimes
string createdSnippets = "Snippets0"
duplicate/o $createdSnippets Responses
variable nTriggers = Dimsize(Triggertimes,0)
make/o/n=(nTriggers) SnippetDurations_frames = Triggertimes_Frame[p+1]-Triggertimes_Frame[p]
make/o/n=(nConditions) SnippetDuration_SingleAxis = SnippetDurations_frames
wavestats/q SnippetDurations_frames
variable SnippetLength = V_min //some snippets are different durations, this takes the minimum.
variable nLoops = nTriggers/nConditions
variable snippet_duration_pp = Triggertimes_Frame[nConditions]-Triggertimes_Frame[0]//again, in frames


//make/o/n=(SnippetDuration_SingleAxis, nRois) Average = 0 //wave with single trials

//for (rr=0;rr<nRois;rr+=1)
//	for (ll=0;ll<nConditions;ll+=1)
//		Average[][rr]+=Responses_average[p + SnippetDuration_SingleAxis*ll][rr]
//	endfor
//endfor

//Average/=nConditions

///// CALCULATE DSi AND DSp ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

make/o/n=(SnippetLength, nConditions, nRois) ResponseMaps = 0
make/o/n=(SnippetLength, nConditions, nRois) ResponseMaps_SVD = 0
make/o/n=(nConditions, nRois) Responses_Directions = 0
make/o/n=(nRois) VectorSum = 0
make/o/n=(nRois) pValues = 0

make /o/n=(SnippetLength,nConditions, nLoops,nRois) Responses_Trials = 0
make /o/n=(SnippetLength,nConditions,nRois) Responses_Trials_Average = 0

for (rr=0;rr<nRois;rr+=1) // nRois
	for (ll=0; ll<nLoops; ll+=1)
		make/o/n=(snippet_duration_pp) CurrentWave = Responses[p][ll][rr]
		for (cc=0; cc<nConditions; cc+=1)
			Responses_Trials[][cc][ll][rr] = CurrentWave[p + (Conditions[cc])*(SnippetDuration_SingleAxis[cc])]
			Responses_Trials_Average[][cc][rr]+=CurrentWave[p + (Conditions[cc])*SnippetDuration_SingleAxis[cc]]/nLoops
		endfor
	endfor
endfor

variable CountSignificant = 0

for (rr=0;rr<nRois;rr+=1)
	// get SVD maps
	for (cc=0; cc<nConditions; cc+=1)
		ResponseMaps[][cc][rr] = Responses_average[p + (Conditions[cc])*(SnippetDuration_SingleAxis[cc])][rr]
	endfor
	make/o/n=(SnippetLength, nConditions) CurrentMap = ResponseMaps[p][q][rr]
	MatrixSVD CurrentMap
	wave M_U, M_VT, W_W //returns V already transposed
	duplicate/o CurrentMap CurrentMap_SVD
	
	CurrentMap_SVD[][] = M_U[p][0]*M_VT[0][q]*W_W[0]
	ResponseMaps_SVD[][][rr] = CurrentMap_SVD[p][q]
	
	// get vector sum
	make/o/n=(nConditions) Activation = M_VT[0][p]
	Wavestats/Q Activation
	if (V_max > abs(V_min)) // normalize to max
		Activation/=V_max
	else
		Activation/=abs(V_min)
	endif
	
	Responses_Directions[][rr] = Activation[p]
	
	make/o/n=(nConditions) CurrentX = 0
	make/o/n=(nConditions) CurrentY = 0
	
	for (cc=0; cc<nConditions; cc+=1)
		CurrentX[cc] = Activation[cc]*Conditions_Offset[cc][0]
		CurrentY[cc] = Activation[cc]*Conditions_Offset[cc][1]
	endfor
	variable xSum = sum(CurrentX)
	variable ySum = sum(CurrentY)
	VectorSum[rr] = sqrt(xSum^2 + ySum^2)
	
	// get p value
	make/o/n=(nIterations) CurrentWave = 0
	for (ii=0; ii<nIterations; ii+=1)
		make/o/n=(nConditions) W_random
		make/n=(nConditions)/FREE random=enoise(1)
		makeIndex random, W_random
		
		make/o/n=(nConditions) CurrentX = 0
		make/o/n=(nConditions) CurrentY = 0
	
		for (cc=0; cc<nConditions; cc+=1)
			CurrentX[cc] = Activation[W_random[cc]]*Conditions_Offset[cc][0] //scrambled activation to get vector sum
			CurrentY[cc] = Activation[W_random[cc]]*Conditions_Offset[cc][1]
		endfor
		xSum = sum(CurrentX)
		ySum = sum(CurrentY)
		CurrentWave[ii] = sqrt(xSum^2 + ySum^2)
	endfor
	
	Sort CurrentWave, CurrentWave
	FindLevel/Q CurrentWave, VectorSum[rr]
	if (V_flag == 0)
		pValues[rr] = 1 - (V_levelX/nIterations)
	else
		pValues[rr] = 0
	endif
	
	if (pValues[rr] < Alpha && pValues[rr]>0)
		CountSignificant+=1
	endif
	
	print "Done with ROI ", rr, " of ", nRois
endfor

print "n =", CountSignificant, " ROIs significant direction selective at alpha=",Alpha

getDSmap(alpha)

killwaves CurrentMap, CurrentMap_SVD, M_U, M_VT, W_W, Activation, CurrentX, CurrentY, CurrentWave

end