#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


function getDSmap(alpha)
variable alpha //threshold for p values to include
wave ROIs, pValues, VectorSum


variable nROIs = Dimsize(pValues, 0)
variable nX = Dimsize(ROIs,0)
variable nY = Dimsize(ROIs,1)


make/o/n=(nX,nY) DSmap = -1 //-1 for pixels not included in analysis because low SD and quality
variable rr, xx, yy

for (yy=0;yy<nY;yy+=1)
	for (xx=0;xx<nX;xx+=1)
		if (ROIs[xx][yy]< 0)
			rr = (ROIs[xx][yy]*-1)-1
			if (pValues[rr]<alpha)
				DSmap[xx][yy]=VectorSum[rr] //Vector sum for DS pixels
			else
				DSmap[xx][yy]=0 //zero for high quality pixels that aren't DS
			endif
		endif
	endfor
endfor
	
end