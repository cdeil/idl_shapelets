function inpoly,Px,Py,xl,yl

; Aug 2004 - Written by Alexie Leithauld

;----------------------------------------------------------------
;   inpoly.pro 
;
;   Determines if a point P(px, py) in inside or outside a polygon.
;   The method used is the ray intersection method based on the
;   Jordan curve theorem. Basically, we trace a "ray" (a semi line)
;   from the point P and we calculate the number of edges that it 
;   intersects. If this number is odd, P is inside.
;    
;   (Px,Py)    Coordinates of point
;   xv         List of x coordinates of vertices of polygon
;   yv         List of y coordinates of vertices of polygon
;
;   The x and y lists may be given clockwise or anticlockwise.
;   The algorithm works for both convex and concave polygons.
;----------------------------------------------------------------

;---------------------- Parameters  -----------------------------

N=n_elements(xl)
xv=dblarr(N)
yv=dblarr(N)

for i=0.,N-1 do begin
 xv[i]=xl[i]
 yv[i]=yl[i]
endfor

nc=0                     ;Number of edge crossings
N=n_elements(xv);        ;Number of vertices

; Test input
if (N lt 3) then print,'A polygon must have at least three vertices'
if (n_elements(xv) ne n_elements(yv)) then print, 'Must have same number of X and Y coordinates'

;---------------------- Change coordinate system -----------------

; Place P at the center of the coordinate system.

for i=0,N-1 do begin 
    xv[i]=xv[i]-Px
    yv[i]=yv[i]-Py  
endfor

;---------------------- Calculate crossings ----------------------

; The positive half of the x axis is chosen as the ray
; We must determine how many edges cross the x axis with x>0

for i=0,N-1 do begin

Ax=xv[i]      ; First vertice of edge
Ay=yv[i]

if (i eq N-1) then begin
  Bx=xv[0]
  By=yv[0]
endif else begin
  Bx=xv[i+1]    ; Second vertice of edge
  By=yv[i+1]
endelse

; We define two regions in the plan: R1/ y<0 and R2/ y>=0. Note that
; the case y=0 (vertice on the ray) is included in R2.

if (Ay lt 0) then signA = -1 else signA=+1
if (By lt 0) then signB = -1 else signB=+1

; The edge crosses the ray only if A and B are in different regions.
; If a vertice is only the ray, the number of crossings will still be
; correct.

if( (signA*signB) lt 0) then begin

    ; if Ax>0 and Bx>o then the edge crosses the positive x axis
    if((Ax gt 0) and (Bx gt 0)) then begin
    nc=nc+1

    ; Otherwise (end points are in diagonally opposite quadrants)
    ; we must calculate the intersection
    endif else begin
	x=Ax-(Ay*(Bx-Ax))/(By-Ay)
        if (x gt 0) then nc=nc+1
    endelse
endif
endfor  ;i

;if inside then uneven
;if outside then even

nc=nc mod 2
return,nc


end
