pro plt_ellipse, x, y, a, b, pa, xe, ye, $
                 NPOINTS=npoints,        $
		 FILL=fill,              $
		 NOPLOT=noplot,          $
		 XSCALE=xs,              $
		 YSCALE=ys,              $
                 COLOR=color,            $
		 LINESTYLE=linestyle,    $
		 THICK=thick

; October 96 - written by A. Refregier
;
; PURPOSE: plot one or several ellipse centered at (x,y) and with
; major and minor axes a and b, respectively, and with
; position angle pa
; INPUT: x,y: center coordinates (can be arrays)
;        a,b: major, minor axis
;	 pa: position angle (counter-clockwise from north) (rad)
; OPTIONAL INPUT: npoints: number of points used to draw the circle
;                 (default=30)
;                 fill: fill the ellipse
;                 noplot: supress plotting (just compute ellipse points)
;                 xscale, yscale: scale for the x,y axis; this is useful
;                                 when the coordinate frame is not square
; OUTPUT: circle plotted on the current plot
; OPTIONAL OUTPUT: xe,ye: coordinates of the ellipse points

; find the number of ellipses to draw
nell=n_elements(x)

; set default parameters
if not keyword_set(xs) then xs=1.
if not keyword_set(ys) then ys=1.
 if not keyword_set(linestyle) then linestyle=0

; compute point coordinates and plot each ellipse
for i=0,nell-1 do begin
  if not keyword_set(npoints) then np=replicate(30,nell) else np=npoints
  phi=[xgen(0.,2*!pi,npoints=np(i)-1),0.]
  xe=(a(i)*cos(phi)*cos(pa(i)+!pi*.5)-b(i)*sin(phi)*sin(pa(i)+!pi*.5))*xs+x(i)
  ye=(a(i)*cos(phi)*sin(pa(i)+!pi*.5)+b(i)*sin(phi)*cos(pa(i)+!pi*.5))*ys+y(i)
  if not keyword_set(noplot) then $
    if keyword_set(fill) then begin
      if keyword_set(color) then polyfill,xe,ye,color=color $
         else polyfill,xe,ye
    endif else begin
      if keyword_set(color) then oplot,xe,ye,color=color,linestyle=linestyle,thick=thick $
         else oplot,xe,ye,linestyle=linestyle,thick=thick,psym=-3
    endelse
endfor

end

