pro shapelets_plot_image, image,              $
                          FRAME=frame,        $
                          COLBAR=colbar,      $
                          CBAR=cbar,          $
                          CRANGE=crange,      $
                          CINVERT=cinvert,    $
                          CLOG=clog,          $
                          CTITLE=ctitle,      $
                          INVERSE=inverse,    $
                          INVERT=invert,      $
                          SCALABLE=scalable,  $
                          NOERASE=noerase,    $
                          CSIZE=csize,        $
                          TRUE=true,          $
                          ISOTROPIC=isotropic, _ref_extra = ex

;$Id: shapelets_plot_image.pro, v2$
;
; Copyright © 2005 Richard Massey and Alexandre Refregier.
;
; This file is a part of the Shapelets analysis code.
; www.astro.caltech.edu/~rjm/shapelets/
;
; The Shapelets code is free software; you can redistribute it and/or
; modify it under the terms of the GNU General Public Licence as published
; by the Free Software Foundation; either version 2 of the Licence, or
; (at your option) any later version.
;
; The Shapelets code is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
; GNU General Public Licence for more details.
;
; You should have received a copy of the GNU General Public Licence
; along with the Shapelets code; if not, write to the Free Software
; Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
;
;+
; NAME:
;      SHAPELETS_PLOT_IMAGE
;
; PURPOSE:
;      Draws a 2D pixellated image.
;
; CATEGORY:
;      Plotting.
;
; PURPOSE:
;      Plot an image scaled to the current coordinate frame.
;      Further plotting (such as contours) can be performed over the
;      resulting plot. Optionally, an annotated color bar can be drawn at
;      the top. The 'scalable' keyword must be invoked when outputing to a
;      postscript file.
;
; NOTES:
;      !p.region=0 can be used to reset the window after using this
;      routine
;
; INPUTS:
;      IMAGE  -    2D monochrome image array
;               OR 3D full colour image array (eg [[[r]],[[g]],[[b]]],true=3)
;               OR Shapelets image/pstamp structure
;
; OPTIONAL INPUTS:
;      Accepts usual graphics inputs, including [X,Y]TITLE, CHARSIZE, etc. 
;      CTITLE	- Title above colour bar.
;      CRANGE	- Range of a for colors.
;                 DEFAULT: [min(a),max(a)]
;      CSIZE	- Vertical size of the colour bar [0-1].
;                 DEFAULT: 0.14
;
; KEYWORD PARAMETERS:
;      Accepts usual graphics keywords, including ISOTROPIC, SCALABLE, etc. 
;      ISOTROPIC - Force scales on x and y axes to be identical.
;      CLOG      - Make colour scale logarithmic.
;      CBAR      - Draw color bar. Also works with colbar.
;      CINVERT   - Invert the color coding.
;      FRAME     - Draw a frame with pixel index limits;
;                  to supress frame but keep image within it, set FRAME=-1
;                  to change limits set frame=[x0,x1,y0,y1].
;                  The latter effect can also be achieved with [X,Y]RANGE.
;      NOERASE   - Supress erasing whatever was previously in the window.
;      SCALABLE  - Force the use of scalable pixels.
;                  DEFAULT: nonscalable used with xterm or win device, 
;                           scalable used with postscript device.
;
; OUTPUTS:
;      Plot of the scaled image with, optionally, a coordinate frame
;      and a color bar.
;
; MODIFICATION HISTORY:
;      Mar 06 - CINVERT keyword added by RM.
;      Aug 05 - ISOTROPIC keyword perfected by RM.
;      Apr 05 - Plotting of 1D arrays and logarithmic color scale added by RM.
;      Mar 05 - Input generalised to include shapelet image structures by RM.
;      Jul 04 - Automatic detection of device type (and adoption of scalable
;               pixels where relevant) by RM.
;      Nov 02 - ISOTROPIC keyword added by Richard Massey. Still temperamental!
;      Aug 00 - Plotting over the full plotting region allowed by AR.
;    3 Sep 97 - Active changes for the color bar rangeadded by AR.
;   14 Oct 97 - Plotting of a color bar and flexible frame labels added by AR.
;    7 Jul 97 - Written by Alexandre Refregier.
;-

;
; Memorise how graphics device was initially set up (to be returned to at end)
;
pregion=!p.region
if keyword_set(inverse) or keyword_set(invert) then cinvert=1B

;
; Parse input image data
;
if n_params() lt 1 then message,"Usage: shapelets_plot_image,image"
if shapelets_structure_type(image,message=message,/silent) then begin
  if image.type eq "image" or image.type eq "pstamp" then a=float(image.image) else message,"You need to input an image!"
endif else begin
  a=float(image)
endelse


;
; Hack to get around the fact that tv can't handle 1D arrays
;
n_dimensions=size(a,/N_DIMENSIONS)
if n_dimensions eq 1 then begin
  a=reform(a,n_elements(a),1)#[1,1]
  n_dimensions=2
endif 


;
; Cope with both monochrome and [r,g,b] full colour images
;
if n_dimensions eq 2 then begin
  if keyword_set(xrange) then a=a[xrange[0]:xrange[1],*]
  if keyword_set(yrange) then a=a[*,yrange[0]:yrange[1]]
  if n_elements(uniq(a)) eq 1 then message,"WARNING: All pixels in the image have the same value!",/info
  nx=(size(a,/DIMENSIONS))[0]
  ny=(size(a,/DIMENSIONS))[1]
endif else if n_dimensions eq 3 then begin
  if not keyword_set(true) then begin
    true_new=where(size(a,/DIMENSIONS) eq 3,n_three_elements)
    if n_three_elements ne 1 then $
      message,"Keyword TRUE is required to plot a colour image!"
  endif else true_new=true
  data_dimensions=where((indgen(3)+1) ne true_new)
  nx=(size(a,/DIMENSIONS))[data_dimensions[0]]
  ny=(size(a,/DIMENSIONS))[data_dimensions[1]]
endif else message,"Input variable is "+strtrim(string(n_dimensions),2)+"D!"


;
; Figure out what kind of plot we are making
; (if we're plotting to postscript, the pixels should be scalable)
;
if n_elements(scalable) eq 0 then $
  if strupcase(!d.name) eq "PS" then scalable=1 else scalable=0


;
; Decide colour range
;
if not keyword_set(crange) then begin
  if keyword_set(clog) then begin
    orders_of_magnitude=4
    crange_lower=(alog10(abs(median(a)))-orders_of_magnitude/2)>min(alog10(abs(a)))
    crange_upper=(alog10(abs(median(a)))+orders_of_magnitude/2)<max(alog10(abs(a)))
    crange=[crange_lower,crange_upper]
  endif else crange=[min(a),max(a)]
endif


;
; Convert image to a logarithmic colour scale
;
if keyword_set(clog) then a=alog10(abs(a))

; Map image to color indices
if not keyword_set(crange) then begin
;  if keyword_set(inverse) then b=-a else $
  b=a
  b_ran=[min(b),max(b)]
  b=(b-b_ran[0])/(b_ran[1]-b_ran[0])*(!d.table_size-1)    ; this is equivalent to tvscl
endif else begin
;  if keyword_set(inverse) then begin    
;    b=(1.-((crange[0]>a<crange[1])-crange[0])/(crange[1]-crange[0]))*(!d.table_size-1)
;  endif else begin
    b=((crange[0]>a<crange[1])-crange[0])/(crange[1]-crange[0])*(!d.table_size-1) 
;  endelse
endelse


;
; Draw color bar
;
if keyword_set(colbar) or keyword_set(cbar) then begin
  ; Decide how big it's going to be
  if not keyword_set(csize) then csize=.14
  ; Reserve upper portion of the window for the color bar
  !p.region=[0.,1.-csize,1.,1.]
  ; Draw colour bar
;  shapelets_plot_image,findgen(8>!d.n_colors<512)*(-1)^keyword_set(inverse),$
;    frame=-1,title=ctitle,noerase=noerase
  shapelets_plot_image,findgen(8>!d.n_colors<512),$
    frame=-1,title=ctitle,noerase=noerase,cinvert=cinvert
  ; Draw a box around it
  if keyword_set(clog) then ctickformat='("10!u",F3.1,"!n")'
  plot,[0],[0],/nodata,/noerase,xrange=crange,/xstyle,title=ctitle,$
    xminor=1,yminor=1,yticks=1,ytickname=[" "," "],xtickformat=ctickformat
  ; Leave the rest of the window for the main image
  !p.region=[0.,0.,1.,1.-csize]
  noerase=1
endif


;
; Scale image to the frame size
;
; Set up the coordinate frame
plot,[0],[0],/nodata,xstyle=12,ystyle=12,noerase=noerase
if not keyword_set(scalable) then begin

  ;
  ; Non-scalable pixel case (to be used for xterm device)
  ;
  if keyword_set(frame) then begin
    px=!x.window[0]*!d.x_vsize                       ; position of window in device pixels
    py=!y.window[0]*!d.y_vsize
    sx=abs((!x.window[1]-!x.window[0])*!d.x_vsize+1) ; desired size of image in pixels
    sy=abs((!y.window[1]-!y.window[0])*!d.y_vsize+1) 
  endif else begin
    px=!x.region[0]*!d.x_vsize                       ; position of window in device pixels
    py=!y.region[0]*!d.y_vsize    
    sx=abs((!x.region[1]-!x.region[0])*!d.x_vsize+1) ; desired size of image in pixels
    sy=abs((!y.region[1]-!y.region[0])*!d.y_vsize+1) 
  endelse
  if keyword_set(isotropic) then begin
    sx=((sx/nx)<(sy/ny))*nx
    sy=((sx/nx)<(sy/ny))*ny
    frame_position=[!x.window[0],!y.window[0],$
  	               (float(sx)*float(!d.y_vsize)/float(sy)/float(!d.x_vsize)*$
  	                float(!y.window[1]-!y.window[0]+1/!d.y_vsize)+!x.window[0]-1/!d.x_vsize)<!x.window[1],$
  	               (float(sy)*float(!d.x_vsize)/float(sx)/float(!d.y_vsize)*$
  	                float(!x.window[1]-!x.window[0]+1/!d.x_vsize)+!y.window[0]-1/!d.y_vsize)<!y.window[1]]
  endif else begin
    frame_position=[!x.window[0],!y.window[0],!x.window[1],!y.window[1]]
  endelse
  
  ;
  ; Repixellate image
  ;
  if n_dimensions le 2 then begin
    b=congrid(b,sx,sy)
  endif else begin
    case true of
      1: b=[[[congrid(reform(b[0,*,*],nx,ny),sx,sy)]],$
            [[congrid(reform(b[1,*,*],nx,ny),sx,sy)]],$
            [[congrid(reform(b[2,*,*],nx,ny),sx,sy)]]]
      2: b=[[[congrid(reform(b[*,0,*],nx,ny),sx,sy)]],$
            [[congrid(reform(b[*,1,*],nx,ny),sx,sy)]],$
            [[congrid(reform(b[*,2,*],nx,ny),sx,sy)]]]
      3: b=[[[congrid(reform(b[*,*,0],nx,ny),sx,sy)]],$
            [[congrid(reform(b[*,*,1],nx,ny),sx,sy)]],$
            [[congrid(reform(b[*,*,2],nx,ny),sx,sy)]]]
    endcase
    true_new=3
  endelse

endif else begin

  ;
  ; Scalable pixel case (to be used for ps device)
  ;
  if keyword_set(frame) then begin
    sx=!x.window[1]-!x.window[0] ; desired size of image in pixels
    sy=!y.window[1]-!y.window[0]
    px=!x.window[0]              ; position of window in device pixels
    py=!y.window[0]
  endif else begin
    sx=!x.region[1]-!x.region[0] ; desired size of image in pixels
    sy=!y.region[1]-!y.region[0]
    px=!x.region[0]              ; position of window in device pixels
    py=!y.region[0]
  endelse   
  if keyword_set(isotropic) then begin
    sx=((sx/nx)<(sy/ny))*nx
    sy=((sx/nx)<(sy/ny))*ny
    frame_position=[!x.window[0],!y.window[0],$
        	    (float(sx)/float(sy)*$
        	     float(!y.window[1]-!y.window[0])+!x.window[0])<!x.window[1],$
        	    (float(sy)/float(sx)*$
        	     float(!x.window[1]-!x.window[0])+!y.window[0])<!y.window[1]]
  endif else begin
    frame_position=[!x.window[0],!y.window[0],!x.window[1],!y.window[1]]
  endelse

endelse	


;
; Flip colour scale
;
if keyword_set(cinvert) then b=(!d.table_size-1)-b


;
; Plot image
;
tv,b,px,py,xsize=sx,ysize=sy,true=true_new,$
   device=1-keyword_set(scalable),norm=keyword_set(scalable)


;
; Overlay frame
;
if keyword_set(frame) then $
if not (n_elements(frame) eq 1 and frame[0] eq -1) then begin
  if n_elements(frame) eq 4 then begin
    xr=[frame[0],frame[1]] & yr=[frame[2],frame[3]]
  endif else begin
    xr=[0.,nx] & yr=[0.,ny]
  endelse
  plot,[0],[0],/nodata,/noerase,position=frame_position,$
    /xstyle,/ystyle,XRANGE=xr,YRANGE=yr,_extra= ex
endif


;
; Return graphics device to initial conditions
;
!p.region=pregion

end


