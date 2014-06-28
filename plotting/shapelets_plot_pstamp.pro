pro shapelets_plot_pstamp, pstamp,              $
                           OFRAME=oframe,       $
                           MASK=mask,	        $
                           NOISE=noise,         $
                           TITLE=title,         $
                           LINESTYLE=linestyle, $
			   CRAN=cran,	        $
                           CSIZE=csize

;$Id: shapelets_plot_pstamp.pro, v1.0$
;
; Copyright © 2005 Richard Massey and Alexandre Refregier.
;
; This file is a part of the Shapelets analysis code.
; www.ast.cam.ac.uk/~rjm/shapelets/
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
;     Jan 05 - idl2 compatability ([]s, etc) implemented by Richard Massey
;     Jan 05 - csize keyword added by AR
;     Mar 02 - Written by A. Refregier
;
; PURPOSE: plot the image and sextractor parameters for an object
; INPUT: pstamp: object structure from shapelets_image2pstamp.pro
; OPTIONAL INPUT: oframe: label the plot with the object frame coordinates
;                 noise: plot noise or S/N for noise=1,2 respectively
;                 mask: plot contours of masked out regions
;                 csize: color bar size
; OUTPUT: plot
;-

COMPILE_OPT idl2

; declarations
si=size(pstamp.image)
n_pstamp=[si[1],si[2]]   ; image size

; select function to plot
if not keyword_set(title) then title=""
if not keyword_set(noise) then begin
  f=pstamp.image
endif else begin
  if noise eq 1 then begin    
    f=pstamp.noise
    gd=where(f gt 0.)
    f[gd]=1./sqrt(f[gd])
    title=title+" rms noise"
  endif else begin
    f=pstamp.image*sqrt(pstamp.noise)
    title=title+" S/N"
  endelse
endelse

; plot image and SExtractor parameters
if not keyword_set(oframe) then begin       ; image frame
  shapelets_plot_image,f,/col,cran=cran,csize=csize,$
    frame=[pstamp.im_ran[0],pstamp.im_ran[1]+1,pstamp.im_ran[2],pstamp.im_ran[3]+1],$
    xtitle="!8x!6 [pixels]",ytitle="!8y!6 [pixels]",title=title
  plt_ellipse,pstamp.x,pstamp.y,$
    pstamp.a,pstamp.b,(pstamp.theta+90.)/!radeg,linestyle=linestyle
  oplot,[pstamp.x],[pstamp.y],psym=1
  if keyword_set(mask) then begin
    id=pstamp.seg[uniq(pstamp.seg,sort(pstamp.seg))]  ; find all unique object ids in postage stamp
    n_id=n_elements(id)
    for i=0,n_id-1 do begin
      if id[i] gt 0 then begin
        segi=pstamp.seg
        bad=where(segi ne id[i])
        segi[bad]=0.
        contour,abs(segi),indgen(n_pstamp[0])+pstamp.im_ran[0],$
          indgen(n_pstamp[1])+pstamp.im_ran[2],/over,$
        level=[abs(id[i])],c_charsize=1.5,/follow,c_annotation=strtrim(string(id[i]),2)
        ;c_annotation=strtrim(string(ids),2),c_charsize=1.,c_labels=replicate(1.,n_gdids)
      endif
    endfor
  endif
endif else begin                            ; object frame
  shapelets_plot_image,f,/frame,/col,cran=cran,csize=csize,$
    xtitle="x (pixels)",ytitle="y (pixels)",title=title  
  plt_ellipse,pstamp.xo,pstamp.yo,$
    pstamp.a,pstamp.b,(pstamp.theta+90.)/!radeg,linestyle=linestyle
  oplot,[pstamp.xo],[pstamp.yo],psym=1
  if keyword_set(mask) then $ 
      contour,pstamp.noise,level=[.1],/over
endelse

;; SExtractor min/max pixels
;plt_rectangle,[pstamp.xmin+pstamp.xo-pstamp.x,pstamp.ymin+pstamp.yo-pstamp.y], $
;		[pstamp.xmax+pstamp.xo-pstamp.x,pstamp.ymax+pstamp.yo-pstamp.y],linestyle=2 

end
