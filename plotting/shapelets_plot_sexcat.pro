pro shapelets_plot_sexcat, sexcat, image, _ref_extra = ex

;$Id: shapelets_plot_sexcat.pro, v2$
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
;      SHAPELETS_PLOT_SEXCAT
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Displays an image (if input) and overlays the positions of all
;      objects found on it by SExtractor. If no image structure is given,
;      it just plots the positions of objects.
; 
; INPUTS:
;      SEXCAT - A shapelets "sexcat" type SExtractor catalogue structure.
;
; OPTIONAL INPUTS:
;      IMAGE  - A 2D image array or a shapelets "image" type image structure. 
;
; KEYWORD PARAMETERS:
;      Standard plotting parameters.
;
; OUTPUTS:
;      Draws a plot to the current output device.
; 
; MODIFICATION HISTORY:
;      Apr 02 - Background image incorporated by RM.
;      Dec 02 - Written by Richard Massey.
;-

; Error catching
if not shapelets_structure_type(sexcat,message=message,/silent) then message,message
if sexcat.type ne "sexcat" then message,'You need to input a "sexcat" type SExtractor catalogue!'
if sexcat.n le 1 then message,"Insufficient objects in the SExtractor catalogue to make a meaningful plot!"

; Is there an image to display in the background?
if keyword_set(image) then begin
  if shapelets_structure_type(image,message=message,/silent) and message eq "image" then begin
    im2plot=image.image
  endif else if size(image,/n_dim) eq 2 then begin
    im2plot=image
  endif
endif

; Draw image in the background if available (or just axes, if it is not)
if keyword_set(im2plot) then begin
  shapelets_plot_image, im2plot, _extra= ex, /fr
endif else begin
  plot,[0,0],/nodata,xrange=[min(sexcat.x[*,0]),max(sexcat.x[*,0])],/xstyle,$
                     yrange=[min(sexcat.x[*,1]),max(sexcat.x[*,1])],/ystyle
endelse

; Overplot locations of objects found
usersym, cos(2*!pi*findgen(21)/20), sin(2*!pi*findgen(21)/20)
for i=0L,sexcat.n-1 do $
  oplot,[sexcat.x[i,0],sexcat.x[i,0]],[sexcat.x[i,1],sexcat.x[i,1]],psym=8,$
    symsize=sexcat.fwhm[i]/median(sexcat.fwhm)

end
