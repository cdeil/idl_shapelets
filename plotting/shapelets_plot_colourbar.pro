pro shapelets_plot_colourbar,range,title=title,scalable=scalable,inverse=inverse,$
  csize=csize,clog=clog,_ref_extra = ex

;$Id: shapelets_plot_colourbar.pro, v2$
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
; October 1997 - Written by A. Refregier
;
; PURPOSE: draw an annotated color bar at the top of the plotting
; window. The window parameters are then set to the remainder of
; the window.
; INPUT: range: range of the scale for the color bar
; OPTIONAL INPUT: title: title for the color bar
; 		  scalable: use scalable pixels (useful when
;                           producing postcript files)
;                 inverse: invert the color coding
;                 csize: vertical color bar size (0-1)
; OUTPUT: an annotated color bar is drawn at the top of the window.
; The window parameters are then set to the remainder of the window.
; NOTE: to restore normal window region set !p.region=[0,0,0,0]
;-

; set color bar vertical size
if not keyword_set(csize) then csize=.14

; reserve upper portion of the window for the color bar
!p.region=[0.,1.-csize,1.,1.]

; Draw colour bar
shapelets_plot_image,(findgen(8>!d.n_colors<512))*(-1)^keyword_set(inverse),$
  frame=-1,title=title,_extra= ex

; Draw a box around it
plot,[0],[0],/nodata,/noerase,xrange=range,/xstyle,title=title,$
  xminor=1,yminor=1,yticks=1,ytickname=[" "," "],_extra= ex

; leave the rest of the window for the image
!p.region=[0.,0.,1.,1.-csize]

end
