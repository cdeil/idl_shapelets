pro plt_rectangle,bl,tr,thick=thick,color=color,linestyle=linestyle

;$Id: plt_rectangle.pro, v2$
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
;       PLT_RECTANGLLE
;
; PURPOSE:
;       Draws a rectangle on the current output device.
;
; CATEGORY:
;       Graphics plotting.
;
; INPUTS:
;       bl: [x,y] coordinates of bottom-left corner.
;       tr: [x,y] coordinates of top-right corner.
;
; KEYWORD PARAMETERS:
;       thick: line thickness (as in IDL plot command)
;       color: line colour index (as in IDL plot command)
;       linestyle: line style parameter (as in IDL plot command)
;
; OUTPUTS:
;       Plot to current output device.
;
; EXAMPLE:
;       
;
; PROCEDURES USED:
;       oplot.
; 
; MODIFICATION HISTORY:
;       Dec 01 - Written by Richard Massey 
;-

COMPILE_OPT idl2

oplot,[bl[0],bl[0]],[bl[1],tr[1]],thick=thick,color=color,linestyle=linestyle
oplot,[bl[0],tr[0]],[bl[1],bl[1]],thick=thick,color=color,linestyle=linestyle
oplot,[bl[0],tr[0]],[tr[1],tr[1]],thick=thick,color=color,linestyle=linestyle
oplot,[tr[0],tr[0]],[bl[1],tr[1]],thick=thick,color=color,linestyle=linestyle

end
