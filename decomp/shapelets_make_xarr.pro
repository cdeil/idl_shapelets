pro shapelets_make_xarr, n, x1, x2, X0=x0

;$Id: shapelets_make_xarr.pro, v2$
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
;      SHAPELETS_MAKE_XARR
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Make arrays containing the values of the x (x1) and  y (x2) coordinates.
;      of a grid. This is convenient for evaluating 2D functions (e.g. shapelet
;      basis functions) on a grid.
;
; INPUTS:
;      n  - Integer array [n1,n2]. Array dimensions (number of pixels).
;
; OPTIONAL INPUTS:
;      x0 - Floating point array [xc1,xc2]. Pixel coordinate of the [0.0,0.0]
;           origin. Default: the centre of the array.
;
; KEYWORD PARAMETERS:
;      None.
;
; OUTPUTS:
;      x1 - Floating point [n1*n2] array containing x values (all of the values
;           in any particular column will be identical).
;      x2 - Floating point [n1*n2] array containing y values (all of the values
;           in any particular row will be identical).
;
; MODIFICATION HISTORY:
;      Dec 01 - Debugged by R. Massey
;      Jul 99 - Written by A. Refregier
;-

COMPILE_OPT idl2

; Set reference point to the centre if not specified
if not keyword_set(x0) then x0=.5*float(n)

; Make x arrays
one = replicate(1.,n[1]>1)
vec = findgen(n[0]>1)-x0[0]+.5  ; x(0) to centre and 0.5 to get into middle of pixel
x1  = vec#one

; Make y arrays
one = replicate(1.,n[0]>1)
vec = findgen(n[1]>1)-x0[1]+.5
x2  = one#vec

end
