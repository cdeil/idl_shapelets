function xgen,x1,x2,NPOINTS=npoints,LOGPLOT=logplot

;$Id: xgen.pro, v2$
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
;       XGEN
;
; PURPOSE:
;       Generate a vector containing equally spaced numbers between
;       x1 and x2. This is typically used to generate an x-axis
;       vector in a plot.
;
; CATEGORY:
;       Miscellaneous.
;
; CALLING SEQUENCE:
;       Xarray = XGEN(x1, x2 [,NPOINTS ,/LOGPLOT] )
;
; INPUTS:
;       x1,x2: interval limits
;
; KEYWORD PARAMETERS:
;       NPOINTS: number of points (default=100)
;       LOGPLOT: produce spacing adapted for a log plot
;
; OUTPUTS:
;       Array of floating point numbers.
;
; PROCEDURES USED:
;       findgen - generate arrays of integers (as floating point variables).
;
; MODIFICATION HISTORY:
;       Jul 03 - Header added by Richard Massey
;       Jun 95 - Written by Alexandre Refregier
;-

COMPILE_OPT idl2

if keyword_set(npoints) then np=npoints else np=100
if keyword_set(logplot) $
  then xx=10.^(float(findgen(np))*(alog10(x2)-alog10(x1))/float(np)+alog10(x1)) $
  else xx=float(findgen(np))*(float(x2)-float(x1))/float(np)+float(x1)

return,xx
end
