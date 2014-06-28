function sign, x

;$Id: sign.pro, v2$
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
;       SIGN
;
; PURPOSE:
;       Compute the sign(x) of a variable x, with sign(x)=0,+1,-1, for
;       x=0,>0,<0 respectively. The variable x can have any dimension.
;
; CATEGORY:
;       Mathematics.
;
; INPUTS:
;       x: scalar or vector
;
; KEYWORD PARAMETERS:
;       None.
;
; OUTPUTS:
;       SIGN(x): integer (either -1, 0 or 1)
;
; EXAMPLE:
;       print,SIGN(12)   ... gives 1
;       print,SIGN(-4.5) ... gives -1
;
; PROCEDURES USED:
;       None.
; 
; MODIFICATION HISTORY:
;       Jul 03 - Header added by RM.
;       Mar 03 - Sped up by Richard Massey, using <> operators rather than
;                a for/next loop.
;       May 99 - Written by Alexandre Refregier.
;-

COMPILE_OPT idl2

n = n_elements(x)
s = intarr(n)

neg=where(x<0) & if ( neg[0] ne -1 ) then s[neg]=-1
pos=where(x>0) & if ( pos[0] ne -1 ) then s[pos]=1

; AR's old version
;for i=0L, n-1 do begin
;   if x(i) gt 0. then s(i) = 1
;   if x(i) lt 0. then s(i) = -1
;endfor

if ( n eq 1 ) then s = s[0]

return, s

end
