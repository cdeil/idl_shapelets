function randome,seed,n,m,l

;$Id: randome.pro, v2$
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

;+
; NAME:
;      RANDOME
;
; CATEGORY:
;      Random number generation.
;
; PURPOSE:
;      Generates a random number (or nxmxl array of) from a rescaled
;      Epanechnikov kernel K(x) = .75(1-x^2) for |x|<1; K(x)=0 elsewhere.
;      See Silverman (1986) page 143 for discussion.
;
; INPUTS:
;      None.
;
; OPTIONAL INPUTS:
;      seed   - Initial seed for random number generation.
;      n,m,l  - Integers, describing the dimensions of a 1D, 2D or 3D array
;               which should be returned (full of random numbers). If these
;               are omitted, a single random number is returned.
;
; KEYWORD PARAMETERS:
;      None.
;
; OUTPUT: 
;      Function returns a random number, or an array of random numbers.
;
; EXAMPLE: 
;      print,randome()
; 
; PROCEDURES USED: 
;      None. 
; 
; MODIFICATION HISTORY:
;      May 02 - Written by Richard Massey
;-

; Determine random seed and the number of random numbers requested.
switch N_Params() of
  0: seed=fix(strmid(systime(0),17,2)+strmid(systime(0),14,2))
  1: n=1
  2: m=1
  3: l=1
  4: break
  else: message,'Need to rewrite randome.pro to handle this many dimensions'
endswitch

; Carry out procedure described in Silverman (1986) page 143
v_one=randomu(seed,n,m,l)*2.-1.
v_two=randomu(seed,n,m,l)*2.-1.
v_thr=randomu(seed,n,m,l)*2.-1.

randval = v_thr + $
 ((abs(v_thr) gt abs(v_two)) and (abs(v_thr) gt abs(v_one))) * (v_two - v_thr)

; Return random number(s).
return,randval

end
