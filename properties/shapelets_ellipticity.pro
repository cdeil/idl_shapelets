function shapelets_ellipticity, decomp,            $
                                N_MAX=n_max,       $
                                GAUSSIAN=gaussian, $
                                FLUX=flux,         $
                                COMPLEX=complex,   $
                                ETHETA=etheta

;$Id: shapelets_ellipticity.pro, v2$
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
;      SHAPELETS_ELLIPTICITY
;
; PURPOSE:
;      Calculates the ellipticity a shapelet model.
;
; EXAMPLE USE:
;       result=shapelets_ellipticity(decomp) & print,result
;
; INPUTS:
;       decomp - One Cartesian shapelet decomposition
;
; OPTIONAL INPUTS: 
;        n_max - Works out the sum only to n=nmax.
;
; KEYWORD PARAMETERS:
;    /GAUSSIAN - Force calculation to only n_max=2, to compute the Gaussian
;                 weighted ellipticity, with a Gaussian of size beta.
;     /COMPLEX - Output an ellipticity in complex form (see below). DEFAULT.
;        /E1E2 - Output an ellipticity as a two component array [e1,e2].
;      /ETHETA - Output an ellipticity in [modulus,argument] form 
;                (the angle is give in radians a/c/w from the x-axis).
;
; OUTPUTS:
;  ELLIPTICITY - The ellipticity of the object. The default format is a complex
;                number e_1 + i e_2, a convenient notation commonly used in the
;                field of weak gravitational lensing. Positive e1 (e2) denotes
;                elongation along the x-axis (the line y=x). Negative e1 (e2) 
;                denotes elongation along the y-axis (the line y=-x). Only a 
;                circle has zero e1 and e2, and therefore zero |e|.
;                This format can be changed using optional keywords.              
;
; OPTIONAL OUTPUTS:
;         FLUX - Total object flux is optionally returned in this variable.
;
; MODIFICATION HISTORY:
;       Mar 05 - Written by Richard Massey.
;-

; Allow for various options
if keyword_set(gaussian) then n_max_calc=2 else if keyword_set(n_max) then n_max_calc=n_max

; Calculate the ellipticity
quadrupole=shapelets_quadrupole(decomp,ellipticity=ellipticity,$
                                n_max=n_max_calc,flux=flux)

; Rearrange output into the desired format
if keyword_set(complex) then begin
  ellipticity=ellipticity
endif else if keyword_set(e1e2) then begin
  ellipticity=[float(ellipticity),imaginary(ellipticity)]
endif else if keyword_set(etheta) then begin
  e=sqrt(float(ellipticity)^2+imaginary(ellipticity)^2)
  theta=atan(imaginary(ellipticity),float(ellipticity))/2.
  ellipticity=[e,theta]
endif

; Tell the world
return,ellipticity

end
