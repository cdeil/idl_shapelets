function shapelets_image_moments,image,printit=printit

;$Id: shapelets_image_moments.pro, v2$
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
;     SHAPELETS_IMAGE_MOMENTS
;
; PURPOSE:
;     Compute the zeroth (flux) and first (centroid) moment of an image.
;
; CATEGORY:
;     Image analysis.
;
; INPUTS:
;     image - 2D image array.
;
; OUTPUTS:
;     mom structure - {f: flux, xc: centroid, r: rms radius}
;
; MODIFICATION HISTORY:
;     Jan 05 - Quadrupole and ellipticity added by Richard Massey
;     Jul 00 - Written by A. Refregier
;-

; make x-arrays
nf=size(image,/dimensions)
shapelets_make_xarr,nf,x1,x2,x0=[0.,0.]

; compute the flux
flux=total(image)

; compute the centroid
x1c=total(image*x1)/flux
x2c=total(image*x2)/flux

; compute the rms radius
r=sqrt(total(((x1-x1c)^2+(x2-x2c)^2)*image)/flux)

; compute the quadrupole moments
q11=total((x1-x1c)^2*image)/flux
q22=total((x2-x2c)^2*image)/flux
q12=total((x1-x1c)*(x2-x2c)*image)/flux

; compute ellipticity and other quantities
e1=(q11-q22)/(q11+q22)   ; ellipticity   [1]
e2=2.*q12/(q11+q22)
a=sqrt(q11+q22+sqrt((q11-q22)^2+4.*q12^2))  ; principal axes [pixels]
b=sqrt(q11+q22-sqrt((q11-q22)^2+4.*q12^2))
alpha=.5*atan(2.*q12,q11-q22)   ; position angle, counter-clockwise from x-axis
                                ;  [rad]

; print out results
if keyword_set(printit) then begin
  print,'flux:',fl
  print,'centroid:',x1c,x2c
  print,'rms radius r:',r
  print,'ellipticity: e1,e2:',e1,e2
  print,'major, minor axis:',a,b
  print,'position angle (deg):',alpha*!radeg
endif

; store in structure
mom={flux:flux,xc:[x1c,x2c],r:r,q11:q11,q22:q22,q12:q12,$
      e1:e1,e2:e2,a:a,b:b,alpha:alpha}

return,mom

end
