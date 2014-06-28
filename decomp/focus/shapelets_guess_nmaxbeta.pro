function shapelets_guess_nmaxbeta, PSTAMP, $
                                   PSF=psf,  $
                                   THETA_MIN_GEOM=theta_min_geom,  $
                                   THETA_MAX_GEOM=theta_max_geom,  $
                                   PRINTIT=printit,          $
                                   CENTRE_GUESS=centre_guess

;$Id: shapelets_guess_nmaxbeta.pro, v2$
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
;      SHAPELETS_GUESS_NMAXBETA
;
; PURPOSE:
;      Guess the optimal beta and nmax for an object.
;
; CATEGORY:
;      Shapelets.
;
; INPUTS:
;      pstamp - IDL "pstamp" structure from shapelets_extract_pstamp.pro.
;
; OPTIONAL INPUTS:
;      th_min_geom - minimum scale for geometrical constraint.
;                    (default: 0.5 pixels)
;
; KEYWORD PARAMETERS:
;      /PRINTIT - if this is set, the guesses are written to stdout.
;
; OUTPUTS:
;      guess - IDL structure containing guessed n_max, beta etc.
;
; MODIFICATION HISTORY:
;      Apr 05 - Th_min/max constraints included in beta guess by RM.
;      Nov 03 - Beta guess made more sensible when deconvolving a PSF by RM.
;      Oct 03 - Beta guess changed from theta_min/max to a calibrated
;               formula based upon SExtractor's FWHM measurement by RM.
;      May 02 - Written by A. Refregier
;-

COMPILE_OPT idl2

; Compute geometrical limits for th_max and th_min
th_minmax_geom=shapelets_geometric_constraints(pstamp, $
                                               psf=psf, $
                                               theta_min_geom=theta_min_geom, $
                                               theta_max_geom=theta_max_geom, $
                                               centre_guess=centre_guess)
if not keyword_set(centre_guess) then centre_guess=[pstamp.xo,pstamp.yo]
th_min_geom=th_minmax_geom[0]
th_max_geom=th_minmax_geom[1]
th_min=th_min_geom
th_max=(2.1*pstamp.a)<th_max_geom

;mom=shapelets_image_moments(pstamp.image)
;th_max=mom.r*1.
;th_min=1.
;beta_guess=sqrt(th_max*th_min)
;n_guess=round(th_max/th_min)
;if not keyword_set(x0) then x0=mom.xc
;print,'initial guess: w, n,x0:',beta_guess,n_guess,x0

; Compute n_max for this guess
n_max=round(th_max/th_min)
if (n_max mod 2) then n_max=n_max+1

; Compute beta for this guess
;beta=sqrt(th_min*th_max)              ; This would definitely be allowed by the geometrical constraints
beta=(pstamp.fwhm/2.86)>1.             ; Constant determined by investigating results using rule above
;beta=(pstamp.fwhm/4.70960)>1.          ; Constant guessed assuming objects are Gaussians (with FWHM 2.3548)
beta_max=th_max/sqrt(n_max+1)          ;
beta_min=th_min*sqrt(n_max+1)<beta_max ;
beta=beta_min>beta<beta_max            ; Now just have to check that this is allowed!
if keyword_set(psf) then beta=(sqrt(beta^2-psf.beta^2))>(beta/3.)

; Print guesses to screen if requested
if keyword_set(printit) then begin
  message,'guess: x0, y0;',x0[0],x0[1],/info
  message,'guess: beta, nmax:',beta,n_max,/info
  message,'guess: th_min, th_max:',th_min,th_max,/info
  message,'geometry: th_min, th_max:',th_min_geom,th_max_geom,/info
endif

; Store in structure
guess={beta:beta,$
       nmax:n_max,$
       centre_guess:centre_guess,$
       th_min:th_min,$
       th_max:th_max,$
       th_min_geom:th_min_geom,$
       th_max_geom:th_max_geom}

return,guess

end
