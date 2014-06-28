function shapelets_geometric_constraints, PSTAMP,                        $
                                          PSF=psf,                       $
                                          CENTRE_GUESS=centre_guess,     $
                                          THETA_MIN_GEOM=theta_min_geom, $
                                          THETA_MAX_GEOM=theta_max_geom

;$Id: shapelets_geometric_contstraints.pro, v2$
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
;      SHAPELETS_GEOMETRIC_CONTSTRAINTS
;
; PURPOSE:
;      Determine theta_min and theta_max for an object.
;
; CATEGORY:
;      Shapelets.
;
; INPUTS:
;      PSTAMP - IDL "pstamp" structure containing a single object from an image.
;
; OPTIONAL INPUTS:
;      THETA_MIN_GEOM - minimum scale for geometrical constraint.
;                    (default: 0.2 pixels)
;      CENTRE_GUESS          - centre of shapelet basis functions.
;
; KEYWORD PARAMETERS:
;      None.
;
; OUTPUTS:
;      THETA_MIN_GEOM - minimum scale for geometrical constraint.
;      THETA_MAX_GEOM - maximum scale for geometrical constraint.
;      
; MODIFICATION HISTORY:
;      Jul 05 - Correction for PSF (de-)convolution included by RM.
;      Jan 05 - Written by Richard Massey.
;-

COMPILE_OPT idl2

; Guess centroid using sextractor centroid
if keyword_set(centre_guess) then centre_guess_local=centre_guess else centre_guess_local=[pstamp.xo,pstamp.yo]

; Guess theta_min_geom as the rms of a square pixel (0.166)
if keyword_set(theta_min_geom) then theta_min_geom_local=theta_min_geom else theta_min_geom_local=0.2

; Guess theta_max_geom
if not keyword_set(theta_max_geom_local) then begin
  n_pixels=[pstamp.im_ran[1]-pstamp.im_ran[0],pstamp.im_ran[3]-pstamp.im_ran[2]]
  case strupcase(pstamp.geometry) of
    "SQUARE": theta_max_geom_local=min([centre_guess_local[0],n_pixels[0]-centre_guess_local[0],centre_guess_local[1],n_pixels[1]-centre_guess_local[1],theta_min_geom_local])
    "CIRCLE": begin
                radius=min(n_pixels)/2.
                centre=n_pixels/2.
                to_centre=centre_guess_local-centre
                dist_to_centre=sqrt(to_centre[0]^2+to_centre[1]^2)
                theta_max_geom_local=(radius-dist_to_centre)>theta_min_geom_local	      
              end
    else: message,"Postage stamp geometry not recognised!"
  endcase
endif


; Correct for PSF (de-)convolution
if keyword_set(psf) then begin
  theta_min_psf=psf.beta/sqrt(psf.n_max+1)
  theta_max_psf=psf.beta*sqrt(psf.n_max+1)
  theta_min_geom_local=sqrt(((theta_min_geom_local)^2-(theta_min_psf)^2)>0)
  theta_max_geom_local=sqrt(((theta_max_geom_local)^2-(theta_max_psf)^2)>theta_min_geom_local^2)
endif

return,[theta_min_geom_local,theta_max_geom_local]

end
