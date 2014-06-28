function shapelets_concatenate_ascii_shapecats,shapecat1,shapecat2

;$Id: shapelets_concatenate_ascii_shapecats.pro, v2$
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
;      SHAPELETS_CONCATENATE_ASCII__SHAPECATS
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Concatenates two (shapecat structure) catalogues.
;
; INPUTS:
;      shapecat1 - Shapelet catalogue structure. 
;      shapecat2 - Shapelet catalogue structure. 
;
; OPTIONAL INPUTS:
;      none.
;
; KEYWORD PARAMETERS:
;      none.
;
; OUTPUTS:
;      shapecat  - Shapelet catalogue structure, containing objects from both
;                  of the input catalogues.
;
; NOTES: 
;      The two catalogues need to have the same n_max, and be both in polar or
;      both in Cartesian shapelets.
;
; MODIFICATION HISTORY:
;      Apr 04 - Obfuscated by RM when changed convention to write catalogues
;               to disc as IDL save files
;      Mar 04 - Written by R.Massey
;-

COMPILE_OPT OBSOLETE

; Check to see the two catalogues are compatible
if (shapecat1.type ne shapecat2.type) then message,$
  "The two structures are not of the same type and are therefore not compatible!"
if (shapecat1.maxn_max ne shapecat2.maxn_max) then message,$
  "The two catalogues are not to the same n_max and are therefore not compatible!"
if ((shapecat1.polar+shapecat2.polar) mod 2 ne 0) then message,$
  "One catalogue contains Cartesian shapelets, the other polar shapelets. They are therefore not compatible!"

; Compile a concatenated list of variables
name=shapecat1.name+"_and_"+shapecat2.name           	 
type=shapecat1.type           	 
n=shapecat1.n+shapecat2.n
maxn_max=shapecat1.maxn_max         	 
polar=shapecat1.polar
seeing=(shapecat1.seeing+shapecat2.seeing)/2.
x=[shapecat1.x,shapecat2.x]
beta=[shapecat1.beta,shapecat2.beta]
n_max=[shapecat1.n_max,shapecat2.n_max]
n_coeffs=shapecat1.n_coeffs
coeffs=[shapecat1.coeffs,shapecat2.coeffs]
error=[shapecat1.error,shapecat2.error]
flux=[shapecat1.flux,shapecat2.flux]
flux_error=[shapecat1.flux_error,shapecat2.flux_error]
snr=[shapecat1.snr,shapecat2.snr]
mag=[shapecat1.mag,shapecat2.mag]     
centroid=[shapecat1.centroid,shapecat2.centroid]
centroid_error=[shapecat1.centroid_error,shapecat2.centroid_error]
rsquared=[shapecat1.rsquared,shapecat2.rsquared]
rsquared_error=[shapecat1.rsquared_error,shapecat2.rsquared_error]
ellipticity=[shapecat1.ellipticity,shapecat2.ellipticity]
ellipticity_error=[shapecat1.ellipticity_error,shapecat2.ellipticity_error]
quadrupole=[shapecat1.quadrupole,shapecat2.quadrupole]
quadrupole_error=[shapecat1.quadrupole_error,shapecat2.quadrupole_error]
flag=[shapecat1.flag,shapecat2.flag]
chisq=[shapecat1.chisq,shapecat2.chisq]
sexid=[shapecat1.sexid,shapecat2.sexid]
sexfwhm=[shapecat1.sexfwhm,shapecat2.sexfwhm]
sexflux=[shapecat1.sexflux,shapecat2.sexflux]
sexmag=[shapecat1.sexmag,shapecat2.sexmag]
sexclass=[shapecat1.sexclass,shapecat2.sexclass]
sexa=[shapecat1.sexa,shapecat2.sexa]
sexb=[shapecat1.sexb,shapecat2.sexb]
sextheta=[shapecat1.sextheta,shapecat2.sextheta]
sexe1=[shapecat1.sexe1,shapecat2.sexe1]
sexe2=[shapecat1.sexe2,shapecat2.sexe2]
sexx=[shapecat1.sexx,shapecat2.sexx]
sexy=[shapecat1.sexy,shapecat2.sexy]
sexarea=[shapecat1.sexarea,shapecat2.sexarea]
objx_max=[shapecat1.objx_max,shapecat2.objx_max]
objx_min=[shapecat1.objx_min,shapecat2.objx_min]
objy_min=[shapecat1.objy_min,shapecat2.objy_min]
objy_max=[shapecat1.objy_max,shapecat2.objy_max]


; Store in structure
shapecat={name:name,                    	  $
     	type:type,		      		  $
     	n:n, 						  $
     	maxn_max:maxn_max,  			  $
     	polar:polar,					  $
     	seeing:seeing, 				  $
     	x:x, 						  $
     	beta:beta,					  $
     	n_max:n_max,					  $
     	n_coeffs:n_coeffs,  			  $
     	coeffs:coeffs, 				  $
     	error:error,					  $
     	flux:flux,					  $
     	flux_error:flux_error,			  $
     	snr:snr,  					  $
     	mag:mag,  					  $
     	centroid:centroid,  			  $
     	centroid_error:centroid_error,	  $
     	rsquared:rsquared,  			  $
     	rsquared_error:rsquared_error,	  $
     	ellipticity:ellipticity, 		  $
     	ellipticity_error:ellipticity_error, $
     	quadrupole:quadrupole,			  $
     	quadrupole_error:quadrupole_error,   $
     	flag:flag,					  $
     	chisq:chisq,					  $
     	sexid:sexid,					  $
     	sexfwhm:sexfwhm,				  $
     	sexflux:sexflux,				  $
     	sexmag:sexmag, 				  $
     	sexclass:sexclass,  			  $
     	sexa:sexa,					  $
     	sexb:sexb,					  $
     	sextheta:sextheta,  			  $
     	sexe1:sexe1,					  $
     	sexe2:sexe2,					  $
     	sexx:sexx,					  $
     	sexy:sexy,					  $
     	sexarea:sexarea,				  $
     	objx_max:objx_max,  			  $
     	objx_min:objx_min,  			  $
     	objy_min:objy_min,  			  $
     	objy_max:objy_max}
	
return,shapecat

end
