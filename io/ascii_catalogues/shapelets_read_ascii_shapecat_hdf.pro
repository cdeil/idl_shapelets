pro shapelets_read_ascii_shapecat_hdf,shapecat,       $
                                      polar=polar,    $
					             moments=moments,$
                                      filt=filt,      $
					             deconv=deconv,  $
					             n_max=n_max,    $
                                      parity=parity,  $
					             resize=resize,  $
                                      v1=v1

;$Id: shapelets_read_ascii_shapecat_hdf.pro, v2$
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
;       SHAPELETS_READ_ASCII_SHAPECAT_HDF
;
; PURPOSE:
;       Reads in and concatenates shapelet catalogues containing all of the
;       objects from both Hubble Deep Fields.
;
; INPUTS:
;       None.
;
; OPTIONAL INPUTS:
;       FILT    - HDF filter e.g. 606 or 814. Default: 814
;       N_MAX   - Truncate catalogue at this shapelet parameter.
;
; KEYWORD PARAMETERS:
;       POLAR   - Convert to polar shapelet coefficents automatically on load,
;                 in reduced {r,theta} form (cf shapelets_polar_reduce.pro).
;       MOMENTS - Calculate object shape moments (astrometry, photometry etc)
;                 automatically on load.
;       DECONV  - Read PSF-deconvolved catalogue (if available)
;       PARITY  - Rotate all objects automatically on load so that their
;                 ellipticity is aligned with the x axis then flip, if
;                 necessary, so that all their isophotes begin to twist in a
;                 clockwise sense.
;       RESIZE  - Rescale so that all shapelets use the same beta.
;       V1      - Load shapelet coefficients as a_ij/a_00, as in version 1
;                 of shapelet code. Default is to store just a_ij.
;
; OUTPUTS:
;       shapecat- name of catalogue (shapecat structure) returned in memory.
;
; MODIFICATION HISTORY:
;       Mar 04 - Split into shapelets_concatenate_shapecats by RM.
;       Jan 02 - Written by Richard Massey.
;-

COMPILE_OPT OBSOLETE

; Set up filenames
if n_params() eq 0 then message,'USE: shapelets_read_shapecat_hdf,shapecat !'
if keyword_set(filt) then filter='f'+strtrim(filt,2) else filter='f814'
in_dir=shapelets_paths(2)


; Read in catalogues
message,"Reading in shape catalogue of HDF north",/info
if keyword_set(polar)   then message,"Conversion to polars may take a few minutes",/info
if keyword_set(moments) then message,"Calculation of objects' moments may take a few minutes",/info
if keyword_set(deconv)  then filter2=filter+"_deconv" else filter2=filter
shapelets_read_ascii_shapecat,cat_n,'HDF-N_'+filter,$
  moments=moments,polar=polar,/silent,parity=parity,resize=resize,n_max=n_max,v1=v1
message,'Reading in shape catalogue of HDF south',/info
shapelets_read_ascii_shapecat,cat_s,'HDF-S_'+filter,$
  moments=moments,polar=polar,/silent,parity=parity,resize=resize,n_max=n_max,v1=v1


; Concatenate catalogues
shapecat=shapelets_concatenate_ascii_shapecats(cat_n,cat_s)
shapecat.name="HDFs"
if not keyword_set(moments) then moments=0B
shapecat=create_struct(shapecat,"moments",moments,"sextractor",1B,"coeffs_error",shapecat.error)
message,'HDF shapelet catalogues successfully in memory',/info


end
