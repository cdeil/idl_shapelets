;$Id: shapelets_read_shapecat.pro, v2$
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
; ************************************************************************
; ************************************************************************
;
; NAME:
;       SHAPELETS_SHAPECAT_MOMENTS
;
; PURPOSE:
;       Calculates the basic moments of galaxies in a shapelet catalgoue,
;       and adds those values to the IDL shapecat structure.
;
; CATEGORY:
;       Shapelet catalogue manipulation.
;
; INPUTS:
;       shapecat - IDL shapecat structure.
;
; OUTPUTS:
;       shapecat - IDL shapecat structure, including additional tags.
;
; MODIFICATION HISTORY:
;       Jul 05 - Full decomp structure used to calculate moments by RM.
;       Mar 05 - N_max keyword added by RM.
;       May 04 - Written by Richard Massey.

pro shapelets_shapecat_moments, shapecat, SILENT=silent

COMPILE_OPT idl2

; Strip catalogue of any old variables to avoid name conflicts and ensure backward compatability
message,'Calculating object moments',/info,noprint=silent
temp_shapecat={name:shapecat.name,type:shapecat.type}
for i=0,n_tags(shapecat)-1 do begin
  n=strupcase((tag_names(shapecat))[i])
  if n ne "NAME" and n ne "TYPE" and n ne "FLUX" and n ne "FLUX_ERROR" $
    and n ne "SNR" and n ne "MAG" and n ne "CENTROID" and n ne "CENTROID_ERROR" $
    and n ne "RSQUARED" and n ne "RSQUARED_ERROR" and n ne "CENTROID" $
    and n ne "CENTROID_ERROR" and n ne "ELLIPTICITY" and n ne "ELLIPTICITY_ERROR" $
    and n ne "QUADRUPOLE" and n ne "QUADRUPOLE_ERROR" and n ne "MOMENTS" then $
    temp_shapecat=create_struct(temp_shapecat,(tag_names(shapecat))[i],shapecat.(i))
endfor
shapecat=temp_shapecat

; Create empty arrays to contain shape moments
shapelets_make_nvec, shapecat.maxn_max, n1, n2, n_coeffs
flux              = fltarr(shapecat.n)
flux_error        = fltarr(shapecat.n)
snr               = fltarr(shapecat.n)
mag               = fltarr(shapecat.n)
centroid          = fltarr(shapecat.n,2)
centroid_error    = fltarr(shapecat.n,2)
rsquared          = fltarr(shapecat.n)
rsquared_error    = fltarr(shapecat.n)
quadrupole        = fltarr(shapecat.n,2,2)
quadrupole_error  = fltarr(shapecat.n,2,2)
ellipticity       = fltarr(shapecat.n,2)
ellipticity_error = fltarr(shapecat.n,2)

; Loop over each object in the catalogue
for i=0,shapecat.n-1 do begin
  ; Extract a (temporary) decomp structure
  decomp=shapelets_shapecat2decomp(shapecat,i)

  ; Calculate desired object moments
  temp_centroid = shapelets_centroid(decomp,/CARTESIAN,/ERROR,/IMAGE_COORDS)
  Q = shapelets_quadrupole(decomp,/CARTESIAN,/ERROR,flux=temp_flux,   $
        rsquared=temp_rsquared,ellipticity=temp_ellipticity,n_max=decomp.n_max)

  ; Store moments in global variables    
  flux[i]                 = temp_flux[0]
  flux_error[i]           = temp_flux[1]
  snr[i]                  = (flux_error[i] eq 0)?!values.f_nan:flux[i]/flux_error[i]
  mag[i]                  = -2.5*alog10(flux[i]>1e-10)
  centroid[i,*]           = temp_centroid[*,0]
  centroid_error[i,*]     = temp_centroid[*,1]
  rsquared[i]             = temp_rsquared[0]
  rsquared_error[i]       = temp_rsquared[1]
  quadrupole[i,*,*]       = Q[*,*,0]
  quadrupole_error[i,*,*] = Q[*,*,1]
  ellipticity[i,*]        = [float(temp_ellipticity[0]),imaginary(temp_ellipticity[0])]
  ellipticity_error[i,*]  = [float(temp_ellipticity[1]),imaginary(temp_ellipticity[1])]

endfor

; Add variables containing moments to shapelet catalogue
shapecat=create_struct(shapecat,"moments",1B,                             $
         "flux",flux,"flux_error",flux_error,                             $
         "snr",snr,"mag",mag,                                             $
         "centroid",centroid,"centroid_error",centroid_error,             $
         "rsquared",rsquared,"rsquared_error",rsquared_error,             $
         "quadrupole",quadrupole,"quadrupole_error",quadrupole_error,     $
         "ellipticity",ellipticity,"ellipticity_error",ellipticity_error)

end


; ************************************************************************
; ************************************************************************
;
;+
; NAME:
;      SHAPELETS_READ_SHAPECAT
;
; PURPOSE:
;      Reads in a shapelet coefficient catalogue of objects, created by
;      shex.pro, then stores them in a shapecat IDL structure. This may
;      then be converted into decomp structures using 
;      shapelets_shapecat2decomp.pro.
;
; CATEGORY:
;      Data I/O.
;
; INPUTS:
;      SHAPECAT  - name of structure where data will be stored.
;      FILENAME  - filename to be read, without path or extension.
;
; OPTIONAL INPUTS
;      N_MAX       - truncates all objects to only n_max coefficients.
;
; KEYWORD PARAMETERS:
;      SILENT      - operates silently.
;      CARTESIAN   - converts all objects to Cartesian shapelet coefs on load.
;      POLAR       - convert all objects to polar shapelet coefficients.
;      MOMENTS     - calculates shapelet unweighted quadrupole moments on load.
;                    This also gives flux and size.
;
; OUTPUTS:
;      SHAPECAT    - Shape information is stored in this structure.
;
; OPTIONAL OUTPUTS:
;      DESCRIPTION - String containing pertinent information.
;
; MEANING OF FIRST FLAG (from shapelets_image2obj.pro):
;      0: OK
;      1: did not normalise noise locally (is object filling postage stamp?)
;      2: nearby object (beware overlapping isophotes!)
;      3: both of above
;      4: near edge or a saturated star
;      5: HUGE object (SExtractor gaves whole image occasionally)
;      NB can only be one of these. If duplicated, a higher number takes priority.
;
; MEANINGS OF SECOND FLAG (from shapelets_focus.pro):
;      0: OK
;      1: has bounced into theta_min or max wall at some point
;      2: nmax_max reached - shapelets incompletely represent object
;      3: chi2>1 but converged (in n_max search) by flatness limit - still incomplete
;      4: did not converge
;      5: centroid off edge
;      6: fatal error
;      NB it can only be one of these. If duplicated, highest number has priority.
;
; MODIFICATION HISTORY:
;      Jul 05 - DESCRIPTION keyword added by RM.
;      May 05 - FULL_PATH keyword added by RM.
;      Mar 05 - GAUSSIAN keyword added by RM.
;      Aug 04 - Reduce option during conversion to polar shapelets added by RM.
;      May 04 - Written by Richard Massey
;-

pro shapelets_read_shapecat, shapecat, filename,  $
                             POLAR=polar,         $
                             CARTESIAN=cartesian, $
                             N_MAX=n_max,         $
                             MOMENTS=moments,     $
                             FULL_PATH=full_path, $
                             PARITY=parity,       $
                             SILENT=silent,       $
                             VERBOSE=verbose,     $
                             DESCRIPTION=decription

COMPILE_OPT idl2

; Set defaults
if not keyword_set(polar) then polar=0B
if keyword_set(cartesian) then polar=0B
if keyword_set(gaussian) then n_max_moments=2
if keyword_set(resize) then moments=1B
if not keyword_set(filename) then message,"You must specify a file name!"

; Strip file extension, then (re-)add it
file=filename ; Don't disrupt input variable
dot=strpos(strlowcase(file),".shapecat")
if dot ne -1 then file=strmid(file,0,dot)
file=file+".shapecat"

; Read in catalogue
if not keyword_set(full_path) then file=shapelets_paths(2,silent=silent)+file
message,"Reading in "+file,/info,noprint=silent
if float(strmid(!version.release,0,3)) ge 6.1 then begin
  restore,file,VERBOSE=verbose,DESCRIPTION=description
endif else begin
  restore,file,VERBOSE=verbose
endelse
if keyword_set(description) then print,description
if size(structure,/TYPE) gt 0 then begin
  shapecat=temporary(structure)
endif else if size(shapecat,/TYPE) eq 0 then begin
  message,"Shapecat was not found within the file!"
endif

; Check that the catalogue is really what we think it is
if not shapelets_structure_type(shapecat,message=message) then message,message
if shapecat.type ne "shapecat" then message,"Restored structure is a "+shapecat.type+", not a shapecat!"

; Adjust n_max of catalogue if desired
if keyword_set(n_max) then shapelets_extend_nmax, $
                         shapecat, n_max-shapecat.maxn_max, silent=silent

; Calculate shape moments if desired
if keyword_set(moments) then begin
  ; Check not already calculated
  already_calculated=0B
  if tag_exist(shapecat,"moments") then if shapecat.moments then already_calculated=1B
  if already_calculated then begin
    message,'Object moments already calculated in save file',/info,noprint=silent
  endif else shapelets_shapecat_moments, shapecat, silent=silent;, n_max=n_moments
endif

; Convert to polar coefficients if desired
if keyword_set(polar) then begin
  shapelets_polar_convert, shapecat, /POLAR, SILENT=silent
endif else if keyword_set(cartesian) and shapecat.polar eq 1 then begin
  shapelets_polar_convert, shapecat, /CARTESIAN, SILENT=silent
endif

; Flip and rotate objects so that they all have the same parity and orentation
if keyword_set(parity) then message,"Need to recode shapelets_reorient.pro!"

; Add in tags for forward compatibility
if not tag_exist(shapecat,"shear_estimates") then shapecat=create_struct(shapecat,"shear_estimates",0B)
if not tag_exist(shapecat,"morphology") then shapecat=create_struct(shapecat,"morphology",0B)

end



