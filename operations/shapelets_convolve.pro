pro shapelets_convolve, structure_a, structure_b, $
                        GAMMA=gamma,	          $
                        N_MAX_G=n_max_g,          $
                        SILENT=silent,            $
                        NOHISTORY=nohistory,      $
			EXTEND=extend,            $
		        MAINTAIN=maintain

;$Id: shapelets_convolve.pro, v2$
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
;      SHAPELETS_CONVOLVE
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Convolves one image with another (e.g. a PSF).
;
; INPUTS:
;      STRUCTURE_A - A shapelet decomp or shapecat structure.
;      STRUCTURE_B - A shapelet decomp structure representing the second image,
;                 or the PSF model (the PSF can be either structure since the
;                 operation is symmetric).
;
; OPTIONAL INPUTS:
;      GAMMA    - Gamma to use for convolved object.
;      N_MAX_G  - N_max to use for convolved object.
;                 NB: Must set both of these or neither (default choice
;                 specified in shapelets_convolution_matrix.pro conserves
;                 information content on al scales).
;
; KEYWORD PARAMETERS:
;      None.
;
; EXMAPLE:
;      decomp_convolved=shapelets_convolve(decomp,psf)
;
; OUTPUTS:
;      Returns a shapelet decomp structure representing the convolved image.
;
; TO DO:
;      Would fit better into the overall scheme of shapelet operation routines
;      if it were a function rather than a procedure.
;
; MODIFICATION HISTORY:
;      Jul 05 - Changed from function to procedure to match other ops by RM.
;      Jul 05 - Shapecats and polar decomp structs accepted as an input by RM.
;      Jun 05 - Unnecessary calls to shapelets_extend_nmax eliminated by RM.
;      Apr 05 - Default choice of gamma left up to subroutines by RM.
;      Jul 04 - Bug fixed (new object didn't get new beta) by RM.
;      Jun 04 - History and errors-on-coefficients tags included by RM.
;      Mar 02 - Written by Richard Massey.
;-

COMPILE_OPT idl2

; Parse inputs
; Check object/catalogue
if not shapelets_structure_type(structure_a,message=message) then message,message
if tag_exist(structure_a,"polar") then if structure_a.polar then shapelets_polar_convert,structure_a,/CARTESIAN,/SILENT
if keyword_set(extend) then shapelets_extend_nmax,structure_a,extend
; Check PSF
if not shapelets_structure_type(structure_b,message=message) then message,message
if structure_b.type ne "decomp" then message,"Second argument must be a decomp structure!"
if tag_exist(structure_b,"polar") then if structure_b.polar then shapelets_polar_convert,structure_b,/CARTESIAN,/SILENT
; Convolved scale size and truncation order
if keyword_set(gamma) then gamma_local=gamma
if keyword_set(n_max_g) then n_max_g_local=n_max_g
                        
; Perform convolution
if structure_a.type eq "decomp" then begin

  ; Obtain convolution matrix (paperI eqn54)
  alpha   = structure_a.beta
  n_max_a = structure_a.n_max
  P_nm    = shapelets_convolution_matrix(structure_b,alpha,n_max_a,gamma_local,n_max_g_local)

  ; Preform convolution
  f_m     = structure_a.coeffs
  h_n     = P_nm#f_m

  ; Might as well operate on errors at the same time 
  ; (assuming these transform in exactly the same way?)
  h_n_error=P_nm#structure_a.coeffs_error
  if n_elements(structure_a.coeffs_error) eq structure_a.n_coeffs $
    then h_error=P_nm#structure_a.coeffs_error $
    else message,"Can't yet handle full covariance error matrix!"
  
  ; Update the structure for output
  if (n_max_g_local ne n_max_a) then $				          ;
    shapelets_extend_nmax,structure_a,n_max_g_local-n_max_a,silent=silent ; Set its new n_max
  structure_a.coeffs	   = h_n				          ; Copy the new shapelet coefficients into the structure (keeping flux constant when we change beta in a moment))
  structure_a.coeffs_error = h_n_error				          ; Copy the error on the new shapelet coefficients
  structure_a.beta 	   = gamma_local 			          ; Set its new scale size
  
  ; Decide number of pixels to use for output
  structure_a.n_pixels=structure_a.n_pixels>structure_b.n_pixels
  structure_a.x=structure_a.x>structure_b.x

endif else if structure_a.type eq "shapecat" then begin
  
  ; Loop over all objects in a shapecat
  n_objects=structure_a.n
  for i=0,n_objects-1 do begin
    ; Extract a (temporary) decomp structure
    decomp=shapelets_shapecat2decomp(structure_a,i)
    ; Perform convolution
    shapelets_convolve,decomp,structure_b,/NOHISTORY,gamma=gamma[i],n_max_g=n_max_g[i],SILENT=silent
    ; Append the convolved object to the end of the catalogue
    shapelets_add,structure_a,decomp,/NOHISTORY
  endfor
  ; Remove original objects from the catalogue
  shapelets_split,structure_a,lindgen(n_objects)+n_objects,/NOHISTORY
  
endif else message,"Cannot convolve "+structure_a.type+" structures!"

; Add operation to object's history record
if not keyword_set(nohistory) then shapelets_update_history,structure_a,"Convolved with "+structure_b.name

end
