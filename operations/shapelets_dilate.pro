pro shapelets_dilate, decomp, scaling,      $
                      SILENT=silent,        $
                      EXTEND=extend ,       $
                      RADIUS=radius,        $
                      AREA=area,            $
                      BETA=beta,            $
                      CARTESIAN=cartesian,  $
                      POLAR=polar,          $
                      ORDER=order,          $
                      FLUX=flux,            $
                      NOHISTORY=nohistory,  $
		                  MAINTAIN=maintain,    $
                      SURFACE_BRIGHTNESS=surface_brightness

;$Id: shapelets_dilate.pro, v2$
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
;      SHAPELETS_DILATE
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Dilates (changes the size of) an object. There are four methods
;      available to perform this operation:
;      (1) Simply rescale beta, the shapelet size parameter.
;      (2) To first order, in Cartesian shapelets, using a&a^ ladder operators.
;      (3) Convert to polar shapelets and use ladder operators, but again only
;          to first order. Then convert back to Cartesian coefficients.
;      (4) An alternative method uses the C_nml formalism for convolution.
;          Convolution with a delta function (ie a PSF with beta -> 0) leaves
;          the object unchanged, but we can now choose a new beta (in
;          shapelets_convolution_matrix.pro, the notation is to change alpha to
;          gamma). This is a general case, no longer limited to first order.
;      Option four is the default; the others are selected using switches.
;
; CAUTION:
;      In general, the default options are suitable in the field of weak
;      gravitational lensing. Care should be taken in other applications that
;      the correct flags are added (especially /FLUX) to produce the desired
;      result.
;
; INPUTS:
;      DECOMP - A Cartesian shapelet decomposition of one object.
;      SCALING- Factor by which the size of the object should change.
;               For the interpretation of this factor, see the options
;               available in "Keyword Parameters" below.
;
; OPTIONAL INPUTS:
;      EXTEND - Final result is up to this n_max. Can be used to extend the
;               coefficient array to take new, higher-order coefficients.
;
; KEYWORD PARAMETERS:
;      RADIUS - Normal scaling factor is ignored; radius scales by factor kappa.
;      AREA   - Normal scaling factor is ignored; area scales by factor kappa.
;               DEFAULT: Scaling is treated as the fractional change in radius,
;                        which is equal to the "convergence" kappa in the
;                        notation of weak gravitational lensing:
;                               +ve Object is enlarged
;                               -ve Object shrinks
;                                0  Object stays the same
;                        The area changes by 1/(1-kappa)^2 =approx= 1+2kappa.
;      FLUX   - Conserve flux during the mapping.
;      SURFACE- Conserve surface brightness during the mapping (DEFAULT).
;      BETA   - Simply rescale the shapelet size parameter.
;      CART   - Perform first order ladder operations on Cartesian shapelets.
;      POLAR  - Perform ladder operations in polar shapelet space.
;               DEFAULT: To general order using C_nml convolution formalism.
;
; NOTES:
;      The RADIUS and AREA options are mutually exclusive. So are the FLUX and
;      SURFACE options, and the BETA, CART and POLAR options. If more than one
;      of any of these combinations is specified, the first takes priority.
;
; OUTPUTS:
;      DECOMP - the Cartesian shapelet decomp strucutre is returned, rescaled.
;
; EXAMPLE USE:
;      A change of scale size beta (without affecting the object) can be
;      implemented via
;      IDL> shapelets_dilate,e,-delta_beta,/flux
;      IDL> shapelets_dilate,e,delta_beta,/flux,/beta
;
; MODIFICATION HISTORY:
;      Sep 05 - Integer division problem fixed with RADIUS keyword by RM.
;      Jun 05 - N_MAX keyword dropped in favour of usual EXTEND keyword by RM.
;      Apr 05 - Higher order transformation via exponenitation added by RM.
;      Aug 04 - Conservation of surface brightness changed to default by RM.
;      Jul 04 - Bug fixed in Cartesian shapelet implementation by RM.
;      Jun 04 - More intuitive options for user interface added by RM.
;      Feb 02 - Written by Richard Massey
;-

COMPILE_OPT idl2

;
; Maintain backwards compatibility
;
if not shapelets_structure_type(decomp,message=message) then message,message
polar_input=decomp.polar
if not keyword_set(order) then order=1

;
; Perform dilatation to higher than first order (only of use with /CART or /POLAR)
;
if fix(order) gt 1 then begin
  if (keyword_set(cartesian) + keyword_set(polar)) mod 2 then begin
    polar=decomp.polar
    cartesian=1-decomp.polar
  endif else shapelets_polar_convert,decomp,cartesian=cartesian,polar=polar,/SILENT
  shapelets_exponentiate_operations, "shapelets_dilate",                      $
    decomp, scaling, order, silent=silent,                                    $
    flux=flux, surface_brightness=surface_brightness, radius=radius,area=area,$
    cartesian=cartesian, polar=polar, extend=extend
  if keyword_set(maintain) then $
    shapelets_polar_convert,decomp,cartesian=1-polar_input,polar=polar_input,/SILENT
endif else begin
  
  ;
  ; Get scale factor right for various input styles
  ;
  if keyword_set(radius) then begin
    ; Treat input kappa as if we want to scale the radius of an object by
    ; that amount. Adjust stored kappa so that it is compatible with WL code.
    kappa=(scaling[0]-1.)/(scaling[0])
  endif else if keyword_set(area) then begin
    ; Treat input kappa as if we want to scale the area of an object by that
    ; amount. Adjust stored kappa so that it is compatible with WL code.
    kappa=1-1./sqrt(scaling[0])
    ;kappa=(scaling-1)/(scaling+1) ;linear equivalent
  endif else kappa=scaling[0]
  
  ;
  ; Increase n_max to contain new coefficients that might be created
  ;
  if keyword_set(extend) then shapelets_extend_nmax,decomp,extend
  
  ;
  ; Perform dilation
  ;
  if keyword_set(beta) then begin

    ; Perform dilation by simply rescaling beta
    decomp.beta=decomp.beta/(1-kappa)
    ; The enlargement also made it brighter
    if keyword_set(flux) then begin
      ; Conserve total flux
      decomp.coeffs=decomp.coeffs*(1-kappa)
    endif else begin
      ; Conserve surface brightness
      decomp.coeffs=decomp.coeffs/(1-kappa)
    endelse
    
  endif else if keyword_set(cartesian) then begin
  
    ; Obtain initial Cartesian shapelet coefficients
    shapelets_polar_convert,decomp,/CARTESIAN,/SILENT
    
    ; Perform dilation using Cartesian step operators
    n1=decomp.n1
    n2=decomp.n2
    if keyword_set(flux) then begin
      ; Conserve total flux
      new_coeffs=decomp.coeffs*(1.-kappa)
    endif else begin
      ; Conserve surface brightness
      new_coeffs=decomp.coeffs*(1.+kappa)
    endelse
  
    for i=0,decomp.n_coeffs-1 do begin
  
      ; Expand leftwards in Cartesian shapelet coefficient space
      j=where(n1 eq n1[i]+2 and n2 eq n2[i])
      if j[0] ne -1 then new_coeffs[i]=new_coeffs[i]-kappa/2.*sqrt((n1[i]+1)*(n1[i]+2))*decomp.coeffs[j[0]]
  
      ; Expand downwards in Cartesian shapelet coefficient space
      j=where(n1 eq n1[i] and n2 eq n2[i]+2)
      if j[0] ne -1 then new_coeffs[i]=new_coeffs[i]-kappa/2.*sqrt((n2[i]+1)*(n2[i]+2))*decomp.coeffs[j[0]]
  
      ; Expand rightwards in Cartesian shapelet coefficient space
      j=where(n1 eq n1[i]-2 and n2 eq n2[i])
      if j[0] ne -1 then new_coeffs[i]=new_coeffs[i]+kappa/2.*sqrt(n1[i]*(n1[i]-1))*decomp.coeffs[j[0]]
  
      ; Expand upwards in Cartesian shapelet coefficient space
      j=where(n1 eq n1[i] and n2 eq n2[i]-2)
      if j[0] ne -1 then new_coeffs[i]=new_coeffs[i]+kappa/2.*sqrt(n2[i]*(n2[i]-1))*decomp.coeffs[j[0]]
  
    endfor
    decomp.coeffs=new_coeffs
     
    ; Recover input format
    if keyword_set(maintain) then if polar_input then $
      shapelets_polar_convert,decomp,/POLAR,/SILENT
  
  endif else if keyword_set(polar) then begin
      
    ; Obtain polar shapelet coefficients
    shapelets_polar_convert,decomp,/POLAR,/SILENT
    shapelets_make_nvec,decomp.n_max,n,m,/POLAR
    
    ; Perform dilation using polar step operators
    if keyword_set(flux) then begin
      ; Conserve total flux
      new_coeffs=(1-kappa)*decomp.coeffs
    endif else begin
      ; Conserve surface brightness
      new_coeffs=(1+kappa)*decomp.coeffs
    endelse
    for i=0,decomp.n_coeffs-1 do begin
  
      ; Expand leftwards in polar shapelet coefficient space
      j=where(m eq m[i] and n eq n[i]-2)
      if j[0] ne -1 then new_coeffs[i] = new_coeffs[i] + $
  	( kappa/2.*sqrt((n[i]-m[i])*(n[i]+m[i]))*decomp.coeffs[j[0]] )
  
      ; Expand rightwards in polar shapelet coefficient space
      j=where(m eq m[i] and n eq n[i]+2)
      if j[0] ne -1 then new_coeffs[i] = new_coeffs[i] - $
  	( kappa/2.*sqrt((n[i]-m[i]+2)*(n[i]+m[i]+2))*decomp.coeffs[j[0]] )
  
    endfor
    decomp.coeffs=new_coeffs
    
    ; Convert back to Cartesian shapelet coefficients
    if keyword_set(maintain) then if not polar_input then $
      shapelets_polar_convert,decomp,/CARTESIAN,/SILENT
  
  endif else begin
  
    ; Perform dilation via a convolution with a delta-function
    alpha = decomp.beta        ; original beta
    beta  = alpha*1e-10        ; tiny! PSF is going to be a delta function
    gamma = alpha*(1.-kappa)   ; new beta
  
    ; Convolve with a delta function
    shapelets_read_psf,psf,beta,/SILENT
    shapelets_convolve,decomp,psf,gamma=gamma,n_max_g=decomp.n_max,silent=silent,maintain=maintain

    ; Enlarge
    decomp.beta=decomp.beta/(1.-kappa)

    ; The enlargement has just made it brighter, but by the wrong factor.
    if keyword_set(flux) then begin
      ; Conserve total flux
      decomp.coeffs=decomp.coeffs*(1-kappa)
    endif else begin
      ; Conserve surface brightness
      decomp.coeffs=decomp.coeffs/(1-kappa)
    endelse

  endelse
endelse

; Add operation to object's history record
if not keyword_set(nohistory) then shapelets_update_history,decomp,"Dilated by "+strtrim(scaling[0],2)

end
