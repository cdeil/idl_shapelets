function shapelets_quadrupole, decomp,              $
                               POLAR=polar,         $
                               CARTESIAN=cartesian, $
                               N_MAX=n_max,         $
                               ERROR=error,         $
                               FLUX=flux,           $
                               RSQUARED=rsquared,   $
                               ELLIPTICITY=ellipticity

;$Id: shapelets_quadrupole.pro, v2$
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
;      SHAPELETS_QUADRUPOLE
;
; PURPOSE:
;      Calculates the unweighted quadrupole moments of a shapelet model.
;
; EXAMPLE USE:
;       result=shapelets_quadrupole(decomp) & print,result
;
; INPUTS:
;       decomp - A Cartesian shapelet decomposition
;
; OPTIONAL INPUTS: 
;        n_max - Works out the sum only to n=nmax.
;
; KEYWORD PARAMETERS:
;        /CART - Perform operation in Cartesian shapelet space (DEFAULT).
;       /POLAR - Perform operation in polar shapelet space.
;       /ERROR - Total object flux is returned in this variable.
;
; OUTPUTS:
;   quadrupole - A 2x2 matrix containing [Q11,Q12] 
;                                        [Q21,Q22]
;                If /ERROR is set, it becomes [Q11,Q12,Q11_error,Q12_error]
;                                             [Q21,Q22,Q21_error,Q22_error]
;
; OPTIONAL OUTPUTS:
;         flux - Total object flux is optionally returned in this variable.
;                If /ERROR is set, it becomes [flux,flux_error]
;     rsquared - Object's rms size is optionally returned in this variable.
;                If /ERROR is set, it becomes [rsquared,rsquared_error]
;  ellipticity - Object's unweighted ellipticity is optionally returned in
;                this variable, using complex {e1,e2} notation.
;                If /ERROR is set, it becomes [{e1,e2},{e1_error,e1_error}]
;
; NOTES:
;       Errors will not be quite correct if we try to perform the calculation
;       in polar shapelet space but do not have the full covariance matrix
;       of Cartesian shapelet coefficients necessary to transform from one to
;       the other. Are all factors of 2 correct in error calculations?
;
; MODIFICATION HISTORY:
;      Sep 05 - Input value of n_max preserved by RM.
;      Oct 03 - Error calculation added by RM.
;      Mar 03 - Polar shapelet version written by Richard Massey.
;        2002 - Cartesian shapelet version written by Tzu-Ching Chang.
;-

COMPILE_OPT idl2

;
; SET DEFAULTS
;
; Maintain backwards compatibility
if not shapelets_structure_type(decomp,message=message) then message,message
; Decide whether default method should be in Cartesian or polar shapelet space
if ((keyword_set(cartesian)+keyword_set(polar)) mod 2) eq 0 then cartesian=1-decomp.polar
; Truncate calculation?
if n_elements(n_max) eq 0 then n_max_local=decomp.n_max else n_max_local=n_max

; If object has no quadrupole m=2 coefficients, calculation is easy!
if n_max_local lt 2 then begin
  Q=[[0.,0.],[0.,0.]]
  ellipticity=complex(0.,0.)
  if keyword_set(error) then begin
    Q=[[[Q]],[[0.,0.],[0.,0.]]]
    ellipticity=[ellipticity,complex(0.,0.)]
  endif
  if arg_present(rsquared) then begin
    rsquared=shapelets_rsquared(decomp,polar=polar,cart=cart,error=error,flux=flux)
  endif else if arg_present(rsquared) then begin
    flux=shapelets_flux(decomp,polar=polar,cart=cart,error=error)
  endif
  return,Q
endif


; Really do need to calculate the moments...
; Remember input variable, so that it can be restored later
decomp_input=decomp
rootpi=sqrt(!pi)
if keyword_set(cartesian) then begin

  ;
  ; WORK IN CARTESIAN SHAPELET SPACE  
  ;
  ; Convert coefficients to Cartesian shapelet format if not already like that
  shapelets_polar_convert,decomp,/CARTESIAN,/SILENT
  shapelets_make_nvec,n_max_local,n1,n2,n_coeffs

  ; Initialise variables to contain coefficient weights
  fweight=fltarr(n_coeffs)
  Q11weight=fltarr(n_coeffs)
  Q12weight=fltarr(n_coeffs)
  Q22weight=fltarr(n_coeffs)
  
  ; Calculate weights
  for i=0,n_coeffs-1 do begin
    if n1[i] mod 2 eq 0 and n2[i] mod 2 eq 0 then begin
      fweight[i]  =2.^(.5*(2.-n1[i]-n2[i]))*                             $
                   sqrt(factorial(n1[i])*factorial(n2[i]))/              $
                   factorial(n1[i]/2)/factorial(n2[i]/2)
      Q11weight[i]=2.^(.5*(2.-n1[i]-n2[i]))*(1.+2.*n1[i])*               $
                   sqrt(factorial(n1[i])*factorial(n2[i]))/              $
                   factorial(n1[i]/2)/factorial(n2[i]/2)
      Q22weight[i]=2.^(.5*(2.-n1[i]-n2[i]))*(1.+2.*n2[i])*               $
                   sqrt(factorial(n1[i])*factorial(n2[i]))/              $
                   factorial(n1[i]/2)/factorial(n2[i]/2)
    endif else if n1[i] mod 2 eq 1 and n2[i] mod 2 eq 1 then begin
      Q12weight[i]=2.^(.5*(2.-n1[i]-n2[i]))*sqrt((n1[i]+1.)*(n2[i]+1.))* $
                   sqrt(factorial(n1[i]+1)*factorial(n2[i]+1))/          $
                   factorial((n1[i]+1)/2)/factorial((n2[i]+1)/2)
    endif
  endfor
  
  flux= rootpi * decomp.beta   * total(fweight*decomp.coeffs)
  Q11 = rootpi * decomp.beta^3 * total(Q11weight*decomp.coeffs)
  Q12 = rootpi * decomp.beta^3 * total(Q12weight*decomp.coeffs)
  Q22 = rootpi * decomp.beta^3 * total(Q22weight*decomp.coeffs)
  Q   = [[Q11,Q12],[Q12,Q22]]
  
  ; Might as well calculate some related quantities very quickly
  ellipticity=complex(Q11-Q22,2*Q12)/(Q11+Q22)
  rsquared=(Q11+Q22)/flux
  
  if keyword_set(error) then begin
    flux_error = rootpi * decomp.beta   * sqrt(total((fweight*decomp.coeffs_error)^2))
    Q11_error  = rootpi * decomp.beta^3 * sqrt(total((Q11weight*decomp.coeffs_error)^2))
    Q12_error  = rootpi * decomp.beta^3 * sqrt(total((Q12weight*decomp.coeffs_error)^2))
    Q22_error  = rootpi * decomp.beta^3 * sqrt(total((Q22weight*decomp.coeffs_error)^2))
    ; Rather than use Q11_error^2+Q22_error^2, they have covariance, so:
    trace_error= rootpi * decomp.beta^3 * sqrt(total(((Q11weight+Q22weight)*decomp.coeffs_error)^2)) 
    ;trace_error= sqrt(Q11_error^2+Q22_error^2)/sqrt(2.) ; gives the same thing

    rsquared_error=abs(rsquared)*$
                   sqrt( (trace_error^2)/(Q11+Q22)^2  $
                        +(flux_error/flux)^2  )
    ; SPH July 05
    e1_error=abs(float(ellipticity)) * trace_error *  $
                   sqrt( 1./(Q11+Q22)^2 + 1./(Q11-Q22)^2 )
                   ; Not strictly true, but other part is highly unstable
                   ;  (and this agrees with the calculation in polar space)
;    e1_error=abs(float(ellipticity)) * trace_error *  $
;                   sqrt( 1./(Q11+Q22)^2  )    ; + 1./(Q11-Q22)^2 )
;                   ; Not strictly true, but other part is highly unstable
;                   ;  (and this agrees with the calculation in polar space)
    e2_error=abs(imaginary(ellipticity))*   $
                   sqrt( (Q12_error/Q12)^2  $
                        +(trace_error^2)/(Q11+Q22)^2 )
    ellipticity_error=complex(e1_error,e2_error)
    ; Concatenate variables to return
    flux=[flux,flux_error]
    rsquared=[rsquared,rsquared_error]
    ellipticity=[ellipticity,ellipticity_error]
    Q=[[[Q]],[[Q11_error,Q12_error],[Q12_error,Q22_error]]]
  endif

endif else begin
  
  ;
  ; WORK IN POLAR SHAPELET SPACE 
  ;
  ; Convert coefficients to polar shapelet format if not already like that
  shapelets_polar_convert,decomp,/POLAR,/SILENT


  ; BUG? m is then populated witho all zeros!....

  ; Get n and m vectors of shapelet coefficients
  shapelets_make_nvec,n_max_local,n,m,/POLAR
  i0=where(m eq 0,n_m0)
  i2=where(m eq 2,n_m2)

  ; Work out various quantities in unweighted ellipticities
  factor = 4*rootpi*decomp.beta^3
  flux   = 2*rootpi*decomp.beta*total((float(decomp.coeffs))[i0])
  ;if keyword_set(gaussian) then begin
  ;  ; Don't need to bother! If working out Gaussian-weighted ellipticities, just set n_max=2
  ;  numerator   = decomp.coeffs[where(m eq 2 and n eq 2)]
  ;  denominator = decomp.coeffs[where(m eq 0 and n eq 2)]
  ;endif else begin
    ; Equation (59) in Shapelets: III. (demoninator of e is FR^2, the trace of Q)
    numerator   = factor*total(( decomp.coeffs*sqrt((n)*(n+2)) )[i2])
    denominator = factor*total((float(decomp.coeffs)*(n+1))[i0])
  ;endelse
  ellipticity = numerator/denominator
  rsquared    = denominator/flux
  ; Extract quadrupole moments
  Q11=(float(numerator)+denominator)/2.
  Q22=(denominator-float(numerator))/2.
  Q12=imaginary(numerator)/2.
  Q=[[Q11,Q12],[Q12,Q22]]


  if keyword_set(error) then begin
    if size(decomp.coeffs_error,/n_dimensions) ne 1 then message,$
      "Warning: error calculation using covariance matrix is not yet written!"
    error_p    = decomp.coeffs_error
    error_pr   = float(error_p)
    error_pi   = imaginary(error_p)
    flux_error = rootpi*decomp.beta*2*sqrt(total((error_pr[i0])^2))
    
    numerator_variancer = factor^2*total((error_pr^2*(n*(n+2)))[i2])
    numerator_variancei = factor^2*total((error_pi^2*(n*(n+2)))[i2])
    denominator_variance= factor^2*total(((error_pr*(n+1))^2)[i0])
    
    Q11_error = sqrt(numerator_variancer+denominator_variance)/2.
    Q22_error = Q11_error ; Can only assume that these are the same, but  it's not necessarily true!
    Q12_error = sqrt(numerator_variancei)/2.

  ; This would, strictly speaking, be the best way of calculating e_errors, 
  ; given the data we have calculated so far:
   e1_error  = abs(float(ellipticity))*$
               sqrt( numerator_variancer/(float(numerator))^2  $
                    +denominator_variance/denominator^2  )
   e2_error  = abs(imaginary(ellipticity))*$
               sqrt( numerator_variancei/(imaginary(numerator))^2 $
                    +denominator_variance/denominator^2 )
  ; However, the numerator involves the cancelling of small numbers and
  ; division by some very different powers to give an average result. The
  ; next method is only approximate, but numerically more stable for e1_error:
  ;  e1_error=abs(float(ellipticity))*$
  ;                 sqrt( (Q11_error^2+Q22_error^2)/(Q11-Q22)^2  $
  ;                      +(Q11_error^2+Q22_error^2)/(Q11+Q22)^2 )
  ;  e2_error=abs(imaginary(ellipticity))*$
  ;                 sqrt( (Q12_error/Q12)^2  $
  ;                      +(Q11_error^2+Q22_error^2)/(Q11+Q22)^2 )
  ; NB: it is equivalent to averaging the variances of the numerator and the
  ; denominator, and may be written in terms of them as:
  ; e1_error  = abs(float(ellipticity))*$
  ;             sqrt( 0.5*(numerator_variancer+denominator_variance)/(float(numerator))^2  $
  ;                  +0.5*(numerator_variancer+denominator_variance)/denominator^2 )
  ; e2_error  = abs(imaginary(ellipticity))*$
  ;             sqrt( (numerator_variancei)/(imaginary(numerator))^2 $
  ;                  +0.5*(numerator_variancer+denominator_variance)/denominator^2 )
    ellipticity_error=complex(e1_error,e2_error)

  ; Similarly for R^2, we have two methods:
    rsquared_error=abs(rsquared)*$
                   sqrt( denominator_variance/(denominator^2) $
                        +(flux_error/flux)^2 )
  ; However, in this case, we only need the denominator, which is the trace
  ; of Q and much larger than the numerator, which is the difference.
  ; It is indeed numerically stable, and agrees with the method in
  ; shapelets_rsquared.pro. For consistency, we will therefore use this 
  ; method. However, note that it gives a ~15% different answer to
  ; the alternative below. NO! I'VE ONLY HAD A ROUGH GUESS AT Q11_error
  ; AND Q22_error, ASSUMING THEY ARE THE SAME - BUT THEY'RE NOT!:
  ; rsquared_error=abs(rsquared)*$
  ;                sqrt( (Q11_error^2+Q22_error^2)/((Q11+Q22)^2)  $
  ;                     +(flux_error/flux)^2  )

    ; Concatenate moments and their errors, to return as one array
    flux=[flux,flux_error]
    rsquared=[rsquared,rsquared_error]
    ellipticity=[ellipticity,ellipticity_error]
    Q=[[[Q]],[[Q11_error,Q12_error],[Q12_error,Q22_error]]]
   
  endif

endelse

; Restore input variable
decomp=decomp_input

return,Q

end
