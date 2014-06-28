;$Id: shapelets_focus_beta.pro, v2$
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
; ***********************************************************************
; ***********************************************************************
;
; NAME:
;      SHAPELETS_AMOEBA_CHISQ
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE: 
;      Compute Chi^2(beta) in the format required for the amoeba
;      minimisation procudure shapelets_1d_amoeba.pro. This version also 
;      finds the optimal centroid, centre_guess.
;
; INPUTS: 
;      BETA         - Shapelet scale size beta
; 
; COMMON VARIABLES:
;      PSTAMP       - Object structure from shapelets_extract_pstamp.pro
;      CENTRE_GUESS - Position of the centre of the shapelet basis functions
;      N_MAX        - Value of n_max to use
;      FLAG         - Flag for decomposition
;
; OPTIONAL INPUTS:
;      None.
;
; KEYWORD PARAMETERS:
;      None.
;
; OUTPUTS:
;      Reduced chi^2(beta) is returned.
;      Updates position of the centre of the shapelet basis functions
;      in the common block. Also possibly alters the flag.
;
; MODIFICATION HISTORY:
;      Jan 06 - Fake decomps made fully, via shapelets_create_decomp.pro by RM
;      Sep 05 - POLAR and DIAMOND options added by RM
;      Aug 05 - Upper geometic constraints changed by RM because Tmax too small
;      Apr 05 - Bug in flagging spotted by Stephane Paulin-Henriksson
;      Apr 04 - Sky background fitting implemented by RM
;      Apr 02 - PSF deconvolution implemented by Richard Massey
;      Apr 02 - Written by Alexandre Refregier


function shapelets_amoeba_chisq, BETA, N_MAX, $
                                 PSTAMP=pstamp,$
                                 PSF=psf,$
                                 SKY=sky,$
                                 DECOMP=decomp,$
                                 FOCUS=focus,$
                                 HIST=hist,$
                                 CENTRE_GUESS=centre_guess,$
                                 RECENTRE=recentre,$
                                 GAUSSIAN_RECENTRING=gaussian_recentring,$
                                 SILENT=silent,$
                                 FULL_FOCUS=full_focus,$
                                 POLAR=polar,$
                                 DIAMOND=diamond

COMPILE_OPT idl2, HIDDEN


; Recentre the basis functions on the previous decomposition's centre of light...
new_centre=decomp.x
if keyword_set(recentre) then $
  if keyword_set(decomp) then $
    if finite(decomp.chisq[1]) then $
      new_centre=shapelets_centroid(decomp,gaussian=gaussian_recentring)



; ...unless this will make it extend outside the postage stamp!
theta_minmax_geom=shapelets_geometric_constraints(pstamp, psf=psf, $
		   theta_min_geom=focus.theta_min_geom, centre_guess=new_centre)
if theta_minmax_geom[1] le theta_minmax_geom[0] then begin
  ; Geometric constraints have crossed
  message,"Centroid has fallen off the edge of the postage stamp",/info,noprint=silent
  focus.flag=focus.flag>10
  decomp=shapelets_create_decomp(n_max)
  decomp.beta=beta
  decomp.chisq=replicate(!values.f_infinity,2)
  decomp.x=new_centre
endif else if beta*(n_max+1)^(0.59) gt theta_minmax_geom[1] then begin
  ; Hits the upper geometric constraint
  message,"Centroid wandering towards the edge of the postage stamp",/info,noprint=silent
  focus.flag=focus.flag>9
  decomp=shapelets_create_decomp(n_max)
  decomp.beta=beta
  decomp.chisq=replicate(!values.f_infinity,2)
  decomp.x=new_centre
endif else if beta/sqrt(n_max+1) lt (theta_minmax_geom[0]>0) then begin
  ; Hits the lower geometric constraint
  message,"Small basis functions",/info,noprint=silent
  focus.flag=focus.flag>1
  decomp=shapelets_create_decomp(n_max)
  decomp.beta=beta
  decomp.chisq=replicate(!values.f_infinity,2)
  decomp.x=new_centre
endif else begin


  ; Decompose image into shapelets with the new parameters and calculate its chi^2 value
  centre_guess=new_centre
  decomp=shapelets_decomp(pstamp,beta,n_max,recomp=recomp,centre=centre_guess,$
         psf=psf,sky=sky,polar=polar,diamond=diamond)
  
  ; Warn of possible singular matrix
  if keyword_set(psf) then if round(sqrt( (decomp.beta^2*(decomp.n_max+1.)+psf.beta^2*(psf.n_max+1.)) / $
		      (decomp.beta^2*(psf.n_max+1.)+psf.beta^2*(decomp.n_max+1.)) * $
		      ((decomp.n_max+1.)*(psf.n_max+1.))  )-1)>0 lt decomp.n_max then focus.flag=focus.flag>2

 
;  print,'AMOEBA: theta',theta_min,theta_max
;  print,'AMOEBA: theta_minmax',theta_minmax_geom
;  print,'AMOEBA: n_pixels',pstamp.n_pixels
;  print,'AMOEBA: xc',xc


endelse

; Update the history records
if keyword_set(full_focus) then shapelets_update_focus, focus, decomp, silent=silent else $
  message,"n_max, beta, chi^2, x_c:"+$
  string(fix(decomp.n_max))+string(decomp.beta)+string(decomp.chisq[1])+$
  string(decomp.x[0])+string(decomp.x[1]),/info,noprint=silent

; Tell the world
return, decomp.chisq[1]

end





; ***********************************************************************
; ***********************************************************************
;+
; NAME:
;       SHAPELETS_FOCUS_BETA
;
; PURPOSE:
;       Find the optimal beta and centroid centre_guess for the decomposition of an
;       image by minimising Chi^2. This is done using a 1-dimensional Amoeba
;       search. The maximum shapelet order n_max is assumed to be known.
;
; CATEGORY:
;       Shapelets.
;
; INPUTS:
;       PSTAMP          - IDL postage stamp structure extracted from a large 
;                         image via shapelets_sexcat2pstamp.pro.
;       N_MAX           - Shapelet truncation order for decomposition.
;
; OPTIONAL INPUTS:
;       BETA_GUESS      - Starting point for shapelet scale size beta.
;       CENTRE_GUESS    - Starting point for centre of shapelet basis functions.
;       AMOEBA_SCALE    - Controls the size of the amoeba's steps.
;       MAX_ITERATIONS  - Maximum number of iterations.
;       PSF             - A shapelet decomp structure representing the local PSF.
;                         The shapelet model of the galaxy will be deconvolved
;                         from this.
;       FOCUS           - Structure containing the prior history of the search.
;                         If this is specified, search parameters contained in 
;                         it also override those specified by hand (see below).
;       BETA_TOLERANCE  - How accurate we want the final beta to be.
;       THETA_MIN_GEOM  - Minimum scale of oscillations that data could possibly
;                         contain.
;
; KEYWORD PARAMETERS:
;       SILENT          - Suppress output to screen.
;       SKY         	  - 1: fit sky background with a constant value around object.
;                         2: fit sky gradient with a plane around object.
;                         DEFAULT: no sky subtraction.
;       NON1            - a_01 and a_10 are forced to be zero, Kuijken-like.
;       FULL_FOCUS      - Record every single attempt made at decomposition during
;                         iteration. Makes things a little slower, but plots nicer. 
;
; OUTPUTS:
;       DECOMP          - IDL shapelet structure containing a shapelet
;                         decomposition of pstamp, using the optimum scale size.
;
; OPTIONAL OUTPUTS:
;       FOCUS           - Structure containing the (updated) history of the search.
;       RECOMP          - Repixellated version of the optimum decomposition.
;
; TO DO:
;       Cut the nonessential calls to shapelets_update_focus using FULL_FOCUS.
;
; MODIFICATION HISTORY:
;       Sep 05 - POLAR and DIAMOND options added by RM
;       Jul 05 - Significant load of bug fixes and clean up by RM.
;       Jan 05 - Optimisation & bugs fixed in centroid iteration by RM
;       May 04 - Sky background fitting added by RM
;       Apr 02 - PSF deconvolution added by Richard Massey
;       Mar 02 - Written by Alexandre Refregier
;-

function shapelets_focus_beta, PSTAMP, N_MAX, FOCUS,                    $
                               RECOMP=recomp,			                      $
                               CENTRE_GUESS=centre_guess,	              $
                               BETA_GUESS=beta_guess,		                $
                               BETA_TOLERANCE=beta_tolerance,	          $
                               MAX_ITERATIONS=max_iterations,	          $
                               THETA_MIN_GEOM=theta_min_geom,	          $
                               AMOEBA_SCALE=amoeba_scale,	              $
                               FULL_FOCUS=full_focus,		                $
                               GAUSSIAN_RECENTRING=gaussian_recentring, $
                               PSF=psf, 			                          $
                               SKY=sky, 			                          $
                               POLAR=polar,                             $
                               DIAMOND=diamond,                         $
                               SILENT=silent

COMPILE_OPT idl2

; Initial declarations
if not keyword_set(max_iterations) then max_iterations=50
message,"Finding beta that minimises chi^2, and simultaneously recentering the decomposition",/info,noprint=silent


; Initialise a focus structure to contain the history of this search
if not keyword_set(focus) then $
  focus=shapelets_create_focus(NAME=pstamp.name,                      $
                               COMMENT="during beta/centroid search", $
                               BETA_TOLERANCE=beta_tolerance,         $
                               THETA_MIN_GEOM=theta_min_geom,         $
                               CHISQ_TARGET=chisq_target,             $
                               CHISQ_TOLERANCE=chisq_tolerance,       $
                               CHISQ_FLATNESS=chisq_flatness)


; Guess starting points for beta and centre_guess if not provided
; (Could be more intelligent here, since n_max is a required input).
if not keyword_set(theta_min_geom) then theta_min_geom=focus.theta_min_geom
if focus.n_iterations gt 0 then begin
  if not keyword_set(beta_guess) then beta_guess=focus.beta
  if not keyword_set(centre_guess) then centre_guess=focus.x
endif else begin
  guess=shapelets_guess_nmaxbeta(pstamp, psf=psf, theta_min_geom=focus.theta_min_geom)
  if not keyword_set(beta_guess) then begin
    beta_guess=guess.beta
    message,'Guessing beta of '+strtrim(string(beta_guess),2),/info,noprint=silent
  endif
  if not keyword_set(centre_guess) then begin
    centre_guess=guess.centre_guess
    message,'Guessing centroid of ('+strtrim(string(centre_guess[0]),2)+','+$
          strtrim(string(centre_guess[1]),2)+')',/info,noprint=silent
  endif 
endelse
if not keyword_set(beta_tolerance) then beta_tolerance=focus.beta_tolerance


; Have a first attempt at decomposition, using the guesses
decomp=shapelets_decomp(pstamp,beta_guess,n_max,recomp=recomp,$
       centre=centre_guess,psf=psf,sky=sky,polar=polar,diamond=diamond)


; Record this first attempt in the focus structure
shapelets_update_focus, focus, decomp, silent=silent
best_decomp=decomp


; Recentre quickly, before placing all of our feet
; (intentionally updates centre_guess variable, even though it has wider scope)
centre_guess=shapelets_centroid(decomp)


; Set scale of steps along beta direction during the amoeba search
if not keyword_set(amoeba_scale) then amoeba_scale=4.3
beta_max=min([centre_guess[0],pstamp.im_ran[1]-pstamp.im_ran[0]-centre_guess[0],$ 
              centre_guess[1],pstamp.im_ran[3]-pstamp.im_ran[2]-centre_guess[1]])$  ; th_max
		  /sqrt(n_max+1)
beta_min=(focus.theta_min_geom*sqrt(n_max+1))
direction=sign(-2.*beta_guess+beta_max+beta_min) ; But which is it closer to? 
scale=direction*beta_guess/amoeba_scale          ; Go in opposite direction!


; Hold centre of basis functions steady while we place our feet and thus
;  initialise the function at the simplex points
simplex_feet=beta_guess+[0.,scale]
chisq_feet=fltarr(2)
for i=0,1 do begin
  chisq_feet[i]=shapelets_amoeba_chisq(simplex_feet[i],n_max,RECENTRE=0,full_focus=full_focus,decomp=decomp,pstamp=pstamp,psf=psf,sky=sky,silent=silent,focus=focus,hist=hist,centre_guess=centre_guess)
  if chisq_feet[i] lt best_decomp.chisq[1] then best_decomp=decomp
endfor
ncalls=2


; Iterate by moving the simplex
converged=0B
while ncalls lt max_iterations and not converged do begin

  ; Sort simplex
  ilo=(chisq_feet[0] gt chisq_feet[1])
  ihi=1-ilo

  ; First reflect about lowest point
  simplex_try=simplex_feet[ilo]-1.*(simplex_feet[ihi]-simplex_feet[ilo])
  chisq_try=shapelets_amoeba_chisq(simplex_try,n_max,/RECENTRE,GAUSSIAN_RECENTRING=gaussian_recentring,full_focus=full_focus,pstamp=pstamp,psf=psf,sky=sky,silent=silent,decomp=decomp,focus=focus,hist=hist,centre_guess=centre_guess)
  ncalls+=1
  if chisq_try lt chisq_feet[ihi] then begin
    chisq_feet[ihi]=chisq_try
    simplex_feet[ihi]=simplex_try
  endif

  ; Is new point lower than lowest point?
  if chisq_try le chisq_feet[ilo] then begin  ; try further expansion if so
    best_decomp=decomp
    ;simplex_try=simplex_feet[ilo]+2*(simplex_feet[ihi]-simplex_feet[ilo])
    simplex_try=simplex_feet[ilo]+1.95*(simplex_feet[ihi]-simplex_feet[ilo])
    chisq_try=shapelets_amoeba_chisq(simplex_try,n_max,/RECENTRE,GAUSSIAN_RECENTRING=gaussian_recentring,full_focus=full_focus,pstamp=pstamp,psf=psf,sky=sky,silent=silent,decomp=decomp,focus=focus,hist=hist,centre_guess=centre_guess)
    ncalls+=1
    if chisq_try lt chisq_feet[ihi] then begin
      best_decomp=decomp
      chisq_feet[ihi]=chisq_try
      simplex_feet[ihi]=simplex_try
    endif

  ; If not, try contraction about lowest point
  endif else begin
    chisq_save=chisq_feet[ihi]
    ;simplex_try=simplex_feet[ilo]+0.5*(simplex_feet[ihi]-simplex_feet[ilo])
    simplex_try=simplex_feet[ilo]+0.55*(simplex_feet[ihi]-simplex_feet[ilo])
    chisq_try=shapelets_amoeba_chisq(simplex_try,n_max,/RECENTRE,GAUSSIAN_RECENTRING=gaussian_recentring,full_focus=full_focus,pstamp=pstamp,psf=psf,sky=sky,silent=silent,decomp=decomp,focus=focus,hist=hist,centre_guess=centre_guess)
    ncalls+=1
    if chisq_try lt chisq_feet[ihi] then begin
      if chisq_try lt chisq_feet[ilo] then best_decomp=decomp
      chisq_feet[ihi]=chisq_try
      simplex_feet[ihi]=simplex_try
    endif else begin   ; if latest point is worse then contract about lowest point
      message,"Have found a local maximum",/info,noprint=silent
      ;simplex_feet[ilo]=simplex_feet[ilo]-.5*(simplex_feet[ilo]-simplex_feet[ihi])
      ;chisq_feet[ilo]=chisq_try
      simplex_feet[ihi]=simplex_try
      chisq_feet[ihi]=chisq_try
    endelse
  endelse

  ; Test for tolerance
  d=total(abs(chisq_feet))
  if d ne 0 then delta_chisq=abs(chisq_feet[ihi]-chisq_feet[ilo])*2./d else delta_chisq=focus.chisq_abs_target/2.
  if delta_chisq lt focus.chisq_flatness*focus.chisq_target then begin  ; done?
    converged=1B
  endif else if abs((simplex_feet[ihi]-simplex_feet[ilo])/simplex_feet[ihi]) lt focus.beta_tolerance then begin
    message,"Iterations in beta smaller than tolerance",/info,noprint=silent
    ;focus.flag=focus.flag>4
    converged=1B
  endif
  
endwhile

; In case the convergence was not achieved with max_iterations steps
if not converged then begin
  message,"Amoeba failed to converge within "+strtrim(string(max_iterations),2)+" steps",$
          /info,noprint=silent
  focus.flag=focus.flag>6
endif


; Try to take advantage of subsequent recentering. If this makes the
; decomposition better, keep it. If we can't recover as good a fit, discard it.
; Bear in mind that, at a new centroid, this decomposition could be worse!
if centre_guess[0] ne best_decomp.x[0] or centre_guess[1] ne best_decomp.x[1] then begin
  decomp=shapelets_decomp(pstamp,best_decomp.beta,n_max,recomp=recomp,$
    centre=centre_guess,psf=psf,polar=polar,diamond=diamond)
  if keyword_set(full_focus) then shapelets_update_focus, focus, decomp, silent=silent else $
    message,"  n_max, beta, chi^2, x_c:"+$
    string(fix(decomp.n_max))+string(decomp.beta)+string(decomp.chisq[1])+$
    string(decomp.x[0])+string(decomp.x[1]),/info,noprint=silent
  if decomp.chisq[1] le decomp.chisq[1] then best_decomp=decomp
endif


; Keep the decomp with the lowest chi squared
decomp=best_decomp


; Make sure history records are up to date with at least the minimum of information
shapelets_update_focus, focus, decomp, silent=silent


; Tell the world
return,decomp

end
