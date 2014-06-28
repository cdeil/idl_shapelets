;$Id: shapelets_focus_nmax.pro, v2$
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
;       SHAPELETS_FOCUS_NMAX
;
; PURPOSE:
;       Find the optimal n_max for the decomposition of an image by exploring
;       different values until Chi^2 is equal to a set value Chi2_0.
;       The shapelet scale beta is assumed to be known.
;
; CATEGORY:
;       Shapelets.
;
; INPUTS:
;       PSTAMP          - IDL postage stamp structure extracted from a large 
;                         image via shapelets_sexcat2pstamp.pro.
;       BETA            - Shapelet scale for decomposition.
;
; OPTIONAL INPUTS:
;       N_MAX_RANGE     - Range of n_max values to explore.
;       CENTRE_GUESS    - Starting point for centre of shapelet basis functions.
;       PSF             - A shapelet decomp structure representing the local PSF.
;                         The shapelet model of the galaxy will be deconvolved
;                         from this.
;       FOCUS           - Structure containing the prior history of the search.
;                         If this is specified, search parameters contained in 
;                         it also override those specified by hand (see below).
;       CHISQ_TARGET	  - Ideal value of reduced chi^2 for the residual image.
;       CHISQ_TOLERANCE - Acceptable accuracy for chi^2, in units of the rms.
;       CHISQ_FLATNESS  - Minimum difference in chi^2 between two decompositions
;                         with n_max differing by two to trigger the flatness
;                         constraint in shapelets_focus_nmax
;       THETA_MIN_GEOM  - Minimum scale of oscillations that data could possibly
;                         contain.
;
; KEYWORD PARAMETERS:
;       SILENT          - Operate silently.
;       VERBOSE         - Operate noisily.
;       SKY             - 1: fit sky background with a constant value around object.
;                         2: fit sky gradient with a plane around object.
;                         DEFAULT: no sky subtraction.
;       FULL_FOCUS      - Record every single attempt made at decomposition during
;                         iteration. Makes things a little slower, but plots nicer. 
;
; OUTPUTS:
;       DECOMP          - IDL shapelet structure containing a shapelet decomposition
;                         of pstamp, to sufficient order to achieve the target chi sqaured
;
; OPTIONAL OUTPUTS:
;       FOCUS           - Structure containing the (updated) history of the search.
;       RECOMP          - Repixellated version of the optimum decomposition.
;
; MODIFICATION HISTORY:
;       Jan 06 - Fake decomps made fully, via shapelets_create_decomp.pro by RM
;       Sep 05 - POLAR and DIAMOND options added by RM
;       Aug 05 - Upper geometic constraints changed by RM because Tmax no good
;       Jul 05 - Significant load of bug fixes and clean up by RM
;       Apr 05 - Iteration stopped if chi^2 is getting worse by RM
;       Apr 05 - Iteration forced to stop at n_max Stephane Paulin-Henriksson
;       Jan 05 - Tidied up to avoid awkward history structure by RM 
;       Jan 05 - Beta adjustments to follow geometric constraints by RM 
;       Dec 04 - Bug fixed in final search around nearby values of n_max by RM
;       May 04 - Sky background fitting added by RM
;       Apr 02 - PSF deconvolution added by Richard Massey
;       Mar 02 - Written by Alexandre Refregier
;-

function shapelets_focus_nmax, PSTAMP, BETA, FOCUS,                 $
                               RECOMP=recomp, 	          	        $
                               N_MAX_RANGE=n_max_range,	            $
                               CENTRE_GUESS=centre_guess,	          $
                               CHISQ_TARGET=chisq_target,	          $
                               CHISQ_TOLERANCE=chisq_tolerance,     $
                               CHISQ_FLATNESS=chisq_flatness,	      $
                               THETA_MIN_GEOM=theta_min_geom,	      $
                               BETA_TOLERANCE=beta_tolerance,	      $
                               FULL_FOCUS=full_focus,               $
                               DIAMOND=diamond,                     $
                               POLAR=polar,                         $
                               PSF=psf,                             $
                               SKY=sky,                             $
                               SILENT=silent

COMPILE_OPT idl2

; Initialise focus structure to contain the history of this search
if not keyword_set(focus) then $
  focus=shapelets_create_focus(NAME=pstamp.name,                    $
                               COMMENT="during n_max optimisation", $
                               BETA_TOLERANCE=beta_tolerance,       $
                               THETA_MIN_GEOM=theta_min_geom,       $
                               CHISQ_TARGET=chisq_target,           $
                               CHISQ_FLATNESS=chisq_flatness,       $
                               CHISQ_TOLERANCE=chisq_tolerance)


; Determine starting points for iteration
if not keyword_set(n_max_range) then n_max_range=[2,20]
if max(n_max_range mod 2) then message,"Odd choice of an odd value in the n_max limits",/info,noprint=silent
if not keyword_set(centre_guess) then begin
  if focus.n_iterations gt 0 then centre_guess=focus.x else begin
    guess=shapelets_guess_nmaxbeta(pstamp, psf=psf, theta_min_geom=focus.theta_min_geom)
    centre_guess=guess.centre
    message,'Guessing centroid of ('+strtrim(string(centre_guess[0]),2)+','+$
          strtrim(string(centre_guess[1]),2)+')',/info,noprint=silent
  endelse
endif


; Compute allowed range of n_max given the size of the image and the pixel size
n_max_range=n_max_range[sort(n_max_range)]
theta_minmax_geom=shapelets_geometric_constraints(pstamp, psf=psf, theta_min_geom=focus.theta_min_geom, centre_guess=centre_guess)
if theta_minmax_geom[1] le (theta_minmax_geom[0]>0) then begin
  message,"Postage stamp not big enough!",/info,noprint=silent
  focus.flag=focus.flag>10
endif
n_max_geom=fix(min([(theta_minmax_geom[1]/beta)^1.69-1.,(beta/theta_minmax_geom[0])^2-1.]))
n_max_range_soft=[n_max_range[0],(n_max_range[1]<n_max_geom)>n_max_range[0]]


; Determine goal
; rms of chi^2 distribution is sqrt{2n}; rms of reduced chi^2 distribution is sqrt{2/n}
focus.chisq_abs_target=focus.chisq_target+sqrt(2./n_elements(pstamp.image))*focus.chisq_tolerance
message,"Finding n_max to achieve chi^2 better than "+strtrim(focus.chisq_abs_target,2),/info,noprint=silent


; Initialise some focus.flags to track the progress of our iteration
converged=0B
flat=0B
out=0B
bounce=0B
worsening=0B


; Start by trying a decomposition at the lowest considered value of n_max
n_max=n_max_range_soft[0]
;n_max=(n_max_range_soft[0]-(n_max_range_soft[0] mod 2))>0
decomp=shapelets_decomp(pstamp,beta,n_max,recomp=recomp,centre=centre_guess,$
  psf=psf,sky=sky,polar=polar,diamond=diamond)
old_chisq=decomp.chisq[1]
if keyword_set(full_focus) then shapelets_update_focus, focus, decomp, silent=silent else $
  message,"  n_max, beta, chi^2, x_c:"+$
  string(decomp.n_max)+string(decomp.beta)+string(decomp.chisq[1])+$
  string(decomp.x[0])+string(decomp.x[1]),/info,noprint=silent

; Gradually increase n_max 
while not converged and not flat and not out and not worsening do begin
    
  ; Store previous decomposition
  old_n_max=n_max
  old_chisq=decomp.chisq[1]
  old_decomp=decomp
  old_recomp=recomp
  
  ; Increment n_max
  n_max=n_max+2
  
  if n_max gt n_max_range[1] then begin   ; Stop at hard (user-input) limit
    
    out=1B
    ;print,'a'     
  
  endif else begin
    
    ; Bounce off geometric constraints
    beta_min=theta_minmax_geom[0]*sqrt(n_max+1)
    beta_max=theta_minmax_geom[1]/sqrt(n_max+1)

    if beta_min gt beta_max then begin
      
      out=1B
      ;print,'b'     
   
    endif else begin
      if (n_max gt n_max_range_soft[1]) or (beta lt beta_min) or (beta gt beta_max) then begin  ; bounce obliquely off th_min_geom/max_geom: nudge beta 
   	bounce=1B
   	beta=beta_min>beta<beta_max
      endif
      
      ; Try decomposition at the new n_max  
      decomp=shapelets_decomp(pstamp,beta,n_max,recomp=recomp,centre=centre_guess,$
   	psf=psf,sky=sky,polar=polar,diamond=diamond)
      chisq=decomp.chisq[1]
      if keyword_set(full_focus) then shapelets_update_focus, focus, decomp, silent=silent else $
        message,"  n_max, beta, chi^2, x_c:"+$
        string(decomp.n_max)+string(decomp.beta)+string(decomp.chisq[1])+$
        string(decomp.x[0])+string(decomp.x[1]),/info,noprint=silent

      ; Warn of possible singular matrix
      if keyword_set(psf) then if round(sqrt( (decomp.beta^2*(decomp.n_max+1.)+psf.beta^2*(psf.n_max+1.)) / $
	    	          (decomp.beta^2*(psf.n_max+1.)+psf.beta^2*(decomp.n_max+1.)) * $
		         ((decomp.n_max+1.)*(psf.n_max+1.))  )-1)>0 lt decomp.n_max then focus.flag=focus.flag>2

      ; Tests for convergence
      if chisq lt (focus.chisq_abs_target) then begin
   	converged=1B  ; Bingo!
      endif else if abs(old_chisq-chisq) lt focus.chisq_flatness then begin
   	flat=1B
      endif else if chisq ge old_chisq then begin
   	worsening=1B  ; We are getting worse (somethimes happens with small pstamps)
        message,"Decomposition is getting worse...",/info,noprint=silent
      endif
  
    endelse
  endelse

endwhile

     
; Have we bounced off the geometrical constraints at any point?
if bounce then focus.flag=focus.flag>1  


; Has the object not been able to be fully modelled as shapelets?
if worsening then focus.flag=focus.flag>5


; Decide which of the neighbouring values is optimal
if out then begin

  ; We reached an n_max limit, either from user-set vlaue or from overlapping
  ; geometric constraints (only likely with very small postage stamps)
  message,'n_max reached maximum',/info,noprint=silent

;  ; Update focus history
;  focus.flag=focus.flag>8
;  print,"n_max",n_max
;  decomp={beta:beta,n_max:n_max,chisq:replicate(!values.f_infinity,2),x:centre_guess}
;  if keyword_set(full_focus) then shapelets_update_focus, focus, decomp, silent=silent else $
;  message,"  n_max, beta, chi^2, x_c:"+$
;    string(fix(decomp.n_max))+string(decomp.beta)+string(decomp.chisq[1])+$
;    string(decomp.x[0])+string(decomp.x[1]),/info,noprint=silent


  ; Try intermediate n_max increase
  if keyword_set(diamond) then begin
    n_max_opt=old_n_max
    ;decomp=old_decomp
    ;recomp=old_recomp
  endif else begin
    n_max=old_n_max+1
    beta_min=theta_minmax_geom[0]*sqrt(n_max+1)
    beta_max=theta_minmax_geom[1]/sqrt(n_max+1)
    beta=beta_min>beta<beta_max
    if (beta_min gt beta_max) or n_max gt n_max_range[1] or n_max gt n_max_range_soft[1] then begin
      chisq=!values.f_infinity
      decomp=shapelets_create_decomp(n_max)
      decomp.beta=beta
      decomp.chisq=replicate(chisq,2)
      decomp.x=centre_guess
    endif else begin
      decomp=shapelets_decomp(pstamp,beta,n_max,recomp=recomp,centre=centre_guess,$
  	         psf=psf,sky=sky,polar=polar,diamond=diamond)
      chisq=decomp.chisq[1]
    endelse
    if keyword_set(full_focus) then shapelets_update_focus, focus, decomp, silent=silent else $
      message,"  n_max, beta, chi^2, x_c:"+$
      string(fix(decomp.n_max))+string(decomp.beta)+string(decomp.chisq[1])+$
      string(decomp.x[0])+string(decomp.x[1]),/info,noprint=silent
  
    ; Is this still bad?
    if chisq gt old_chisq then begin
      n_max_opt=old_n_max
      decomp=old_decomp
      recomp=old_recomp
      if keyword_set(full_focus) then shapelets_update_focus, focus, decomp, silent=silent else $
    	message,"  n_max, beta, chi^2, x_c:"+$
    	string(fix(decomp.n_max))+string(decomp.beta)+string(decomp.chisq[1])+$
    	string(decomp.x[0])+string(decomp.x[1]),/info,noprint=silent
    endif else n_max_opt=n_max
  endelse

endif else if flat then begin   ; did we stop because of flatness? keep last value
  message,"Chi^2 flatness criterion reached",/info,noprint=silent
  focus.flag=focus.flag>3
  n_max_opt=n_max
endif else begin     ; no: then pick best value from 2/3 nearest neighbors
  if keyword_set(diamond) then begin
    chisq_three=[old_chisq,!values.f_infinity,chisq]
  endif else begin
    decomp_mid=shapelets_decomp(pstamp,beta,n_max-1,recomp=recomp_mid,centre=centre_guess,$
      psf=psf,sky=sky,polar=polar,diamond=diamond)
    if keyword_set(full_focus) then shapelets_update_focus, focus, decomp_mid, silent=silent else $
      message,"  n_max, beta, chi^2, x_c:"+$
      string(fix(decomp_mid.n_max))+string(decomp_mid.beta)+string(decomp_mid.chisq[1])+$
      string(decomp_mid.x[0])+string(decomp_mid.x[1]),/info,noprint=silent
    chisq_three=[old_chisq,decomp_mid.chisq[1],chisq]
  endelse
  n_max_three=[n_max-2,n_max-1,n_max]
  opt=where(chisq_three lt focus.chisq_abs_target,n_opt)
  if n_opt gt 0 then opt=opt[0] else dummy=min(chisq_three,opt)
  n_max_opt=n_max_three[opt]
  case opt of
    0: begin
         decomp=old_decomp
         recomp=old_recomp
         shapelets_update_focus, focus, decomp, silent=silent
       end
    1: begin
         decomp=decomp_mid
         recomp=recomp_mid
       end
    2: shapelets_update_focus, focus, decomp, silent=silent
    else: stop
  endcase
endelse

; Keep focus structre up to date
shapelets_update_focus, focus, decomp, silent=silent
if focus.n_max ne n_max_opt or decomp.n_max ne n_max_opt then begin
  message,"Error in decomp coefficient assignement",/info,noprint=silent
  focus.flag=focus.flag>10
endif

; Tell the world
return,decomp

end




