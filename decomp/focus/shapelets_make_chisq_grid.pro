;$Id: shapelets_make_chisq_grid.pro, v1$
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
;      SHAPELETS_MAKE_CHISQ_GRID
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Compute the chi^2 difference between the input and the reconstructed
;      image on a grid of values of beta and n_max.
;
; INPUTS:
;      PSTAMP - object structure from shapelets_extract_pstamp.pro
;
; OPTIONAL INPUTS:
;      PLOTIT: 1: plot chi^2 array; 2: also plot reconstructed images
;      N_MAX_RANGE: n_max range to investigate (default: educated guess)
;      BETA_RANGE: beta range to investigate (default: educated guess)
;      N_BETA: number of beta values to consider (default=10)
;      X0: centroid position (default: guess it from image moments)
;
; KEYWORD PARAMETERS:
;      SILENT: don't print results
;
; OUTPUTS:
;      CHISQ_GRID: result structure
;
; MODIFICATION HISTORY:
;      Sep 05 - POLAR and DIAMOND options added by RM.
;      Jul 05 - Updated by Richard Massey.
;      Apr 02 - Modified by AR to use least-square fitting rather than
;               linear decomposition
;      Aug 00 - Written by A. Refregier
;-

function shapelets_make_chisq_grid, PSTAMP,                  $
                                    N_MAX_RANGE=n_max_range, $
                                    BETA_RANGE=beta_range,   $
                                    FOCUS=focus,             $
                                    SHAPECAT=shapecat,       $
                                    N_BETA=n_beta,           $
                                    X0=x0,                   $
                                    PSF=psf,                 $
                                    SKY=sky,                 $
                                    PLOTIT=plotit,           $
                                    DIAMOND=diamond,         $
                                    POLAR=polar,             $
                                    SILENT=silent

COMPILE_OPT idl2

; Run automatic focus routine
; (This will always be useful, since it calculates a guess, and also properly 
; recentres the object)
message,"Running automatic focus",/info,noprint=silent
decomp_focus=shapelets_focus(PSTAMP,                  $
                             FOCUS=focus,             $
                             RECOMP=recomp_focus,     $
                             N_MAX_RANGE=n_max_range, $
                             PSF=psf,                 $
                             /FULL_FOCUS)


; Select range of n_max and beta to explore fully
if not keyword_set(n_max_range) then n_max_range=[0,n_guess*2<20]  ; range of n_max values
n_n_max=1+n_max_range[1]-n_max_range[0]
n_max=lindgen(n_n_max)+n_max_range[0]
if not keyword_set(beta_range) then beta_range=[beta_guess/2.,beta_guess*2.]
if not keyword_set(n_beta) then n_beta=10
case n_beta of
  1: beta=mean(float(beta_range))
  2: beta=float(beta_range)
  else: begin
          beta=xgen(beta_range[0],beta_range[1],np=n_beta)
          beta=beta[0:n_beta-1]+(beta[1]-beta[0])/2.
        end
endcase


; Set plotting options
if keyword_set(plotit) then begin
  device,get_screen_size=screen_size
  !p.multi=[0,n_n_max,n_beta,0,0]
  window,0,xsize=screen_size[0]*0.75,ysize=screen_size[1]*0.75,$
           title="Shapelet decompositions"
endif	   


; Decompose the image for different values of beta and max and compute chi^2
chisq_grid=fltarr(n_n_max,n_beta)
;gamma_grid=beta##replicate(1,n_n_max)
;n_max_gamma_grid=n_max#replicate(1,n_beta)
chisq_min=!values.f_infinity
for i=n_beta-1,0L,-1 do begin
  message,"Trying beta="+strtrim(beta[i],2),/info,noprint=silent
  for j=0L,n_n_max-1 do begin
    decomp=shapelets_decomp(pstamp.image,beta[i],n_max[j],recomp=recomp,centre=focus.x,/over,/ls,$
      noise=pstamp.noise,psf=psf,sky=sky,polar=polar,diamond=diamond);,gamma=gamma,n_max_gamma=n_max_gamma
;    if keyword_set(psf) then begin
;      gamma_grid[j,i]=gamma
;      n_max_gamma_grid[j,i]=n_max_gamma
;    endif
    chisq_grid[j,i]=decomp.chisq[1]
    if chisq_grid[j,i] lt chisq_min then begin
      chisq_min=chisq_grid[j,i]
      beta_best_fit=beta[i]
      n_max_best_fit=n_max[j]
    endif
    if keyword_set(plotit) then shapelets_plot,recomp
    if arg_present(shapecat) then begin
      if i eq (n_beta-1) and j le 1 then shapecat=shapelets_create_shapecat()
      shapelets_add,shapecat,decomp
    endif
  endfor
endfor


; Report best fit values
message,"best fit: beta, n_max: "+string(beta_best_fit)+string(n_max_best_fit),/info,noprint=silent
message,"theta_min,theta_max:"+string(beta_best_fit/sqrt(n_max_best_fit+1))+string(beta_best_fit*sqrt(n_max_best_fit+1)),/info,noprint=silent

theta_minmax_geom=shapelets_geometric_constraints(pstamp, psf=psf, centre_guess=focus.x)
theta_min_constraint=theta_minmax_geom[0]*sqrt(n_max+1)
theta_max_constraint=theta_minmax_geom[1]/sqrt(n_max+1)


; Store best fit values in structure
chisq_grid={name:pstamp.name,               	         $
            type:"chisq_grid",              	         $
            chisq:chisq_grid,               	         $
            n_max:n_max,                    	         $
            beta:beta,                      	         $
            x:focus.x,                      	         $
            chisq_target:focus.chisq_target,           $
            chisq_tolerance:focus.chisq_tolerance,     $
            chisq_flatness:focus.chisq_flatness,       $
            theta_min_geom:focus.theta_min_geom,       $
            theta_min_constraint:theta_min_constraint, $
            theta_max_constraint:theta_max_constraint, $
            beta_tolerance:focus.beta_tolerance,       $
            beta_guess:focus.beta_guess,               $
            nmax_guess:focus.n_max_guess,              $
            n_max_best_fit:n_max_best_fit,             $
            beta_best_fit:beta_best_fit}

; 
if arg_present(focus) then begin
  
endif

return,chisq_grid

end


