;$Id: shapelets_create_focus.pro, v1$
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
; NAME:
;      SHAPELETS_UPDATE_FOCUS
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      
;
; INPUTS:
;      
;
; OPTIONAL INPUTS:
;      None.
;
; KEYWORD PARAMETERS:
;      None.
;
; OUTPUTS:
;      
;
; MODIFICATION HISTORY:
;      Jul 05 - Written by Richard Massey.
;-

pro shapelets_update_focus, focus, decomp, history=history, silent=silent

COMPILE_OPT idl2

; Parse input
if not keyword_set(focus) then focus=shapelets_create_focus()
;if not shapelets_structure_type(focus,message=message) then message,message
;if focus.type ne "focus" then message,"Not a focus structure!"
if not keyword_set(decomp) then message,"Nothing supplied with which to update the focus structure!"
;if not shapelets_structure_type(decomp,message=message) then message,message
;if decomp.type ne "decomp" then message,"Not a decomp structure!"

; Update latest values in focus structure
focus.beta=decomp.beta
focus.n_max=decomp.n_max
focus.x=decomp.x
focus.chisq=decomp.chisq[1]

if focus.n_iterations eq 0 then begin
  ; Add parameter guesses to new structures - inteded for use after first decomposition
  focus.x_guess=decomp.x
  focus.n_max_guess=decomp.n_max
  focus.beta_guess=decomp.beta
  ; Initialise the history of decomposition parameters
  focus.beta_history=decomp.beta
  focus.x_history=transpose(decomp.x)
  focus.n_max_history=decomp.n_max
  focus.chisq_history=decomp.chisq[1]
  if keyword_set(history) then focus.history=history
endif else begin
  ; Append to the history of decomposition parameters
  focus_new={name:focus.name,type:"focus"}
  names=tag_names(focus)
  for i=0,n_tags(focus)-1 do begin
    case strupcase(names[i]) of
      "NAME":
      "TYPE":
      "HISTORY":       if keyword_set(history) then focus_new=create_struct(focus_new,"history",[focus.history,history]) else focus_new=create_struct(focus_new,"history",focus.history)
      "BETA_HISTORY":  focus_new=create_struct(focus_new,"beta_history",[focus.beta_history,decomp.beta])
      "N_MAX_HISTORY": focus_new=create_struct(focus_new,"n_max_history",[focus.n_max_history,decomp.n_max])
      "X_HISTORY":     focus_new=create_struct(focus_new,"x_history",[focus.x_history,transpose(decomp.x)])
      "CHISQ_HISTORY": focus_new=create_struct(focus_new,"chisq_history",[focus.chisq_history,decomp.chisq[1]])
      else: focus_new=create_struct(focus_new,names[i],focus.(i))
    endcase
  endfor
  focus=focus_new
endelse
focus.n_iterations+=1

; Display update
message,"n_max, beta, chi^2, x_c:"+$
    string(fix(focus.n_max))+string(focus.beta)+string(focus.chisq)+$
    string(focus.x[0])+string(focus.x[1]),/info,noprint=silent

end

; ***********************************************************************
; ***********************************************************************
;
;+
; NAME:
;      SHAPELETS_CREATE_FOCUS
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Initiliase a band new decomp structure.
;      Works like an IDL __DEFINE procedure, but for an anonymous structure
;      (which we need because they will shrink or expand to accomodate
;      varying amounts of data).
;
; INPUTS:
;      N_MAX  - The desired truncation order of the shapelet expansion.
;
; OPTIONAL INPUTS:
;      None.
;
; KEYWORD PARAMETERS:
;      POLAR  - Whether or not the expansion should be in polar (set) or
;               Cartesian (unset) shapelet coefficients.
;
; OUTPUTS:
;      FOCUS - Returns the new focus structure.
;
; MODIFICATION HISTORY:
;      Jul 05 - Written by Richard Massey.
;-

function shapelets_create_focus, NAME=name,                      $
                                 COMMENT=comment,                $
                                 CHISQ_TARGET=chisq_target,      $
                                 CHISQ_TOLERANCE=chisq_tolerance,$
                                 CHISQ_FLATNESS=chisq_flatness,  $
                                 BETA_TOLERANCE=beta_tolerance,  $
                                 THETA_MIN_GEOM=theta_min_geom

COMPILE_OPT idl2

; Parse input
if not keyword_set(name)    then name=""
if not keyword_set(comment) then comment=""

; Initialise history
systime=systime()
date=strmid(systime,0,4)+strmid(systime,8,3)+strmid(systime,4,4)+strmid(systime,20,4)
time=strmid(systime,11,5)
history="Created at "+time+" on "+date+" "+comment+"."

; Set default numerical values
if n_elements(chisq_target) eq 0    then chisq_target_local=1.0    else chisq_target_local    = chisq_target      
if n_elements(chisq_tolerance) eq 0 then chisq_tolerance_local=1.0 else chisq_tolerance_local = chisq_tolerance
if n_elements(chisq_flatness) eq 0  then chisq_flatness_local=0.01 else chisq_flatness_local  = chisq_flatness  
if n_elements(beta_tolerance) eq 0  then beta_tolerance_local=1e-3 else beta_tolerance_local  = beta_tolerance  
if n_elements(theta_min_geom) eq 0  then theta_min_geom_local=0.2  else theta_min_geom_local  = theta_min_geom  

; List flag definitions
flag_interpret=["OK",							                         	                    $ ; 0
                "Bounced off geometrical constraints",                                           	                    $ ; 1
                "Entered the regime where the least-squares fitting matrix may be singular",     	                    $ ; 2
                "Converged by flatness limit - shapelets may incompletely represent object",     	                    $ ; 3
                "Not used",                                             $ ; 4
                "Chi^2 not monotonic in n_max - shapelets may incompletely represent object",                               $ ; 5
                "Amoeba dithered about and did not converge to target chi^2 - shapelets may incompletely represent object", $ ; 6
                "Focus iteration did not converge to target chisq - shapelets may incompletely represent object",           $ ; 7
                "Maximum n_max reached - shapelets may incompletely represent object",           	                    $ ; 8
                "Centroid wandered, pushing basis functions off the edge of the postage stamp",                             $ ; 9
                "Fatal crash during focus routine"]                                                                           ;10
flag_interpret_mini=["OK","Bounced","Singular matrix","Flatness","N/A","Worsening (!8n!X!imax!n)","Amoeba dithered","Focus didn't coverge","Max !8n!X!imax!n","Centroid wandered","Fatal crash"]

; Create actual structure
focus={name:name,            		           $
       type:"focus",         		           $
       history:history,      		           $
       beta:0.,              		           $
       n_max:0,              		           $
       x:fltarr(2),         		           $
       chisq:!values.f_infinity,                   $
       chisq_target:chisq_target_local,            $
       chisq_tolerance:chisq_tolerance_local,      $
       chisq_flatness:chisq_flatness_local,        $
       chisq_abs_target:float(chisq_target_local), $
       beta_tolerance:beta_tolerance_local,        $
       theta_min_geom:theta_min_geom_local,        $
       n_max_guess:0,           	           $
       beta_guess:0.,           	           $
       x_guess:fltarr(2),      		           $
       n_iterations:0,          	           $
       beta_history:0.,         	           $
       n_max_history:0,         	           $
       x_history:fltarr(1,2),  		           $
       chisq_history:0.,        	           $
       flag:0,                  	           $
       flag_interpret:flag_interpret,              $
       flag_interpret_mini:flag_interpret_mini}

; Return the shiny new focus structure
return,focus

end
