pro shapelets_plot_chisq_grid, chisq_grid, focus, cran=cran, clog=clog,psf=psf

;$Id: shapelets_plot_chisq_grid.pro, v2$
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
; August 2000 - Written by A. Refregier
; PURPOSE: plot Chi^2 on the beta-nmax grid
; INPUT: chisq_grid: Chi^2 structure from shapelets_make_chisq_grid.pro
; OUTPUT: plot
;-

; Compute beta range
n_beta=n_elements(chisq_grid.beta)
beta_min=chisq_grid.beta(0)
beta_max=chisq_grid.beta(n_beta-1)
delta_beta=chisq_grid.beta(1)-chisq_grid.beta(0)

; Plot chi^2 vs beta and n_max
if keyword_set(clog) then title="!8log!6 (!7v!6!e2!n)" else title="!7v!6!e2!n"
shapelets_plot_image,chisq_grid.chisq,/col,$
  frame=[min(chisq_grid.n_max)-.5,max(chisq_grid.n_max)+1-.5,beta_min-delta_beta/2,beta_max-delta_beta/2],$
  xtitle="!8n!6!imax!n",ytitle="!7b!6",title=title, $
  cran=cran, clog=clog
;oplot,[chisq_grid.nmax_g],[chisq_grid.beta_g],psym=6
;oplot,[chisq_grid.nmax_bf],[chisq_grid.beta_bf],psym=1

; Overlay contours
contour,chisq_grid.chisq,chisq_grid.n_max,chisq_grid.beta,$
    level=chisq_grid.chisq_target+chisq_grid.chisq_tolerance*[-2.,-1.,0.,1.,2.],$
    /over,c_lines=[1,2,0,2,1]

; Plot geometrical constraints
print,chisq_grid.theta_min_constraint
oplot,chisq_grid.n_max,chisq_grid.theta_min_constraint,psym=-3,linestyle=3
oplot,chisq_grid.n_max,chisq_grid.theta_max_constraint,psym=-3,linestyle=3
;oplot,chisq_grid.n_max,chisq_grid.theta_min_geom/sqrt(chisq_grid.n_max+1),psym=-3,linestyle=3

; Plot singular matrix constraint
if keyword_set(psf) then begin
  N_alpha=chisq_grid.n_max+1
  N_beta=psf.n_max+1.
  N_gamma=N_alpha-0.5
  singular_alpha=psf.beta*sqrt((  N_alpha*(N_beta^2-N_gamma^2)/$
                                 (N_beta*(N_gamma^2-N_alpha^2))    )>0)
  print,singular_alpha
  oplot,chisq_grid.n_max,singular_alpha,psym=-3
endif

if keyword_set(focus) then begin
  if focus.n_iterations gt 1 then begin
    oplot,focus.n_max_history,focus.beta_history,psym=-3
    oplot,focus.n_max_history,focus.beta_history,psym=1
  endif
  usersym, cos(2*!pi*findgen(21)/20), sin(2*!pi*findgen(21)/20)
  oplot,[focus.n_max_guess<max(chisq_grid.n_max)],[focus.beta_guess],psym=8,symsize=3
  oplot,[focus.n_max],[focus.beta],psym=7,symsize=3
  xyouts,.95*max(chisq_grid.n_max),.95*max(chisq_grid.beta),/data,charsize=charsize,align=1,$
       "!7b!3="+string(focus.beta,format="(f6.2)")+", !8n!3!dmax!n="+strtrim(string(focus.n_max),2)+$
       "!C!7v!3!e2!n="+strtrim(string(focus.chisq),2)+$;", snr="+strtrim(string(snr),2)+$
       "!Cfocus: "+focus.flag_interpret_mini[focus.flag];+" ("+strtrim(focus.flag,2)+")"

endif

end


