;$Id: shapelets_select_stars.pro, v1$
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
;      SHAPELETS_SELECT_STARS_PLOT
;
; CATEGORY:
;      Subroutine of shapelets_select_stars.pro
;
; PURPOSE:
;      Draws or refreshes the basic plot used to select objects.
;
; INPUTS:
;      SEXCAT
;
; OPTIONAL INPUTS:
;      TITLE
;      MAG_RANGE
;      MU_MAX_RANGE
;      HAPPY_POLY_X
;      HAPPY_POLY_Y
;      CANCEL_POLY_X
;      CANCEL_POLY_Y
;
; KEYWORD PARAMETERS:
;      None.
;
; OUTPUTS:
;      Draws a plot to STDOUT, which is hopefully the screen.
;
; MODIFICATION HISTORY:
;      May 05 - Written by Richard Massey.


pro shapelets_select_stars_plot, sexcat,$
                                 TITLE=title,$
                                 MAG_RANGE=mag_range,$
                                 MU_MAX_RANGE=mu_max_range,$
                                 HAPPY_POLY_X=happy_poly_x,$
                                 HAPPY_POLY_Y=happy_poly_y,$
                                 CANCEL_POLY_X=cancel_poly_x,$
                                 CANCEL_POLY_Y=cancel_poly_y

COMPILE_OPT idl2, HIDDEN

; Draw the points in a mag vs mu_max plane
plot,sexcat.mag,sexcat.mu_max,psym=3,$
     xran=mag_range,/xstyle,yran=mu_max_range,/ystyle,$
     title=title,xtitle="!6!8I!6 band magnitude",ytitle="!7l!6!dmaximum!n"

; Draw the happy box
if keyword_set(happy_poly_x) and keyword_set(happy_poly_y) then begin
  oplot,[happy_poly_x,happy_poly_x[0]],[happy_poly_y,happy_poly_y[0]],psym=-3
  xyouts,mean(happy_poly_x),mean(happy_poly_y),align=0.5,"!6Happy?  :-)"
endif

; Draw the cancel box
if keyword_set(cancel_poly_x) and keyword_set(cancel_poly_y) then begin
  oplot,[cancel_poly_x,cancel_poly_x[0]],[cancel_poly_y,cancel_poly_y[0]],psym=-3
  xyouts,mean(cancel_poly_x),mean(cancel_poly_y),align=0.5,"!6Clear  :-("
endif

end

;***************************************************************************
;***************************************************************************
;
;+
; NAME:
;      SHAPELETS_SELECT_STARS
;
; CATEGORY:
;      Shapelets.
;
; PURPOSE:
;      Interactively select various different object types.
;      Follow the instructions in the plot title and on the command line.
;      Right-click anywhere (or click in "cancel" box) to undo.
;
; INPUTS:
;      SEXCAT - A shapelets SExtractor catalogue structure.
;
; OPTIONAL INPUTS:
;      None.
;
; KEYWORD PARAMETERS:
;      None.
;
; OUTPUTS:
;      Write .iselect and _sat.reg files to disc
;
; OPTIONAL OUTPUTS:
;      ISELECT - Structure containing the selected parameters, a list
;                of saturated objects and a list of unsaturated stars.
;
; MODIFICATION HISTORY:
;      May 05 - Written by Richard Massey.
;-

pro shapelets_select_stars, sexcat, iselect, DS9REGION=ds9region

COMPILE_OPT idl2

; Plotting parameters
loadct,2,/SILENT
saturated_colour=100 ; Red
star_colour=30       ; Green
galaxy_colour=170    ; Purple
mag_range=[16<min(sexcat.mag<29,max=max_mag),27>max_mag]
mu_max_range=[19<min(sexcat.mu_max,max=max_mu_max),26>max_mu_max]
mag_range=[min(sexcat.mag<29,max=max_mag)-1,max_mag+1]
mu_max_range=[min(sexcat.mu_max,max=max_mu_max)-1,max_mu_max+1]
happy_poly_x=[mag_range[0]+1,mag_range[0]+2,mag_range[0]+2,mag_range[0]+1]
happy_poly_y=[mu_max_range[1]-2,mu_max_range[1]-2,mu_max_range[1]-1,mu_max_range[1]-1]
cancel_poly_x=happy_poly_x
cancel_poly_y=happy_poly_y-1.2

; Draw mag vs mu_max plane
device,get_screen_size=screen_size
window,0,title="Interactively selecting stars in "+sexcat.name,$
  xsize=fix(screen_size[0]), ysize=fix(screen_size[1]*0.75);,$
select_saturation:
shapelets_select_stars_plot, sexcat,                       $
  title="Select saturation point",                         $
  mag_range=mag_range, mu_max_range=mu_max_range,          $
  happy_poly_x=happy_poly_x, happy_poly_y=happy_poly_y,    $
  cancel_poly_x=cancel_poly_x, cancel_poly_y=cancel_poly_y

; Select saturation point
print,"Select saturation point"
happy=0B
cursor,junk,mu_max_sat,/down
saturated=where(sexcat.mu_max lt mu_max_sat,n_saturated)

; Refresh plot
outline_stellar_locus:
shapelets_select_stars_plot,sexcat,$
  title="Outline stellar locus",$
  mag_range=mag_range,mu_max_range=mu_max_range,$
  happy_poly_x=happy_poly_x,happy_poly_y=happy_poly_y,$
  cancel_poly_x=cancel_poly_x,cancel_poly_y=cancel_poly_y
oplot,mag_range,replicate(mu_max_sat,2),psym=-3
if n_saturated gt 0 then oplot,sexcat.mag[saturated],sexcat.mu_max[saturated],psym=3,symsize=2,color=saturated_colour

; Select first point in 
print,"Outline stellar locus"
cursor,mag_stellar,mu_max_stellar,/down
if inpoly(mag_stellar,mu_max_stellar,happy_poly_x,happy_poly_y) then goto,outline_stellar_locus
if !mouse.button ge 2 or inpoly(mag_stellar,mu_max_stellar,cancel_poly_x,cancel_poly_y) then goto,select_saturation
mouse_time=!mouse.time
happy=0B
n_clicks=0
while not happy do begin
  cursor,mag_stellar_next,mu_max_stellar_next,/down
  if n_clicks ge 2 and (!mouse.time-mouse_time lt 200 or inpoly(mag_stellar_next,mu_max_stellar_next,happy_poly_x,happy_poly_y)) then begin
    happy=1B
  endif else if !mouse.button ge 2 or inpoly(mag_stellar_next,mu_max_stellar_next,cancel_poly_x,cancel_poly_y) then begin
    goto,outline_stellar_locus
  endif else begin
    n_clicks+=1
    mag_stellar=[mag_stellar,mag_stellar_next]
    mu_max_stellar=[mu_max_stellar,mu_max_stellar_next]
    oplot,mag_stellar,mu_max_stellar,color=star_colour,psym=-3
  endelse
  mouse_time=!mouse.time
endwhile

; Select stars
stars=inside(sexcat.mag,sexcat.mu_max,mag_stellar,mu_max_stellar,/INDEX)
if stars[0] gt -1 then stars=stars[where(sexcat.mu_max[stars] ge mu_max_sat,n_stars)]
print,"Selected "+strtrim(string(n_stars),2)+" stars"

; Refresh plot
shapelets_select_stars_plot,sexcat,$
  title='Click "Happy?" box again to finish',$
  mag_range=mag_range,mu_max_range=mu_max_range,$
  happy_poly_x=happy_poly_x,happy_poly_y=happy_poly_y,$
  cancel_poly_x=cancel_poly_x,cancel_poly_y=cancel_poly_y
oplot,mag_range,replicate(mu_max_sat,2),psym=-3
if n_saturated gt 0 then oplot,sexcat.mag[saturated],sexcat.mu_max[saturated],psym=3,color=saturated_colour
oplot,[mag_stellar,mag_stellar[0]],[mu_max_stellar,mu_max_stellar[0]],psym=-3
if n_stars gt 0 then oplot,sexcat.mag[stars],sexcat.mu_max[stars],psym=3,color=star_colour

; Ask user to confirm that everything is OK
print,'Click "Happy?" box again to finish'
happy=0B
while not happy do begin
  cursor,mag_position,mu_max_position,/down
  if !mouse.button ge 2 or inpoly(mag_position,mu_max_position,cancel_poly_x,cancel_poly_y) then goto,outline_stellar_locus
  if inpoly(mag_position,mu_max_position,happy_poly_x,happy_poly_y) then happy=1B 
endwhile
wdelete,0

; Assume any objects left over to be galaxies
stars_and_saturated=([stars,saturated])[sort([stars,saturated])]
n_galaxies=sexcat.n-n_elements(stars_and_saturated)
if n_galaxies eq 0 then galaxies=-1 else begin
  for i=0,sexcat.n-1 do begin
    if min(abs(stars_and_saturated-i)) gt 0 then begin
      if n_elements(galaxies) eq 0 then galaxies=i else galaxies=[galaxies,i]
    endif
  endfor
endelse
print,"Left with "+strtrim(string(n_galaxies),2)+" galaxies"

; Save selection to disc
iselect={name:sexcat.name,        $
         type:"iselect",          $
         saturated:saturated,     $
         stars:stars,             $
         galaxies:galaxies,       $
         n_saturated:n_saturated, $
         n_stars:n_stars,         $
         n_galaxies:n_galaxies,   $
         mu_max_sat:mu_max_sat,   $
         mag_stellar:mag_stellar, $
         mu_max_stellar:mu_max_stellar}
save,iselect,filename=shapelets_paths(3,/SILENT)+(strsplit(sexcat.name,/EXTRACT))[0]+".iselect"

; Make a region (.reg) file to view saturated objects in DS9
if keyword_set(ds9region) then begin
  regfile=shapelets_paths(3,/SILENT)+(strsplit(sexcat.name,/EXTRACT))[0]+"_saturated.reg"
  openw,lun,regfile,/get_lun
  for i=0,n_stars-1 do begin
    printf,lun,"point "+$
      strtrim(sexcat.x[stars[i],0]+1.5,2)+" "+$
      strtrim(sexcat.x[stars[i],1]+1.0,2)+" "+$
             "# point=diamond color=red"
  endfor
  close,lun
endif

print,"Finished!"

end
