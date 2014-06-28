pro shapelets_read_sexcat, SEXCAT, FILENAME,    $
	                         FULL_PATH=full_path, $
                           ASCII=ascii,         $
                           COMMENT=comment,     $
                           UNIT=unit,           $
                           XBUGFIXED=xbugfixed, $
                           SILENT=silent

;$Id: shapelets_read_sexcat.pro, v2$
;
; Copyright © 2005 Richard Massey and Alexandre Refregier.
;
; This file is a part of the Shapelets analysis code.
; www.ast.cam.ac.uk/~rjm/shapelets/
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
;      SHAPELETS_READ_SEXCAT
;
; PURPOSE:
;      Reads in a SExtractor catalogue.
;
; INPUTS:
;      FILENAME  - Name of input image
;
; KEYWORD PARAMETERS:
;      SILENT    - Operates silently.
;      ASCII     - Assume catalogue is in (shapelets version 1&2) ASCII format.
;      FULL_PATH - Accepts the input filename as is (can be absolute or
;                  relative). Default behaviour is to prepend a path name
;                  with a relevant path obtained from from shapelets_paths.pro.
;      XBUGFIXED - SExtractor seems to have an intermittent bug that increases
;                  x coordinates by one in some images. Since this happens more
;                  often than not in my experience, the default behaviour of
;                  this routine is to automatically subtract 1. Setting this
;                  flag prevents that from happening. Note that SExtractor's
;                  pixel indexing system is already different from the IDL
;                  default, because it labels the bottom-left corner of pixels
;                  rather than their centre.
;
; OUTPUTS:
;      SEXCAT    - SExtractor catalogue structure.
;      UNIT      - Units of each field in the catalogue.
;      COMMENT   - Strings containing other pertinent information about each
;                  field in the catalogue.
;
; NOTES:
;      This routine assumes that any ascii format SExtractor catalogues were
;      created using the configuration file shex_ascii.param, which defines
;      each of the columns to be:
;      x
;      y
;      x_peak
;      y_peak
;      x_min
;      y_min
;      x_max
;      y_max
;      a     }
;      b     } ellipticity
;      theta }
;      class_star
;      fwhm
;      magnitude
;      flux
;      area
;      ID number
;      DO NOT CHANGE THAT CONFIURATION! (It is also hard-wired into the image 
;      simulation package). If you wish to use different SExtractor outputs,
;      use the (now default) fits format output, which can handle any columns.
;      In this case, the columns can be selected by editing config/shex.param.
;
; MODIFICATION HISTORY:
;      Nov 05 - XBUGFIXED keyword added for Barnaby Rowe.
;      May 05 - Stripping of file extensions from input implemented by RM.
;      Mar 05 - Recognition of different file formats improved by RM.
;      Mar 05 - Comments and units optionally returned with catalogue by RM.
;      Jan 05 - SExtractor's bug with ALL x coordinates compensated for by RM.
;      Aug 04 - File format changed from ascii to (fits) binary by RM.
;      Jun 04 - Default file extension changed from .sex to .sexcat by RM.
;      Mar 03 - Simlified by AR to only read catalogue and not run SExtractor.
;      Sep 03 - Actual execution of SExtractor moved out to sex.pro by RM.
;      May 03 - Ellipticity e,e1,e2 notation added to structure by RM.
;      Mar 03 - Filepaths generalised and routine renamed by RM.
;      Apr 02 - shapelets_read_sexcat.pro adapted from shex.pro by RM.
;      Dec 01 - First version of shex.pro written by Richard Massey.
;-

COMPILE_OPT idl2
ON_ERROR,2


;
;
; SET DEFAULTS
;
;
if not keyword_set(silent)    then silent=0B  ; Operate silently?
if not keyword_set(ascii)     then ascii=0B   ; Read in ascii version?
if not keyword_set(filename)  then message,"You must specify a file name!"
file=filename                                 ; Don't disrupt input variable
if not keyword_set(full_path) then file=shapelets_paths(2,silent=silent)+file
; Strip file extension
;dot=strpos(strlowcase(file),".")
;if dot ne -1 then file=strmid(file,0,dot)
   
       

;
;
; DETERMINE FORMAT OF SEXTRACTOR CATALOGUE, BASED ON ITS FILE EXTENSION
;
;
; Ensure backwards compatability
if ascii then begin
  if file_test(file+".sex") then begin
    extension=".sex"
  endif else if file_test(file+".sxt") then begin
    extension=".sxt"
  endif else if file_test(file+".sexcat") then begin
    extension=".sexcat"
    message,"File "+file+".sex (ASCII format) does not exist! Perhaps you want "+filename+".sexcat (FITS format). To obtain this, run sex.pro without the /ASCII option."
  endif else message,"File "+file+".sex does not exist!"
endif else begin
  if file_test(file+".sexcat") then begin
    ; Expected behaviour for shapelets versions 3 and above
    extension=".sexcat"
  endif else if file_test(file+".sex") then begin
    ; Ensure backwards compatability with shapelets v1 and v2, plus output
    ; from the image simulation package
    message,"WARNING: Catalogue only exists in ASCII format!",/info,noprint=silent
    extension=".sex"
    ascii=1B
  endif else if file_test(file+".sxt") then begin
    ; Deal with Jason's more prudish setups!
    extension=".sxt"
  endif else begin
    ; File doesn't seem to exist in any of our incarnations
    message,"File "+file+".sexcat does not exist!"
  endelse
endelse



;
;
; READ IN SExtractor CATALOGUE:
;
;
if ascii then begin
  
  ; Read from disc in ASCII format
  message,'Reading in ASCII format SExtractor catalogue from '+filename+extension,$
          /info,noprint=silent
  readcol,file+extension,x,y,xpeak,ypeak,xmin,ymin,xmax,ymax,a,b,theta,$
          class,fwhm,mag,flux,area,id,$
          format="F,F,I,I,I,I,I,I,F,F,F,F,F,F,F,I,L",/SILENT
  
  ; Store in structure
  sexcat={name:filename,           $
          type:"sexcat",           $
          n:n_elements(x),         $
          x:[[x],[y]],             $
          xpeak:[[xpeak],[ypeak]], $
          xmin:[[xmin],[ymin]],    $
          xmax:[[xmax],[ymax]],    $
          a:a,                     $
          b:b,                     $
          theta:theta,             $
          class:class,             $
          fwhm:fwhm,               $
          mag:mag,                 $
          flux:flux,               $
          area:area,               $
          id:id }

  ; Remember column units and relevant information
  comment={n:"Number of objects in catalogue",                      $
           x:"Object position [pixel]",                             $
           xpeak:"Coordinate of the brightest pixel [pixel]",       $
           xmin:"Minimum coordinate among detected pixels [pixel]", $
           xmax:"Maximum coordinate among detected pixels [pixel]", $
           id:"Running object number",                              $
           class:"S/G classifier output",                           $
           a:"Profile RMS along major axis [pixel]",                $
           b:"Profile RMS along minor axis [pixel]",                $
           theta:"Position angle (CCW/x) [deg]",                    $
           flux:"Best of FLUX_AUTO and FLUX_ISOCOR [count]",        $
           mag:"Best of MAG_AUTO and MAG_ISOCOR [mag]",             $
           fwhm:"FWHM assuming a gaussian core [pixel]",            $
           area:"Isophotal area above Analysis threshold [pixel**2]" }

  ; Remember column units and relevant information
  unit={n:"",          $
        x:"pixel",     $
        xpeak:"pixel", $
        xmin:"pixel",  $
        xmax:"pixel",  $
        id:"",         $
        class:"",      $
        a:"pixel",     $
        b:"pixel",     $
        theta:"deg",   $
        flux:"count",  $
        mag:"mag",     $
        fwhm:"pixel",  $
        area:"pixel**2" }
 
endif else begin

  ; Read from disc in FITS format
  message,'Reading in FITS format SExtractor catalogue from '+$
          filename+extension,/info,noprint=silent
  catalogue = mrdfits(file+extension, 2, header, /SILENT)

  ; Create a list of the units
  variable_name_list = strtrim( sxpar( header, "TTYPE*", COUNT=count, COMMENT=comment_list ),2)

  ; Initialise output structures
  sexcat={name:filename,type:"sexcat",n:n_elements(catalogue)}
  comment={n:"Number of objects in catalogue"}
  unit={n:""}

  ; Recover my own meta-data
  ;param_file = strtrim( sxpar( header, "PARAMFL", COUNT=n_paramfile, COMMENT=comment_list ),2)
  
  ; Fill up SExtractor catalogue structure, one field at a time.
  for i=0,count-1 do begin

    ; Obtain name of variable, plus its units and any pertinent comments
    variable_name=variable_name_list[i]
    variable_comment=strtrim(comment_list[i],2)
    variable_unit=sxpar(header, "TUNIT"+strtrim(i+1,2), count=unit_count)
    if unit_count ne 0 then begin
      variable_unit=strtrim(variable_unit[0],2)
      variable_comment=variable_comment + " [" + variable_unit + "]"
    endif else begin
      variable_unit="" 
    endelse

    ; Match 1D vectors of x and y variables (e.g. positions) into 2D arrays, and update comments accordingly
    yequiv=where(variable_name_list eq "Y"+strmid(variable_name,1),n_yequiv)
    if strupcase(strmid(variable_name,0,1)) eq "X" and n_yequiv gt 0 then begin
      variable=[[catalogue.(i)],[catalogue.(yequiv)]]
      position=strpos(variable_comment," along x") & if position ne -1 then variable_comment=strmid(variable_comment,0,position)+strmid(variable_comment,position+8)
      position=strpos(variable_comment,"x-coordinate ") & if position ne -1 then variable_comment=strmid(variable_comment,0,position)+"coordinate"+strmid(variable_comment,position+12)
      variable_comment=strupcase(strmid(variable_comment,0,1))+strmid(variable_comment,1)
    endif else variable=catalogue.(i)

    ; Rename a few variables for backwards compatibility
    case variable_name of
      "CLASS_STAR": variable_name="CLASS"
      "ISOAREA_IMAGE": variable_name="AREA"
      "NUMBER": variable_name="ID"
      "MAG_BEST": variable_name="MAG"
      "FLUX_BEST": variable_name="FLUX"
      "ALPHA_J2000": variable_name="RA"
      "DELTA_J2000": variable_name="DEC"
      else:
    endcase

    ; Strip label "_IMAGE" from all variable names. This is assumed by default!
    position=strpos(variable_name,"_IMAGE") & if position ne -1 then variable_name=strmid(variable_name,0,position)

    ; Insert variable (plus comments and units) into relevant output structures
    if strmid(variable_name,0,1) ne "Y" then begin
      sexcat=create_struct(sexcat,variable_name,variable)
      comment=create_struct(comment,variable_name,variable_comment)
      unit=create_struct(unit,variable_name,variable_unit)
    endif
  
  endfor

endelse


;
;
; CALCULATE ELLIPTICITY IN NICER FORM
;
;
if tag_exist(sexcat,"a") and tag_exist(sexcat,"b") then begin
  e=(sexcat.a^2-sexcat.b^2)/(sexcat.a^2+sexcat.b^2)
  sexcat=create_struct(sexcat,"e",e)
  unit=create_struct(unit,"e","")
  comment=create_struct(comment,"e", "Ellipticity from ratio of major/minor axes")
  if tag_exist(sexcat,"theta") then begin
    e1=e*cos(sexcat.theta*!pi/90)
    e2=e*sin(sexcat.theta*!pi/90)
    sexcat=create_struct(sexcat,"e1",e1,"e2",e2)
    unit=create_struct(unit,"e1","","e2","")
    comment=create_struct(comment,"e1","Ellipticity component 1 (weak lensing notation)",$
                                  "e2","Ellipticity component 2 (weak lensing notation)")
  endif
endif


;
;
; CONVERT SEXTRACTOR INDICES, STARTING AT 1, TO IDL INDICES, STARTING AT 0
;
;
if tag_exist(sexcat,"id") then sexcat.id=sexcat.id-1L
if tag_exist(sexcat,"x") then begin
  sexcat.x[*,0]=sexcat.x[*,0]-1.5 ; SExtractor definitely has a bug with x coords...
  sexcat.x[*,1]=sexcat.x[*,1]-0.5 ; ...or its conventions are utterly horrible!
  if keyword_set(xbugfixed) then sexcat.x[*,0]=sexcat.x[*,0]+1
endif 
if tag_exist(sexcat,"xpeak") then begin
  sexcat.xpeak[*,0]=fix(sexcat.xpeak[*,0]-2)
  sexcat.xpeak[*,1]=fix(sexcat.xpeak[*,1]-1)
  if keyword_set(xbugfixed) then sexcat.xpeak[*,0]=sexcat.xpeak[*,0]+1
endif 
if tag_exist(sexcat,"xmin") then begin
  sexcat.xmin[*,0]=fix(sexcat.xmin[*,0]-2)
  sexcat.xmin[*,1]=fix(sexcat.xmin[*,1]-1)
  if keyword_set(xbugfixed) then sexcat.x[*,0]=sexcat.xmin[*,0]+1
endif 
if tag_exist(sexcat,"xmax") then begin
  sexcat.xmax[*,0]=fix(sexcat.xmax[*,0]-2)
  sexcat.xmax[*,1]=fix(sexcat.xmax[*,1]-1)
  if keyword_set(xbugfixed) then sexcat.x[*,0]=sexcat.xmax[*,0]+1
endif 
;if tag_exist(sexcat,"area") then sexcat.area=fix(sexcat.area)


;
;
; TAKE A GUESS AT THE AVERAGE SEEING IN THE IMAGE
;
;
if tag_exist(sexcat,"class") and tag_exist(sexcat,"fwhm") then begin
  if ( max(sexcat.class) gt 0.8 ) then begin
    seeing=median([sexcat.fwhm[where(sexcat.class gt 0.8)]])
  endif else seeing=median([sexcat.fwhm])
  sexcat=create_struct(sexcat,"seeing",seeing)
  comment=create_struct(comment,"seeing","Guesstimate of seeing from size of stars [pixel]")
  unit=create_struct(unit,"seeing","pixel")
endif


;
;
; PRINT OUT A FRIENDLY MESSAGE THAT ALL HAS GONE WELL
;
;
message,'Found '+strtrim(sexcat.n,2)+' objects.',/info,noprint=silent


end
