pro shapelets_read_sexcat, SEXCAT, FILENAME,  $
	                      FULLPATH=fullpath, $
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
;      FILENAME - Name of input image
;
; KEYWORD PARAMETERS:
;      SILENT   - Operates silently.
;      FULLPATH - Accepts the input filename as is (can be absolute or
;                 relative). Default behaviour is to prepend a path name
;                 from shapelets_paths.pro.
;
; OUTPUTS:
;      SEXCAT   - SExtractor catalogue structure
;
; NOTES:
;      This routine assumes that the SExtractor catalogue was created using
;      the configuration file shex.param, which defines each of the columns
;      to be:
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
;
; MODIFICATION HISTORY:
;      Jun 04 - Default file extension changed from .sex to .sexcat by RM.
;      Mar 03 - Simlified by AR to only read catalogue and not run SExtractor.
;      Sep 03 - Actual execution of SExtractor moved out to sex.pro by RM.
;      May 03 - Ellipticity e,e1,e2 notation added to structure by RM.
;      Mar 03 - Filepaths generalised and routine renamed by RM.
;      Apr 02 - shapelets_read_sexcat.pro adapted from shex.pro by RM.
;      Dec 01 - First version of shex.pro written by Richard Massey.
;-

COMPILE_OPT idl2, OBSOLETE

;
;
; SET DEFAULTS
;
;
if not keyword_set(silent)    then silent=0              ; Operate silently?
if not keyword_set(filename)  then message,"You must specify a file name!"
file=filename                                            ; Don't disrupt input variables
if not keyword_set(fullpath)  then file=shapelets_paths(3,silent=silent)+file
       

;
;
; ADD RELEVANT PATH NAME AND EXTENSION
;
;
if file_test(file+".sexcat") then begin
  ; Add standard file extension
  extension=".sexcat"
endif else if file_test(file+".sex") then begin
  ; Ensure backwards compatability
  extension=".sex"
endif else if file_test(file+".sxt") then begin
  ; Deal with more prudish setups!
  extension=".sxt"
endif else begin
  ; File doesn't seem to exist in any of our incarnations
  message,"File "+file+".sexcat does not exist!"
endelse


;
;
; READ IN SExtractor CATALOGUE:
;
;
message,'Reading in SExtractor catalogue from '+filename+extension,$
        /info,noprint=silent
readcol,file+extension,x,y,xpeak,ypeak,xmin,ymin,xmax,ymax,a,b,theta,$
        class,fwhm,mag,flux,area,id,$
        format="F,F,I,I,I,I,I,I,F,F,F,F,F,F,F,I,I",/SILENT
; Convert SExtractor indices, starting at 1, to IDL indices, starting at 0
x=x-0.5
y=y-0.5
xpeak=fix(xpeak-1)
ypeak=fix(ypeak-1)
xmin=fix(xmin-1)
xmax=fix(xmax-1)
ymin=fix(ymin-1)
ymax=fix(ymax-1)
area=fix(area)
id=id-1
; UGLY HACK! But SExtractor seems to make a mistake with x co-ordinates?
x=x-1


;
;
; MAKE A GUESS AT THE AVERAGE SEEING IN THE IMAGE
;
;
if ( max(class) gt 0.5 ) then seeing=median(fwhm[where(class gt 0.5)]) $
  else seeing=median(fwhm)


;
;
; CALCULATE ELLIPTICITY IN NICER FORM
;
;
e=(a^2-b^2)/(a^2+b^2)
e1=e*cos(theta*!pi/90)
e2=e*sin(theta*!pi/90)


;
;
; STORE IN STRUCTURE
;
;
sexcat={name:filename,   	   $
        type:"sexcat",   	   $
	   n:n_elements(x), 	   $
        seeing:seeing,   	   $
        x:[[x],[y]],     	   $
        xpeak:[[xpeak],[ypeak]], $
        xmin:[[xmin],[ymin]],    $
        xmax:[[xmax],[ymax]],    $
        id:id,           	   $
        class:class,     	   $
        a:a,             	   $
        b:b,             	   $
        theta:theta,     	   $
        e:e,             	   $
        e1:e1,           	   $
        e2:e2,           	   $
        flux:flux,       	   $
        mag:mag,         	   $
        fwhm:fwhm,       	   $
        area:area }


;
;
; PRINT OUT A FRIENDLY MESSAGE THAT ALL HAS GONE WELL
;
;
message,'Found '+strtrim(n_elements(x),1)+' objects.',/info,noprint=silent


end
