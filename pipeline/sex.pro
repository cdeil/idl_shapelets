pro sex, filename, SEXCAT=sexcat,                     $
                   SILENT=silent,                     $
                   FULL_PATH=full_path,               $
                   ASCII=ascii,                       $
                   TELESCOPE=telescope,               $
                   SEEING=seeing,                     $
                   FILTER=filter,                     $
                   BACKGROUND_LEVEL=background_level, $
                   RMS_BACKGROUND=rms_background,     $
                   SEGMENTATION=segmentation
           
;$Id: sex.pro, v2$
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
;       SEX
;
; PURPOSE:
;       Runs SExtractor (by spawning it as a child process under the shell
;       used to run IDL) to locate objects within an image.
;
; CATEGORY:
;       ROI finder.
;
; CALLING PROCEDURE:
;       sex,"example"
;
; INPUTS:
;       FILENAME - Name of input image
;
; OPTIONAL INPUTS:
;       TELESCOPE - Gets SExtractor parameters right. Default is HDF.
;       FILTER    - SExtractor convolution filter. Default is "gauss_2.0_5x5"
;       SEEING    - Used for star/galaxy separation.
;
; KEYWORD PARAMETERS:
;       /SILENT   - Operates silently.
;       /ASCII    - Writes out an ASCII catalogue, rather than a binary (fits)
;                   table. There is a similar switch in
;                   shapelets_read_sexcat.pro to read in the two file types.
;       /SEGMENTA - Also output a segmentation map.
;       /BACKGROU - Also ouptut an image containing the local background level.
;       /RMS_BACK - Also ouptut an image containing the local noise rms.
;       /FULL_PAT - Assume the input filename is absolute, rather than relative
;                   to some location specified in shapelets_paths.pro.
;
; OUTPUTS:
;       SExtractor is run and outputs various files to disc.
;
; OPTIONAL OUTPUTS:
;       SEXCAT    - SExtractor catalogue data structure, automatically read
;                   back in from disc.
;
; TO DO:
;       Try a two-stage version of SExtractor. Run once to get rough catalogue
;       of stars. Then run again with a new seeing setting. Bertin's rule of
;       thumb for setting the seeing in default.sex file is SEEING=FLUX_RADIUS*2
;       This should improve star/galaxy separation in SExtractor, at a minimal
;       cost incomputing time.
;
;       Add information (to fits headers? to a new file format?) storing the
;       parameter file, the date, the image name, etc. for future reference.
;
; PROCEDURES USED:
;       None. But requires SExtractor to be installed on the computer, with
;       the full path of the "sex" executable stored in shapelets_paths.pro.
;
; MODIFICATION HISTORY:
;       Dec 05 - Bug in the full_path option fixed by RM
;       Oct 05 - Spurious UNIX commands removed by RM and replaced with IDL code
;       Oct 05 - Catalogue written to intended folder by Joel Berge
;       Apr 05 - RMS keyword added by Will high and RM
;       Apr 05 - SEXCAT keyword added by RM
;       Mar 05 - SEEING keyword added by RM
;       Mar 05 - Can change columns in SExtractor fits (but not ASCII) output by RM
;       Jan 05 - Appropriate message displayed after errors spawning sex by RM
;       Jan 05 - Option to write out ASCII/FITS format catalogues added by RM
;       Jul 04 - The segmentation file is now called _seg.fits by Alexandre Refregier
;       May 04 - Convention that weight maps stored as _weight.fits adopted by RM
;       Sep 03 - sex.pro adapted from shapelets_read_sexcat.pro by RM
;       May 03 - Ellipticity e,e1,e2 notation added to structure by RM
;       Mar 03 - Filepaths generalised and routine renamed shapelets_read_sexcat by RM
;       Apr 02 - rd_scat.pro adapted from shapelets_extractor.pro by RM
;       Dec 01 - First version of shapelets_extractor.pro written by Richard Massey
;-

COMPILE_OPT idl2
ON_ERROR,2


; Set defaults
if not keyword_set(silent)    then silent=0              ; Operate silently?
if not keyword_set(telescope) then telescope="default"	 ; Gets SEx params right
if not keyword_set(filter)    then filter="gauss_2.0_5x5"; SExtractor detection filter
if not keyword_set(filename)  then message,"You must specify a file name!"
if keyword_set(full_path) then begin
  ; Strip file extension
  input_file=strmid(filename,0,strpos(strlowcase(filename),".fit"))
  output_file=strmid(filename,0,strpos(strlowcase(filename),".fit"))
endif else begin
  ; Add relevant path name
  input_file=shapelets_paths(3,silent=silent)+filename
  output_file=shapelets_paths(2,silent=silent)+filename
endelse


; Check to see if the image file exists
if ( not file_test(input_file+".fits") ) then $
	message,"Image "+input_file+".fits not available!"


; Print intent to screen
message,"Locating objects with SExtractor",/info,noprint=silent


; Decide format and file extension for backwards compatability
if keyword_set(ASCII) then begin
  format="ASCII_HEAD"
  extension="sex"
  paramfile="shex_ascii.param"
endif else begin
  format="FITS_LDAC"  
  extension="sexcat"
  paramfile="shex.param"
endelse


; Form shell command
sex_dir=shapelets_paths(7,/silent)
command_line=shapelets_paths(8,/silent)+" "+input_file+".fits"     +$
  " -c "              +sex_dir+telescope+".sex"              +$
  " -STARNNW_NAME "   +sex_dir+"default.nnw"                 +$
  " -FILTER_NAME "    +sex_dir+filter+".conv"                +$
  " -PARAMETERS_NAME "+sex_dir+paramfile                     +$
  " -CATALOG_TYPE "   +format                                +$
  " -CATALOG_NAME "   +output_file+"."+extension
if keyword_set(segmentation) or keyword_set(background_level) or keyword_set(rms_background) then begin
  checkimage_type=" -CHECKIMAGE_TYPE "
  checkimage_name=" -CHECKIMAGE_NAME "
endif else begin
 checkimage_type=" "
 checkimage_name=" "
endelse
if keyword_set(segmentation) then begin
  checkimage_type=checkimage_type+'"SEGMENTATION",'
  checkimage_name=checkimage_name+'"'+output_file+'_seg.fits",'
endif
if keyword_set(background_level) then begin
  checkimage_type=checkimage_type+'"-OBJECTS",'
  checkimage_name=checkimage_name+'"'+output_file+'_back.fits",'
endif
if keyword_set(rms_background) then begin
  checkimage_type=checkimage_type+'"BACKGROUND_RMS",'
  checkimage_name=checkimage_name+'"'+output_file+'_rms.fits",'
endif
checkimage_type=strmid(checkimage_type,0,strlen(checkimage_type)-1)
checkimage_name=strmid(checkimage_name,0,strlen(checkimage_name)-1)
command_line=command_line+checkimage_type+" "+checkimage_name
if keyword_set(silent) then $
  command_line=command_line+" -VERBOSE_TYPE QUIET "
if ( file_test(input_file+"_weight.fits") ) then $
  command_line=command_line             +$
  " -WEIGHT_TYPE MAP_WEIGHT "           +$
  " -WEIGHT_IMAGE "+input_file+"_weight.fits"
if keyword_set(seeing) then begin
    message,'using input SEEING ='+strtrim(seeing),/continue
    command_line=command_line        +$
  " -SEEING_FWHM "+strtrim(string(seeing),2)+" "
endif


; Run SExtractor
spawn,command_line,exit_status=exit_status


;; Add some useful information to the header
;if not keyword_set(ASCII) then begin
;  catalogue = mrdfits(file+"."+extension, 2, header, /SILENT)  
;  paramfile='shex.param'
;  param_file_string="PARAMFL = '"+paramfile+"'"
;  param_file_string=param_file_string+strmid(header[n_elements(header)-1],strlen(param_file_string))
;  insert_here=where(strupcase(strmid(header,0,7)) eq "EXTNAME",insert)
;  if insert gt 0 then header=[header[0:insert_here-1],param_file_string,header[insert_here:*]]
;  ; This next line doesn't write it out in the correct format!
;  mwrfits, catalogue, file+"b."+extension, header
;endif


; Report success or failure
if exit_status eq 0 then begin
  message,"SExtractor catalogue written to "+output_file+"."+extension,/info,noprint=silent
endif else begin
  message,"SExtractor failed to run"
endelse


; Read SExtractor catalogue into memory if desired
if arg_present(sexcat) then shapelets_read_sexcat,sexcat,filename,     $
                                                  full_path=full_path, $
                                                  silent=silent


end
