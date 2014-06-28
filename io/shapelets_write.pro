;$Id: shapelets_write.pro, v2$
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
; ************************************************************************
; ************************************************************************
;+
; NAME:
;      SHAPELETS_WRITE
;
; PURPOSE:
;      Writes a shapelet coefficient catalogue of objects to disc.
;
; CATEGORY:
;      Data I/O.
;
; INPUTS:
;      STRUCTURE   - IDL shapecat or decomp structure to be stored.
;      FILENAME    - Desired name of saved file, without path or extension.
;
; OPTIONAL INPUTS
;      DESCRIPTION - String containing pertinent information.
;
; KEYWORD PARAMTERS:
;      SILENT      - Operates silently.
;      FULL_PATH   - Accepts the input filename as is (can be absolute or
;                    relative). Default behaviour is to prepend a path name
;                    with a relevant path obtained from from shapelets_paths.pro.
;
; OUTPUTS:
;      Shape information is written to disc.
;
; MODIFICATION HISTORY:
;      Jan 06 - Various image and FITS header options introduced by RM.
;      Jul 05 - Possibility of using DESCRIPTION keyword added by RM. 
;      Jul 05 - Default "IDL" filename adopted by RM.
;      Jul 05 - Renamed to shapelets_write by RM.
;      Jun 05 - Made able to cope with decomp structures by RM.
;      May 05 - FULL_PATH option added by RM.
;      May 04 - Shapelets_write_shapecat written by Richard Massey.
;-

pro shapelets_write, structure, filename, 	    $
                     FULL_PATH=full_path, 	    $
                     SILENT=silent,       	    $
                     IMAGE=image,               $
                     MASK=mask,                 $
                     SEGMENTATION=segmentation, $
                     NOISE=noise,               $
                     _ref_extra = ex

COMPILE_OPT idl2

if not keyword_set(filename) then filename="idl"
if not shapelets_structure_type(structure,message=message) then message,message
filepath=""
case strupcase(message) of
  "SHAPECAT": begin
                if not keyword_set(full_path) then filepath=shapelets_paths(2,silent=silent)
                if structure.polar then message,"Maybe should convert catalogue to Cartesian shapelets!",/info
                message,"Writing shapelet catalogue to "+filepath+filename+"."+structure.type,$
                        /info,noprint=silent
                save,structure,filename=filepath+filename+"."+structure.type,_extra= ex
	      end
  "DECOMP":   begin
                if not keyword_set(full_path) then filepath=shapelets_paths(2,silent=silent)
                message,"Writing decomp structure to "+filepath+filename+"."+structure.type+". Restore using RESTORE.",$
                       /info,noprint=silent
                save,structure,filename=filepath+filename+"."+structure.type,_extra= ex
	      end
  "IMAGE":    begin
                if keyword_set(full_path) then begin
                  fileext=".fits"
                endif else begin
                  filepath=shapelets_paths(3,silent=silent)
                  fileext=".fits"
                endelse
                n_written=0
                if keyword_set(segmentation) then begin
                  extension="_seg"
                  header=structure.header
                  sxaddpar,header,"DATAMAX",max(structure.seg,min=datamin)
                  sxaddpar,header,"DATAMIN",datamin
                  sxaddpar,header_noise,"OBJECT","Segmentaion map"
                  sxaddpar,header_noise,"BUNIT",0B
                  message,"Writing segmentation image to "+filepath+filename+extension+fileext,/info,noprint=silent
                  file_delete,filepath+filename+extension+fileext,/QUIET
                  writefits,filepath+filename+extension+fileext,structure.seg,header
                  n_written=n_written+1
                endif
                if keyword_set(noise) then begin
                  extension="_weight"
                  header=structure.header
                  sxaddpar,header,"DATAMAX",max(structure.noise,min=datamin)
                  sxaddpar,header,"DATAMIN",datamin
                  sxaddpar,header_noise,"OBJECT","Inverse variance map"
                  sxaddpar,header_noise,"BUNIT","counts^(-2)sec^2"
                  message,"Writing inverse variance weight map to "+filepath+filename+extension+fileext,/info,noprint=silent
                  file_delete,filepath+filename+extension+fileext,/QUIET
                  writefits,filepath+filename+extension+fileext,structure.noise,header
                  n_written=n_written+1
	              endif
                if keyword_set(mask) then begin
                  extension="_mask"
                  header=structure.header
                  sxaddpar,header,"DATAMAX",max(structure.mask,min=datamin)
                  sxaddpar,header,"DATAMIN",datamin
                  sxaddpar,header_noise,"OBJECT","Image mask"
                  sxaddpar,header_noise,"BUNIT",0B
                  message,"Writing image mask to "+filepath+filename+extension+fileext,/info,noprint=silent
                  file_delete,filepath+filename+extension+fileext,/QUIET
                  writefits,filepath+filename+extension+fileext,structure.mask,header
                  n_written=n_written+1
                endif
                if keyword_set(image) or n_written eq 0 then begin
                  message,"Writing image to "+filepath+filename+fileext,/info,noprint=silent
                  file_delete,filepath+filename+fileext,/QUIET
                  writefits,filepath+filename+fileext,structure.image,structure.header
                endif
	      end
  else: message,"Cannot save "+structure.type+"s with this routine!"
endcase

end
