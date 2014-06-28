function shapelets_paths, n, SILENT=silent

;$Id: shapelets_paths.pro, v2$
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
;       SHAPELTS_PATHS
;
; PURPOSE:
;       This stores the locations of data on a locally accesible hard disc.
;       Please update this file so that all the strings point to the correct
;       locations!
;
; CATEGORY:
;       Disc management.
;
; CALLING PROCEDURE:
;       filepath=shapelets_paths(1,/SILENT)
;
; INPUTS:
;       n - integer scalar. Which path?
;
; KEYWORD PARAMETERS:
;       [/SILENT] - If this flag is set, the routine operates silently.
;
; OUTPUTS:
;       A string containing a filepath.
;
; EXAMPLE:
;       Find the location where data is stored on the local disc.
;
;       filepath=shapelets_paths(1) & print,filepath
;
;       In this case (1), it is the base directory where data such as images,
;       shapelet/SExtractor catalogues or PSF images are stored.
;
; NOTES:
;       I have created this function so that the shapelet routines can be
;       quickly set up to run on a new computer, and only this file needs to
;       be altered.
;
; PROCEDURES USED:
;       Often calls itself recursively.
;
; MODIFICATION HISTORY:
;       Apr 05 - Compatability with other OSs improved by RM.
;       Jul 03 - Finalised for version 1 of shapelet code release by RM.
;       Jan 01 - Written by Richard Massey
;-

message,"Shapelet export version of path names!",/info,noprint=silent

case n of

    0: begin  ; Directory containing the unpacked shapelets IDL code
         case !version.os_family of
         "Windows": out="C:\Documents and Settings\Rich Massey\My Documents\IDL routines\shapelets\"
              else: out=filepath("",root_dir=expand_path("~"),subdirectory=["idl","shapelets"])
         endcase
       end

    1: $ ; Generic base data directory
       out=filepath("",root_dir=expand_path("~"),subdirectory="data")

    2: $ ; Directory containing object catalogues (*.sex and *.shapelets)
       out=shapelets_paths(1,/silent)

    3: $ ; Directory containing images (*.fits)
       out=shapelets_paths(1,/silent)

    4: $ ; Directory containing PSF images
       out=filepath("",root_dir=shapelets_paths(1,/silent),subdirectory=["PSF"])

    5: $ ; Directory in which to put diagnostic plots and other files
       out=filepath("",root_dir=shapelets_paths(1,/silent),subdirectory=["diagnostics"])
    
    7: $ ; SExtractor directory
       out=filepath("",root_dir=shapelets_paths(0,/silent),subdirectory=["pipeline","config"])

    8: $ ; SExtractor executable
       out=filepath("sex",root_dir="/",subdirectory=["usr","local","optical","sextractor","bin"])
	  
 else: message,"Please set up shapelets_paths.pro to include the location of this file!"

endcase

return,out

end
