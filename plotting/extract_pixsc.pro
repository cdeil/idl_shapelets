function extract_pixsc,telescope


; Written October 2005 by Joel Berge
;
;

; Purpose
;  extract the value of PIXEL_SCALE from SExtractor configuration
;  files
;
; Inputs
;  none
;
; Optional input
;  telescope - telescope used for the image
;
; Output
;  value of PIXEL_SCALE



spawn,"gawk -F'PIXEL_SCALE' '{print  $2}' "+shapelets_paths(7)+telescope+".sex",aa

fcc=float(strmid(aa,0,8))

gd=where(fcc ne 0)
pix_sc=fcc(gd)
pix_sc=pix_sc(0)

return,pix_sc

end
