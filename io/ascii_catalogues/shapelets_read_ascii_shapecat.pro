;$Id: shapelets_read_ascii_shapecat.pro, v2$
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
;
; NAME:
;      SHAPELETS_READ_ASCII_SHAPEHEAD
;
; PURPOSE:
;      Reads just the header information in shapelet catalogues.
;
; CATEGORY:
;      A component of shapelets_read_shapecat.
;
; INPUTS:
;      filename - complete filename to be read, including path and extension.
;
; OUTPUTS:
;      Returns a structure containing some details about the catalogue.
;
; EXAMPLE CATALOGUE HEADER FORMATTING:
;# Shapelet coefficients of objects in example_image.fits:
;#   num objs       out of      n_max       n_coeffs       seeing  n_pixels(x)  n_pixels(y)  pixel scale
;#         28           29          5             21     5.310000          136          136     0.000000
;#
;#         x0           y0         beta     SEx fwhm    magnitude        n_max         Flag       SEx ID  Chi squared    SEx class        SEx A        SEx B    SEx theta       SEx e1       SEx e2     SEx flux        SEx x        SEx y     SEx area Object:x_min        x_max        y_min        y_max
;# n1       0            0            1            0            1            2            0            1            2            3            0            1            2            3            4            0            1            2            3            4            5
;# n2       0            1            0            2            1            0            3            2            1            0            4            3            2            1            0            5            4            3            2            1            0
;
; MODIFICATION HISTORY:
;      Oct 03 - Extra information added to headers by RM
;      Aug 01 - Written by Richard Massey

function shapelets_read_ascii_shapehead,filename

COMPILE_OPT idl2, OBSOLETE, HIDDEN

; Initialise local variables
junk=' '
n=0L
n_max=0
n_coeffs=0
n_total=0L
seeing=0.
textheader=strarr(1)

; Open shapelet catalogue
openr,lun_h,filename,/get_lun

; Read in file name
readf,lun_h,textheader
namestart=strpos(textheader,"in")+3
nameend  =strpos(textheader,".fit")
name=strmid(textheader,namestart,nameend-namestart)

; Skip second line ("#   num objs       out of     n_coeffs     n_coeffs       seeing  n_pixels(x)  n_pixels(y)  pixel scale")
readf,lun_h,textheader

; Read in global data
readf,lun_h,format='(A1,I11,3(I13),F13.6,2(I13),F13.6)',$
      junk,n,n_total,n_max,n_coeffs,seeing,n_pix_x,n_pix_y,pixel_scale

; Skip next two lines("#" and "#         x0           y0         beta     SEx fwhm      SEx mag        n_max         Flag       SEx ID  Chi squared    SEx class        SEx A        SEx B    SEx theta       SEx e1       SEx e2     SEx flux        SEx x        SEx y     SEx area Object:x_min        x_max        y_min        y_max")
readf,lun_h,textheader
readf,lun_h,textheader

; Read in n_1 values for each coefficient
thirdline=" "
readf,lun_h,thirdline
thirdline="    "+strmid(thirdline,4)+" "
n1=intarr(n_coeffs)
for i=0,n_coeffs-1 do n1[i]=fix(strmid(thirdline,(13*i),(13*(i+1))-1))

; Read in n_2 values for each coefficient
fourthline=" "
readf,lun_h,fourthline
fourthline="    "+strmid(fourthline,4)+" "
n2=intarr(n_coeffs)
for i=0,n_coeffs-1 do n2[i]=fix(strmid(fourthline,(13*i),(13*(i+1))-1))

; Close shapelet catalogue
close,lun_h
free_lun,lun_h

; Store information in an IDL structure and return
return,{name:name[0],            $
        n:n,                     $
        n_max:n_max,             $
        n_coeffs:n_coeffs,       $
        n1:n1,                   $
        n2:n2,                   $
        n_total:n_total,         $
        seeing:seeing,           $
        n_pix:[n_pix_x,n_pix_y], $
        pixel_scale:pixel_scale }

end



; ************************************************************************
; ************************************************************************
;
;+
; NAME:
;      SHAPELETS_READ_ASCII_SHAPECAT
;
; PURPOSE:
;      Reads in a shapelet coefficient catalogue of objects, created by
;      shapelets_extractor.pro, then stores them in a shapecat IDL structure. This may
;      then be converted into decomp structures using 
;      shapelets_shapecat2decomp.pro.
;
; CATEGORY:
;      Data I/O.
;
; INPUTS:
;      cat      - name of structure where data will be stored.
;      filename - filename to be read, without path or extension.
;
; OPTIONAL INPUTS
;      obj      - read in only object #obj. Cat will now have only 1 
;                 dimension.
;      n_max    - truncates all objects to only n_max coefficients.
;
; KEYWORD PARAMTERS:
;      /POLAR   - converts all objects to polar shapelet coefs on load.
;                 Convention is to have f_00 but then f_ij/f_00 for the rest.
;      /MOMENTS - calculates shapelet unweighted quadrupole moments on load.
;                 This also gives flux and size.
;      /PARITY  - reorientates to horizontal (based on phase{a_22}) and
;                 flips parity of those with right(?)-handedness (based on 
;                 phase{a_42}-phase{a_22}) on load.
;                 NOT YET FULLY TESTED AND PROBABLY NOT WORKING PROPERLY!
;      /RESIZE  - rescales objects on load to have same beta
;                 NOT REALLY VERY USEFUL AND BLOODY SLOW!
;                 a quicker alternative is suggested in a comment below
;      /V1      - Load shapelet coefficients as a_ij/a_00, as in version 1
;                 of shapelet code. Default is to store just a_ij.
;
; OUTPUTS:
;      cat      - Shape information is stored in this structure.
;      header   - Optionally returns all information from the catalogue's
;                 header into this variable (see shapelets_read_shapehead
;                 routine above).
;
; MEANING OF FIRST FLAG (from shapelets_image2obj.pro):
;      0: OK
;      1: did not normalise noise locally (is object filling postage stamp?)
;      2: nearby object (beware overlapping isophotes!)
;      3: both of above
;      4: near edge or a saturated star
;      5: HUGE object (SExtractor gaves whole image occasionally)
;      NB can only be one of these. If duplicated, a higher number takes priority.
;
; MEANINGS OF SECOND FLAG (from shapelets_focus.pro):
;      0: OK
;      1: has bounced into theta_min or max wall at some point
;      2: nmax_max reached - shapelets incompletely represent object
;      3: chi2>1 but converged (in n_max search) by flatness limit - still incomplete
;      4: did not converge
;      5: centroid off edge
;      6: fatal error
;      NB it can only be one of these. If duplicated, highest number has priority.
;
; MODIFICATION HISTORY:
;      Apr 04 - Routine obfuscated as shapelet catalogues now stored on disc in IDL save format
;      Mar 03 - RM added inclusion of ellipticity/quadrupole moments calculation
;      Apr 02 - RM speeded up conversion to polars by looking only at used coefficients
;      Apr 02 - RM added inclusion of size/magnitude calculations within the structure
;      Feb 02 - RM added conversion to polar shapelets
;      Dec 01 - RM added option to select just one object
;      Aug 01 - Written by R.Masey
;-


pro shapelets_read_ascii_shapecat, cat , filename, obj=obj, header_info=header_info, $
                  n_max=n_max, polar=polar, moments=moments, $
                  parity=parity, resize=resize, silent=silent, v1=v1

COMPILE_OPT OBSOLETE

; Set defaults
if not keyword_set(polar) then polar=0
if keyword_set(resize) then moments=1
if n_params() lt 2 then filename='HDF-N_f814'

; Read in header
if n_elements(header_info) eq 0 then header_info=shapelets_read_ascii_shapehead(shapelets_paths(2,/silent)+filename+'.shape')
n_coeffs = header_info.n_coeffs
n   = header_info.n
message,'Reading in catalogue from '+shapelets_paths(2,silent=silent)+filename+'.shape',/info,noprint=silent

; Open catalogue file for reading
openr,lun_h,shapelets_paths(2,/silent)+filename+'.shape',/get_lun

; Skip header
junk=strarr(1) & for i=1,6 do readf,lun_h,junk


if ( keyword_set(obj) ) then begin

  ;
  ; EXTRACT ONLY ONE OBJECT FROM CATALOGUE. SKIP THROUGH ALL THE REST.
  ;
  
  n=1
  for i=0,3*(obj-1) do readf,lun_h,junk

  ; See general extraction below for comments
  object_data  = strarr(1)
  coefficients = strarr(1)
  errors       = strarr(1)
  readf,lun_h,object_data
  readf,lun_h,coefficients
  readf,lun_h,errors

  x0      = float([strmid(object_data,0,12),strmid(object_data,13,12)])
  beta    = float(strmid(object_data,26,12)) & beta=beta(0)
  fwhm    = float(strmid(object_data,39,12)) & fwhm=fwhm(0)
  mag     = float(strmid(object_data,52,12)) &  mag=mag(0)
  nmax    =   fix(strmid(object_data,65,12)) & nmax=nmax(0)
  if keyword_set(n_max) then nmax=nmax<n_max
  a_used  =(nmax+1)*(nmax+2)/2
  flag    = fix(strmid(object_data,78,12))   & flag=[flag(0)/10,flag(0) mod 10]
  sexid   = fix(strmid(object_data,91,12))   & sexid=sexid(0)
  chisq   = float(strmid(object_data,104,12))& chisq=chisq(0)
  sexclass= (float(strmid(object_data,117,12)))[0]
  sexa    = (float(strmid(object_data,130,12)))[0]
  sexb    = (float(strmid(object_data,143,12)))[0]
  sexth   = (float(strmid(object_data,156,12)))[0]
  sexe1   = (float(strmid(object_data,169,12)))[0]
  sexe2   = (float(strmid(object_data,182,12)))[0]
  sexflux = (float(strmid(object_data,195,12)))[0]
  sexx    = (float(strmid(object_data,208,12)))[0]
  sexy    = (float(strmid(object_data,221,12)))[0]
  sexarea = (fix(strmid(object_data,234,12)))[0]
  objxmin = (fix(strmid(object_data,247,12)))[0]
  objxmax = (fix(strmid(object_data,260,12)))[0]
  objymin = (fix(strmid(object_data,273,12)))[0]
  objymax = (fix(strmid(object_data,286)))[0]

  ;help,x0,beta,fwhm,mag,nmax

  coefficients = coefficients+" "
  a=fltarr(header_info.n_coeffs)
  for i=0,a_used-1 do a(i)=float(strmid(coefficients,(13*i),12))

  errors = errors+" "
  e=fltarr(header_info.n_coeffs)
  for i=0,a_used-1 do e(i)=float(strmid(errors,(13*i),12))

  ; CHANGE STORAGE CONVENTION (APRIL 2004) TO KEEP TRUE VALUES OF ALL COEFFICIENTS
  if not keyword_set(v1) then begin
    a[1:a_used-1]=a[1:a_used-1]*a[0]
    e[1:a_used-1]=e[1:a_used-1]*e[0]
  endif
 
  ; Calculate size and magnitude from the coefficients themselves
  if keyword_set(moments) then begin
    if nmax ge 2 then begin ; quadrupole moments only nonzero if nmax>2
      oquad = shapelets_quadrupole({a:a(0:a_used-1), n_coeffs:a_used, n1:header_info.n1, $
              n2:header_info.n2, beta:beta, n_max:nmax},flux=fluxtemp) ; minimal decomp structure
      oell  = [oquad[0,0]-oquad[1,1],2*oquad[0,1]]/(oquad[0,0]+oquad[1,1])
      oflux = fluxtemp
      orsq  = (oquad[0,0]+oquad[1,1])/fluxtemp
      omag  = -2.5*alog10(fluxtemp) ; converting flux to magnitude-type
     
    endif else if nmax(obj) gt 0 then begin ; can still calculate flux
      orsq  = shapelets_rsquared({a:a(0:a_used-1), n_coeffs:a_used, n1:header_info.n1, $
                   n2:header_info.n2, beta:beta, n_max:nmax},flux=oflux)
      omag  = 5.*alog10(oflux)
      oquad = [[0.,0.],[0.,0.]]
      oell  = 0.

    endif
  endif else begin
   orsq=0. & omag=0. & oflux=0. & oell=0. & oquad=[[0.,0.],[0.,0.]]
  endelse
  
  ; Haven't yet coded up the calculation of flux errors for 1 object at a time!
  orsq_e=0. & oflux_e=0. & osnr=oflux/oflux_e
  
  ; Need to calculate centroid here, too!
  ocentre=[0.,0.]

  ; Convert to polar shapelets if requested
  if keyword_set(polar) then a=shapelets_polar_reduce(shapelets_c2p(a[0:a_used-1]));(0:a_used-1)))



endif else begin

  ;
  ; LOAD ENTIRE CATALOGUE
  ;
  if keyword_set(polar) then message,'Conversion to polars may take a few minutes',/info,noprint=silent

  ; Create variables to store information
  x0           = fltarr(header_info.n,2)
  beta         = fltarr(header_info.n)
  fwhm         = fltarr(header_info.n)
  mag          = fltarr(header_info.n)
  nmax         = intarr(header_info.n)
  oflux        = fltarr(header_info.n)
  oflux_e      = fltarr(header_info.n)
  osnr         = fltarr(header_info.n)
  orsq         = fltarr(header_info.n)
  orsq_e       = fltarr(header_info.n)
  oell         = fltarr(header_info.n,2)
  oell_e       = fltarr(header_info.n,2)
  oquad        = fltarr(header_info.n,2,2)
  oquad_e      = fltarr(header_info.n,2,2)
  ocentre      = fltarr(header_info.n,2)
  ocentre_e    = fltarr(header_info.n,2)
  omag         = fltarr(header_info.n)
  flag         = intarr(header_info.n,2)
  chisq        = fltarr(header_info.n)
  sexid        = intarr(header_info.n)
  sexclass     = fltarr(header_info.n)
  sexa         = fltarr(header_info.n)
  sexb         = fltarr(header_info.n)
  sexth        = fltarr(header_info.n)
  sexe1        = fltarr(header_info.n)
  sexe2        = fltarr(header_info.n)
  sexflux      = fltarr(header_info.n)
  sexx         = fltarr(header_info.n)
  sexy         = fltarr(header_info.n)
  sexarea      = intarr(header_info.n)
  objxmax      = intarr(header_info.n)
  objxmin      = intarr(header_info.n)
  objymin      = intarr(header_info.n)
  objymax      = intarr(header_info.n)
  object_data  = strarr(1)
  coefficients = strarr(1)
  errors       = strarr(1)
  n1=header_info.n1 & n2=header_info.n2
  readf,lun_h,junk
  if keyword_set(resize) then begin
    sizemag = 1
    new_n_max = header_info.n_max+5
    new_beta  = 5.
    a    = fltarr(header_info.n,(new_n_max+1)*(new_n_max+2)/2)
    e    = fltarr(header_info.n,(new_n_max+1)*(new_n_max+2)/2)
  endif else begin
    a    = fltarr(header_info.n,header_info.n_coeffs)
    e    = fltarr(header_info.n,header_info.n_coeffs)
  endelse

  ; Calculate the first coefficient which will represent a phase
  first_phase=(header_info.n_coeffs+header_info.n_max/2+1)/2
  ph_fact=1

  ; Work out Cartesian to Polar shapelet conversion matrix in advance (once only!)
  ;if keyword_set(polar) then junk=shapelets_c2p(fltarr(header_info.n_coeffs),matrix=c2pmatrix)
 
  ; Read in each object. The 3 lines correspond to: general info, coefficients, errors
  for obj=0,n-1 do begin

    ;print,obj+1;," /",n
    
    readf,lun_h,object_data & if strcmp(object_data(0),'# End',5) then break
    readf,lun_h,coefficients
    readf,lun_h,errors

    x0[obj,*]    = float([strmid(object_data,0,12),strmid(object_data,13,12)])
    beta[obj]    = float(strmid(object_data,26,12))
    fwhm[obj]    = float(strmid(object_data,39,12))
    mag[obj]     = float(strmid(object_data,52,12))
    nmax[obj]    = fix(strmid(object_data,65,12))
    flagtemp     = fix(strmid(object_data,78,12))
    flag[obj,0]  = flagtemp[0]/10 & flag[obj,1] = flagtemp[0] mod 10
    chisq[obj]   = float(strmid(object_data,104,12))
    sexid[obj]   = fix(strmid(object_data,91,12))
    sexclass[obj]= float(strmid(object_data,117,12))
    sexa[obj]    = float(strmid(object_data,130,12))
    sexb[obj]    = float(strmid(object_data,143,12))
    sexth[obj]   = float(strmid(object_data,156,12))
    sexe1[obj]   = float(strmid(object_data,169,12))
    sexe2[obj]   = float(strmid(object_data,182,12))
    sexflux[obj] = float(strmid(object_data,195,12))
    sexx[obj]    = float(strmid(object_data,208,12))
    sexy[obj]    = float(strmid(object_data,221,12))
    sexarea[obj] = fix(strmid(object_data,234,12))
    objxmin[obj] = fix(strmid(object_data,247,12))
    objxmax[obj] = fix(strmid(object_data,260,12))
    objymin[obj] = fix(strmid(object_data,273,12))
    objymax[obj] = fix(strmid(object_data,286))

    if keyword_set(n_max) then nmax(obj)=nmax(obj)<n_max

    a_used    = (nmax(obj)+1) * (nmax(obj)+2)/2 ; Don`t include zero coeffs for speed

    ; Read in shapelet coefficients for that object
    for i=0,a_used-1 do a(obj,i)=float(strmid(coefficients,(13*i),12))
    for i=0,a_used-1 do e(obj,i)=float(strmid(errors,(13*i),12))

    ; CHANGE STORAGE CONVENTION (APRIL 2004) TO KEEP TRUE VALUES OF ALL COEFFICIENTS
    if not keyword_set(v1) then begin
      a[obj,1:a_used-1]=a[obj,1:a_used-1]*a[obj,0]
      e[obj,1:a_used-1]=e[obj,1:a_used-1]*e[obj,0]
    endif
 
    ;; Just convert to polar shapelets if requested
    ;atemp=shapelets_polar_reduce(shapelets_c2p(reform(a(obj,0:a_used-1)),matrix=matrix))
    ;n_n0=nmax(obj)/2+1       ; number of coefs where n=0
    ;n_phases=(a_used-n_n0)/2 ; number of coefs representing complex phases
    ;n_moduli=(a_used+n_n0)/2 ; number of coefs representing complex moduli
    ;a(obj,0:n_moduli-1)=atemp(0:n_moduli-1)
    ;a(obj,first_phase:first_phase+n_phases-1)=atemp(n_moduli:n_moduli+n_phases-1)
    ;endif

    ; Calculate unweighted quadrupole moments from the coefficients 
    ; This gives size, magnitude and ellipticity
    if keyword_set(moments) then begin
       
      ; Whatever happens, we will need a fake decomp structure for this object
      decomptemp = {coeffs:[a(obj,0),a(obj,0)*reform(a(obj,1:a_used-1))], $
     			error:[e(obj,0),e(obj,0)*reform(e(obj,1:a_used-1))], $
     			coeffs_error:[e(obj,0),e(obj,0)*reform(e(obj,1:a_used-1))], $
     			x:x0[obj,*], $
     			n_coeffs:a_used, n1:n1, n2:n2, beta:beta(obj), n_max:nmax(obj)}
      
      if nmax[obj] ge 2 then begin ; quadrupole moments only nonzero if nmax>2
       
        Q 		    = shapelets_quadrupole(decomptemp,/cartesian,/error,$
     				   flux=fluxtemp,rsquared=rsqtemp,ellipticity=elltemp)
        oquad[obj,*,*]	= Q[*,*,0]
        oquad_e[obj,*,*] = Q[*,*,0]
        oell[obj,*] 	= [float(elltemp[0]),imaginary(elltemp[0])];[Q[0,0]-Q[1,1],2*Q[0,1]]/(Q[0,0]+Q[1,1])
        oell_e[obj,*]	= [float(elltemp[1]),imaginary(elltemp[1])]
        oflux[obj]  	= fluxtemp[0]
        oflux_e[obj]	= fluxtemp[1]
        osnr[obj]  	 = oflux[obj]/oflux_e[obj]
        omag[obj]		= -2.5*alog10(fluxtemp[0]) ; convert flux to magnitude
        orsq[obj]		= rsqtemp[0]
        orsq_e[obj] 	= rsqtemp[1]
        cnttemp		= shapelets_centroid(decomptemp,/cartesian,/error)
        ocentre[obj,*]	= cnttemp[*,0]
        ocentre_e[obj,*] = cnttemp[*,1]

      endif else if nmax[obj] eq 1 then begin ; can still calculate centroid and flux
         
        matrix=0
        orsq[obj]		= shapelets_rsquared(decomptemp,flux=fluxtemp,matrix=matrix)
        oflux[obj]  	= fluxtemp[0]
        omag[obj]		= -2.5*alog10(fluxtemp[0])
        cntrd  		= shapelets_centroid(decomptemp,matrix=matrix,/error)
        ocentre[obj,*]	= cntrd[0:1]
        ocentre_e[obj,*] = cntrd[2:3]
      	   
      endif else if nmax[obj] eq 0 then begin ; can still calculate flux
         
        orsq[obj]  = shapelets_rsquared(decomptemp,flux=fluxtemp)
        oflux[obj] = fluxtemp
        omag[obj]  = -2.5*alog10(fluxtemp) ; converting flux to magnitude-type
     
      endif	    
    endif

  
;   ; Get all objects to have the same beta
;   if keyword_set(resize) then begin
;      decomptemp  = {a:[a[obj,0],a[obj,0]*reform(a(obj,1:a_used-1))], n_coeffs:a_used, n1:n1, n2:n2, w:beta(obj), $
;      nmax:nmax(obj), x0:[15,15], nf:[30,30], over:1, integrate:1, error:1}
;      ; Resize - or could do by repeating 1st order until rsq same as before
;      shapelets_dilate,decomptemp,(beta(obj)/new_beta),n_max=new_n_max  ; would be quicker to set n_max=n_max+4 then increase up to new_nmax
;      a(obj,0:*)  = decomptemp.a
;      a_used	   = (new_n_max+1)*(new_n_max+2)/2
;      nmax(obj)   = new_n_max
;      decomptemp.w=new_beta
;      ;print,beta(obj),new_beta,orsq(obj),rsquared(decomptemp)
;      beta(obj)   = new_beta
;   endif



    ; Do parity flip and rotation alignment
    if keyword_set(parity) then begin
    
      ; First convert to polars
      ; long-winded but works
      ;a(obj,*)=transpose(shapelets_polar_reduce(shapelets_c2p(transpose(a(obj,*)))))
      ; quicker but need to rearrange phase coefs afterwards
      atemp=shapelets_c2p(reform(a[obj,*]),/REDUCE,matrix=matrix)
      ;atemp=shapelets_polar_reduce(shapelets_c2p(reform(a(obj,0:a_used-1)),matrix=matrix))

      n_n0=nmax(obj)/2+1       ; number of coefs where n=0
      n_phases=(a_used-n_n0)/2 ; number of coefs representing complex phases
      n_moduli=(a_used+n_n0)/2 ; number of coefs representing complex moduli

      ; Now look at orientation
      if nmax(obj) ge 2 then begin
	   orientation=atemp(n_moduli+1)
      endif else orientation=0
	 
      ; Now look at parity
      if nmax(obj) ge 4 then begin
        ph_fact = sign( sin(( atemp(n_moduli+1) - atemp(n_moduli+4) ))  )
      endif else ph_fact=1

      a(obj,0:n_moduli-1)=atemp(0:n_moduli-1)
      a(obj,first_phase:first_phase+n_phases-1)=ph_fact*(atemp(n_moduli:n_moduli+n_phases-1)-orientation)
      
	 ; Convert back to Cartesians
      if not keyword_set(polar) then begin
	   atemp=shapelets_p2c(shapelets_polar_expand(reform(a[obj,0:a_used-1])))
	   a[obj,0:a_used-1]=atemp
      endif
    
    endif else if keyword_set(polar) then begin

      a[obj,*]=shapelets_c2p(a[obj,*],/REDUCE,matrix=c2p_matrix)
 
    endif

  endfor

endelse

if keyword_set(n_max) then maxn_max=header_info.n_max<nmax else maxn_max=header_info.n_max

if keyword_set(v1) then type="shapecat_v1" else type="shapecat"

;Store in a structure
cat={name:header_info.name,      $
     type:type,                  $
     n:n,                        $
     maxn_max:maxn_max,          $
     polar:byte(polar),          $
     seeing:header_info.seeing,  $
     x:x0,          	        $
     beta:beta,         	        $
     n_max:nmax,     	        $
     n_coeffs:n_coeffs,          $
     coeffs:a,            	   $
	coeffs_error:e,             $
	error:e,        	        $
     flux:oflux,     	        $
     flux_error:oflux_e,         $
     snr:osnr,                   $
     mag:omag,                   $
     centroid:ocentre,           $
     centroid_error:ocentre_e,   $
     rsquared:orsq,              $
     rsquared_error:orsq_e,      $
     ellipticity:oell,           $
     ellipticity_error:oell_e,   $
     quadrupole:oquad,           $
     quadrupole_error:oquad_e,   $
     flag:flag,                  $
     chisq:chisq,                $
     sexid:sexid,                $
     sexfwhm:fwhm,               $
     sexflux:sexflux,            $
     sexmag:mag,                 $
     sexclass:sexclass,          $
     sexa:sexa,                  $
     sexb:sexb,                  $
     sextheta:sexth,             $
     sexe1:sexe1,                $
     sexe2:sexe2,                $
     sexx:sexx,                  $
     sexy:sexy,                  $
     sexarea:sexarea,            $
     objx_max:objxmax,           $
     objx_min:objxmin,           $
     objy_min:objymin,           $
     objy_max:objymax}
	

; Close up and finish
close,lun_h
free_lun,lun_h
message,'Catalogue successfully in memory',/info,noprint=silent

end



