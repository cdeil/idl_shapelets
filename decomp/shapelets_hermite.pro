;$Id: shapelets_hermite.pro, v2$
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


function shapelets_hermite_coeffs, n
;
; NAME:
;       SHAPELETS_HERMITE_COEFFS
;
; CATEGORY:
;       A component of shapelets_hermite.pro.
;
; PURPOSE:
;       Compute the polynomial coefficients for Hn(x), the 1D Hermite
;       polynomial of order n
;
; INPUTS:
;       n - order of the Hermite polynomial
;
; KEYWORD PARAMETERS:
;       None.
;
; OUTPUTS:
;       h - vector of coefficients ci for Hn(x)=Sum_i ci*x^i, i=0,..,n
;
; MODIFICATION HISTORY:
;       Jul 05 - Floating point overflows beyond n=20 fixed by RM.
;       Sep 03 - Ccombined with shapelets_hermite.pro by Richard Massey.
;       Jul 99 - Written by A. Refregier
;

COMPILE_OPT idl2, HIDDEN

if n le 20 then begin

  coeffs = lon64arr(n+1)
  for i=0L,fix(n/2) do $
    coeffs[n-2*i] = (-1)^i * 2ULL^(n-2*i) * ( factorial(n,/UL64) / ( factorial(i,/UL64) * factorial(n-2*i,/UL64) ))

endif else begin

  ; n! doesn't fit into a UL64 integer if n>20, so use floating points instead
  coeffs = dblarr(n+1)
  for i=0L,fix(n/2) do $
    coeffs[n-2*i] = (-1)^i * 2.d0^(n-2*i) * factorial(double(n)) / ( factorial(double(i)) * factorial(double(n-2*i)) )

endelse

return, coeffs

end






function shapelets_hermite, n, x
;+
; NAME:
;       SHAPELETS_HERMITE
;
; CATEGORY:
;       Mathematical functions.
;
; PURPOSE:
;       Compute the 1D Hermite polynomial Hn(x) of order n.
;       Faster than the astlib HERMITE function because it stores a
;       pre-compiled set of low-order hermite functions.
;
; CALLING PROCEDURE:
;       result=shapelets_hermite(n,x)
;
; INPUTS:
;       n - order of the Hermite polynomial
;       x - can be a scalar, vector or array(coordinate grid)
;
; OUTPUTS:
;       Hn(x) - Hermite polynomial order n evaluated at position x
;
; PROCEDURES USED:
;       shapelets_hermite_coeffs.
;
; MODIFICATION HISTORY:
;       Jul 05 - Sped up by RM using IDL's built-in poly function.
;       Oct 03 - Number of calculations done by default increased by R. Massey
;       Jul 99 - Written by A. Refregier
;-

COMPILE_OPT idl2

nx = n_elements(x)
x=double(x)

case n of
   0: begin
        sx = size(x)
        case sx[0] of
           0: h = 1.
           1: h = replicate(1., sx[1])
           2: h = replicate(1., sx[1], sx[2])
        else: message, 'x must be a scalar, a vector or an array'
        endcase
      end
;   1: h = 2.*x
;   2: h = 4.*x^2-2.
;   3: h = 8*x^3-12.*x
;   4: h = 16*x^4-48.*x^2+12.
;   5: h = 32*x^5-160.*x^3+120.*x
;   6: h = 64.*x^6-480.*x^4+720.*x^2-120.
;   7: h = 128.d0*x^7-1344.d0*x^5+3360.d0*x^3-1680*x
;   8: h = 256.d0*x^8-3584.d0*x^6+13440.d0*x^4-13440.d0*x^2+1680.d0
;   9: h = 512.d0*x^9-9216.d0*x^7+48348.d0*x^5-80640.d0*x^3+30240.d0*x
;  10: h = 1024.d0*x^10-23040.d0*x^8+161280.d0*x^6-403200.d0*x^4+302400.d0*x^2-30240.d0
;  11: h = -665280.d0*x       + 2217600.d0*x^3         - 1774080.d0*x^5        + 506880.d0*x^7         -  56320.d0*x^9         + 2048.d0*x^11
;  12: h =  665280.d0         - 7983360.d0*x^2         + 13305600.d0*x^4       - 7096320.d0*x^6        + 1520640.d0*x^8        - 135168.d0*x^10        + 4096.d0*x^12
;  13: h =  17297280.d0*x     - 69189120.d0*x^3        + 69189120.d0*x^5       - 26357760.d0*x^7       + 4392960.d0*x^9        - 319488.d0*x^11        + 8192.d0*x^13
;  14: h = -17297280.d0       + 242161920.d0*x^2       - 484323840.d0*x^4      + 322882560.d0*x^6      - 92252160.d0*x^8       + 12300288.d0*x^10      - 745472.d0*x^12       + 16384.d0*x^14
;  15: h = -518918400.d0*x    + 2421619200.d0*x^3      - 2905943040.d0*x^5     + 1383782400.d0*x^7     - 307507200.d0*x^9      + 33546240.d0*x^11      - 1720320.d0*x^13      + 32768.d0*x^15
;  16: h =  518918400.d0      - 8302694400.d0*x^2      + 19372953600.d0*x^4    - 15498362880.d0*x^6    + 5535129600.d0*x^8     - 984023040.d0*x^10     + 89456640.d0*x^12     - 3932160.d0*x^14     + 65536.d0*x^16
;  17: h =  17643225600.d0*x  - 94097203200.d0*x^3     + 131736084480.d0*x^5   - 75277762560.d0*x^7    + 20910489600.d0*x^9    - 3041525760.d0*x^11    + 233963520.d0*x^13    - 8912896.d0*x^15     + 131072.d0*x^17
;  18: h = -17643225600.d0    + 317578060800.d0*x^2    - 846874828800.d0*x^4   + 790416506880.d0*x^6   - 338749931520.d0*x^8   + 75277762560.d0*x^10   - 9124577280.d0*x^12   + 601620480.d0*x^14   - 20054016.d0*x^16   + 262144.d0*x^18
;  19: h = -670442572800.d0*x + 4022655436800.d0*x^3   - 6436248698880.d0*x^5  + 4290832465920.d0*x^7  - 1430277488640.d0*x^9  + 260050452480.d0*x^11  - 26671841280.d0*x^13  + 1524105216.d0*x^15  - 44826624.d0*x^17   + 524288.d0*x^19
;  20: h =  670442572800.d0   - 13408851456000.d0*x^2  + 40226554368000.d0*x^4 - 42908324659200.d0*x^6 + 21454162329600.d0*x^8 - 5721109954560.d0*x^10 + 866834841600.d0*x^12 - 76205260800.d0*x^14 + 3810263040.d0*x^16 - 99614720.d0*x^18 + 1048576.d0*x^20
;else: begin
;        c = shapelets_hermite_coeffs(n)
;        if size(x,/n_dimensions) eq 0 then h = total(c*x^findgen(n+1)) else begin
;           h = c[0]+c[1]*x   ; Get h to have correct dimensions
;           for i=2, n do h=h+c[i]*x^i
;        endelse
;      end
   1: h = poly(x,[0,2])
   2: h = poly(x,[-2,0,4])
   3: h = poly(x,[0,-12,0,8])
   4: h = poly(x,[12,0,-48,0,16])         
   5: h = poly(x,[0,120,0,-160,0,32])     
   6: h = poly(x,[-120,0,720,0,-480,0,64])
   7: h = poly(x,[0,-1680,0,3360,0,-1344,0,128])
   8: h = poly(x,[1680,0,-13440,0,13440,0,-3584,0,256])
   9: h = poly(x,[0,30240,0,-80640,0,48348,0,-9216,0,512])
  10: h = poly(x,[-30240,0,302400,0,-403200,0,161280,0,-23040,0,1024])
  11: h = poly(x,[0,-665280,0,2217600,0,-1774080,0,506880,0,-56320,0,2048])
  12: h = poly(x,[665280,0,-7983360,0,13305600,0,-7096320,0,1520640,0,-135168,0,4096])
  13: h = poly(x,[0,17297280,0,-69189120,0,69189120,0,-26357760,0,4392960,0,-319488,0,8192])
  14: h = poly(x,[-17297280,0,242161920,0,-484323840,0,322882560,0,-92252160,0,12300288,0,-745472,0,16384])
  15: h = poly(x,[0,-518918400,0,2421619200,0,-2905943040,0,1383782400,0,-307507200,0,33546240,0,-1720320,0,32768])
  16: h = poly(x,[518918400,0,-8302694400,0,19372953600,0,-15498362880,0,5535129600,0,-984023040,0,89456640,0,-3932160,0,65536])
  17: h = poly(x,[0,17643225600,0,-94097203200,0,131736084480,0,-75277762560,0,20910489600,0,-3041525760,0,233963520,0,-8912896,0,131072])
  18: h = poly(x,[-17643225600,0,317578060800,0,-846874828800,0,790416506880,0,-338749931520,0,75277762560,0,-9124577280,0,601620480,0,-20054016,0,262144])
  19: h = poly(x,[0,-670442572800,0,4022655436800,0,-6436248698880,0,4290832465920,0,-1430277488640,0,260050452480,0,-26671841280,0,1524105216,0,-44826624,0,524288])
  20: h = poly(x,[670442572800,0,-13408851456000,0,40226554368000,0,-42908324659200,0,21454162329600,0,-5721109954560,0,866834841600,0,-76205260800,0,3810263040,0,-99614720,0,1048576])
else: h = poly(x,shapelets_hermite_coeffs(n))
endcase

return, h

end
