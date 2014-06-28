; $Id: gamma.pro,v 1.10 2001/05/04 21:01:57 ali Exp $
;
; Copyright (c) 1995-2001, Research Systems, Inc.  All rights reserved.
;       Unauthorized reproduction prohibited.

;+
; NAME:
;       GAMMA
; PURPOSE:
;	The introduction of the GAMMA function as a built in system routine
;	in IDL 4.0 caused inconvenience to customers with existing code
;	in which GAMMA had been used as a variable, because IDL system
;	routines have precedence over variable names. To minimize this
;	problem, RSI has renamed GAMMA to NR_GAMMA.
;
;	This wrapper serves to make NR_GAMMA available under the name
;	GAMMA as documented in the IDL Reference Manual. However, since
;	IDL library routines have lower precedence than variables, programs
;	that use GAMMA as a variable name will work as before.
;
;	See the documentation for GAMMA in the IDL Reference manual for
;	details on arguments, keywords, and results.
;
;
; MODIFICATION HISTORY:
;	3 July 1995, AB, RSI.
;	AB, 5/4/2001, Switch from using _EXTRA to _STRICT_EXTRA, so that
;		incorrect keywords will cause issue proper messages to
;		be issued instead of being silently ignored.
;-

function gamma, x, _REF_EXTRA=extra
  ON_ERROR, 2
  return, NR_GAMMA(x, _strict_extra=extra)
end
