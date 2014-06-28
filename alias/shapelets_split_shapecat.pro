function shapelets_split_shapecat,shapecat_in,_ref_extra = ex
COMPILE_OPT IDL2, HIDDEN, OBSOLETE
shapecat_out=shapecat_in
shapelets_split,shapecat_out,selected,_extra= ex
return,shapecat_out
end
