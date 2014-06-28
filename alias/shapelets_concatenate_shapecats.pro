function shapelets_concatenate_shapecats, shapecat1, shapecat2
COMPILE_OPT idl2, OBSOLETE
shapecat=shapecat1
shapelets_add, shapecat, shapecat2
return, shapecat
end
