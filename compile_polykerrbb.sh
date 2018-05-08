
hmake clean

xspec <<END

initpackage package_polykerrbb lmodel.dat `pwd`

lmod package_polykerrbb `pwd`

END
