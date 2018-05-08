# polykerrbb
Modification to the XSPEC model kerrbb, adding polynomial interpolation of table values.

Files included:
Source code:
kerrbbpoly.cxx
runkbbpoly.f

XSPEC model definition:
lmodel.dat

Compilation script:
compile_polykerrbb.sh


Installation:
1. Copy all four files to an appropriate directory - this directory should replace <path> in the following steps.
2. Run <path>/compile_polykerrbb.sh
  This may produce a conversion and many unused variable warnings, but should end with:

"
Local model library has been built from model definition and code files in:
/data/djkb2/kerrbb/public

XSPEC12>
XSPEC12>lmod package_polykerrbb /data/djkb2/kerrbb/public
Model package package_polykerrbb successfully loaded.
XSPEC12>
XSPEC12>
 XSPEC: quit
"

Model loading:
a. for an individual XSPEC session:
XSPEC12>lmod package_polykerrbb <path>

b. automatically for all future sessions:
add the line from step a. to your xspec.rc file (which may need to be created, typically in ~/.xspec).
