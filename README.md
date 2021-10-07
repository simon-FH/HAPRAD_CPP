HAPRAD_CPP
==========

DONT FORGET TO COPY DICT FILE INTO SLIB/!!!!

All the steps to make Pi Plus Analysis (Multilpicity Ratio and Braodening) are detailed here.

To compile the library please make sure you have installed make, imake (xutils-dev in Ubuntu), cernlib and ROOT

First it is necessary to set the enviorment variables

setenv GITPIPLUS PiPlusDirectory
setenv CLASTOOL ${GITPIPLUS}/ClasToolPack/ClasTool
setenv OS_NAME Linux
setenv CLAS_PACK ${GITPIPLUS}/ClasToolPack/CLAS_PACK
setenv CLAS_LIB ${GITPIPLUS}/ClasToolPack/CLAS_LIB
setenv ANALYSER ${GITPIPLUS}/ClasToolPack/Analyser/analysis_lib
setenv MYSQL_INCLUDE $MYSQLINC
setenv MYSQL_LIB $MYSQLIB
setenv LD_LIBRARY_PATH ${CLASTOOL}/slib/${OS_NAME}:${LD_LIBRARY_PATH}
setenv LD_LIBRARY_PATH ${ANALYSER}/slib:$${LD_LIBRARY_PATH}

To compile HAPRAD_CPP it is necessary to have libMathMore in root, to do this libgsl need to be installed and in configure set --enable-gsl-shared (make sure configure output is "Checking whether to build libMathMore ... yes")

cernlib
pdflib804 (libpdflib804-2-dev)
phtools (libphtools2-dev)

## Authors

H. Hakobyan, R. Oyarzun and S. Mancilla.
