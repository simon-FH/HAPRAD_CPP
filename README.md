# PiPlusAnalysis

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

To compile HAPRAD_CPP it is necessary to have libMathMore in root, to do this libgsl need to be installed and
in configure set --enable-gsl-shared
(make sure configure output is "Checking whether to build libMathMore ... yes")
* cernlib
* pdflib804 (libpdflib804-2-dev)
* phtools (libphtools2-dev)

## Instructions

1. Primero que todo, deben clonar el repositorio de la pagina de Orlando:

~ $ git clone http://github.com/orsosa/PiPlusAnalysis

2. Luego, generar la variable de entorno editando el archivo .bashrc agregando

export HAPRAD_CPP=/CARPETADEUSUARIOODONDEHAYANCLONADO/PiPlusAnalysis/RCFactor/HAPRAD_CPP

(No olvidar que tienen que cerrar la terminal y abrirla de nuevo, o escribir source .bashrc)

3. Despues, compilar cada programa ejecutando el Makefile de cada carpeta

cd PiPlusAnalysis
cd RCFactor
cd HAPRAD_CPP
make
cd ../Gen5DimData
make
cd ../PhiHist
make
cd ../GenRC
make

4. Si todo salio bien, correr cada programa (pueden empezar desde el segundo programa usando mi linea de instruccion)

  ((Otra advertencia: pueden ver el modelo de archivo de datos con mis archivos en: /user/b/borquez/PiData/))

* El primer programa es Gen5DimData:

  cd Gen5DimData/

  La instruccion para ejecutarlo tiene el siguiente formato:

  ./5dim_acc_corr Q2_MIN Q2_MAX N_Q2 XB_MIN XB_MAX N_XB NU_MIN NU_MAX N_NU ZH_MIN ZH_MAX N_ZH PT_MIN PT_MAX N_PT PHI_MIN PHI_MAX N_PHI /LOCATIONOFDATA/ fDataExt fSimuExt nSimul elecExt pionExt XF_CUT Metal

  En mi caso:

  ./5dim_acc_corr 1 4 6 0.1 0.55 5 2.2 4.2 10 0 1 10 0 1 5 -180 180 12 /user/b/borquez/PiData/ _data _simul 1 DOESNTMATTER .root 0 C

  (El problema es que en mi carpeta no tengo las simulaciones, ojala tengan todo junto)

  Finalmente, se deberia generar un programa llamado C_5_dim_dist.root (en mi caso) en la carpeta del programa.

* El segundo programa es PhiHist:

  cd PhiHist/

  La instruccion para ejecutarlo tiene el siguiente formato:

  ./phihist /ubicaciondeMETAL_5_dim_dist.root/ PHI_MIN PHI_MAX N_PHI Metal(OPT)

  En mi caso:

  ./phihist /user/b/borquez/PiPlusAnalysis/RCFactor/Gen5DimData/ -180 180 12 C

  Finalmente, se deberia generar un programa llamado Cnewphihist.root en la carpeta del programa.

* El ultimo programa es GenRC:

  cd GenRC/

  La instruccion para ejecutarlo tiene el siguiente formato:

  ./genRC Q2_MIN Q2_MAX N_Q2 XB_MIN XB_MAX N_XB ZH_MIN ZH_MAX N_ZH PT_MIN PT_MAX N_PT PHI_MIN PHI_MAX N_PHI /ubicaciondeMetalnewphihist.root/ Metal

  En mi caso:

  ./genRC 1 4 6 0.1 0.55 5 0 1 10 0 1 5 -180 180 12 /user/b/borquez/PiPlusAnalysis/RCFactor/PhiHist/ C

  Finalmente, se deberia generar un programa llamado RCFactorC.txt en la carpeta del programa.
