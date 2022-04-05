#!/bin/sh

EFIT=../EFIT/eqdsk.py
EXEC=../ENV/exec_dirs.py
FCHK=../GFIT/fit_checkfile.py
FIT=../GFIT/fittool.py
GEFIT=../EFIT/get_efit.py
MDS=../MDS/MDS.py
PROG=../ENV/progress.py
REFIT=../EFIT/rtefit.py
BATCH=../ENV/batch_run.py
CHTOO=../CHEASE/ch_tool.py
NODE=../ENV/nodelist.py
NAME=../ENV/read_namelist.py
GMSE=../GEFIT/gefit_mse.py
GTOO=../GEFIT/gefit_tool.py
GECE=../GFIT/gfit_ece_multi.py
KNOT=../ENV/knots_tool3.py
MMDS=../MDS/multi_mds.py
AEFIT=../EFIT/aeqdsk.py
CSAPS=../CSAPS/csaps/

cd bin
chmod 777 gefit gfit
rm -f aeqdsk bs2k dena denaf rdena rdenaf reqdsk fgefit gped pchease ercal gjastab readp
ln -s ../GEFIT/aeqdsk.py aeqdsk
ln -s ../BS2K/bs2kstar bs2k
ln -s ../DENA/dena.py dena
ln -s ../DENA/denaf.py denaf
ln -s ../DENA/rdena.py rdena
ln -s ../DENA/rdenaf.py rdenaf
ln -s ../EFIT/eqdsk.py reqdsk
ln -s ../FGEFIT/fgefit.py fgefit
ln -s ../EPED/gui_eped.py gped
ln -s ../CHEASE/chease pchease
ln -s ../ERGEN/ercal.py ercal
ln -s ../JASTAB/gui_jastab.py gjastab
ln -s ../GFIT/read_pfile.py readp
cd ../

cd BS2K
chmod 777 bs2kstar
rm -f exec_dirs.py
ln -s $EXEC
cd ../

cd CHEASE
chmod 777 chease
rm -f eqdsk.py exec_dirs.py
ln -s $EXEC
ln -s $EFIT
cd ../

cd DENA
chmod 777 denaf.py dena.py rdenaf.py rdena.py
rm -f eqdsk.py exec_dirs.py fit_checkfile.py fittool.py get_efit.py MDS.py progress.py rtefit.py
rm -rf csaps
ln -s $EFIT 
ln -s $EXEC 
ln -s $FCHK 
ln -s $FIT 
ln -s $GEFIT 
ln -s $MDS 
ln -s $PROG 
ln -s $REFIT
ln -s $CSAPS
cd ../

cd EFIT
chmod 777 aeqdsk.py rtefit.py efit*
rm -f exec_dirs.py
ln -s $EXEC
cd ../

cd EPED
chmod 777 eped_gfit.py eped.py gui_eped.py
rm -f batch_run.py ch_tool.py eqdsk.py exec_dirs.py fit_checkfile.py nodelist.py progress.py read_namelist.py
ln -s $BATCH
ln -s $CHTOO
ln -s $EFIT
ln -s $EXEC
ln -s $FCHK
ln -s $NODE
ln -s $PROG
ln -s $NAME
cd ../

cd ERGEN
chmod 777 ercal.py
rm -f exec_dirs.py eqdsk.py
ln -s $EFIT
ln -s $EXEC
cd ../

cd FGEFIT
chmod 777 fgefit.py
rm -f batch_run.py ch_tool.py eqdsk.py exec_dirs.py gefit_mse.py gefit_tool.py get_efit.py gfit_ece_multi.py 
rm -f knots_tool3.py MDS.py multi_mds.py progress.py
ln -s $BATCH
ln -s $CHTOO
ln -s $EFIT
ln -s $EXEC
ln -s $GMSE
ln -s $GTOO
ln -s $GEFIT
ln -s $GECE
ln -s $KNOT
ln -s $MDS
ln -s $MMDS
ln -s $PROG
cd ../

cd GEFIT
chmod 777 efit_rmp.py gefit.py
rm -f aeqdsk.py batch_run.py ch_tool.py eqdsk.py exec_dirs.py get_efit.py knots_tool3.py MDS.py multi_mds.py
ln -s $AEFIT
ln -s $BATCH
ln -s $CHTOO
ln -s $EFIT
ln -s $EXEC
ln -s $GEFIT
ln -s $KNOT
ln -s $MDS
ln -s $MMDS
cd ../

cd GFIT
chmod 777 gfitp.py guifit.py read_pfile.py TS_NE_dummy.dat
rm -f batch_run.py ch_tool.py csaps eqdsk.py exec_dirs.py gefit_mse.py gefit_tool.py get_efit.py 
rm -f knots_tool3.py MDS.py multi_mds.py progress.py
rm -rf csaps
ln -s $BATCH
ln -s $CHTOO
ln -s $CSAPS
ln -s $EFIT
ln -s $EXEC
ln -s $GMSE
ln -s $GTOO
ln -s $GEFIT
ln -s $KNOT
ln -s $MDS
ln -s $MMDS
ln -s $PROG
cd ../

cd JASTAB
chmod 777 gui_jastab.py jadiag.py japlot.py
rm -f batch_run.py ch_tool.py eqdsk.py exec_dirs.py nodelist.py progress.py read_namelist.py
ln -s $BATCH
ln -s $CHTOO
ln -s $EFIT
ln -s $EXEC
ln -s $NODE
ln -s $PROG
ln -s $NAME
cd ../

cd MDS
chmod 777 gui_mds.py gui_mds.py_old plot_da.py
rm -f aeqdsk.py eqdsk.py exec_dirs.py get_efit.py little_gui5.py progress.py
ln -s $AEFIT
ln -s $EFIT
ln -s $EXEC
ln -s $GEFIT
ln -s gui_mds.py little_gui5.py
ln -s $PROG
cd ../

cd NUBEAM
chmod 777 mdescr_A123B123.dat nubeam_run.py sconfig_A123B123.dat
rm -f ch_tool.py eqdsk.py exec_dirs.py
ln -s $CHTOO
ln -s $EFIT
ln -s $EXEC
cd ../

cd PEDSCAN
chmod 777 gui_pedscan.py ped_scanner
rm -f batch_run.py ch_tool.py eqdsk.py exec_dirs.py fit_checkfile.py fittool.py nodelist.py progress.py read_namelist.py
rm -rf csaps
ln -s $BATCH
ln -s $CHTOO
ln -s $EFIT
ln -s $EXEC
ln -s $FCHK
ln -s $FIT
ln -s $NODE
ln -s $PROG
ln -s $NAME
ln -s $CSAPS
cd ../
