##---machine inform
#job scheduler
scheduler = 'pbs'       #sge
qsub_exec = '/opt/pbs/bin/qsub'
qstat_exec = '/opt/pbs/bin/qstat'
qdel_exec = '/opt/pbs/bin/qdel'
#node list option
node_machine = 'ukstar' #fusma
node_force=True         #True
node_init='compute'     #'none'
node_default='compute'  #old_group.q@node10
#efit links
efit_source_dir = '/home/users/efit'
efit_source_dir = '/EFIT'
efit_address = '172.17.250.23'
#efit years
shotk = dict()
years = ['2015','2016','2017','2018','2019','2020','2021','2022']
for i in years: shotk[i] = dict()
shotk['2015']['shot'] = range(12289,14389)
shotk['2016']['shot'] = range(14954,17364)
shotk['2017']['shot'] = range(17857,19393)
shotk['2018']['shot'] = range(19815,21760)
shotk['2019']['shot'] = range(21761,24080)
shotk['2020']['shot'] = range(24180,27400)
shotk['2021']['shot'] = range(27401,30450)
shotk['2022']['shot'] = range(30451,36000)
for year in years:
      for efit_no in range(1,6):
              shotk[year][efit_no] = '/EFIT_RUN/EFITDATA_%s/EFIT%02i'%(year,efit_no)
#mse favorable channels.
mse_good_ch = dict()
for i in years: mse_good_ch[i] = []
mse_good_ch['2015'] = [0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0]
mse_good_ch['2016'] = [1,1,1,0,1, 1,1,1,1,1, 1,1,0,0,0, 0,0,0,0,0, 0,0,0,0,0]
mse_good_ch['2017'] = [0,1,1,1,1, 1,1,1,1,0, 0,1,1,1,0, 0,0,1,1,0, 0,0,0,0,0]
mse_good_ch['2018'] = [0,0,1,1,1, 1,1,1,1,1, 1,1,1,1,0, 0,0,1,0,0, 0,0,0,0,0]
mse_good_ch['2019'] = [0,0,1,1,1, 0,1,0,0,0, 1,1,1,0,1, 0,0,1,1,0, 0,0,0,0,0]
mse_good_ch['2020'] = [1,1,1,0,1, 1,1,1,1,1, 1,1,0,0,0, 0,0,1,0,0, 0,0,0,0,0]
mse_good_ch['2021'] = [1,1,1,0,1, 1,1,1,1,1, 1,1,0,0,0, 0,0,1,0,0, 0,0,0,0,0]
mse_good_ch['2022'] = [1,1,1,0,1, 1,1,1,1,1, 1,1,0,0,0, 0,0,1,0,0, 0,0,0,0,0]

##---Executables
#chease
chease_exec  = '/usr/local/analysis/EQUI/CHEASE/bin/chease'
#helena
helena_exec  = '/usr/local/analysis/EQUI/HELENA/bin/hel13'
#mishka
mis_exec     = '/usr/local/analysis/STAB/MISHKA/bin/mishka1fast_'
#elite
elite_exec   = '/usr/local/analysis/STAB/ELITE/bin/elite'
elite_dir    = '/usr/local/analysis/STAB/ELITE/bin'
#gzip
gzip_dir     = 'gzip'
#python
python2_exec = '/usr/local/anaconda2/bin/python2'
python3_exec = '/usr/local/anaconda3/bin/python3'
pythonc_exec = python3_exec

##---python scripts
#python_home  = '/home/ksk911211/PYTHON/KEFIT_py/'
python_home  = '/usr/local/analysis/KEFIT/'
#chease
chease_dir    = python_home+'/CHEASE/chease'
#nubeam
nubeam_dir    = python_home+'/NUBEAM/nubeam.py'
nubeam_dir2   = python_home+'/NUBEAM/nubeam_run.py'
Mfile         = python_home+'/NUBEAM/mdescr_A123B123.dat'
Sfile         = python_home+'/NUBEAM/sconfig_A123B123.dat'
Ifile         = python_home+'/NUBEAM/nubeam_init.dat'
stepfile      = python_home+'/NUBEAM/nubeam_step.dat'
plasma_state_test_exec='/usr/local/analysis/KEFIT/NTCC/LINUX/test/plasma_state_test'
nubeam_comp_exec      ='/usr/local/analysis/KEFIT/NTCC/LINUX/test/mpi_nubeam_comp_exec'
adasdir               ='/usr/local/analysis/KEFIT/NTCC/adas310_fortran_driver'
preactdir             ='/usr/local/analysis/KEFIT/NTCC/preact'
mpirun                ='/usr/local/mpich/bin/mpirun'
#eped
stab_dir      = python_home+'/EPED/eped_stab.py'
gped_dir      = python_home+'/EPED/gui_eped.py'
gfit2_dir     = python_home+'/EPED/eped_gfit.py'
eped_dir      = python_home+'/EPED/eped.py'
#pedscan
pedscan_dir   = python_home+'/PEDSCAN/gui_pedscan.py'
pedscane_dir  = python_home+'/PEDSCAN/ped_scanner'
pedscane_dir2 = python_home+'/PEDSCAN/pedscan.py'
pedstab_dir   = python_home+'/PEDSCAN/pedstab_chease.py'
#jastbility
jastab_dir    = python_home+'/JASTAB/jadiag.py'
jastabc_dir   = python_home+'/JASTAB/ja_stab_chease.py'
jastabh_dir   = python_home+'/JASTAB/ja_stab_helena.py'
japlot_dir    = python_home+'/JASTAB/japlot.py'
#mds
mds_dir       = python_home+'/MDS/gui_mds.py'
mds_dir2      = python_home+'/MDS/'
mds_tci       = python_home+'/MDS/tci.py'
mds_ref       = python_home+'/MDS/reflec.py'
mds_over      = python_home+'/MDS/mds_overview.py'
mds_ts        = python_home+'/MDS/ts5.py'
mds_ces       = python_home+'/MDS/ces5.py'
mds_lit       = python_home+'/MDS/little_gui5.py'
mds_da        = python_home+'/MDS/plot_da.py'
mse_corr      = ''
mse_dir       = ''
#fittings
gfitp_dir     = python_home+'/GFIT/gfitp.py'
gfit_dir      = python_home+'/GFIT/guifit.py'
dummy_dir     = python_home+'/GFIT/TS_NE_dummy.dat'
rdena_db_dir  = '/home/ksk911211/DENA/DBs'
#efit
efit_dir      = python_home+'/EFIT'
kindata_dir   = python_home+'/EFIT/kindata'
efit_rmp      = python_home+'/GEFIT/efit_rmp.py'
mpraw_file    = python_home+'/EFIT/rmpcomp_for_rtefit.dat'
mpraw_dat     = python_home+'/EFIT/mp_coils_comp.pickle'
#dena
dena_dir      = python_home+'/DENA/dena.py'

##---etc (author, version)
author = dict()
author['eped']        = '||             developed by S.K.Kim & PLARE(SNU)               ||'
author['gfit']        = '||             developed by S.K.Kim & PLARE(SNU)               ||'
author['chease']      = '||      developed by S.K.Kim, C.Y.Lee, B.S.Kim & PLARE(SNU)    ||'
author['jatool']      = '||         developed by S.K.Kim, C.B.Lim & PLARE(SNU)          ||'
author['pedscanner']  = '||             developed by S.K.Kim & PLARE(SNU)               ||'
author['fgefit']      = '||               Developed by S.K.Kim & PLARE(SNU)             ||'
author['gfitp']       = 'Developed by SNU (S.K.Kim)'
author['gefit']       = 'Developed by SNU (S.K.Kim,Y.Lee,C.Lee,B.Kim and Y.S.Na)'
author['gefit2']      = 'Supported by NFRI(H.S.Kim and L.Terzolo)'
author['jatool2']     = 'Developed by S.K.Kim & PLARE(SNU)'
author['eped2']       = 'Developed by S.K.Kim & PLARE(SNU)'
author['bs2kstar']    = 'Developed by S.K.Kim & PLARE(SNU)'

comment = dict()
comment['gefit']      = 'Bug report: sk42@princeton.edu'
comment['jatool']     = 'Bug report: sk42@princeton.edu'
comment['gfit']       = '||                Bug report: sk42@princeton.edu               ||'
comment['eped']       = 'Bug report: sk42@princeton.edu'
comment['fgefit']     = '||                Bug report: sk42@princeton.edu               ||'

version = dict()
version['eped']       = '1.0'
version['gfit']       = '3.3' 
version['chease']     = '2.2' 
version['jatool']     = '1.1' 
version['pedscanner'] = '1.0'
version['gfitp']      = '2.0' 
version['gefit']      = '3.1' 
version['bs2kstar']   = '1.0'
version['fgefit']     = '1.0'
version['dena']       = '1.2' 
version['rdena']      = '1.0'
version['mds']        = '2.0'
