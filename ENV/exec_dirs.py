##---machine inform
machine = 'nkstar' #nkstar, ukstar_root, ukstar

##---Experimental setup
#efit years
shotk = {}
years = []
for i in range(2007,2027): years.append('%i'%i)
for i in years: shotk[i] = dict()
shotk['2007']['shot'] = range(  635, 1284)
shotk['2008']['shot'] = range(  635, 1284)
shotk['2009']['shot'] = range( 1284, 2457)
shotk['2010']['shot'] = range( 2457, 4469)
shotk['2011']['shot'] = range( 4469, 5569)
shotk['2012']['shot'] = range( 7232, 8355)
shotk['2013']['shot'] = range( 8355, 9428)
shotk['2014']['shot'] = range(10616,11725)
shotk['2015']['shot'] = range(13302,14408)
shotk['2016']['shot'] = range(16277,17377)
shotk['2017']['shot'] = range(18370,19397)
shotk['2018']['shot'] = range(20648,21759)
shotk['2019']['shot'] = range(21801,24082)
shotk['2020']['shot'] = range(24082,27401)
shotk['2021']['shot'] = range(27401,30446)
shotk['2022']['shot'] = range(30446,32769)
shotk['2023']['shot'] = range(32769,34837)
shotk['2024']['shot'] = range(34925,37900)
shotk['2025']['shot'] = range(37925,39900)
shotk['2026']['shot'] = range(39925,40900)
for year in years:
    for efit_no in range(1,6):
        shotk[year][efit_no] = '/EFIT_RUN/EFITDATA_%s/EFIT%02i'%(year,efit_no)
#mse favorable channels.
mse_good_ch = {}
for i in years: mse_good_ch[i] = []
for i in range(2007,2016):
    mse_good_ch['%i'%i] = [0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0]
mse_good_ch['2016'] = [1,1,1,0,1, 1,1,1,1,1, 1,1,0,0,0, 0,0,0,0,0, 0,0,0,0,0]
mse_good_ch['2017'] = [0,1,1,1,1, 1,1,1,1,0, 0,1,1,1,0, 0,0,1,1,0, 0,0,0,0,0]
mse_good_ch['2018'] = [0,0,1,1,1, 1,1,1,1,1, 1,1,1,1,0, 0,0,1,0,0, 0,0,0,0,0]
mse_good_ch['2019'] = [0,0,1,1,1, 0,1,0,0,0, 1,1,1,0,1, 0,0,1,1,0, 0,0,0,0,0]
mse_good_ch['2020'] = [1,1,1,0,1, 1,1,1,1,1, 1,1,0,0,0, 0,0,1,0,0, 0,0,0,0,0]
mse_good_ch['2021'] = [1,1,1,0,1, 1,1,1,1,1, 1,1,0,0,0, 0,0,1,0,0, 0,0,0,0,0]
mse_good_ch['2022'] = [1,1,1,0,1, 1,1,1,1,1, 1,1,0,0,0, 0,0,1,0,0, 0,0,0,0,0]
mse_good_ch['2023'] = [0,1,1,1,1, 0,1,1,1,1, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0]
mse_good_ch['2024'] = [0,0,1,1,1, 0,1,1,1,1, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0]
mse_good_ch['2025'] = [0,0,1,1,1, 0,1,1,1,1, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0]
mse_good_ch['2026'] = [0,0,1,1,1, 0,1,1,1,1, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0]

##--Diagnostic exception
ts_location_core = {}
ts_location_edge = {}
ts_location_core[2019] = [1806.0,1826.0,1848.0,1871.0,1894.0,1917.0,1942.0,1966.0,1991.0,2016.0,2041.0,2068.0,2093.0,2120.0]
ts_location_edge[2019] = [2124.0,2137.0,2143.0,2149.0,2156.0,2162.0,2177.0,2191.0,2202.0,2216.0,2229.0,2242.0,2257.0,2271.0,2285.0,2297.0,2311.0]
ts_location_core[2024] = [1790.6,1813.6,1836.1,1857.8,1881.9,1905.8,1927.6,1952.4,1976.9,2001.8,2027.7,2053.3,2080.0,2105.6]
ts_location_edge[2024] = [2115.6,2124.3,2135.1,2143.8,2154.0,2162.1,2170.9,2179.7,2190.1,2199.7,2212.4,2221.3,2231.1,2252.2,2273.5,2295.7,2305.7]
ts_location_core[2025] = [1785.2,1806.7,1828.2,1850.8,1873.3,1896.2,1919.6,1944.0,1967.6,1992.8,2017.8,2043.2,2069.6,2094.9]
ts_location_edge[2025] = [2096.3,2104.8,2113.4,2122.7,2131.3,2140.7,2149.4,2158.1,2167.6,2177.9,2197.0,2217.0,2237.2,2258.3,2279.5]


ces_location = {}
ces_location[2011] = [1.795,1.800,1.850,1.900,1.950,2.000,2.050,2.100,2.140,2.160,2.170,2.180,2.190,2.200,2.205,2.210,2.215,2.220,2.225,2.230,2.235,2.240,2.245,2.250,2.255,2.265,2.275,2.280,2.285,2.290,2.295,2.300]
ces_location[2012] = [1.800,1.850,1.900,1.950,2.000,2.050,2.100,2.150,2.170,2.180,2.190,2.200,2.205,2.210,2.215,2.220,2.225,2.230,2.235,2.240,2.245,2.250,2.255,2.260,2.265,2.270,2.275,2.280,2.285,2.290,2.295,2.300]
ces_location[2013] = [1.795,1.850,1.900,1.950,2.000,2.050,2.100,2.150,2.170,2.180,2.190,2.200,2.205,2.210,2.215,2.220,2.225,2.230,2.235,2.240,2.245,2.250,2.255,2.260,2.265,2.270,2.275,2.280,2.285,2.290,2.295,2.300]
ces_location[2014] = [1.795,1.850,1.900,1.950,2.000,2.050,2.100,2.150,2.170,2.180,2.190,2.200,2.205,2.210,2.215,2.220,2.225,2.230,2.235,2.240,2.245,2.250,2.255,2.260,2.265,2.270,2.275,2.280,2.285,2.290,2.295,2.300]
ces_location[2015] = [1.801,1.822,1.843,1.874,1.895,1.945,1.995,2.016,2.047,2.078,2.099,2.125,2.150,2.171,2.192,2.203,2.213,2.223,2.228,2.233,2.238,2.243,2.248,2.253,2.259,2.264,2.269,2.273,2.280,2.286,2.291,2.296]

if 'ukstar' in machine:
    python_home     = '/usr/local/analysis/KEFIT/'
    mds_address     = 'nkstar.kstar.kfe.re.kr:8005' #mds_address
    efit_source_dir = '/EFIT/'                      #efit_source 
    scheduler       = 'pbs'                         #scheduler
    qsub_exec       = '/opt/pbs/bin/qsub'           
    qstat_exec      = '/opt/pbs/bin/qstat'          
    qdel_exec       = '/opt/pbs/bin/qdel'
    node_machine    = machine                       #node list option
    node_force      = True
    node_init       = 'workq'
    node_default    = 'workq'
    chease_exec     = '/usr/local/analysis/EQUIL/CHEASE/bin/chease'
    #helena
    helena_exec     = '/usr/local/analysis/EQUIL/HELENA/bin/hel13'
    #mishka
    mis_exec        = '/usr/local/analysis/STAB/MISHKA/bin/mishka1fast_'
    #elite
    elite_exec      = '/usr/local/analysis/STAB/ELITE/bin/elite'
    elite_dir       = '/usr/local/analysis/STAB/ELITE/bin'
    #gzip
    gzip_dir        = 'gzip'
    #python
    python2_exec    = '/usr/bin/python2'
    python3_exec    = '/usr/local/analysis/Python_env/gefit_env/miniconda3-py38/bin/python3'
    ##---Default DBs
    rdena_db_dir    = '/UKSTAR_HOME/ksk911211/DENA/DBs'    
    #nubeam
    plasma_state_test_exec='/usr/local/analysis/NUBEAM/ntcc/LINUX/test/plasma_state_test'
    nubeam_comp_exec      ='/usr/local/analysis/NUBEAM/ntcc/LINUX/test/mpi_nubeam_comp_exec'
    adasdir               ='/usr/local/analysis/NUBEAM/ntcc/LINUX/adas'
    preactdir             ='/usr/local/analysis/NUBEAM/ntcc/LINUX/preact'
    mpirun                ='/usr/local/mpich/bin/mpirun'

else:
    python_home     = '/home/users/ksk911211/PYTHON/KEFIT_PACKAGE/'
    mds_address     = 'nkstar.kstar.kfe.re.kr:8005'
    efit_source_dir = '/EFIT/'
    scheduler       = ''
    qsub_exec       = ''           
    qstat_exec      = ''          
    qdel_exec       = ''
    node_machine    = '' 
    node_force      = True
    node_init       = ''
    node_default    = ''
    chease_exec     = ''
    helena_exec     = ''
    #mishka
    mis_exec        = ''
    #elite
    elite_exec      = ''
    elite_dir       = ''
    #gzip
    gzip_dir        = 'gzip'
    #python
    python2_exec    = '/usr/bin/python2'
    python3_exec    = '/usr/bin/python3'
    ##---Default DBs
    rdena_db_dir    = ''
    #nubeam
    plasma_state_test_exec=''
    nubeam_comp_exec      =''
    adasdir               =''
    preactdir             =''
    mpirun                =''

if machine == 'ukstar':
    python_home     = '/UKSTAR_HOME/ksk911211/PYTHON/KEFIT_PACKAGE/'
    chease_exec     = '/UKSTAR_HOME/ksk911211/CODE/EQUIL/CHEASE/bin/chease'
    helena_exec     = '/UKSTAR_HOME/ksk911211/CODE/EQUIL/HELENA/bin/hel13'
    mis_exec        = '/UKSTAR_HOME/ksk911211/CODE/STAB/MISHKA/bin/mishka1fast_'
    elite_exec      = '/UKSTAR_HOME/ksk911211/CODE/STAB/ELITE/bin/elite'
    elite_dir       = '/UKSTAR_HOME/ksk911211/CODE/STAB/ELITE/bin'

##---python scripts
pythonc_exec  = python3_exec
#gefit
gefit_exec2   = python_home+'/GEFIT/gefit.py'
gefit_exec3   = python_home+'/GEFIT/gefit.py'
#gfit
gfit_exec2    = python_home+'/GFIT/guifit.py'
gfit_exec3    = python_home+'/GFIT/guifit.py'
#chease
chease_dir    = python_home+'/CHEASE/chease'
#nubeam
nubeam_dir    = python_home+'/NUBEAM/nubeam.py'
nubeam_dir2   = python_home+'/NUBEAM/nubeam_run.py'
nubeam_config = python_home+'/NUBEAM/configs/'
Mfile         = python_home+'/NUBEAM/mdescr_A123B123.dat'
Sfile         = python_home+'/NUBEAM/sconfig_A123B123.dat'
Ifile         = python_home+'/NUBEAM/nubeam_init.dat'
stepfile      = python_home+'/NUBEAM/nubeam_step.dat'
#infos
popup_dir     = python_home+'/ENV/popup.py'
gefit_info    = python_home+'/INFO/gefit.txt'
gfit_info     = python_home+'/INFO/gfit.txt'
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
version['gfit']       = '3.5'
version['chease']     = '2.2' 
version['jatool']     = '1.1' 
version['pedscanner'] = '1.0'
version['gfitp']      = '2.0' 
version['gefit']      = '3.3' 
version['bs2kstar']   = '1.0'
version['fgefit']     = '1.0'
version['dena']       = '1.3' 
version['rdena']      = '1.0'
version['mds']        = '2.1'
