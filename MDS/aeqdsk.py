#!/usr/local/miniconda3/bin/python3
import os,sys
import numpy as np
import fortranformat

f1040 = fortranformat.FortranRecordReader('1x,4e16.9')
f1041 = fortranformat.FortranRecordReader('1x,4i5')
f1060 = fortranformat.FortranRecordReader('A1,f8.3,9x,i5,11x,i5,1x,a3,1x,i3,1x,i3,1x,a3,1x,2i5')

def read_f1040(input_str):
    try:
        return f1040.read(input_str)
    except Exception:
        return [0.0] * 4


def _read_afile(filename):

    adat = dict()
    with open(filename, 'r') as f:
        AFILE = iter(f.readlines())

    adat['__header__'] = ''
    k = 0
    for line in AFILE:
        if line[0] == '*':
            break
        else:
            adat['__header__'] += line
            if k == 1:
                adat['shot'] = int(line[:7])
            k += 1
    try:

        (
            dummy,
            adat['time'],
            adat['jflag'],
            adat['lflag'],
            adat['limloc'],
            adat['mco2v'],
            adat['mco2r'],
            adat['qmflag'],
            adat['nlold'],
            adat['nlnew'],
        ) = f1060.read(line)

    except Exception:

        (
            adat['jflag'],
            adat['lflag'],
            adat['limloc'],
            adat['mco2v'],
            adat['mco2r'],
            adat['qmflag'],
            adat['nlold'],
            adat['nlnew'],
        ) = line.split()[1:]
        print('bad time in aEQDSK file: %s' % (filename))
        adat['time'] = 0
        for k in ['jflag', 'lflag', 'mco2v', 'mco2r', 'nlold', 'nlnew']:
            adat[k] = int(adat[k])

    for k in ['nlold', 'nlnew']:
        if adat[k] is None:
            adat[k] = 0

    adat['rseps'] = np.zeros(2)
    adat['zseps'] = np.zeros(2)

    # fmt: off
    adat['tsaisq'],   adat['rcencm'],   adat['bcentr'],   adat['pasmat']   = read_f1040(next(AFILE))
    adat['cpasma'],   adat['rout'],     adat['zout'],     adat['aout']     = read_f1040(next(AFILE))
    adat['eout'],     adat['doutu'],    adat['doutl'],    adat['vout']     = read_f1040(next(AFILE))
    adat['rcurrt'],   adat['zcurrt'],   adat['qsta'],     adat['betat']    = read_f1040(next(AFILE))
    adat['betap'],    adat['ali'],      adat['oleft'],    adat['oright']   = read_f1040(next(AFILE))
    adat['otop'],     adat['obott'],    adat['qpsib'],    adat['vertn']    = read_f1040(next(AFILE))
    # fmt: on

    for arr, adim in (('rco2v', 'mco2v'), ('dco2v', 'mco2v'), ('rco2r', 'mco2r'), ('dco2r', 'mco2r')):
        tmp = []
        for k in range(int(np.ceil(adat[adim] / 4.0))):
            tmp.extend(read_f1040(next(AFILE)))
        adat[arr] = tmp[: adat[adim]]

    # fmt: off
    adat['shearb'],   adat['bpolav'],   adat['s1'],       adat['s2']       = read_f1040(next(AFILE))
    adat['s3'],       adat['qout'],     adat['olefs'],    adat['orighs']   = read_f1040(next(AFILE))
    adat['otops'],    adat['sibdry'],   adat['areao'],    adat['wplasm']   = read_f1040(next(AFILE))
    adat['terror'],   adat['elongm'],   adat['qqmagx'],   adat['cdflux']   = read_f1040(next(AFILE))
    adat['alpha'],    adat['rttt'],     adat['psiref'],   adat['xndnt']    = read_f1040(next(AFILE))
    adat['rseps'][0], adat['zseps'][0], adat['rseps'][1], adat['zseps'][1] = read_f1040(next(AFILE))
    adat['sepexp'], adat['obots'], adat['btaxp'], adat['btaxv'] = read_f1040(next(AFILE))
    adat['aaq1'], adat['aaq2'], adat['aaq3'], adat['seplim'] = read_f1040(next(AFILE))
    adat['rmagx'], adat['zmagx'], adat['simagx'], adat['taumhd'] = read_f1040(next(AFILE))
    adat['betapd'], adat['betatd'], adat['wplasmd'], adat['diamag'] = read_f1040(next(AFILE))
    adat['vloopt'], adat['taudia'], adat['qmerci'], adat['tavem'] = read_f1040(next(AFILE))

    adat['nsilop0'],  adat['magpri0'],  adat['nfcoil0'],  adat['nesum0']   = f1041.read(next(AFILE))
    # fmt: on

    tmp = []
    for k in range(int(np.ceil((adat['nsilop0'] + adat['magpri0']) / 4.0))):
        tmp.extend(read_f1040(next(AFILE)))
    adat['csilop'] = tmp[: adat['nsilop0']]
    adat['cmpr2'] = tmp[adat['nsilop0'] : (adat['nsilop0'] + adat['magpri0'])]

    adat['ccbrsp'] = []
    for k in range(int(np.ceil((adat['nfcoil0']) / 4.0))):
        adat['ccbrsp'].extend(read_f1040(next(AFILE)))

    adat['eccurt'] = []
    for k in range(int(np.ceil((adat['nesum0']) / 4.0))):
        adat['eccurt'].extend(read_f1040(next(AFILE)))

    try:
        # fmt: off
        # Machine specific ?
        adat['pbinj'],  adat['rvsin'],   adat['zvsin'],   adat['rvsout']  = read_f1040(next(AFILE))
        adat['zvsout'], adat['vsurfa'],  adat['wpdot'],   adat['wbdot']   = read_f1040(next(AFILE))
        adat['slantu'], adat['slantl'],  adat['zuperts'], adat['chipre']  = read_f1040(next(AFILE))
        adat['cjor95'], adat['pp95'],    adat['ssep'],    adat['yyy2']    = read_f1040(next(AFILE))
        adat['xnnc'],   adat['cprof'],   adat['oring'],   adat['cjor0']   = read_f1040(next(AFILE))
        adat['fexpan'], adat['qqmin'],   adat['chigamt'], adat['ssi01']   = read_f1040(next(AFILE))
        adat['fexpvs'], adat['sepnose'], adat['ssi95'],   adat['rqqmin']  = read_f1040(next(AFILE))
        adat['cjor99'], adat['cj1ave'],  adat['rmidin'],  adat['rmidout'] = read_f1040(next(AFILE))
        adat['psurfa'], adat['peak'],    adat['dminux'],  adat['dminlx']  = read_f1040(next(AFILE))
        adat['dolubaf'],adat['dolubafm'],adat['diludom'], adat['diludomm']= read_f1040(next(AFILE))
        adat['ratsol'], adat['rvsiu'],   adat['zvsiu'],   adat['rvsid']   = read_f1040(next(AFILE))
        adat['zvsid'],  adat['rvsou'],   adat['zvsou'],   adat['rvsod']   = read_f1040(next(AFILE))
        adat['zvsod'],  adat['condno'],  adat['psin32'],  adat['psin21']  = read_f1040(next(AFILE))
        adat['rq32in'], adat['rq21top'], adat['chilibt'], adat['ali3']    = read_f1040(next(AFILE))
        adat['xbetapr'],adat['tflux'],   adat['tchimls'], adat['twagap']  = read_f1040(next(AFILE))
        # fmt: on
    except StopIteration as _excp:
        pass

    # anything extra go in the footer
    adat['__footer__'] = ''
    try:
        for line in AFILE:
            adat['__footer__'] += line
    except Exception:
        pass

    # add betaN calculation to a-file
    # if it is not a vacuum shot
    if adat['cpasma'] != 0.0:
        i = adat['cpasma'] / 1e6
        a = adat['aout'] / 100.0
        bt = adat['bcentr'] * adat['rcencm'] / adat['rout']
        i_n = i / a / bt
        adat['betan'] = abs(adat['betat'] / i_n)

    # lists into arrays
    for var in adat:
        if isinstance(adat[var], list):
            adat[var] = np.array([_f for _f in adat[var] if _f is not None])

    # remove NaN from aEQDSK file to allow saving
    for k in adat:
        if isinstance(adat[k], np.ndarray) and np.any(np.isnan(adat[k])):
            adat[k][np.isnan(adat[k])] = 0
            printe('%s array is NaN in aEQDSK file: %s' % (k, filename))
        elif isinstance(adat[k],float) and np.any(np.isnan(adat[k])):
            adat[k] = 0
            printe('%s entry is NaN in aEQDSK file: %s' % (k, filename))

    return adat

if __name__ == "__main__":

    import aeqdsk

    nfile = sys.argv[1]
    if not os.path.isfile(nfile): print('>>> No input file'); exit()
    adat = _read_afile(nfile)
    for k in adat.keys():
      print('>>> %s %s'%(k,adat[k]))
