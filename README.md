## KEFIT_PACK
Python package for KSTAR profile and equilibrium recontruction.<br />
Official Version 3.

## Installation
This python package requires following libraries.<br />
Find in github<br />
-lmfit [Copyright 2021, Matthew Newville, Till Stensitzki, Renee Otten, and others.]<br />
-csaps [MIT License Copyright (c) 2017 Eugene Prilepin]<br />
<br />
It also requires source codes.<br />
-MISHKA/ELITE[official ver]<br />
-CHEASE/HELENA[customized ver]<br />
-NUBEAM[official ver]<br />
<br />
All enviromental variables are stored in<br />
-ENV/exec_dirs.py<br />
<br />
Execute install.sh to produce links for all executables.<br />

## Usage
For KSTAR plasma, it is dedicated to reproducing Kinetic-profile/EFIT and ideal pedestal stability analysis.

## Authors and Support
SangKyeun Kim: sk42@princeton.edu<br />
HyunSeok Kim: hskim0618@kfe.re.kr<br />
BoSeong Kim: bobokim@snu.ac.kr<br />
Changyoung Lee: leecyid@snu.ac.kr<br />
SeongMoo Yang: syang@pppl.gov<br />

## Roadmap
Extending/optimizing the capabilities for KSTAR diagnostics and stability runs.

## Project status
Before KSTAR 2020 campaign (V.2)<br />
Updated for KSTAR 2021 campaign (V.3)<br />
Updated for KSTAR 2022 campaign (V.3.1)<br />
<br />
@BS2K<br />
-Ver 1.0: BS current generator<br />
<br />
@CEHASE<br />
-Ver 1.0: Initial python chease wrapper with current models.<br />
-Ver 2.0: NUBEAM routine is added.<br />
-Ver 2.1: Reproduce/output neoclassical variables<br />
-Ver 2.2: 1D zeff and carbon density profiles<br />
<br />
@DENA<br />
-Ver 1.0: Density fit using interferometers <br />
-Ver 1.1: General EFIT schemes are launched <br />
-Ver 1.2: Improved scaler/ and pre-run reading subroutines<br />
<br />
@EPED<br />
-Ver 1.0: Gui-based EPED<br />
<br />
@ERGEN<br />
-Ver 1.0: Er generator<br />
<br />
@FGEFIT<br />
-Ver 1.0: Forced brutal kEFIT tool<br />
<br />
@GEFIT<br />
-Ver 1.0: Gui-based EFIT tool<br />
-Ver 2.0: Synthetic MSE + RMP compensation added.<br />
-Ver 3.0: Improved MDS loader + general EFIT schemes<br />
-Ver 3.1: Pfile generator<br />
-Ver 3.2: Add Update popup <br />
<br />
@GFIT<br />
-Ver 1.0: Single channel based fitting tool<br />
-Ver 2.0: Extended fitting functions<br />
-Ver 3.0: Multi channel based fitting<br />
-Ver 3.1: DENA merged<br />
-Ver 3.2: New smooth spline scheme using CSPAS<br />
-Ver 3.3: TS CORE/EDGE calibration factor<br />
<br />
@JASTAB<br />
-Ver 1.0: Gui-based ideal pedestal stability tool<br />
-Ver 1.1: Jedge is added.<br />
<br />
@PEDSCAN<br />
-Ver 1.0: Gui-based pedestal modifier<br />
<br />
@GMDS<br />
-Ver 1.0: Gui-based MDS loader (beta, GEFIT2.0 compatible)<br />
-Ver 2.0: Gui-based MDS loader (GEFIT 3.0 compatiable), supporting EFIT01-05/CES/TS/TCI/MSE/RMP compensation
<br />

