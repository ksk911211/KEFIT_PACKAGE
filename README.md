## KEFIT_PACK
Official Version 3.1

## Installation
This python package requires following libraries.
Find in github
-lmfit
-csaps

It also requires source codes.
-MISHKA/ELITE[official ver]
-CHEASE/HELENA[customized ver]
-NUBEAM[official ver]

All enviromental variables are stored in
-ENV/exec_dirs.py

Execute install.sh to produce links for all executables.

## Usage
For KSTAR plasma, dedicated to reproducing Kinetic-profile/EFIT and ideal pedestal stability analysis.

## Authors and Support
SangKyeun Kim: sk42@princeton.edu
HyunSeok Kim: hskim0618@kfe.re.kr
BoSeong Kim: bobokim@snu.ac.kr
Changyoung Lee: leecyid@snu.ac.kr
SeongMoo Yang: syang@pppl.gov

## Roadmap
Extending/optimizing the capabilities for KSTAR diagnostics and stability runs.

## Project status
Updated from KSTAR 2022 campaign

#BS2K
-Ver 1.0: BS current generator

#CEHASE
-Ver 1.0: Initial python chease wrapper with current models.
-Ver 2.0: NUBEAM routine is added.
-Ver 2.1: Reproduce/output neoclassical variables

#DENA
-Ver 1.0: Density fit using interferometers 
-Ver 1.1: General EFIT schemes are launched 
-Ver 1.2: Improved scaler/ and pre-run reading subroutines

#EPED
-Ver 1.0: Gui-based EPED

#ERGEN
-Ver 1.0: Er generator

#FGEFIT
-Ver 1.0: Forced brutal kEFIT tool

#GEFIT
-Ver 1.0: Gui-based EFIT tool
-Ver 2.0: Synthetic MSE + RMP compensation added.
-Ver 3.0: Improved MDS loader + general EFIT schemes
-Ver 3.1: Pfile generator

#GFIT
-Ver 1.0: Single channel based fitting tool
-Ver 2.0: Extended fitting functions
-Ver 3.0: Multi channel based fitting
-Ver 3.1: DENA merged
-Ver 3.2: New smooth spline scheme using CSPAS

#JASTAB
-Ver 1.0: Gui-based ideal pedestal stability tool
-Ver 1.1: Jedge is added.

#PEDSCAN
-Ver 1.0: Gui-based pedestal modifier
