
# SISS
Simple Inelastic Scattering Simulator is an event generator for inelastic electron-proton scattering. It uses equivalent-radiator method to generate events for electron-proton scattering.


## Requirements

It's based on C++ and uses Root framework ( from CERN ).
To run these macros, you need to have Root installed ( can be found at https://root.cern/install/all_releases/ ).\
A C++ compiler is also needed to compile and run C++ files.


## Usage

There are three types of simulations that can be run from this project.

### 1.  Scattering at an angle

To simulate scattering at a fixed angle, run `simulateScattering` macro and pass the energy of incident elector ( in eV ) and the number of events to be generated. The value of the angle ( in radians ) must be set in `constants.C` file as `EXACT_DETECTOR_ANGLE`

For example, to generate 1000000 events of scattering for electron incident with 2.5 GeV energy, one can run in shell:

``` $ root -l 'simulateScattering.C(2.5e9,1000000)' ```

### 2.  Scattering in an angular range

To simulate scattering between two angles, run `simulateScatteringRange` macro and pass the energy of incident elector ( in eV ) and the number of events to be generated. The lower and the upper values of the angle ( in radians ) can be set in `constants.C` file as `MINIMUM_DETECTOR_ANGLE` and `MAXIMUM_DETECTOR_ANGLE` respectively.

For example, to generate 1000000 events of scattering for electron incident with 2.5 GeV energy, one can run in shell:

``` $ root -l 'simulateScatteringRange.C(2.5e9,1000000)' ```

### 3.  Scattering using Maximon's formula

To use Maximon's formula ( from [5] ), run `maximon` macro and pass the energy of incident elector ( in eV ) and the number of events to be generated. The value of the angle ( in radians ) is set in `constants.C` file as `EXACT_DETECTOR_ANGLE`

For example, to generate 1000000 events of scattering for electron incident with 2.5 GeV energy, one can run in shell:

``` $ root -l 'maximon.C(2.5e9,1000000)' ```

## Output

The output root files are saved as a file whose location is given by `DESTINATION_FILE` macro defined in `constants.C` file. The histogram of lepton energy ( E_l ) is divided into `nbins_E_l` bins with energy range from `xmin_E_l` to `xmax_E_l` MeV, which are also defined in `constants.C.
These commons abbreviations are used frequently in the plots:

E_l, theta_l = Energy, polar angle of lepton\
E_p, theta_p = Energy, polar angle of proton\
E_g = Energy of photon


## References
[1] Ent, R. ; Filippone, B. W. ; Makins, N. C.R. ; Milner, R. G. ; O'Neill, T. G. ; Wasson, D. A. / Radiative corrections for (e,e′p) reactions at GeV energies. In: Physical Review C - Nuclear Physics. 2001 ; Vol. 64, No. 5. pp. 546101-5461023.

[2]  R. A. Montgomery,  J. R. M. Annand,  D. Dutta, C. E. Keppel,  P. King,B. Wojtsekhowski, and J. Zhang.  Proposed measurement of tagged deepinelastic scattering in Hall A of Jefferson lab.AIP  Conf.  Proc.

[3] Povh, B. (2015). Particles and nuclei: An introduction to the physical concepts. Heidelberg, Germany: Springer.

[4] L.M. Mo and Y.S. Tsai, Rev. Mod. Phys. 41, 205 (1969).

[5]  L. C. Maximon and J. A. Tjon.  Radiative corrections to electron-protonscattering.Physical Review C, 62(5), Oct 2000.  ISSN 1089-490X.  doi:  10.1103/physrevc.62.054320

[6] A.V. Gramolin, V.S. Fadin, A.L. Feldman, et al. A new event generator for the
elastic scattering of charged leptons on protons. arXiv:1401.2959.
