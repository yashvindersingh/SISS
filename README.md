# SISS
Simple Inelastic Scattering Simulator ( ... under progress ... )

This is an event generator for inelastic electron-proton scattering. It uses modified-equivalent radiator method to generate events for scattering.
It's based on C,C++ and uses root classes.

To run a simulation, run simulateScattering macro and pass the energy of incident elector ( in eV ) and the number of events to be generated.
For example, to generate 1000000 events of scattering for electron incident with 2.5 GeV energy, one can run in shell:

``` $ root -l 'simulateScattering.C(2.5e9,1000000)' ```

The output root files are saved as a file whose location is given by DESTINATION_FILE macro defined in constants.C file.
These commons abbreviations are used frequently in the code:
E_l, theta_l = Energy, polar angle of lepton
E_p, theta_p = Energy, polar angle of proton
E_g = Energy of photon

Generation of 10 million events takes around 10-15 seconds on a standard home computer.

## References
[1] Povh, B. (2015). Particles and nuclei: An introduction to the physical concepts. Heidelberg, Germany: Springer.

[2] A.V. Gramolin, V.S. Fadin, A.L. Feldman, et al. A new event generator for the
elastic scattering of charged leptons on protons. arXiv:1401.2959.

[3] Ent, R. ; Filippone, B. W. ; Makins, N. C.R. ; Milner, R. G. ; O'Neill, T. G. ; Wasson, D. A. / Radiative corrections for (e,e′p) reactions at GeV energies. In: Physical Review C - Nuclear Physics. 2001 ; Vol. 64, No. 5. pp. 546101-5461023.

[4] L.M. Mo and Y.S. Tsai, Rev. Mod. Phys. 41, 205 (1969).
