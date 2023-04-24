# For this project, the following software versions were used:
Julia==v1.8
FFTW==v1.6.0
HDF5==v0.16.14
Optim==v1.7.4
ProgressMeter==v1.7.2
ImageFiltering==v0.7.3

# Externally supplied data
1. ../Share/ABidentities_blobel2021_chr2_35Mb_60Mb.npy contains AB identities derived from the data in "Zhang, H., Lam, J., Zhang, D. et al. CTCF and transcription influence chromatin structure re-configuration after mitosis. Nat Commun 12, 5157 (2021). https://doi.org/10.1038/s41467-021-25418-5".
2. ../Share contains simulation data that we analyze

# Scripts
1. get_mechanical_properties.jl fits the mechanical properties of chains.
2. get_activity_profiles-reference.jl runs costly precomputations for getting the activity profiles in a way that makes the results comparable.
3. analyze_activity_profiles.jl fits and analyzes the activity profiles.
4. predict_msd_discrete_polymer_dynamics.jl predicts the subdiffusion dynamics of loci on a heterogeneous active chain.