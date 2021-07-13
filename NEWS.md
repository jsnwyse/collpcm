# collpcm 1.2

## Minor changes

- Added checks for bad values of latent space dimension argument in `collpcm.fit()`
- Removed erroneous reshaping of the latent position array samples in `collpcm.fit()`; this code had no impact as its effect was undone in a later reshaping step
- `plot` function for `collpcm` class changed to automatically include pie charts of membership

## Larger changes

- Rewrite of the `collpcm.control()` function for better clarity 
- `collpcm.fit()` now has a native implementation of Procrustes matching which was previously sourced from `vegan` package
- A vignette has been added to this version


