
# Representating spectra

There are currently two core types for representing spectra:

* `SampledSpectrum` represents a continuous spectrum sampled at even intervals in some
  wavelength range.

* `SingleLine` represents a single very narrow line, essentially a $\delta$-distribution
  that is zero everywhere, but has a finite integral over wavelength.

`SampledSpectrum` should be used for all surface materials, as well as light sources with
continuous spectra. `SingleLine` can be used for narrow light sources like lasers. 

The advantage of using `SingleLine` is that it makes certain computations faster. When a
`SingleLine` light source is reflected by a `SampledSpectrum` surface, the result is a
`SingleLine`. Computing the effects of `SingleLine` light (energy reaching a camera, or
radiation pressure) is much faster than for a `SampledSpectrum` which must be integrated
and possibly convolved with some kind of spectral response function.


## SampledSpectrum

```@docs 
SampledSpectrum
```

Note that `SampledSpectrum` uses integers for the wavelength range.

When used for most light sources, the physical units of the `values` array are assumed to
be W / m$^{-2}$ sr$^{-1}$ nm$^{-1}$. The exception is parallel light, which has an
infinitely narrow angular distribution and therefore units of W / m$^{-2}$ nm$^{-1}$.

For surfaces, the values represent reflectivity (wavelength-dependent Bond albedo) and
should have values between 0 and 1.


## SingleLine

```@docs 
SingleLine
```

When used for most light sources, the physical units of the `value` array are assumed to be
W / m$^{-2}$ sr$^{-1}$.
