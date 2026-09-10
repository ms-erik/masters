# Five-pair strict electronic ARROW

This package increases confinement using five thin high-k InGaAs
antiresonant films on each side of a lower-k InGaAlAs central guide.

The active layers remain propagating at 150 meV. Increasing the potential
to the full 510 meV InAlAs offset would make the layers evanescent and turn
the device back into a resonant-tunneling structure.

## Run

```matlab
run_multipair_electronic_arrow_demo
```

## Nominal material model

- high-k film/access: In0.53Ga0.47As, V = 0 meV
- buffer: In0.53Al0.08Ga0.39As, design V about 85 meV
- guide: In0.53Al0.10Ga0.37As, design V about 106 meV
- target energy: 150 meV
- five high-k films per side

The InGaAlAs band edges and masses are currently interpolation-based.
Before fixing an MBE recipe, use a calibrated quaternary band-lineup model,
HRXRD composition measurements, and test quantum wells.

## Fabrication interpretation

More reflector films reduce leakage but increase:

- total stack thickness;
- accumulated thickness error;
- sensitivity to interface roughness;
- number of composition transitions;
- risk of alloy disorder and phase detuning.

This package is therefore a trade-off study, not an assertion that more
periods are always better.

Reference for improved InGaAlAs/InP band interpolation:

Y. S. Yong et al., Applied Physics B 99, 477-486 (2010),
DOI: 10.1007/s00340-010-4000-3.
