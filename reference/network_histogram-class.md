# Network histogram class hierarchy

All fitted objects returned by
[`nethist`](https://enigmasong.github.io/nethist/reference/nethist.md),
[`multinethist`](https://enigmasong.github.io/nethist/reference/multinethist.md),
and
[`hnethist`](https://enigmasong.github.io/nethist/reference/hnethist.md)
inherit from the common base class `network_histogram`:

## Details

- `nethist` fits have class `c("nethist", "network_histogram")`.

- `multinethist` fits have class
  `c("multinethist", "network_histogram")`.

- `hnethist` fits have class
  `c("hnethist", "nethist", "network_histogram")`.

Every `network_histogram` object carries `cluster`, `thetahat`,
`rho_hat`, `normalized_LL`, `method`, and `h`. Methods that only depend
on those fields are defined once on the base class, so
[`logLik.network_histogram`](https://enigmasong.github.io/nethist/reference/logLik.network_histogram.md)
and
[`summary.network_histogram`](https://enigmasong.github.io/nethist/reference/summary.network_histogram.md)
are shared by all three. Methods whose behaviour depends on the layer
dimension, such as
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`plot3d()`](https://enigmasong.github.io/nethist/reference/plot3d.md),
and [`fitted()`](https://rdrr.io/r/stats/fitted.values.html), remain
specific to the subclasses.

## See also

[`logLik.network_histogram`](https://enigmasong.github.io/nethist/reference/logLik.network_histogram.md),
[`summary.network_histogram`](https://enigmasong.github.io/nethist/reference/summary.network_histogram.md),
[`print.nethist`](https://enigmasong.github.io/nethist/reference/print.nethist.md)
