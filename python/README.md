# Fabrica for MTEX for Python 1.0.0

`fabrica.py` ports the field-direction and specimen-to-geographic rotation parts
of the MATLAB Fabrica scripts. It was exercised with the locally installed
MTEX Python package `mtex==1.0.0`. This package has its own version number,
separate from the MATLAB MTEX release number.

Run from this directory, or add it to `PYTHONPATH`:

```python
from mtex import plottingConvention, xvector, zvector
from fabrica import SDT2or, SDP2or, TP2vec, V2TP, makeGeoSample

# +X east, +Y north, +Z out of the screen.
pc = plottingConvention(zvector, xvector)

fab_or, fab_rot, strike_v, pole_v, line_v = SDT2or(
    strike=0, dip=60, trend=90, convention=pc
)
print(V2TP(pole_v, pc))  # approximately (270, 30)
print(V2TP(line_v, pc))  # approximately (90, 60)

# When rake is measured instead of a lineation trend:
fab_or, fab_rot, strike_v, pole_v, line_v = SDP2or(0, 60, 90, pc)

# If an EBSD map is present, name a phase to retain its rotated orientations.
# sample = makeGeoSample(0, 60, 90, ebsd, "phase name", convention=pc)
# rotated_map = sample["geo"]["rotebsd"][0]
```

Angles are degrees. Strike follows the right-hand rule; positive plunge is
downward. `TP2vec` accepts signed plunge. `V2TP` returns the downward axial
representative, so an upward vector may come back with the opposite trend.
Trend is `NaN` for a vertical vector. An up-dip trend supplied to `SDT2or`
also becomes its opposite down-going lineation. For a vertical foliation, a
lineation trend must be parallel to strike. `makeGeoSample` retains full
precision for solved plunge and reports rake on 0–180°, so an obtuse angle
from strike is not folded into an acute one. In its `geo` dictionary, `trend`
keeps the measured input while `lineTrend` records the plotted downward
lineation. For strike 0°, dip 60°, and measured trend 270°, `lineTrend` is 90°.

All functions accept an explicit MTEX `plottingConvention`. The default in this
Python file is the original Fabrica identity convention (+X east, +Y north,
+Z out). MATLAB Fabrica now follows `plottingConvention.default`, which is a
session setting. To match a MATLAB session, pass its convention explicitly in
Python. For example, `plottingConvention(-zvector, xvector)` gives +X east,
-Y north, and +Z into the screen. The direction functions use the convention's
east/north/into-screen axes, so all four X cardinal alignments work with either
Z screen sense. `SDT2or` and `SDP2or` preserve the MATLAB specimen defaults:
lineation +X and foliation pole -Y when +Z is out, with both signs reversed
when +Z is into the screen.

`geoFab(strike, dip, trend, specimen_lineation, specimen_pole, convention)`
accepts explicit perpendicular specimen axes. `FandL2Or(dip_direction, dip,
rake, convention)` retains its separate historical specimen frame, +Z for the
foliation pole and +X for lineation. Its field vectors use the same corrected
right-hand-rule geometry as `SDP2or`.

`rotTex2Geo` left-composes a fabric rotation with MTEX orientations.
`rotEBSD2Geo` rotates EBSD orientations with `keepXY=True`, preserving map
positions. `rotODFGeo` rotates an MTEX `SO3Fun` on its specimen side. A plain
MTEX `rotation` or Fabrica `orientation` can be passed to these helpers;
`rotODFGeo` also accepts a dictionary returned by `makeGeoSample`. The source
map and orientations are not mutated. `makeGeoSample` returns nested
`sample["spec"]` and `sample["geo"]` dictionaries rather than MATLAB structs.
The convention is retained at `sample["geo"]["how2plot"]` and
`sample["spec"]["how2plot"]`. It can be given as the `convention` keyword
or as a MATLAB-style positional optional argument.

For pole figures and direction plots, use MTEX Python's native `plot`,
`annotate`, and `plotPF` with the same `how2plot=pc` and lower hemisphere.
For example, after the `SDT2or` call above:

```python
from mtex import annotate, hold, newMtexFigure, plot

fig = newMtexFigure()
plot(pole_v, plane=True, lower=True, how2plot=pc)
hold(True)
annotate(pole_v, lower=True, how2plot=pc,
         marker="s", markerSize=12, markerFaceColor="k")
annotate(line_v, lower=True, how2plot=pc,
         marker="o", markerSize=12, markerFaceColor="k")
hold(False)
```

The lower-hemisphere great-circle trace is plotted without `antipodal=True`;
in the installed MTEX Python version that flag also draws a closing chord.

The [MATLAB MTEX 6.1.1 edition](../matlab/mtex61/README.md) includes separate
plotting examples around the same field geometry.

To make the geographic invariance visible, run
`plot_geographic_convention_proof.py`. Its default is the **synthetic**
`SYNTHETIC DOWN-DIP` row in this directory's `field_fabric_synthetic.csv`: right-hand-rule
strike 30°, dip 60° SE, lineation trend 120°, plunge 60°. The script writes a
PNG, PDF, and JSON check report in this `python/` directory. The eight lower-hemisphere
panels cover +Z/-Z out of screen with +X/+Y/-X/-Y east. Each panel plots the
same geographic foliation great circle, foliation pole square, and plunging
lineation circle. All eight specimen line/pole axis pairs are illustrative
perpendicular alignments. The plot uses MTEX Python vectors and rotations
with a direct equal-area Schmidt projection, matching MTEX's default spherical
projection. It checks the supplied plunge against the plane-constrained solve,
the rotation of both specimen axes, lineation-on-plane geometry, and identical
geographic screen coordinates across all eight panels.

Activate an environment with MTEX Python, NumPy, and Matplotlib installed,
then run these commands from the repository root:

```sh
python python/plot_geographic_convention_proof.py
```

To use a field row in your own CSV, pass `--csv my_field_data.csv`,
`--sample MY_ROW_NAME`, and `--output /tmp/my_geographic_proof.png`. The explicit
output prevents a measured field plot from replacing the bundled synthetic
figure. The CSV needs `name,strike,dip,dipDir,trend,plunge` columns, with
azimuths in degrees clockwise from north and positive plunge downward. `strike`
is right-hand-rule. External field data are optional and are not distributed
with Fabrica.

Run the behavioral and spherical-grid matrix checks in the installed MTEX
Python environment:

```sh
python -m unittest discover -s python -p 'test_*.py' -v
```

`test_fabrica_matrix.py` uses 54 spherical specimen line directions, six
perpendicular pole angles per line, 180 valid field strike/dip/rake cases,
and eight plotting conventions (466,560 mappings). An independent numeric
east/north/down oracle checks MTEX directions, source-axis rotation,
pole/plane orthogonality, trend/plunge, native MTEX equal-area projection,
and coincidence of geographic screen coordinates. Thirty exactly vertical
foliation cases with an oblique rake are recorded as skipped because trend
alone cannot recover their plunge. A deliberately wrong plotting convention
must fail the direction metric. Run the file directly to save the detailed
worst-case report:

```sh
python python/test_fabrica_matrix.py
```

The MATLAB and Python matrix reports include complete log-scale histograms
and quantiles for the mapped specimen line/pole direction residual and the
geographic equal-area XY residual. Both projection residuals use a unit-radius
horizon. After generating both reports, render the two-example geographic
figure and standalone residual-density plots with:

```sh
python python/plot_geographic_matrix_figure.py
```

The figure's first stereonet uses the synthetic CSV row above by default;
the second is a separate synthetic oblique grid case. Supply `--csv`,
`--sample`, and `--output` to replace the first case and keep the bundled
figure synthetic. The script requires both matrix reports, generated by the
MATLAB and Python matrix tests.
