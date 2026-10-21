"""Fabrica field directions and specimen-to-geographic rotations for MTEX Python.

Angles are degrees. Trend is clockwise from geographic north; positive plunge
points downward. A plotting convention describes how east, north, and out of the
screen lie in MTEX coordinates. The default here is the legacy Fabrica identity
convention: +X east, +Y north, +Z out of the screen. Pass a convention explicitly
when the specimen or plot uses another layout.
"""

from __future__ import annotations

import numpy as np
from mtex import EBSD, SO3Fun, crystalFrame, degree, orientation, plottingConvention, rotation, vector3d, xvector, yvector, zvector


def _convention(convention: plottingConvention | None) -> plottingConvention:
    return plottingConvention(zvector, xvector) if convention is None else plottingConvention(convention)


def TP2vec(trend, plunge, convention: plottingConvention | None = None) -> vector3d:
    """Convert geographic trend/plunge to an MTEX vector, retaining signed plunge.

    `trend` and `plunge` may be broadcastable NumPy arrays. A positive plunge is
    into the screen under the supplied geographic plotting convention.
    """
    pc = _convention(convention)
    if np.any(np.abs(np.asarray(plunge, dtype=float)) > 90):
        raise ValueError("Plunge must be between -90 and 90 degrees")
    t, p = np.broadcast_arrays(np.deg2rad(trend), np.deg2rad(plunge))
    east = np.asarray(pc.east.data)
    north = np.asarray(pc.north.data)
    down = np.asarray(pc.intoScreen.data)
    coords = (
        (np.cos(p) * np.sin(t))[..., None] * east
        + (np.cos(p) * np.cos(t))[..., None] * north
        + np.sin(p)[..., None] * down
    )
    return vector3d(coords, isNormalized=True)


def V2TP(v: vector3d, convention: plottingConvention | None = None):
    """Return (trend, plunge) in degrees for the downward axial representative.

    Upward vectors are reversed before measurement. Trend is NaN for a vertical
    vector, where it has no geometric meaning. Array shape follows `v`.
    """
    pc = _convention(convention)
    coords = np.asarray(v.data, dtype=float)
    length = np.linalg.norm(coords, axis=-1)
    if not np.all(np.isfinite(coords)) or np.any(length == 0):
        raise ValueError("V2TP needs nonzero finite vectors")
    down = np.asarray(pc.intoScreen.data)
    flip = np.sum(coords * down, axis=-1) < -1e-10 * length
    coords = np.where(flip[..., None], -coords, coords)
    east_component = np.sum(coords * np.asarray(pc.east.data), axis=-1)
    north_component = np.sum(coords * np.asarray(pc.north.data), axis=-1)
    down_component = np.sum(coords * down, axis=-1)
    east_component = np.where(np.abs(east_component) < 1e-12 * length, 0, east_component)
    north_component = np.where(np.abs(north_component) < 1e-12 * length, 0, north_component)
    down_component = np.where(np.abs(down_component) < 1e-12 * length, 0, down_component)
    horizontal = np.hypot(east_component, north_component)
    trend = np.mod(np.rad2deg(np.arctan2(east_component, north_component)), 360)
    vertical = horizontal < 1e-10 * length
    trend = np.where(vertical, np.nan, trend)
    plunge = np.rad2deg(np.arctan2(down_component, horizontal))
    plunge = np.where(vertical, 90, plunge)
    return trend, plunge


def _lower_hemisphere(v: vector3d, pc: plottingConvention) -> vector3d:
    coords = np.asarray(v.data)
    up = np.asarray(pc.outOfScreen.data)
    flip = np.sum(coords * up, axis=-1) > 1e-10 * np.linalg.norm(coords, axis=-1)
    return vector3d(np.where(flip[..., None], -coords, coords))


def _specimen_axes(pc: plottingConvention):
    zsign = float(np.asarray(pc.outOfScreen.data) @ np.asarray(zvector.data))
    if abs(abs(zsign) - 1) > 1e-8:
        raise ValueError("Fabrica specimen defaults require +Z out of or into the screen")
    return (xvector, -yvector) if zsign > 0 else (-xvector, yvector)


def _fabric_rotation(
    lineation: vector3d,
    pole: vector3d,
    specimen_lineation: vector3d,
    specimen_pole: vector3d,
):
    a = np.asarray(specimen_lineation.data)
    b = np.asarray(specimen_pole.data)
    norms = np.linalg.norm(a, axis=-1) * np.linalg.norm(b, axis=-1)
    if np.any(norms == 0) or np.any(np.abs(np.sum(a * b, axis=-1)) > 1e-8 * norms):
        raise ValueError("specimen lineation and pole must be nonzero and perpendicular")
    fab_rot = rotation.map(specimen_lineation, lineation, specimen_pole, pole)
    fab_or = orientation(fab_rot, crystalFrame("mmm"))
    return fab_or, fab_rot


def geoFab(
    strike,
    dip,
    trend,
    specimen_lineation: vector3d,
    specimen_pole: vector3d,
    convention: plottingConvention | None = None,
):
    """Map explicitly supplied specimen lineation and pole axes into the field.

    The foliation uses right-hand-rule strike/dip. The lineation lies in that
    plane and has the supplied horizontal trend; its plunge is solved from the
    plane. `specimen_lineation` and `specimen_pole` are MTEX vectors in the
    original section's frame and must be perpendicular.
    Returns `(fab_or, fab_rot, strike_v, pole_v, lineation_v)` as in MATLAB.
    """
    pc = _convention(convention)
    dip = np.asarray(dip, dtype=float)
    if np.any(~np.isfinite(dip)) or np.any((dip < 0) | (dip > 90)):
        raise ValueError("Dip must be between 0 and 90 degrees")
    strike_v = TP2vec(strike, 0, pc)
    trend_v = TP2vec(trend, 0, pc)
    pole_v = TP2vec(np.asarray(strike) - 90, 90 - dip, pc)
    pole_down = np.sum(np.asarray(pole_v.data) * np.asarray(pc.intoScreen.data), axis=-1)
    horizontal_dot = np.sum(np.asarray(trend_v.data) * np.asarray(pole_v.data), axis=-1)
    vertical_plane = np.abs(pole_down) < 1e-10
    if np.any(vertical_plane & (np.abs(horizontal_dot) > 1e-10)):
        raise ValueError("this trend is not in the vertical foliation plane")
    vertical_part = np.zeros(np.broadcast_shapes(np.shape(pole_down), np.shape(horizontal_dot)))
    np.divide(-horizontal_dot, pole_down, out=vertical_part, where=~vertical_plane)
    raw_xyz = np.asarray(trend_v.data) + vertical_part[..., None] * np.asarray(pc.intoScreen.data)
    lineation_v = _lower_hemisphere(vector3d(raw_xyz / np.linalg.norm(raw_xyz, axis=-1)[..., None]), pc)
    pole_v = _lower_hemisphere(pole_v, pc)
    fab_or, fab_rot = _fabric_rotation(lineation_v, pole_v, specimen_lineation, specimen_pole)
    return fab_or, fab_rot, strike_v, pole_v, lineation_v


def SDT2or(strike, dip, trend, convention: plottingConvention | None = None):
    """Fabric from right-hand strike, dip, and lineation trend.

    Uses specimen lineation/pole +X/-Y when +Z is out of screen, and -X/+Y
    when +Z is into screen, matching the MATLAB Fabrica defaults.
    """
    pc = _convention(convention)
    specimen_lineation, specimen_pole = _specimen_axes(pc)
    return geoFab(strike, dip, trend, specimen_lineation, specimen_pole, pc)


def SDP2or(strike, dip, pitch, convention: plottingConvention | None = None):
    """Map specimen lineation/pole using right-hand strike, dip, and pitch.

    Positive pitch rotates from strike toward down-dip about the downward pole.
    Returns `(fab_or, fab_rot, strike_v, pole_v, lineation_v)` as in MATLAB.
    """
    pc = _convention(convention)
    dip = np.asarray(dip, dtype=float)
    if np.any(~np.isfinite(dip)) or np.any((dip < 0) | (dip > 90)):
        raise ValueError("Dip must be between 0 and 90 degrees")
    strike_v = TP2vec(strike, 0, pc)
    pole_v = TP2vec(np.asarray(strike) - 90, 90 - dip, pc)
    lineation_v = rotation.byAxisAngle(pole_v, np.asarray(pitch) * degree) * strike_v
    pole_v = _lower_hemisphere(pole_v, pc)
    lineation_v = _lower_hemisphere(lineation_v, pc)
    specimen_lineation, specimen_pole = _specimen_axes(pc)
    fab_or, fab_rot = _fabric_rotation(lineation_v, pole_v, specimen_lineation, specimen_pole)
    return fab_or, fab_rot, strike_v, pole_v, lineation_v


def FandL2Or(dip_direction, dip, rake, convention: plottingConvention | None = None):
    """Fabric orientation from dip direction, dip, and right-hand rake.

    This historical entry point uses a distinct specimen frame: +Z is the
    foliation pole and +X the lineation. The geographic vectors use the same
    right-hand-rule field geometry as `SDP2or`. Returns `(fab_or, strike_v,
    pole_v, lineation_v)` as in the MATLAB function.
    """
    _, _, strike_v, pole_v, lineation_v = SDP2or(
        np.asarray(dip_direction) - 90, dip, rake, convention
    )
    fab_or, _ = _fabric_rotation(lineation_v, pole_v, xvector, zvector)
    return fab_or, strike_v, pole_v, lineation_v


def _plain_rotation(fabric_rotation: rotation) -> rotation:
    if not isinstance(fabric_rotation, rotation):
        raise TypeError("fabric rotation must be an MTEX rotation or orientation")
    return rotation(fabric_rotation)


def rotTex2Geo(o: orientation, fabric_rotation: rotation) -> orientation:
    """Left-compose an MTEX orientation with the specimen-to-field rotation."""
    if not isinstance(o, orientation):
        raise TypeError("rotTex2Geo needs MTEX orientation data")
    return _plain_rotation(fabric_rotation) * o


def rotEBSD2Geo(ebsd: EBSD, fabric_rotation: rotation) -> EBSD:
    """Rotate EBSD orientations while retaining the original pixel positions."""
    if not isinstance(ebsd, EBSD):
        raise TypeError("rotEBSD2Geo needs an MTEX EBSD map")
    return ebsd.rotate(_plain_rotation(fabric_rotation), keepXY=True)


def rotODFGeo(odf: SO3Fun, fabric_rotation_or_sample) -> SO3Fun:
    """Turn an ODF on its specimen side by a fabric rotation or sample dict."""
    if not isinstance(odf, SO3Fun):
        raise TypeError("rotODFGeo needs an MTEX SO3Fun ODF")
    source = fabric_rotation_or_sample
    if isinstance(source, dict):
        source = source["geo"]["fabRot"]
    return odf.rotate(_plain_rotation(source), side="left")


def makeGeoSample(strike, dip, trend, *items, convention: plottingConvention | None = None):
    """Collect a field fabric, optional EBSD maps, and named phase orientations.

    The MATLAB-like call is `makeGeoSample(strike, dip, trend, ebsd, 'phase')`.
    Returns nested Python dictionaries `sample['spec']` and `sample['geo']`.
    Every source map is kept under `spec['ebsd']`; `geo['rotebsd']` contains
    orientation-rotated maps with the same XY positions. Named phase results
    come from the first source map, matching the original MATLAB helper.
    """
    conventions = [item for item in items if isinstance(item, plottingConvention)]
    if len(conventions) > 1 or (conventions and convention is not None):
        raise ValueError("supply only one plotting convention")
    pc = _convention(conventions[0] if conventions else convention)
    fab_or, fab_rot, strike_v, pole_v, lineation_v = SDT2or(strike, dip, trend, pc)
    specimen_lineation, specimen_pole = _specimen_axes(pc)
    line_trend, plunge = V2TP(lineation_v, pc)
    rake = np.rad2deg(lineation_v.angle(strike_v, antipodal=False, symmetry=False))
    sample = {
        "spec": {"lin": specimen_lineation, "pole": specimen_pole, "how2plot": pc},
        "geo": {
            "strike": strike,
            "dip": dip,
            "trend": trend,
            "lineTrend": line_trend,
            "plunge": plunge,
            "rake": rake,
            "pole": pole_v,
            "lin": lineation_v,
            "fabOr": fab_or,
            "fabRot": fab_rot,
            "how2plot": pc,
        },
    }
    maps = [item for item in items if isinstance(item, EBSD)]
    phases = [item for item in items if isinstance(item, str)]
    if len(maps) + len(phases) + len(conventions) != len(items):
        raise TypeError("optional items must be MTEX EBSD maps, phase names, or a plotting convention")
    if phases and not maps:
        raise ValueError("phase names require an EBSD map")
    if maps:
        sample["spec"]["ebsd"] = maps
        sample["geo"]["rotebsd"] = [rotEBSD2Geo(ebsd, fab_rot) for ebsd in maps]
        for phase in phases:
            sample["geo"][phase] = rotTex2Geo(maps[0][phase].orientations, fab_rot)
    return sample


__all__ = [
    "TP2vec", "V2TP", "geoFab", "SDT2or", "SDP2or", "FandL2Or", "rotTex2Geo",
    "rotEBSD2Geo", "rotODFGeo", "makeGeoSample",
]
