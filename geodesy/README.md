# geodesy — implements geodetic computations and map projections on the WGS-84 reference ellipsoid

[![Go Reference](https://pkg.go.dev/badge/github.com/malcolmston/algebra/geodesy.svg)](https://pkg.go.dev/github.com/malcolmston/algebra/geodesy)

Package geodesy implements geodetic computations and map projections on the
WGS-84 reference ellipsoid (and other configurable ellipsoids) using only the
Go standard library.

- Great-circle (spherical) calculations: Haversine, spherical law of cosines
and equirectangular distances, initial and final bearings, midpoints,
intermediate/interpolated points, destination points, cross-track and
along-track distances, path intersection, maximum latitude of a great circle,
and spherical polygon/triangle area. - Rhumb-line (loxodrome) distance,
bearing, midpoint and destination. - Ellipsoidal geodesics via Vincenty's
direct and inverse formulae, including geodesic interpolation and waypoint
generation. - Map projections: spherical and Web (EPSG:3857) Mercator, the
Krüger-series transverse Mercator, and the Universal Transverse Mercator
(UTM) grid, with forward and inverse transforms. - Military Grid Reference
System (MGRS) encoding and decoding. - Earth-Centred Earth-Fixed (ECEF) <->
geodetic conversions and local East-North-Up (ENU) tangent-plane transforms. -
Angle helpers: degree/radian conversion, normalisation, and
degrees-minutes-seconds (DMS) parsing and formatting.

## Install

```sh
go get github.com/malcolmston/algebra@v0.8.0
```

```go
import "github.com/malcolmston/algebra/geodesy"
```

## Exported surface

### Functions

| Function | What it does |
| --- | --- |
| `func AlongTrackDistance(p, start, end LatLon) float64` | AlongTrackDistance is AlongTrackDistanceR on the WGS-84 mean sphere. |
| `func AlongTrackDistanceR(p, start, end LatLon, radius float64) float64` | AlongTrackDistanceR returns the distance in metres from the start point to the closest point on the great circle start->end to the foot of the… |
| `func AngularDistance(p1, p2 LatLon) float64` | AngularDistance returns the great-circle angular separation between two points, in radians, using the Haversine formula. |
| `func AuthalicLatitude(latDeg float64, e Ellipsoid) float64` | AuthalicLatitude returns the authalic (equal-area) latitude ξ (degrees) at the given geodetic latitude on the ellipsoid. |
| `func BackAzimuth(bearing float64) float64` | BackAzimuth returns the reverse bearing (degrees, [0,360)) of a forward bearing. |
| `func BearingToCompass16(bearing float64) string` | BearingToCompass16 returns the 16-point compass abbreviation (N, NNE, NE, ...) nearest to the given bearing (degrees). |
| `func BearingToCompass32(bearing float64) string` | BearingToCompass32 returns the 32-point compass abbreviation nearest to the given bearing (degrees). |
| `func BearingToCompass8(bearing float64) string` | BearingToCompass8 returns the 8-point compass abbreviation (N, NE, E, ...) nearest to the given bearing (degrees). |
| `func ChordDistance(p1, p2 LatLon, e Ellipsoid) float64` | ChordDistance returns the straight-line distance in metres through the Earth between two geodetic surface points on the given ellipsoid (height 0). |
| `func CompassToBearing(name string) (float64, error)` | CompassToBearing returns the central bearing (degrees) of the named compass point (case-insensitive, for example "NE" or "SSW"). |
| `func ConformalLatitude(latDeg float64, e Ellipsoid) float64` | ConformalLatitude returns the conformal latitude χ (degrees) at the given geodetic latitude on the ellipsoid. |
| `func CrossTrackDistance(p, start, end LatLon) float64` | CrossTrackDistance is CrossTrackDistanceR on the WGS-84 mean sphere. |
| `func CrossTrackDistanceR(p, start, end LatLon, radius float64) float64` | CrossTrackDistanceR returns the signed distance in metres of point p from the great circle defined by the start and end points, on a sphere of the… |
| `func DMSToDegrees(d, m, s float64) float64` | DMSToDegrees converts a degrees/minutes/seconds triple to decimal degrees. |
| `func DegToGrad(deg float64) float64` | DegToGrad converts degrees to gradians (400 gradians = 360 degrees). |
| `func DegToRad(deg float64) float64` | DegToRad converts an angle from degrees to radians. |
| `func DegreesToDMS(deg float64) (d int, m int, s float64)` | DegreesToDMS decomposes decimal degrees into a signed degree component and non-negative minutes and seconds. |
| `func ECEFToGeodetic(c ECEF, e Ellipsoid) (lat, lon, height float64)` | ECEFToGeodetic converts ECEF coordinates to geodetic latitude/longitude (degrees) and ellipsoidal height (metres) on the given ellipsoid, using… |
| `func ENUToGeodetic(l ENU, refLat, refLon, refHeight float64, e Ellipsoid) (lat, lon, height float64)` | ENUToGeodetic converts local East-North-Up coordinates back to a geodetic position (degrees latitude/longitude, metres height). |
| `func EllipsoidalMercatorForward(lat, lon float64, e Ellipsoid) (x, y float64)` | EllipsoidalMercatorForward projects geodetic latitude/longitude (degrees) to the conformal ellipsoidal Mercator projection (metres) on the given… |
| `func EllipsoidalMercatorInverse(x, y float64, e Ellipsoid) (lat, lon float64)` | EllipsoidalMercatorInverse recovers geodetic latitude/longitude (degrees) from ellipsoidal Mercator coordinates (metres) on the given ellipsoid. |
| `func EquirectangularDistance(p1, p2 LatLon, radius float64) float64` | EquirectangularDistance returns an approximate distance in metres using the equirectangular (flat-Earth) projection. |
| `func FeetToMeters(ft float64) float64` | FeetToMeters converts feet to metres. |
| `func FinalBearing(p1, p2 LatLon) float64` | FinalBearing returns the bearing at which the great circle from p1 arrives at p2, in degrees clockwise from north in [0,360). |
| `func FinalBearingRad(p1, p2 LatLon) float64` | FinalBearingRad returns the bearing (radians, [0,2π)) at which the great circle from p1 arrives at p2. |
| `func FormatDMS(deg float64, pos, neg byte, secPrec int) string` | FormatDMS renders decimal degrees as a degrees/minutes/seconds string with the given hemisphere characters for positive and negative values (for… |
| `func GeocentricLatitude(latDeg float64, e Ellipsoid) float64` | GeocentricLatitude returns the geocentric latitude (degrees) of a point on the ellipsoid surface at the given geodetic latitude. |
| `func GeocentricToGeodeticLatitude(geocentricDeg float64, e Ellipsoid) float64` | GeocentricToGeodeticLatitude returns the geodetic latitude (degrees) for a point on the ellipsoid surface at the given geocentric latitude. |
| `func GradToDeg(grad float64) float64` | GradToDeg converts gradians to degrees. |
| `func GradToRad(g float64) float64` | GradToRad converts gradians to radians. |
| `func GridConvergence(lat, lon, lon0, k0, falseEasting, falseNorthing float64, e Ellipsoid) float64` | GridConvergence returns the meridian convergence angle (degrees) of a transverse Mercator projection at the given geodetic point: the angle from grid… |
| `func HaversineDistance(p1, p2 LatLon) float64` | HaversineDistance returns the great-circle distance in metres on a sphere of the WGS-84 mean radius. |
| `func HaversineDistanceR(p1, p2 LatLon, radius float64) float64` | HaversineDistanceR returns the great-circle distance between two points on a sphere of the given radius, using the Haversine formula. |
| `func InitialBearing(p1, p2 LatLon) float64` | InitialBearing returns the initial bearing from p1 to p2 along the great circle, in degrees clockwise from north in [0,360). |
| `func InitialBearingRad(p1, p2 LatLon) float64` | InitialBearingRad returns the initial (forward azimuth) bearing from p1 to p2 along the great circle, in radians measured clockwise from north in… |
| `func IsometricLatitude(latDeg float64, e Ellipsoid) float64` | IsometricLatitude returns the isometric latitude ψ (dimensionless) at the given geodetic latitude (degrees) on the ellipsoid. |
| `func KilometersToMeters(km float64) float64` | KilometersToMeters converts kilometres to metres. |
| `func LatLonToMGRS(p LatLon, precision int) (string, error)` | LatLonToMGRS encodes a WGS-84 geodetic position as an MGRS reference string with the given precision (1..5). |
| `func LengthOfDegreeOfLatitude(latDeg float64, e Ellipsoid) float64` | LengthOfDegreeOfLatitude returns the north-south distance (metres) spanned by one degree of latitude at the given latitude (degrees) on the ellipsoid. |
| `func LengthOfDegreeOfLongitude(latDeg float64, e Ellipsoid) float64` | LengthOfDegreeOfLongitude returns the east-west distance (metres) spanned by one degree of longitude at the given latitude (degrees) on the ellipsoid. |
| `func LonLatToTile(lat, lon float64, zoom int) (tx, ty int)` | LonLatToTile returns the slippy-map tile column and row (x, y) containing the given WGS-84 latitude/longitude (degrees) at the given zoom level,… |
| `func MGRSLatBand(lat float64) byte` | MGRSLatBand returns the MGRS/UTM latitude band letter for the given latitude (degrees). |
| `func MaxLatitude(p LatLon, bearingDeg float64) float64` | MaxLatitude returns the maximum latitude reached by a great circle that crosses the given point at the given bearing (degrees), via Clairaut's… |
| `func MercatorForward(lat, lon, radius float64) (x, y float64)` | MercatorForward projects a spherical latitude/longitude (degrees) to planar (x, y) coordinates in metres on a sphere of the given radius, with the… |
| `func MercatorInverse(x, y, radius float64) (lat, lon float64)` | MercatorInverse recovers spherical latitude/longitude (degrees) from planar Mercator (x, y) coordinates in metres on a sphere of the given radius. |
| `func MeridianArcBetween(lat1, lat2 float64, e Ellipsoid) float64` | MeridianArcBetween returns the meridian arc distance (metres) between two latitudes (degrees) along the same meridian on the given ellipsoid. |
| `func MeridianArcLength(latDeg float64, e Ellipsoid) float64` | MeridianArcLength returns the meridian arc distance (metres) from the equator to the given geodetic latitude (degrees) on the given ellipsoid. |
| `func MetersToFeet(m float64) float64` | MetersToFeet converts metres to feet. |
| `func MetersToKilometers(m float64) float64` | MetersToKilometers converts metres to kilometres. |
| `func MetersToMiles(m float64) float64` | MetersToMiles converts metres to statute miles. |
| `func MetersToNauticalMiles(m float64) float64` | MetersToNauticalMiles converts metres to nautical miles. |
| `func MilesToMeters(mi float64) float64` | MilesToMeters converts statute miles to metres. |
| `func NauticalMilesToMeters(nm float64) float64` | NauticalMilesToMeters converts nautical miles to metres. |
| `func NormalizeBearing(bearing float64) float64` | NormalizeBearing reduces a compass bearing to the half-open interval [0, 360). |
| `func NormalizeDegrees(deg float64) float64` | NormalizeDegrees reduces an angle to the half-open interval [0, 360). |
| `func NormalizeLatitude(lat float64) float64` | NormalizeLatitude clamps a latitude to the closed interval [-90, 90]. |
| `func NormalizeLongitude(lon float64) float64` | NormalizeLongitude wraps a longitude into the half-open interval [-180, 180); the antimeridian is represented as -180. |
| `func NormalizeRadians(r float64) float64` | NormalizeRadians reduces an angle in radians to the half-open interval [0, 2π). |
| `func ParallelArcLength(latDeg, dLonDeg float64, e Ellipsoid) float64` | ParallelArcLength returns the distance (metres) along a parallel of latitude spanning dLon degrees of longitude at the given latitude on the… |
| `func ParallelChordDistance(latDeg, dLonDeg float64, e Ellipsoid) float64` | ParallelChordDistance returns the true east-west distance (metres) along the parallel at latitude latDeg spanning dLon degrees of longitude on the… |
| `func ParseDMS(s string) (float64, error)` | ParseDMS parses a coordinate written in degrees/minutes/seconds into decimal degrees. |
| `func PointInSphericalPolygon(p LatLon, vertices []LatLon) bool` | PointInSphericalPolygon reports whether p lies inside the spherical polygon defined by vertices (in order, closed automatically), using a winding /… |
| `func PointScaleFactor(lat, lon, lon0, k0, falseEasting, falseNorthing float64, e Ellipsoid) float64` | PointScaleFactor returns the local scale factor of a transverse Mercator projection at the given geodetic point (the ratio of grid distance to true… |
| `func QuarterMeridian(e Ellipsoid) float64` | QuarterMeridian returns the length (metres) of a meridian quadrant (equator to pole) on the given ellipsoid. |
| `func RadToDeg(rad float64) float64` | RadToDeg converts an angle from radians to degrees. |
| `func RadToGrad(r float64) float64` | RadToGrad converts radians to gradians. |
| `func RadiusOfParallel(latDeg float64, e Ellipsoid) float64` | RadiusOfParallel returns the radius (metres) of the circle of latitude (small circle) at the given latitude (degrees), N·cosφ, on the ellipsoid. |
| `func RectifyingLatitude(latDeg float64, e Ellipsoid) float64` | RectifyingLatitude returns the rectifying latitude μ (degrees) at the given geodetic latitude, the latitude scaled so that meridian arc length is… |
| `func ReducedLatitude(latDeg float64, e Ellipsoid) float64` | ReducedLatitude returns the reduced (parametric) latitude β in degrees for a point on the ellipsoid surface at geodetic latitude latDeg. |
| `func RelativeBearing(from, to float64) float64` | RelativeBearing returns the signed angle (degrees, [-180,180)) from bearing from to bearing to, that is the turn required to go from one heading to… |
| `func RhumbBearing(p1, p2 LatLon) float64` | RhumbBearing returns the constant compass bearing of the rhumb line from p1 to p2, in degrees clockwise from north in [0,360). |
| `func RhumbDistance(p1, p2 LatLon) float64` | RhumbDistance is RhumbDistanceR on the WGS-84 mean sphere. |
| `func RhumbDistanceR(p1, p2 LatLon, radius float64) float64` | RhumbDistanceR returns the rhumb-line (loxodrome) distance in metres between two points on a sphere of the given radius. |
| `func SphericalExcess(a, b, c float64) float64` | SphericalExcess returns the spherical excess (E = A + B + C - π) in radians of a spherical triangle with the three given interior angles in radians. |
| `func SphericalLawOfCosinesDistance(p1, p2 LatLon, radius float64) float64` | SphericalLawOfCosinesDistance returns the great-circle distance in metres on a sphere of the given radius, using the spherical law of cosines. |
| `func SphericalPolygonArea(vertices []LatLon, radius float64) float64` | SphericalPolygonArea returns the area in square metres enclosed by the spherical polygon whose vertices are given in order (the polygon is closed… |
| `func SphericalPolygonPerimeter(vertices []LatLon, radius float64) float64` | SphericalPolygonPerimeter returns the great-circle perimeter in metres of the polygon (closed automatically) on a sphere of the given radius. |
| `func SphericalTriangleArea(p1, p2, p3 LatLon, radius float64) float64` | SphericalTriangleArea returns the area in square metres of the spherical triangle with the three given vertices, on a sphere of the given radius. |
| `func TileToLonLat(tx, ty, zoom int) (lat, lon float64)` | TileToLonLat returns the WGS-84 latitude/longitude (degrees) of the north-west corner of the slippy-map tile (tx, ty) at the given zoom level. |
| `func TransverseMercatorForward(lat, lon, lon0, k0, falseEasting, falseNorthing float64, e Ellipsoid) (easting, northing float64)` | TransverseMercatorForward projects geodetic latitude/longitude (degrees) to a transverse Mercator grid on the given ellipsoid, using the Krüger… |
| `func TransverseMercatorInverse(easting, northing, lon0, k0, falseEasting, falseNorthing float64, e Ellipsoid) (lat, lon float64)` | TransverseMercatorInverse recovers geodetic latitude/longitude (degrees) from a transverse Mercator grid coordinate (easting, northing in metres) on… |
| `func UTMCentralMeridian(zone int) float64` | UTMCentralMeridian returns the central meridian (degrees) of the given UTM zone. |
| `func UTMGridConvergence(p LatLon) float64` | UTMGridConvergence returns the grid convergence angle (degrees) of the UTM projection at the given WGS-84 point. |
| `func UTMPointScale(p LatLon) float64` | UTMPointScale returns the local scale factor of the UTM projection at the given WGS-84 point. |
| `func UTMToMGRS(u UTMCoord, lat float64, precision int) (string, error)` | UTMToMGRS encodes a UTM coordinate as an MGRS reference string with the given precision (number of digits per easting/northing, 1..5; 5 = 1 m… |
| `func UTMZone(lat, lon float64) int` | UTMZone returns the UTM zone number (1..60) for the given longitude and latitude (degrees), applying the Norway and Svalbard zone exceptions. |
| `func VincentyDistance(p1, p2 LatLon) (float64, error)` | VincentyDistance returns the geodesic distance in metres between two points on the WGS-84 ellipsoid. |
| `func VincentyFinalBearing(p1, p2 LatLon) (float64, error)` | VincentyFinalBearing returns the final azimuth (degrees) of the geodesic from p1 to p2 on the WGS-84 ellipsoid. |
| `func VincentyInitialBearing(p1, p2 LatLon) (float64, error)` | VincentyInitialBearing returns the initial azimuth (degrees) of the geodesic from p1 to p2 on the WGS-84 ellipsoid. |
| `func WebMercatorForward(lat, lon float64) (x, y float64)` | WebMercatorForward projects WGS-84 latitude/longitude (degrees) to Web Mercator (EPSG:3857 / "Pseudo-Mercator") coordinates in metres, as used by… |
| `func WebMercatorInverse(x, y float64) (lat, lon float64)` | WebMercatorInverse recovers WGS-84 latitude/longitude (degrees) from Web Mercator (EPSG:3857) coordinates in metres. |
| `func WrapPi(r float64) float64` | WrapPi wraps an angle in radians to the half-open interval (-π, π]. |
| `func WrapPi180(d float64) float64` | WrapPi180 wraps an angle in degrees to the half-open interval [-180, 180). |

### Types

| Type | What it is |
| --- | --- |
| `AER` | AER is a local look-angle: azimuth and elevation in degrees and slant range in metres, relative to an observer at a geodetic reference point. |
| `BoundingBox` | BoundingBox is an axis-aligned latitude/longitude rectangle, defined by its south-west (minimum) and north-east (maximum) corners. |
| `DirectResult` | DirectResult holds the solution of the geodesic direct problem: the point reached from a start point along a geodesic, and the azimuth there. |
| `ECEF` | ECEF is a point in the Earth-Centred, Earth-Fixed Cartesian frame, with the origin at the centre of mass, the Z axis toward the reference pole, the X… |
| `ENU` | ENU is a local tangent-plane coordinate: East, North and Up offsets in metres relative to a reference geodetic origin. |
| `Ellipsoid` | Ellipsoid is a biaxial reference ellipsoid defined by its semi-major axis A (metres) and flattening F. |
| `InverseResult` | InverseResult holds the solution of the geodesic inverse problem: the geodesic distance between two points and the forward azimuths at each end. |
| `LatLon` | LatLon is a geodetic position given by latitude and longitude in degrees. |
| `UTMCoord` | UTMCoord is a Universal Transverse Mercator grid coordinate. |

<details>
<summary><code>AER</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ENUToAER(l ENU) AER` | ENUToAER converts a local East-North-Up offset to azimuth/elevation/range look-angles. |
| `func LookAngles(obsLat, obsLon, obsHeight, tgtLat, tgtLon, tgtHeight float64, e Ellipsoid) AER` | LookAngles returns the azimuth/elevation/range from an observer at the given geodetic reference point (degrees, metres height) to a target geodetic… |

</details>

<details>
<summary><code>BoundingBox</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func BoundingBoxAround(center LatLon, radius float64) BoundingBox` | BoundingBoxAround returns a bounding box centred on the given point that encloses a circle of the given radius (metres), computed on the WGS-84 mean… |
| `func BoundingBoxOfPoints(points []LatLon) BoundingBox` | BoundingBoxOfPoints returns the smallest bounding box containing all the given points. |
| `func NewBoundingBox(sw, ne LatLon) BoundingBox` | NewBoundingBox returns the bounding box with the given south-west and north-east corners, normalising the corner ordering. |
| `func (b BoundingBox) Center() LatLon` | Center returns the centre point of the box. |
| `func (b BoundingBox) Contains(p LatLon) bool` | Contains reports whether the point lies within the box (inclusive of edges). |
| `func (b BoundingBox) Corners() [4]LatLon` | Corners returns the four corners of the box in the order south-west, south-east, north-east, north-west. |
| `func (b BoundingBox) Expand(meters float64) BoundingBox` | Expand returns a copy of the box grown outward by the given amount of metres on all sides (WGS-84 mean sphere). |
| `func (b BoundingBox) HeightMeters() float64` | HeightMeters returns the north-south extent of the box in metres on the WGS-84 mean sphere. |
| `func (b BoundingBox) Intersects(o BoundingBox) bool` | Intersects reports whether two boxes overlap (sharing an edge counts). |
| `func (b BoundingBox) NorthEast() LatLon` | NorthEast returns the north-east corner of the box. |
| `func (b BoundingBox) SouthWest() LatLon` | SouthWest returns the south-west corner of the box. |
| `func (b BoundingBox) Union(o BoundingBox) BoundingBox` | Union returns the smallest box containing both b and o. |
| `func (b BoundingBox) WidthMeters() float64` | WidthMeters returns the east-west extent of the box in metres at its centre latitude on the WGS-84 mean sphere. |

</details>

<details>
<summary><code>DirectResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func VincentyDirect(start LatLon, bearingDeg, distance float64, e Ellipsoid) DirectResult` | VincentyDirect solves the direct geodesic problem on the given ellipsoid: from a start point, travel distance metres along the initial azimuth… |

</details>

<details>
<summary><code>ECEF</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func ENUToECEF(l ENU, refLat, refLon, refHeight float64, e Ellipsoid) ECEF` | ENUToECEF converts a local East-North-Up offset back to an ECEF point, relative to a geodetic reference origin on the given ellipsoid. |
| `func GeodeticToECEF(lat, lon, height float64, e Ellipsoid) ECEF` | GeodeticToECEF converts geodetic latitude/longitude (degrees) and ellipsoidal height (metres) to ECEF coordinates on the given ellipsoid. |
| `func (c ECEF) Add(o ECEF) ECEF` | Add returns the vector sum of two ECEF points. |
| `func (c ECEF) Cross(o ECEF) ECEF` | Cross returns the cross product of two ECEF vectors. |
| `func (c ECEF) Distance(o ECEF) float64` | Distance returns the straight-line (chord) distance between two ECEF points. |
| `func (c ECEF) Dot(o ECEF) float64` | Dot returns the dot product of two ECEF vectors. |
| `func (c ECEF) Norm() float64` | Norm returns the Euclidean length of the ECEF vector. |
| `func (c ECEF) Normalize() ECEF` | Normalize returns the unit vector in the direction of c (the zero vector is returned unchanged). |
| `func (c ECEF) Scale(s float64) ECEF` | Scale returns the ECEF vector scaled by s. |
| `func (c ECEF) Sub(o ECEF) ECEF` | Sub returns the vector difference c - o. |

</details>

<details>
<summary><code>ENU</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func AERToENU(a AER) ENU` | AERToENU converts azimuth/elevation/range look-angles to a local East-North-Up offset. |
| `func ECEFToENU(c ECEF, refLat, refLon, refHeight float64, e Ellipsoid) ENU` | ECEFToENU converts an ECEF point to local East-North-Up coordinates relative to a geodetic reference origin (degrees, metres) on the given ellipsoid. |
| `func GeodeticToENU(lat, lon, height, refLat, refLon, refHeight float64, e Ellipsoid) ENU` | GeodeticToENU converts a geodetic point (degrees, metres) to local East-North-Up coordinates relative to a geodetic reference origin. |
| `func (l ENU) HorizontalRange() float64` | HorizontalRange returns the horizontal (ground) distance in metres of the ENU offset, ignoring the vertical component. |
| `func (l ENU) SlantRange() float64` | SlantRange returns the full 3-D distance in metres of the ENU offset. |

</details>

<details>
<summary><code>Ellipsoid</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func NewEllipsoid(a, invFlattening float64) Ellipsoid` | NewEllipsoid builds an ellipsoid from a semi-major axis and inverse flattening 1/f. |
| `func (e Ellipsoid) AuthalicRadius() float64` | AuthalicRadius returns the radius of a sphere with the same surface area. |
| `func (e Ellipsoid) FirstEccentricity() float64` | FirstEccentricity returns the first eccentricity e. |
| `func (e Ellipsoid) FirstEccentricitySquared() float64` | FirstEccentricitySquared returns e² = f(2-f). |
| `func (e Ellipsoid) Flattening() float64` | Flattening returns the flattening f = (a-b)/a. |
| `func (e Ellipsoid) GaussianRadius(latDeg float64) float64` | GaussianRadius returns the Gaussian (geometric mean) radius of curvature sqrt(M·N) at the given geodetic latitude (degrees). |
| `func (e Ellipsoid) InverseFlattening() float64` | InverseFlattening returns 1/f, or +Inf for a sphere. |
| `func (e Ellipsoid) IsSphere() bool` | IsSphere reports whether the ellipsoid is a perfect sphere (f == 0). |
| `func (e Ellipsoid) MeanRadius() float64` | MeanRadius returns the arithmetic mean radius R1 = (2a+b)/3. |
| `func (e Ellipsoid) MeridianRadius(latDeg float64) float64` | MeridianRadius returns the meridional radius of curvature M at the given geodetic latitude (degrees). |
| `func (e Ellipsoid) PrimeVerticalRadius(latDeg float64) float64` | PrimeVerticalRadius returns the radius of curvature N in the prime vertical at the given geodetic latitude (degrees). |
| `func (e Ellipsoid) SecondEccentricity() float64` | SecondEccentricity returns the second eccentricity e'. |
| `func (e Ellipsoid) SecondEccentricitySquared() float64` | SecondEccentricitySquared returns e'² = e²/(1-e²). |
| `func (e Ellipsoid) SemiMajorAxis() float64` | SemiMajorAxis returns the equatorial radius a in metres. |
| `func (e Ellipsoid) SemiMinorAxis() float64` | SemiMinorAxis returns the polar radius b = a(1-f) in metres. |
| `func (e Ellipsoid) ThirdFlattening() float64` | ThirdFlattening returns n = f/(2-f), used by the Krüger series. |
| `func (e Ellipsoid) VolumetricRadius() float64` | VolumetricRadius returns the radius of a sphere with the same volume, (a²b)^(1/3). |

</details>

<details>
<summary><code>InverseResult</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func VincentyInverse(p1, p2 LatLon, e Ellipsoid) (InverseResult, error)` | VincentyInverse solves the inverse geodesic problem on the given ellipsoid using Vincenty's formula, returning the distance and initial/final… |

</details>

<details>
<summary><code>LatLon</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func Destination(start LatLon, bearingDeg, distance float64) LatLon` | Destination returns the point reached from start along a great circle on the given initial bearing (degrees) for distance metres, on the WGS-84 mean… |
| `func DestinationR(start LatLon, bearingDeg, distance, radius float64) LatLon` | DestinationR returns the point reached by travelling distance metres along a great circle from the start point on the given initial bearing… |
| `func GeodesicInterpolate(p1, p2 LatLon, fraction float64, e Ellipsoid) (LatLon, error)` | GeodesicInterpolate returns the point at the given fraction of the geodesic distance from p1 (fraction 0) to p2 (fraction 1) on the given ellipsoid. |
| `func GeodesicWaypoints(p1, p2 LatLon, n int, e Ellipsoid) ([]LatLon, error)` | GeodesicWaypoints returns n+1 points evenly spaced by geodesic distance from p1 to p2 inclusive on the given ellipsoid. |
| `func IntermediatePoint(p1, p2 LatLon, fraction float64) LatLon` | IntermediatePoint returns the point at the given fraction of the great-circle distance from p1 (fraction 0) to p2 (fraction 1). |
| `func InterpolatePoints(p1, p2 LatLon, n int) []LatLon` | InterpolatePoints returns n+1 points evenly spaced (by fraction of the total great-circle distance) from p1 to p2 inclusive. |
| `func Intersection(p1 LatLon, bearing1 float64, p2 LatLon, bearing2 float64) (LatLon, error)` | Intersection returns the point of intersection of two great-circle paths, each defined by a point and an initial bearing (degrees). |
| `func MGRSToLatLon(s string) (LatLon, error)` | MGRSToLatLon decodes an MGRS reference string into a WGS-84 geodetic position (the south-west corner of the referenced square when the numeric… |
| `func Midpoint(p1, p2 LatLon) LatLon` | Midpoint returns the half-way point along the great circle between p1 and p2. |
| `func NewLatLon(lat, lon float64) LatLon` | NewLatLon constructs a LatLon from latitude and longitude in degrees. |
| `func RhumbDestination(start LatLon, bearingDeg, distance float64) LatLon` | RhumbDestination is RhumbDestinationR on the WGS-84 mean sphere. |
| `func RhumbDestinationR(start LatLon, bearingDeg, distance, radius float64) LatLon` | RhumbDestinationR returns the point reached by travelling distance metres along a rhumb line from start on the given constant bearing (degrees), on a… |
| `func RhumbMidpoint(p1, p2 LatLon) LatLon` | RhumbMidpoint returns the point half-way along the rhumb line between p1 and p2 (loxodromic midpoint). |
| `func UTMToLatLon(u UTMCoord) (LatLon, error)` | UTMToLatLon converts a UTM grid coordinate to a WGS-84 geodetic position. |
| `func UTMToLatLonEllipsoid(u UTMCoord, e Ellipsoid) (LatLon, error)` | UTMToLatLonEllipsoid converts a UTM coordinate to a geodetic position on the given ellipsoid. |
| `func VincentyDestination(start LatLon, bearingDeg, distance float64) LatLon` | VincentyDestination returns the point reached from start along a WGS-84 geodesic on the given initial bearing (degrees) for distance metres. |
| `func (p LatLon) Antipode() LatLon` | Antipode returns the point diametrically opposite on the sphere. |
| `func (p LatLon) CrossTrackDistanceTo(start, end LatLon) float64` | CrossTrackDistanceTo returns the signed distance in metres of p from the great circle defined by start and end (WGS-84 mean sphere). |
| `func (p LatLon) DestinationPoint(bearing, distance float64) LatLon` | DestinationPoint returns the point reached from p along a great circle on the given bearing (degrees) for distance metres on the WGS-84 mean sphere. |
| `func (p LatLon) DistanceTo(q LatLon) float64` | DistanceTo returns the great-circle distance in metres from p to q on the WGS-84 mean sphere (Haversine). |
| `func (p LatLon) Equal(q LatLon, tol float64) bool` | Equal reports whether two positions are within tol degrees on both axes. |
| `func (p LatLon) FinalBearingTo(q LatLon) float64` | FinalBearingTo returns the final great-circle bearing (degrees) from p to q. |
| `func (p LatLon) GeodesicDistanceTo(q LatLon) (float64, error)` | GeodesicDistanceTo returns the geodesic distance in metres from p to q on the WGS-84 ellipsoid (Vincenty inverse). |
| `func (p LatLon) InitialBearingTo(q LatLon) float64` | InitialBearingTo returns the initial great-circle bearing (degrees) from p to q. |
| `func (p LatLon) LatRad() float64` | LatRad returns the latitude in radians. |
| `func (p LatLon) LonRad() float64` | LonRad returns the longitude in radians. |
| `func (p LatLon) MidpointTo(q LatLon) LatLon` | MidpointTo returns the great-circle midpoint between p and q. |
| `func (p LatLon) Normalized() LatLon` | Normalized returns the position with latitude clamped to [-90,90] and longitude wrapped to [-180,180). |
| `func (p LatLon) RhumbBearingTo(q LatLon) float64` | RhumbBearingTo returns the constant rhumb-line bearing (degrees) from p to q. |
| `func (p LatLon) RhumbDistanceTo(q LatLon) float64` | RhumbDistanceTo returns the rhumb-line distance in metres from p to q on the WGS-84 mean sphere. |
| `func (p LatLon) String() string` | String renders the position in signed decimal degrees. |
| `func (p LatLon) ToECEF(height float64) ECEF` | ToECEF converts p at the given ellipsoidal height (metres) to ECEF coordinates on the WGS-84 ellipsoid. |
| `func (p LatLon) ToMGRS(precision int) (string, error)` | ToMGRS encodes p as an MGRS reference string with the given precision. |
| `func (p LatLon) ToUTM() (UTMCoord, error)` | ToUTM converts p to a UTM grid coordinate on WGS-84. |
| `func (p LatLon) Valid() bool` | Valid reports whether the latitude is in [-90,90] and longitude in [-180,180]. |

</details>

<details>
<summary><code>UTMCoord</code> — constructors and methods</summary>

| Signature | What it does |
| --- | --- |
| `func LatLonToUTM(p LatLon) (UTMCoord, error)` | LatLonToUTM converts a WGS-84 geodetic position to a UTM grid coordinate. |
| `func LatLonToUTMEllipsoid(p LatLon, e Ellipsoid) (UTMCoord, error)` | LatLonToUTMEllipsoid converts a geodetic position to a UTM coordinate on the given ellipsoid. |
| `func ParseMGRS(s string) (UTMCoord, error)` | ParseMGRS decodes an MGRS reference string into a UTM coordinate. |
| `func ParseUTM(s string) (UTMCoord, error)` | ParseUTM parses a UTM string in the form produced by UTMCoord.String, for example "31 N 448252 5411932" or "31N 448252 5411932". |
| `func (u UTMCoord) String() string` | String renders a UTM coordinate as "zone hemisphere easting northing", for example "31 N 448252 5411932". |

</details>

### Constants

`MetersPerNauticalMile`, `MetersPerStatuteMile`, `MetersPerFoot`, `MetersPerKilometer`, `EarthRadiusMean`, `EarthRadiusAuthalic`, `EarthRadiusVolumetric`, `EarthRadiusEquatorial`, `EarthRadiusPolar`, `UTMScaleFactor`, `UTMFalseEasting`, `UTMFalseNorthingSouth`, `WebMercatorMaxLatitude`

### Variables

`ErrNoConvergence`, `ErrInvalidCoordinate`, `ErrInvalidMGRS`, `ErrInvalidUTM`, `ErrPolarRegion`, `WGS84`, `GRS80`, `Airy1830`, `International1924`, `Clarke1866`, `Bessel1841`, `SphereWGS84`

Full signatures, doc comments and every runnable example are on
[pkg.go.dev](https://pkg.go.dev/github.com/malcolmston/algebra/geodesy).

## License

MIT, as part of [`github.com/malcolmston/algebra`](..).
An independent re-implementation, not affiliated with or endorsed by the
original project.
