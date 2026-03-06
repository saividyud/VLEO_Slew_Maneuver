# Aerospace Toolbox Public API (MATLAB R2025b)

This list was derived from the installed MATLAB R2025b tree on this machine and cross-checked against:

- `/usr/local/MATLAB/R2025b/toolbox/aero/aero/Contents.m`
- installed public files under `/usr/local/MATLAB/R2025b/toolbox/aero`
- `which(...)` lookups for satellite-scenario APIs

It is intended to capture public Aerospace Toolbox entry points and omit `private`, `internal`, demo, and adapter code.

## Axes Transformations

- `angle2dcm`
- `angle2quat`
- `angle2rod`
- `dcm2alphabeta`
- `dcm2angle`
- `dcm2latlon`
- `dcm2quat`
- `dcm2rod`
- `dcmbody2stability`
- `dcmbody2wind`
- `dcmecef2ned`
- `dcmeci2ecef`
- `ecef2lla`
- `eci2aer`
- `eci2lla`
- `flat2lla`
- `geoc2geod`
- `geod2geoc`
- `lla2ecef`
- `lla2eci`
- `lla2flat`
- `quat2angle`
- `quat2dcm`
- `quat2rod`
- `rod2angle`
- `rod2dcm`
- `rod2quat`

## Environment

- `atmoscira`
- `atmoscoesa`
- `atmoshwm`
- `atmoshwm07`
- `atmosisa`
- `atmoslapse`
- `atmosnonstd`
- `atmosnrlmsise00`
- `atmospalt`
- `geoidheight`
- `geoidegm96`
- `gravitycentrifugal`
- `gravitysphericalharmonic`
- `gravitywgs84`
- `gravityzonal`
- `igrfmagm`
- `igrf11magm`
- `wrldmagm`
- `zonalplanetparams`

## Flight Parameters

- `airspeed`
- `alphabeta`
- `correctairspeed`
- `dpressure`
- `geocradius`
- `machnumber`
- `rrdelta`
- `rrsigma`
- `rrtheta`

## Quaternion Math

- `quatconj`
- `quatdivide`
- `quatexp`
- `quatinv`
- `quatinterp`
- `quatlog`
- `quatmod`
- `quatmultiply`
- `quatnorm`
- `quatnormalize`
- `quatpower`
- `quatrotate`

## Time

- `decyear`
- `deltaUT1`
- `juliandate`
- `leapyear`
- `mjuliandate`
- `tdbjuliandate`

## Unit Conversion

- `convacc`
- `convang`
- `convangacc`
- `convangvel`
- `convdensity`
- `convforce`
- `convlength`
- `convmass`
- `convpres`
- `convtemp`
- `convvel`
- `getunitdata`
- `unitconversion`

## Gas Dynamics

- `flowfanno`
- `flowisentropic`
- `flownormalshock`
- `flowprandtlmeyer`
- `flowrayleigh`

## Celestial Phenomena

- `aeroReadSpaceWeatherData`
- `deltaCIP`
- `earthNutation`
- `fluxSolarAndGeomagnetic`
- `moonLibration`
- `planetEphemeris`
- `polarMotion`

## Utilities and File Reading

- `aeroDataPackage`
- `aeroReadIERSData`
- `datcomimport`

## Flight Instruments

- `uiaeroairspeed`
- `uiaeroaltimeter`
- `uiaeroclimb`
- `uiaeroegt`
- `uiaeroheading`
- `uiaerohorizon`
- `uiaerorpm`
- `uiaeroturn`

## Animation

- `Aero.Animation`
- `Aero.Body`
- `Aero.Camera`
- `Aero.FlightGearAnimation`
- `Aero.Geometry`
- `Aero.Node`
- `Aero.Viewpoint`
- `Aero.VirtualRealityAnimation`
- `CustomReadBodyTSData`
- `doFirstOrderChaseCameraDynamics`
- `fganimation`
- `staticCameraPosition`

## Aircraft

- `aircraftEnvironment`
- `aircraftProperties`
- `Aero.Aircraft.ControlState`
- `Aero.Aircraft.Environment`
- `Aero.Aircraft.Properties`

## Fixed Wing

- `fixedWingAircraft`
- `fixedWingCoefficient`
- `fixedWingState`
- `fixedWingStateCustom`
- `fixedWingSurface`
- `fixedWingThrust`
- `Aero.FixedWing`
- `Aero.FixedWing.Coefficient`
- `Aero.FixedWing.State`
- `Aero.FixedWing.Surface`
- `Aero.FixedWing.Thrust`

## Satellite Scenarios

- `satelliteScenario`
- `satelliteScenarioViewer`
- `satellite`
- `walkerDelta`
- `walkerStar`
- `groundStation`
- `gimbal`
- `conicalSensor`
- `access`
- `fieldOfView`
- `eclipse`

Notes:

- `satelliteScenario`, `satelliteScenarioViewer`, `satellite`, `walkerDelta`, `walkerStar`, and `groundStation` resolve as installed entry points.
- `gimbal`, `conicalSensor`, `access`, `fieldOfView`, and `eclipse` are exposed as method-style APIs on scenario objects rather than simple top-level files.

## Spacecraft

- `propagateOrbit`
- `tleread`
- `ommread`

## Spacecraft Transformations

- `ecef2eci`
- `eci2ecef`
- `ijk2keplerian`
- `keplerian2ijk`
- `siderealTime`
- `greenwichSRT`

## Graphics

- `altitudeEnvelopeContour`
- `boundaryline`
- `shortPeriodCategoryAPlot`
- `shortPeriodCategoryBPlot`
- `shortPeriodCategoryCPlot`
