/*
 * geoobserver.h
 * Created by Serhii Tsyba (sertsy@gmail.com) on 19.07.10.
 *
 * This file contains definitions and routines of some phenomena apparent to the observer on the surface of the Earth.
 *
 * First, this file contains function that compute time of the day of rising and setting of the celestial for the
 * geographical location of the observer.
 *
 * Second, two functions that perform corrections to apparent position of the celestial body are given. These are
 * correction to body's apparent celestial position due to effect of diurnal parallax and correction to the body's
 * altitude due to the effect of atmospheric refraction.
 *
 * Finally, a function is given to compute parallactic angle of the clestial body. Parallactic angle is the angle
 * between apparent North pole of celestial body and and its zenith point - an uppermost point of the disk at the sky
 * as seen by observer at the given instant. Parallactic angle is not related to parallax and name derives from the
 * word 'parallel'.
 *
 * Note, that functions, given in this file make sense only, when apparent position of the celestial body is used.
 *
 * Also note that functions, given in this file prefer equatorial coordinate system over ecliptic, that is used
 * throughout the rest of the library. And when convertin from ecliptic coordinate system, the value of the obliquity
 * of the ecliptic should be corrected for the effect of nutation.
 */

#ifndef GEOOBSERVER_H
#define GEOOBSERVER_H

#include "calendar.h"
#include "coordinates.h"

#define SEA_LEVEL 0.0                   // defines sea level height constant for convinience
#define STANDARD_TEMPERATURE 283.15     // standard atmospheric temperature at sea level measured in K, ( = 10 °C)
#define STANDARD_PRESSURE 101325.0      // standard atmospheric pressure at sea level measured in Pa ( = 760 ㎜ Hg)

/*
 * Computes paralactic angle, commonly denotes as q of celestial body given time instant (d), goegraphic position (gp)
 * of the observer and body's position (eq) in equatorial coordinate system. Output value is measured in radians.
 *
 * Source: Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 94.
 */
double parallactic_angle(date d, geographic_point gp, equatorial_point ep);

/*
 * Computes rising time, in Universal Time (UT), of the celestial body on a given date (d), given geographic location
 * of the observer (gp), apparent position (ep) of the body at 0ʰ UT on a given date in equatorial coordinates and
 * body's 'standard' altitude (sa), i.e. geometric altitude of the center of the body at the time of apparent rising or
 * setting expressed in radians.
 * The following values can be used as 'standard' altitude (expressed in degrees, convert to radians when using):
 *      -0°.5667            for stars and planets of the Solar System
 *      -0°.8333            for the Sun
 *      -0.7275π - 0°.5667  for the Moon, where π is the Moon's horizontal parallax
 * Function returns rising time, in Universal Time, expressed in hours, of the body or -1 if body does not rise on a
 * given date.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 97.
 */
double rising(date d, geographic_point gp, equatorial_point ep, double sa);

/*
 * Computes setting time, in Universal Time (UT), of the celestial body on a given date (d), given geographic location
 * of the observer (gp), apparent position (ep) of the body at 0ʰ UT on a given date in equatorial coordinates and
 * body's 'standard' altitude (sa), i.e. geometric altitude of the center of the body at the time of apparent rising or
 * setting expressed in radians.
 * The following values can be used as 'standard' altitude (expressed in degrees, convert to radians when using):
 *      -0˚.5667            for stars and planets of the Solar System
 *      -0˚.8333            for the Sun
 *      -0.7275π - 0˚.5667  for the Moon, where π is the Moon's horizontal parallax
 * Function returns setting time, in Universal Time, expressed in hours, of the body or -1 if body does not rset on a
 * given date.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 97.
 */
double setting(date d, geographic_point gp, equatorial_point ep, double sa);

/*
 * Computes the value of apparent displacement of altitude of the celestial body due to effect of atmospheric
 * refraction given observers altitude (a) expressed in radians, atmospheric temperature (t) measured in Kelvins (°K)
 * and atmospheric pressure (p) measured in Pascals (Pa).
 * For adopted standard environmental conditions at sea level, use provided definitions.
 * Output value is measured in radians.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 105.
 */
double atmospheric_refraction(double a, double t, double p);

/*
 * Computes topocentric apparent position of the celestial body from geocentric apparent position (ep). Topocentric
 * apparent position of the body slightly changes due to effect of diurnal parallax and/or geographical location of the
 * observer. Both positions are given in equatorial coordinate systems.
 * Also observer's height above sea level (a) measured in kilometers, observer's geographic location (gp) and
 * equatorial horizontal parallax of the celesteial body (ehp) expressed in radians are required.
 * When observer's hieght above sea level needs to be ignored or unknown, use provided definition for convinience.
 * Output value is the input apparent position of the celestial body, corrected for diurnal parallax.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 263.
 */
equatorial_point diurnal_parallax(date d, geographic_point gp, double a, equatorial_point ep, double ehp);

#endif // GEOOBSERVER_H