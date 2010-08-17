/*
 * earth.h
 * Created by Serhii Tsyba (sertsy@gmail.com) on 24.06.10.
 *
 * This file contains definitions and routines that compute effects of some phenomena produced by the Earth. These
 * phenomena include:
 *      precession - a slow change in the orientation of the Earth's rotational axis due to gravitational attraction
 *                   from the Sun and the Moon;
 *      nutation - a periodic oscillation of the Earth's rotational axis around its mean position, mainly due to
 *                 gravitational attraction from the Moon;
 *      annual aberration - an apparent movement of celestial object due to finite speed of light;
 *      obliquity of the ecliptic - an angle that plane of ecliptic makes with the plane of equator;
 *
 * Effects of these phenomena are needed to be taken into account when computing apparent position of celestial bodies.
 * Computations of obliquity of the ecliptic uses the formula adopted by International Astronomical Union (IAU).
 * Computation of nutation uses numerical method proposed by IAU Theory of Nutation in 1980.
 *
 * All of the above mentioned routines work on ecliptic coordinates.
 *
 * Another function tha might come handy computes geodesic distance between to locations on Earth. The distance is
 * computed along the surface of the planet and assumes geoid to be ellipsoid.
 */

#ifndef EARTH_H
#define EARTH_H

#include "calendar.h"
#include "coordinates.h"

#define EARTH_EQUATORIAL_RADIUS 6378.14         // equatorial radius (a) of the Earth measured in kilometers
#define EARTH_POLAR_RADIUS 6356.755             // polar radius (b) of the Earth measured in kilometers
#define EARTH_FLATTERING 0.00335281             // flatterinf of the Earth (f)
#define EARTH_MERIDIAN_ECCENTRICITY 0.08181922  // eccentricity of the Earth meridian

#define ABERRATION_CONSTANT 20.49552            // value of the constant of aberration at J2000 measured in arcseconds
#define ASTRONOMICAL_UNIT 149597871             // astronomical unit (AU) measured in kilometers

/*
 * Computes geodesic distance (shortest distance measured along Earth's surface) between two locations (gp1 and gp2) on
 * the Earth geoid surface. Result is measured in kilometers.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 81.
 */
double geodesic_distance(geographic_point gp1, geographic_point gp2);

/*
 * Performs reduction of ecliptic coordinates (ep) from one epoch (jd1) to another (jd2) due to effect of precession of
 * Earth's axis. Starting and target epochs are given as Julian dates. Output is the input reduced to a given epoch.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, pp. 128-129.
 */
ecliptic_point precession(ecliptic_point ep, double jd1, double jd2);

/*
 * Computes value of nutation in longitude, commonly denoted as Δψ, for a given date (d). Uses 1980 IAU Theory of
 * Nutation. Output value is measured in radians.
 *
 * Source: P.K. Seidelman. 1980 IAU Theory of Nutation: The Final Report of the IAU Working Group on Nutation.
 *         U.S. Naval Observatory, Nautical Almanac Office, Washington, D.C. 20390, U.S.A., 1981.
 */
double nutation_in_longitude(date d);

/*
 * Computes value of nutation in obliquity, commonly denoted as Δε for a given date (d). Uses 1980 IAU Theory of
 * Nutation. Output value is measured in radians.
 *
 * Source: P.K. Seidelman. 1980 IAU Theory of Nutation: The Final Report of the IAU Working Group on Nutation.
 *         U.S. Naval Observatory, Nautical Almanac Office, Washington, D.C. 20390, U.S.A., 1981.
 */
double nutation_in_obliquity(date d);

/*
 * Computes changes in longitude and latitude of a given ecliptic point (ep) due to effect of annual aberration on a
 * given date (d). Output value is a change in longitude and latitude (Δλ, Δβ) due to effect of annual aberration.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 139.
 */
ecliptic_point aberration(date d, ecliptic_point ep);

/*
 * Computes mean obliquity of the ecliptic (inclination of Earth's rotational axis of the mean equator), commonly
 * denoted as e for a given date (d). Uses the formula adopted by IAU. Output value is measured in radians.
 *
 * Source: J.H. Lieske, T. Lederle, W. Fricke and B. Morando. Expressions for the precession Quantities Based upon
 *         the IAU (1976) System of Astronomical Constants, Astronomy ans Astrophysics, vol. 58, 1977, p. 15.
 */
double obliquity_of_ecliptic(date d);


#endif // EARTH_H