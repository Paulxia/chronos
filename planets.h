/*
 * planets.h
 * Created by Serhii Tsyba (sertsy@gmail.com) on 30.07.10.
 *
 * This file contains various routines related to observations of major planets of the Solar System.
 *
 * File provides two groups of functions: first one used to construct planetary ephemerides and second one useful for
 * observations of the planet.
 *
 * Functions related to planetary ephemerides provide functionality for computing
 *      True (geometric) position of the planet.
 *      Apparent postion of the planet. This function uses previous routine to determine geometric position of the
 *          planet and then performs corrections to acheive apparent position. These corrections include:
 *          - correction for light-time effect, i.e. planet being seen where it was when light left it;
 *          - correction for planetary aberration, i.e. effect of motion of the Earth along with the effect that finite
 *          velocity of light couses the apparent displament of the planet;
 *          - correction for nutations of the rotational axis of the Earth;
 *          Note, that with all the abovementioned correction apparent position of the planet is geocentric, i.e. as
 *          if it was seen from the centre of the Earth. In order to get topocentric position of the planet (i.e. as
 *          seen from a certain location on the Earth globe) one should also perform correction for diurnal parallax
 *          and, if necessary, for the effect of atmospheric refraction.
 *      Distance to the Sun. Result of this function is measured in astronomical units (AU).
 *      Distance to the Earth. Result is measurd in astronomical units (AU) and corresponds to the true distance to the
 *          Earth, not apparent.
 *
 * These functions rely on semi-analytical planetary theory VSOP87 version D developed by P. Bretagnon and G. Francou
 * at Service de Mèchanique Céleste du Bureau des Longitudes, Paris, France. This theory considers perturbed motion
 * of the planets, i.e. being effected by gravitational attraction not only Sun, but each other and thus is more
 * precise than analytical Keplerian approach. This implementation uses full series given in the theory and thus
 * claimed accuracy of the results does not exceed one second of arc for time span since 3000 B.C. to 3000 A.D.
 * (compared to planetary ephemeris DE200 developed by NASA Jet Propulsion Laboratory ?).
 *
 * For more information on theory VSOP87 refer to source material [1].
 *
 * Functions related to observations of the major planets provide functionality for computing
 *      Phase angle of the planet, which is the angle between the light incident to a planet and the light reflected
 *          from it, as seen from the Earth. Namely, this is an angle between Sun, planet and Earth. This angle is
 *          always in range [0, π]. For instance, for Venus and Mercury this angle can reach edges of a given range,
 *          however, for Mars the maxmum value of phase angle is never greater than about π/4.
 *      Illuminated fraction of the disk of a planet, i.e. the percentage of the planet's disk being illuminated by
 *          the Sun as seen from Earth.
 *      Apparent magnitude of the planet.
 *
 * Each function in the library accepts planet as a second argument. Current implementation adopts enumeration used in
 * the implementation of the VSOP87 theory, which is advised to be used for convinience. This enumeration indexes
 * planets in their order of appearence in the Solar System by distance from the Sun. Namely Mercury being 0, Venus
 * being 1, and Neptune being indexed 7.
 *
 * [1] P. Bretagnon and G. Francou. Planetary theories in rectangular and spherical variables. VSOP87 solutions.
 *     Astronomy and Astrophysics, vol. 202, 1988, pp 309-315.
 */

#ifndef PLANETS_H
#define PLANETS_H

#include "calendar.h"
#include "coordinates.h"
#include <vsop87d.h>

/*
 * Computes true (geometrical) geocentric position of a major planet (p) on a given date (d). Output is expressed in
 * ecliptic coordinates referred to the mean ecliptic and equinox of date. Returns (0, 0) if Earth is specified as a
 * target planet.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 209.
 */
ecliptic_point planet_true_position(date d, Planet p);

/*
 * Computes apparent geocentric position of a major planet (p) on a given date (d). Output is expressed in ecliptic
 * coordinates referred to the mean ecliptic and equinox of date. Returns (0, 0) if Earth is specified as a target
 * planet.
 * Note, that apparent position of a planet is geocentric, which means it does not involve correction for diurnal
 * parallax and correction for atmospheric refraction. This corrections should be taken care of separately, if
 * necessary.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 209.
 */
ecliptic_point planet_apparent_position(date d, Planet p);

/*
 * Computes distance from the centre of a major planet (p) to the centre of the Sun on a given date (d). Output is
 * measured in astronomical units (AU).
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 209.
 */
double planet_distance_to_sun(date d, Planet p);

/*
 * Computes distance from the centre of a given major planet (p) to the centre of the Earth on a given date (d). Output
 * is measured in astronomical units (AU). Note, that computed distance is true (geometrical), not apparent.
 * Function returns 0, if the Earth is specified as a target planet.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 209.
 */
double planet_distance_to_earth(date d, Planet p);

/*
 * Computes the phase angle of a major planet (p) on a given date (d). Output is measured in radians and is in range
 * [0, π]. Function returns 0 if the Earth is specified as a target planet.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 267.
 */
double planet_phase_angle(date d, Planet p);

/*
 * Computes persentage of the illuminated fraction of a disk of a given major planet (p) on a given date (d) as seen
 * from Earth. Output is a value in range [0, 1], where 0 mean that planet's disk is not visible at all and 1 meaning
 * that planet's disk is fully illuminated.
 * Function returns -1, if the Earth is specified as a target planet.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 267.
 */
double planet_disk_illuminated_fraction(date d, Planet p);

/*
 * Computes the value of the apparent magnitude of a major planet (p) on a given date (d). The brighter the planet
 * appears, the lower the output value.
 * Since function may output no reasonable value if Earth is given as a target planet, result in that case should be
 * ignored, i.e. it would contain some random value.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 269.
 */
double planet_apparent_magnitude(date d, Planet p);

#endif // PLANETS_H