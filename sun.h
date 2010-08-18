/*
 * sun.h
 * Created by Serhii Tsyba (sertsy@gmail.com) on 26.07.10.
 *
 * This file contains various routine related to observations of the Sun.
 *
 * Three functions are given that allow computations of solar ephemerides:
 *      True (geometric) geocentric position of the Sun.
 *      Apparent positions of the Sun uses previous routine to determine geometric position and then performs
 *          corrections for effects of solar aberration and nutation of rotational axis of the Earth. It should be
 *          noted that this apparent position is geocentric (i.e. as seen from the centre of the Earth), and thus,
 *          requires additional corrections for diurnal parallax and atmospheric refraction if exact topocenric
 *          position (i.e. as seen from the surface of the Earth) is required.
 *      Distance from the Sun to the Earth, computed from centres of the bodies and is measured in astronomical units
 *          (AU).
 *
 * Positions, computed by abovementioned routines are expressed in ecliptic coordinated referred to the mean ecliptic
 * and equinox of date.
 *
 * These functions rely on semi-analytical planetary theory VSOP87 version D developed by P. Bretagnon and G. Francou
 * at Service de Mèchanique Céleste du Bureau des Longitudes, Paris, France. This theory considers perturbed motion
 * of the planets, i.e. being effected by gravitational attraction not only Sun, but each other and thus is more
 * precise than analytical Keplerian approach. This implementation uses full series given in the theory and thus
 * claimed accuracy of the results does not exceed one second of arc for time span since 3000 B.C. to 3000 A.D.
 * (compared to planetary ephemeris DE200 developed by NASA Jet Propulsion Laboratory ?)
 *
 * For more information on theory VSOP87 refer to source material [1]. 
 *
 * Functions that compute dates of solstice and equinoxes should be accurate up to seconds given the accuracy of
 * planetary theory VSOP87. Times of solstices and equinoxes is given in Dynamical Time and should be corrected to
 * Universal Time, if necessary.
 *
 * It should be noted, that names of the equinoxes and solstices are the ones used in Northern hemisphere. Since
 * Southern hemisphere has seasons opposite to the ones in Northern, to eliminate naming confusion, refer to the
 * description of solstice and equinox parameters to see their occurences dates.
 *
 * Also a function is given to compute the value of the equiation of time, commonly denoted as E. Equiation of time is
 * the difference between apparent and mean time, i.e. between the hour angle of the true Sun and the mean Sun.
 *
 * [1] P. Bretagnon and G. Francou. Planetary theories in rectangular and spherical variables. VSOP87 solutions.
 *     Astronomy and Astrophysics, vol. 202, 1988, pp 309-315.
 */

#ifndef SUN_H
#define SUN_H

#include "calendar.h"
#include "coordinates.h"

/*
 * An enumeration provided for convinient use with routine that computes equinoxes.
 */
typedef enum {
    VERNAL_EQUINOX = 0,         // index of vernal equinox (occurs in March)
    AUTUMNAL_EQUINOX = 2        // index of autumnal equinox (occurs in September)
} Equinox;

/*
 * An enumeration provided for convinient use with routine that computes solstices.
 */
typedef enum {
    SUMMER_SOLSTICE = 1,        // index of Summer (Northern) solstice (occurs in June)
    WINTER_SOLSTICE = 3         // index of Winter (Southern) solstice (occurs in December)
} Solstice;

/*
 * Computes true (geometrical) geocentric position of the Sun on a given date (d). Output is expressed in ecliptic
 * coordinates referred to the mean ecliptic and equinox of date.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 154.
 */
ecliptic_point sun_true_position(date d);

/*
 * Computes apparent geocentric position of the Sun on a given date (d). Output is expressed in ecliptic coordinates
 * referred to the mean ecliptic and equinox of date.
 * Note, that apparent position of the Sun is geocentric, which means it does not involve correction for diurnal
 * parallax and correction for atmospheric refraction. This corrections should be taken care of separately, if
 * necessary.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 154.
 */
ecliptic_point sun_apparent_position(date d);

/*
 * Computes distance from the centre of the Sun to the centre of the Earth. Output is measured in astronomical units
 * (AU). Note, that computed distance is true (geometrical), not apparent.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 154.
 */
double sun_distance_to_earth(date d);

/*
 * Computes date of equinox (e) of a given year (y). Output is measured in Dynamical Time.
 * Second argument specifies equinox to be computed; use of provided enumeration is advised.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 168.
 */
date equinox(int y, Equinox e);

/*
 * Computes date of solstice (s) of a given year (y). Output is measured in Dynamical Time.
 * Second argument specifies solstice to be computed; use of provided enumeration is advised.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 168.
 */
date solstice(int y, Solstice s);

/*
 * Computes the value of equiation of time, commonly denoted as E at a given date (d). Output value is measured in
 * hours.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 171.
 */
double solve_equation_of_time(date d);

#endif // SUN_H