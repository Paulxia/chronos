/*
 * coordinates.h
 * Created by Serhii Tsyba (sertsy@gmail.com) on 28.06.10.
 *
 * This file defines data types describing various coordinate systems used in astronomy and geography and routines to
 * convert amongst them. These coordinate systems are:
 *
 * Geographic. Used to describe a point on the surface of an Earth globe. Such point is given by two values:
 *      longitude - an angle measured positively westwards from Greenwich meridian along Earth equator;
 *      latitude - an angle measured positive north of Earth equator and negative south of it;
 *
 * Horizontal. Used to describe a point on the surface of the celestial sphere. Such point is given by two values:
 *      azimuth - an angle measured westwards from South along the plane of local horizon;
 *      elevation (also called altitude) - an angle positive above local horizon and negative below it;
 *
 * Equatorial. Used to describe a point on the surface of the celestial sphere. Such point is given by two values:
 *      right ascention - an angle measured eastwards of the vernal equinox along the plane of the celestial equator;
 *      declination - an angle measured positive north of the celestial equator and negative south of it;
 *
 * Ecliptic. Used to describe a point on the surface of the celestial sphere. Such point is given by two values:
 *      longitude - an angle measured eastwards from vernal equinox along the plane of ecliptic;
 *      latitude - an angle measured positive north of the ecliptic and negative south of it;
 *
 * Note, geography considers geographical longitude measured positive eastwards, but astronomy adopts universal
 * planetographical longitude to be measured westwards. For instance, geographical longitude of Washington D.C., U.S.A.
 * considered to be +77°02', but of Vienna, Austria: -16°23'.
 *
 * In similar manner, navigation considers azimuth measured from North, however, in astronomy azimuth is considered
 * to be measured from South, since hour angles are measured from South.
 *
 * All angles, and thus all of the above mentioned data types, in this library are expressed in radians, unless stated
 * otherwise.
 *
 * Ecliptic and coordinate system points are always referred to the mean equinox of date, unless stated otherwise.
 * To change reference to other epoch, proper correction for precession of the Earth's rotational axis must be made,
 * as well as any other corrections, if necessary.
 *
 * When converting between equatorial and ecliptic coordinate systems, a proper value of the obliquity of the ecliptic
 * must be given. If coordinates being converted represent true (geometric) position, then true value of the obliquity
 * of the ecliptic should be given. When apparent position is represented, then the value of the obliquity of the
 * ecliptic should be corrected for nutation of the Earth's rotational axis.
 *
 * Ecliptic coordinate system is the preffered celestial coordinate system throughout this library.
 */

#ifndef COORDINATES_H
#define COORDINATES_H

#include "calendar.h"

/*
 * Defines a datatype describing geographical coordinates consisting of geographical longitude (L) and latitude (ϕ)
 * measured in radians. Geographical longitude is measured positevely westwards from the Greenwich meridian and
 * negatively eastwards. Geographical latitude is positive in the northern hemisphere and negative in the southern.
 */
typedef struct {
    double longitude;           // L
    double latitude;            // ϕ
} geographic_point;

/*
 * A datatype describing horizontal coordinates consisting of azimuth (A) and elevation (altitude, h) measured in
 * radians. Azimuth is measured westwards from South. Elevation is positive above the horizon and negative below.
 */
typedef struct {
    double azimuth;             // A
    double elevation;           // h
} horizontal_point;

/*
 * A datatype defining equatorial coordinates consisting of right ascension (α) and declination (δ) measured in radians.
 * Decliantion is positive in in the northern part of celestial semisphere and is negative in the southern part.
 */
typedef struct {
    double right_ascension;     // α
    double declination;         // δ
} equatorial_point;

/*
 * A datatype defining ecliptic coordinates consisting of ecliptical (celestial) longitude (λ) and ecliptical
 * (celestial) latitude (β) measured in radians. Ecliptical longitude is measured from the vernal equinox along the
 * ecliptic and ecliptical latitude is positive if north of the ecliptic and negative if south.
 */
typedef struct {
    double longitude;           // λ
    double latitude;            // β
} ecliptic_point;

/*
 * Converts a point in equatorial coordinate system (eqp) to a point in ecliptical coordinate system given the
 * value of the obliquity of the ecliptic (e). If position of the celestial body being considered is apparent, then
 * value of the obliquity of the ecliptic should be corrected for nutation of the Earth's rotational axis.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 88.
 */
ecliptic_point equatorial_to_ecliptic(equatorial_point eqp, double e);

/*
 * Converts a point in ecliptic coordinate (ecp) system to a point in equatorial coordinate system given the value of
 * obliquity of the ecliptic (e). If position of the celestial body being considered is apparent, then value of the
 * obliquity of the ecliptic should be corrected for nutation of the Earth's rotational axis.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 89.
 */
equatorial_point ecliptic_to_equatorial(ecliptic_point ecp, double e);

/*
 * Converts a point in equatorial coordinate system (ep) located at geographical position (gp) to a point in local
 * horizontal coordinate system on a given date (d).
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 89.
 */
horizontal_point equatorial_to_horizontal(date d, geographic_point gp, equatorial_point ep);

/*
 * Converts a point in local horizontal coordinate system (hp) located at geographical position (gp) to a point in
 * equatorial coordinate system on a given date (d).
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 89.
 */
equatorial_point horizontal_to_equatorial(date d, geographic_point gp, horizontal_point hp);

#endif // COORDINTES_H