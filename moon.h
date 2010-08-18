/*
 * moon.h
 * Created by Serhii Tsyba (sertsy@gmail.com) on 03.08.10.
 *
 * This file contains various routines related to observations of the Moon.
 *
 * File provides two groups of functions: the first one may be used to construct lunar ephemerides and the second one
 * is useful for observations of the properties of the Moon.
 *
 * Functions related to lunar ephemerides provide functionality for computing
 *      True (geometric) position of the Moon.
 *      Apparent postion of the Moon. This function uses previous routine to determine geometric position of the
 *          Moon and then performs correction due to effect of nutation of rotational axis of the Earth to acheive
 *          apparent position.
 *          Note, that abovementioned apparent position of the Moon is geocentric, i.e. as if it was seen from the
 *          centre of the Earth. In order to get topocentric position of the planet (i.e. as seen from a certain
 *          location on the Earth globe) one should also perform correction for diurnal parallax and, if necessary, for
 *          the effect of atmospheric refraction. These corrections are very important in the case of the Moon due to
 *          Moon's close location to the Earth. For instance, correction for diurnal parallax may reach as much as one
 *          degree of arc.
 *      Distance to the Earth. Result is measured from centres of the bodies in astronomical units (AU).
 *
 * These functions rely on semi-analytical lunar theory ELP version ELP2000-82B developed by  M. Chapront-Touzé and
 * J. Chapront from Bureau des Longitudes, Paris, France. This theory considers perturbed motion of the Moon, i.e.
 * being effected by gravitational attraction not only Earth, but by the Sun and the planets. This implementation uses
 * full series given in the theory and thus claimed accuracy of the results does not exceed 800 seconds of arc in
 * longitude, 100 seconds of arc in latitude and 100 meters in distance for time span 1900 - 2000 A.D. compared to
 * lunar ephemeris LE51 developed by NASA Jet Propulsion Laboratory.
 *
 * For more information on theory ELP refer to source materials listed below.
 *
 * Functions related to observations of the physical properties of the Moon provide functionality for computing
 *      Phase angle of the Moon, which is the angle between the light incident to the Moon and the light reflected from
 *          it, as seen from the Earth. Namely, this is an angle between the Sun, the Moon and the Earth.
 *      Illuminated fraction of the disk of the Moon, i.e. the percentage of the Moon's disk being illuminated by
 *          the Sun as seen from Earth. This function may be used to determine current phase of the Moon.
 *      Position angle of the briht limb of the Moon. This is the angle that the bright limb of the Moon (part
 *          illuminated by the Sun) is turned measured from northernmost point of the Moon's disk seen as from Earth,
 *          but not true North Pole of the Moon.
 *
 * [1] M. Chapront-Touzé and J. Chapront. ELP 2000-85: a semi-analytical lunar ephemeris adequate for historical times.
 *     Astronomy and Astrophysics, vol. 190, 1988, pp. 342-352.
 * [2] M. Chapront-Touzé and J. Chapront. The lunar ephemeris ELP 2000. Astronomy and Astrophysics, vol. 124, 1983,
 *     pp. 50-62.
 */

#include "calendar.h"
#include "coordinates.h"

/*
 * Computes true (geometrical) geocentric position of the Moon on a given date (d). Output is expressed in ecliptic
 * coordinates referred to the mean ecliptic and equinox of date.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 307.
 */
ecliptic_point moon_true_position(date d);

/*
 * Computes apparent geocentric position of the Moon on a given date (d). Output is expressed in ecliptic coordinates
 * referred to the mean ecliptic and equinox of date.
 * Note, that apparent position of the Moon is geocentric, which means it does not involve correction for diurnal
 * parallax and correction for atmospheric refraction. This corrections should be taken care of separately, since it is
 * very important in case of the Moon.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 307.
 */
ecliptic_point moon_apparent_position(date d);

/*
 * Computes distance from the centre of the Moon to the centre of the Earth on a given date (d). Output is measured in
 * astronomical units (AU).
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 307.
 */
double moon_distance_to_earth(date d);

/*
 * Computes the phase angle of the Moon on a given date (d). Output is measured in radians and is in range [0, π].
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 315.
 */
double moon_phase_angle(date d);

/*
 * Computes persentage of the illuminated fraction of a disk of the Moon on a given date (d) as seen from Earth. Output
 * is a value in range [0, 1], where 0 mean that Moon's disk is not visible at all and 1 meaning that Moons's disk is
 * fully illuminated.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 315.
 */
double moon_disk_illuminated_fraction(date d);

/*
 * Computes position angle of Moon's bright limb on a given date (d), i.e. the angle of the midpoint of the illuminated
 * limb of the Moon reckoned eastward from the North Point of the disk (not from the axis of rotation of the lunar
 * globe). Output is measured in radians.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 316.
 */
double moon_bright_limb_position_angle(date d);