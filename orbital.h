/*
 * orbital.h
 * Created by Serhii Tsyba (sertsy@gmail.com) on 29.07.10.
 *
 * This files provides routines that can compute following mean orbital elements of the major planets of the Solar
 * System:
 *          mean longitude of the planet (λ)
 *          semimajor axis of the orbit (a)
 *          eccentricity of the orbit (e)
 *          mean inclination of the plane of the ecliptic (i)
 *          mean longitude of the ascending node (Ω)
 *          mean longitude of the perhelion (ϖ)
 *
 * Orbital elements are referred to by their denotion given in paranthesis above throughout comments for shorter
 * notation.
 *
 * All mean orbital elements, except a and e are effected by precession of Earth's rotational axis. For this reason
 * two functions are given that compute mean orbital elements in different referrence frames: in one, elemetns are
 * referred to the standard equinox of J2000 and in second one to the mean dynamical ecliptic and equinox of date.
 *
 * Orbital elemetns λ, i, Ω and ϖ are expressed in radians and a is expressed in astronomical units (AU).
 *
 * Routines given in this file use semi-analytic planetary theory VSOP82 developed by P. Bretagnon. Refer to the source
 * material for more information.
 *
 * Source: P. Bretagnon. Théorie du mouvement de l'ensemble des planetès. Solution VSOP82. Astronomy and Astrophisics,
 *         vol. 114, 1982, pp. 277-287.
 *
 *         P. Bretagnon, G. Francou. Planetary thories in rectangular and spherical variables. VSOP87 solution.
 *         Astronomy and Astrophysics, vol. 202, 1988, pp. 309-315.
 */

#ifndef ORBITAL_H
#define ORBITAL_H

#include "calendar.h"
#include <vsop87d.h>

#define TOTAL_ORBITAL_ELEMENTS 6    // total amount of mean orbital elements

/*
 * An enumeration provided for convinient access to mean orbital arguments being computed using routines of this file.
 */
enum {
    PLANET_MEAN_LONGITUDE = 0,      // index of mean longitude of the planet (λ)
    SEMI_MAJOR_AXIS = 1,            // index of semimajor axis of the orbit (a)
    ECCENTRICITY_OF_ORBIT,          // index of eccentricity of the orbit (e)
    ECLIPTIC_INCLINATION,           // index of mean inclination of the plane of the ecliptic (i)
    ASCENDING_NODE_LONGITUDE,       // index of mean longitude of the ascending node (Ω)
    PERHELION_LONGITUDE = 5         // index of mean longitude of the perhelion (ϖ)
};

/*
 * Computes mean orbital arguments of a major (p) planet according to planetary theory VSOP82 given time instant (t) at
 * which elements are to be computed expressed in Julian millenia since the beginning of the epoch J2000.
 * Output is written into a given array (elements) in the following order: λ, a, e, i, Ω, ϖ. Use of provided
 * enumeration is advised for convinience. For the quantities denoted, refer to the top comment.
 * Result orbital elements are reckoned to the standard equinox of J2000.
 * Elements λ, i, Ω and ϖ are expressed in radians; a is expressed in astronomical units (AU).
 *
 * Source: P. Bretagnon. Théorie du mouvement de l'ensemble des planetès. Solution VSOP82. Astronomy and Astrophisics,
 *         vol. 114, 1982, pp. 277-287.
 */
void compute_orbital_elements_of_J2000(date d, Planet p, double elements[]);

/*
 * Computes mean orbital arguments of a major (p) planet according to planetary theory VSOP82 given time instant (t) at
 * which elements are to be computed expressed in Julian millenia since the beginning of the epoch J2000.
 * Output is written into a given array (elements) in the following order: λ, a, e, i, Ω, ϖ. Use of provided
 * enumeration is advised for convinience. For the quantities denoted, refer to top comment.
 * Result orbital elements are reckoned to the mean dynamical ecliptic and equinox of date.
 * Elements λ, i, Ω and ϖ are expressed in radians; a is expressed in astronomical units (AU).
 *
 * Source: P. Bretagnon. Théorie du mouvement de l'ensemble des planetès. Solution VSOP82. Astronomy and Astrophisics,
 *         vol. 114, 1982, pp. 277-287.
 */
void compute_orbital_elements_of_date(date d, Planet p, double elements[]);

#endif // ORBITAL_H