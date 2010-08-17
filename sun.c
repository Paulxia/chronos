/*
 * sun.c
 * Created by Serhii Tsyba (sertsy@gmail.com) on 26.07.10.
 */

#include "sun.h"
#include "earth.h"
#include "orbital.h"
#include "coordinates.h"
#include <vsop87d.h>
#include <math.h>

ecliptic_point sun_true_position(date d)
{
    double t;                   // time interval since the beginning of the epoch J2000 to given date in Julian centuries
    spherical_point sp;         // heliocentric position of the Earth
    double lp;                  // auxiliary computational variable
    ecliptic_point ep;          // result geocentric true position of the Sun
    
    // calculating time since the beginning of the epoch J2000 measured in Julian centuries
    t = (julian_ephemeris_date(d) - J2000) / DAYS_IN_JULIAN_CENTURY;
    
    // calculating heliocentric position of the Earth using semi-analytical planetary theory VSOP87
    // VSOP87D uses time measured in Julian millenia instead of centuries, hence division by ten
    sp = heliocentric_planetary_position(t / 10.0, EARTH);
    
    // since VSOP87D position is reckoned to mean dynamical ecliptic, converting to FK5 system
    lp = sp.longitude - (1.397 * t + 0.00031 * t * t) * M_PI / 180.0;
    sp.longitude -= 0.09033 * M_PI / 180.0 / 3600.0;
    sp.latitude += 0.03916 * (cos(lp) - sin(lp)) * M_PI / 180.0 / 3600.0;
    
    // calculating geocentric position of the Sun, which is inverse heliocentric position of the Earth
    ep.longitude = sp.longitude + M_PI;
    ep.latitude = -sp.latitude;
    
    // shifting longtitude in range [0, 2π]
    ep.longitude = fmod(ep.longitude, 2 * M_PI);
    if (ep.longitude < 0.0)
        ep.longitude += 2 * M_PI;
    
    return ep;
}

ecliptic_point sun_apparent_position(date d)
{
    double t;                   // time interval since the beginning of the epoch J2000 to given date in Julian centuries
    spherical_point sp;         // heliocentric position of the Earth
    double lp;                  // auxiliary computational variable
    ecliptic_point ep;          // result geocentric true position of the Sun
    
    // calculating time since the beginning of the epoch J2000 measured in Julian centuries
    t = (julian_ephemeris_date(d) - J2000) / DAYS_IN_JULIAN_CENTURY;
    
    // calculating heliocentric position of the Earth using semi-analytical planetary theory VSOP87
    // VSOP87D uses time measured in Julian millenia instead of centuries, hence division by ten
    sp = heliocentric_planetary_position(t / 10.0, EARTH);
    
    // since VSOP87D position is reckoned to mean dynamical ecliptic, converting to FK5 system
    lp = sp.longitude - (1.397 * t + 0.00031 * t * t) * M_PI / 180.0;
    sp.longitude -= 0.09033 * M_PI / 180.0 / 3600.0;
    sp.latitude += 0.03916 * (cos(lp) - sin(lp)) * M_PI / 180.0 / 3600.0;
    
    // calculating geocentric position of the Sun, which is inverse heliocentric position of the Earth
    ep.longitude = sp.longitude + M_PI;
    ep.latitude = -sp.latitude;
    
    // performing correction for nutation
    ep.longitude += nutation_in_longitude(d);
    
    // performing correction for aberration
    // in the case of Sun, aberration correction formula simplifies to -k / R, where k is aberration constant and R is
    // Sun - Earth distance expressed in AU
    // aberration constant is given in arcseconds, thus converting to radians, (π = 648000")
    ep.longitude -= (ABERRATION_CONSTANT * M_PI / 648000) / sp.distance;
    
    // shifting longtitude in range [0, 2π]
    ep.longitude = fmod(ep.longitude, 2 * M_PI);
    if (ep.longitude < 0.0)
        ep.longitude += 2 * M_PI;
    
    return ep;
}

double sun_distance_to_earth(date d)
{
    double t;                   // time interval since the beginning of the epoch J2000 to given date in Julian millenia
    spherical_point sp;         // true position of the Sun at a given date
    
    // calculating time since the beginning of the epoch J2000 measured in Julian millenia till a given date 
    t = (julian_ephemeris_date(d) - J2000) / DAYS_IN_JULIAN_MILLENIUM;
    
    // calculating Sun's true position using semi-analytical planetary theory VSOP87
    sp = heliocentric_planetary_position(t, EARTH);
    
    return sp.distance;
}

/*
 * Computes date of the solstice or equinox (k) on a given year (y). A provided enumeration is adviced for specifying
 * equinoxes and solstices.
 * Function uses iterative approach with result accuracy up to seconds.
 */
static date equinox_solstice(int y, int k)
{
    date d;                     // result calendar date of solstice/equinox
    ecliptic_point sp;          // apparent position of the Sun on solstice/equinox
    double jde;                 // Julian ephemeris date of solstice/equinox
    double c;                   // correction to Julian ephemeris date of solstice/equinox on each iteration
    
    // required order of precision of correction to Julian ephemeris date of equinox/solstice that guarantees
    // computational precision up to seconds; since 1ᵈ = 86400ˢ, then in order to have computational accuracy up to
    // seconds, correction must be order of 1 / 86400 = 0.0000115, i.e. order of 10⁻⁷
    const double p = 10e-7;
    
    // initializing with very rough approximation of equinox/solstice date; equinoxes occur around 21 March/September
    // and solstices around 21 June/December
    d.day = 21.0;
    d.month = (k + 1) * 3;
    d.year = y;
    
    // computing starting rough Julian ephemeris date of equinox/solstice
    jde = julian_ephemeris_date(d);
    
    // on soltice/equinox apparent longitude of the Sun (that is including the effects of aberration and nutattion)
    // must be multiple of 90°; to find this instance an iterative approach is used which starts at approximate time
    // of solstice/equinox and computes Sun's apparent position; if needed accuracy is not reached, time instant is
    // corrected using the formula given in source material;
    
    // iteratively computing Julian ephemeris date of equinox/solstice
    do {
        // computing apparent position of the Sun on a given date
        sp = sun_apparent_position(d);
        // computing correction to Julian ephemeris date of solstice/equinox
        c = 58.0 * sin(k * M_PI / 2 - sp.longitude);
        // correcting Julian ephemeris date of the solstice/equinox
        jde += c;
        // computing calendar date from Julian ephemeris date
        d = calendar_date(jde - dynamical_time_difference(d));
    }
    // iteration continues while correction for Julian ephemeris date is higher than 10⁻⁷
    while (c > p);
    
    return d;
}

date equinox(int y, Equinox e)
{
    date d;                     // result date of equinox
    
    // checking whether given value is really equinox and not arbitrary integer
    if (e == VERNAL_EQUINOX || e == AUTUMNAL_EQUINOX)
        // computing given equinox
        d = equinox_solstice(y, (int)e);
    else {
        // creating invalid date structure
        d.day = -1.0;
        d.month = UNKNOWN_MONTH;
        d.year = y;
    }

    return d;
}

date solstice(int y, Solstice s)
{
    date d;                     // result date of solstice
    
    // checking whether given value is really solstice and not arbitrary integer
    if (s == SUMMER_SOLSTICE || s == WINTER_SOLSTICE)
        // computing given solstice
        d = equinox_solstice(y, (int)s);
    else {
        // creating invalid date structure
        d.day = -1.0;
        d.month = UNKNOWN_MONTH;
        d.year = y;
    }
    
    return d;   
}

double solve_equation_of_time(date d)
{
    double eoe[TOTAL_ORBITAL_ELEMENTS];     // mean orbital elements of the Earth
    double l;                   // Sun's mean longitude
    ecliptic_point ap;          // Sun's apparent position in ecliptic coordiantes
    equatorial_point ep;        // Sun's apparent position in equatorial coordinates
    double n;                   // value of nutation in longitude
    double e;                   // obliquity of the ecliptic
    double r;                   // result value of equation of time
    
    // computing mean orbital elements of the earth referred to the equinox of date
    compute_orbital_elements_of_date(d, EARTH, eoe);
    // Sun's geocentric mean longitude is Earth's heliocentric mean longitude plus π
    l = eoe[PLANET_MEAN_LONGITUDE] + M_PI;
    
    // computing nutation in longitude
    n = nutation_in_longitude(d);
    // computing obliquity of the ecliptic corrected for nutation
    e = obliquity_of_ecliptic(d) + n;
    
    // computing apparent position of the Sun
    ap = sun_apparent_position(d);
    // converting apparent position of the Sun from ecliptic coordinates to equatorial
    ep = ecliptic_to_equatorial(ap, e);
    
    // computing value of equation of time measured in radians
    r = l - ep.right_ascension + n * cos(e);
    
    // shifting value in range [0, 2π]
    r = fmod(r, 2 * M_PI);
    if (r < 0.0)
        r += 2 * M_PI;
    
    // converting value of equiation of time from radians to hours (π = 12ʰ)
    r *= 12.0 / M_PI;
    
    return r;
}