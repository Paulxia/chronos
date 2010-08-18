/*
 * moon.c
 * Created by Serhii Tsyba on 03.08.10.
 *
 *
 * TODO: Consider moving from ELP theory version ELP2000-2B to version ELP/MPP02
 */

#include "moon.h"
#include "earth.h"
#include "sun.h"
#include <elp2000-82b/elp2000-82b.h>
#include <math.h>

ecliptic_point moon_true_position(date d)
{
    double t;                   // time in Julian centuries since the beginning of the epoch J2000 till a given date
    spherical_point sp;         // true position of the Moon as computed by semi-analytical theory ELP2000
    ecliptic_point ep;          // result true position of the Moon
    
    // computing time in Julian centuries since the beginning of the epoch J2000 till a given date
    t = (julian_ephemeris_date(d) - J2000) / DAYS_IN_JULIAN_CENTURY;
    
    // computing true position of the moon by means of semi-analytic lunar theory ELP2000
    sp = geocentric_moon_position_of_date(t);
    
    // since theory ELP2000 expresses ecliptic coordinates in arcseconds, converting from arcseconds to radians (π = 648000")
    ep.longitude = sp.longitude * M_PI / 648000;
    ep.latitude = sp.latitude * M_PI / 648000;
    
    // shifting longitude to range [0, 2π)
    ep.longitude = fmod(ep.longitude, 2 * M_PI);
    if (ep.longitude < 0.0)
        ep.longitude += 2 * M_PI;
    
    return ep;
}

ecliptic_point moon_apparent_position(date d)
{
    ecliptic_point ep;          // result apparent position of the Moon
    
    // computing true position of the Moon
    ep = moon_true_position(d);
    
    // correcting true position of the Moon for nutation to get apparent position
    ep.longitude += nutation_in_longitude(d);
    
    // even thought correction for nutation is small, it may still shift longitude out of desired interval [0, 2π);
    // making sure, that longitude to range [0, 2π)
    ep.longitude = fmod(ep.longitude, 2 * M_PI);
    if (ep.longitude < 0.0)
        ep.longitude += 2 * M_PI;
    
    return ep;
}

double moon_distance_to_earth(date d)
{
    double t;                   // time in Julian centuries since the beginning of the epoch J2000 till a given date
    spherical_point sp;         // true position of the Moon as computed by means of semi-analytic theory ELP2000
    double ds;                  // result disatnce from the Earth to the Moon
    
    // computing time in Julian centuries since the beginning of the epoch J2000 till a given date
    t = (julian_ephemeris_date(d) - J2000) / DAYS_IN_JULIAN_CENTURY;
    
    // computing true position of the moon by means of semi-analytic lunar theory ELP2000
    sp = geocentric_moon_position_of_date(t);
    
    // since theory ELP2000 expresses distance in kilometers, converting from kilometers to astronomical units (1 AU = 149597871 km)
    ds = sp.distance / ASTRONOMICAL_UNIT;
    
    return ds;
}

double moon_phase_angle(date d)
{
    ecliptic_point mep;         // geocentric apparent position of the Moon
    ecliptic_point sep;         // geocentric apparent position of the Sun
    double des;                 // distance from the Earth to the Sun
    double dem;                 // distance from the Earth to the Moon
    double p;                   // geocentric elongation of the Moon from the Sun (ψ)
    double i;                   // result phase angle of the Moon
    
    // computing geocentric apparent position of the Moon
    mep = moon_apparent_position(d);
    // computing geocentric apparent position of the Sun
    sep = sun_apparent_position(d);
    
    // computing geocentric elongation of the Moon from the Sun
    p = acos(cos(mep.latitude) * cos(mep.longitude - sep.longitude));
    
    // computing distance from the Earth to the Sun
    des = sun_distance_to_earth(d);
    // computing distance from the Earth to the Moon
    dem = moon_distance_to_earth(d);
    
    // computing phase angle of the Moon
    i = atan2(des * sin(p), dem - des * cos(p));
    
    return i;
}

double moon_disk_illuminated_fraction(date d)
{
    double i;                   // phase angle of the Moon
    double k;                   // result illuminated fraction of the Moon's disk
    
    // computing phase angle of the Moon
    i = moon_phase_angle(d);
    
    // computing illuminated fraction of the disk of the Moon
    k = (1.0 + cos(i)) / 2;
    
    return k;
}

double moon_bright_limb_position_angle(date d)
{
    ecliptic_point sp;          // geocentric apparent position of the Sun
    ecliptic_point mp;          // geocentric apparent position of the Moon
    equatorial_point sep;       // geocentric apparent position of the Sun expressed in equatorial coordinates
    equatorial_point mep;       // geocentric apparent position of the Moon expressed in equatorial coordinates
    double e;                   // true obliquity of the ecliptic
    double x;                   // result position angle of the Moon's bright limb
    
    // computing apparent position of the Sun expressed in ecliptic coordinates
    sp = sun_apparent_position(d);
    // computing apparent position of the Moon expressed in ecliptic coordinates
    mp = moon_apparent_position(d);
    
    // computing obliquity of the ecliptic, corrected for the effect of nutation
    e = obliquity_of_ecliptic(d) + nutation_in_obliquity(d);
    
    // converting position of the Sun from ecliptic coordinates to equatorial
    sep = ecliptic_to_equatorial(sp, e);
    // converting position of the Moon from ecliptic coordinates to equatorial
    mep = ecliptic_to_equatorial(mp, e);
    
    // computing position angle of the bright limb of the Moon
    x = atan2(cos(sep.declination) * sin(sep.right_ascension - mep.right_ascension),
              sin(sep.declination) * cos(mep.declination) - cos(sep.declination) * sin(mep.declination) * cos(sep.right_ascension - mep.right_ascension));
    
    // shifting result value to range [0, 2π)
    if (x < 0.0)
        x += 2 * M_PI;
    
    return x;
}