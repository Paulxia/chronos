/*
 * planets.c
 * Created by Serhii Tsyba (sertsy@gmail.com) on 30.07.10.
 */

#include "planets.h"
#include "orbital.h"
#include "earth.h"
#include "sun.h"
#include <math.h>

ecliptic_point planet_true_position(date d, Planet p)
{
    double t;                   // time in Julian millenia since the beginning of the epoch J2000 till a given date
    spherical_point psp;        // planet's heliocentric true position as computed by VSOP87 theory
    spherical_point esp;        // Earth's heliocentric true position as computed by VSOP87 theory
    double x, y, z;             // rectangular coordinates of planet's geocentric true position
    ecliptic_point ep;          // result planet's geocentric true position
    double lp;                  // variable needed to calculate correction for the reference point of ecliptic coordinates
    
    if (p != EARTH){
        // computing time measured in Julian millenia since the beginning of the epoch J2000 till a given date
        t = (julian_ephemeris_date(d) - J2000) / DAYS_IN_JULIAN_MILLENIUM;
        
        // computing heliocentric position of a given planet
        psp = heliocentric_planetary_position(t, p);
        // computing heliocentric position of of Earth
        esp = heliocentric_planetary_position(t, EARTH);
        
        // computing rectangular coordinates of geocentric position of the planet from heliocentric position of the
        // planet and heliocentric position of the Earth
        x = psp.distance * cos(psp.latitude) * cos(psp.longitude) - esp.distance * cos(esp.latitude) * cos(esp.longitude);
        y = psp.distance * cos(psp.latitude) * sin(psp.longitude) - esp.distance * cos(esp.latitude) * sin(esp.longitude);
        z = psp.distance * sin(psp.latitude) - esp.distance * sin(esp.latitude);
        
        // converting geocentric position of the planet from rectangular coordinates to ecliptic
        ep.longitude = atan2(y, x);
        ep.latitude = atan2(z, sqrt(x * x + y * y));
        
        // since VSOP87 reference frame differs slightly from FK5, computing correction to ecliptic coordinates
        lp = ep.longitude - 1.397 * t * M_PI / 180.0 - 0.00031 * t * t * M_PI / 180.0;
        ep.longitude += (-0.09033 + 0.03916 * (cos(lp) + sin(lp) * tan(ep.latitude))) * M_PI / 64800.0;
        ep.latitude += (0.03916 * (cos(lp) - sin(lp))) * M_PI / 64800.0;
        
        // shifting longitude to range [0, 2π)
        ep.longitude = fmod(ep.longitude, 2 * M_PI);
        if (ep.longitude < 0.0)
            ep.longitude += 2 * M_PI;
    }
    // if Earth is given as a planet to compute geocentric position of, then the only meaningful position may have
    // coordinates λ = 0, β = 0
    else {
        ep.longitude = 0.0;
        ep.latitude = 0.0;
    }
    
    return ep;
}

ecliptic_point planet_apparent_position(date d, Planet p)
{
    double t;                   // time in Julian millenia since the beginning of the epoch J2000 till a given date
    double lt0, lt1;            // values of light-time correction during iterative computation of this effect
    spherical_point psp;        // planet's heliocentric true position as computed by VSOP87 theory
    spherical_point esp;        // Earth's heliocentric true position as computed by VSOP87 theory
    double x, y, z;             // rectangular coordinates of planet's geocentric true position
    double ds;                  // distance from a given planet to the Earth (Δ)
    double lp;                  // variable needed to calculate correction for the reference point of ecliptic coordinates
    ecliptic_point a;           // correction to planet's position due to effect of aberration
    ecliptic_point ep;          // result planet's geocentric apparent position
    
    // required order of precision when computing correction to the light-time effect
    const double pr = 10e-7;
    
    if (p != EARTH) {
        // computing time measured in Julian centuries since the beginning of the epoch J2000
        t = (julian_ephemeris_date(d) - J2000) / DAYS_IN_JULIAN_CENTURY;
        
        // VSOP87 calculations require time to be measured in Julian millenia instead of centuries, hence there is a
        // division by ten in the 
        
        // computing heliocentric position of Earth
        esp = heliocentric_planetary_position(t / 10.0, EARTH);
        
        // A note on light-time correction. An apparent position of the planet at time t is its true position a time
        // t - τ, where τ is the time taken by light to reach the surface of the Earth from the surface of a given
        // planet. This time can be computed by the following formula:
        //
        //                                  τ = 0.0057755183Δ
        //
        // where Δ is the apparent distance between planet and the Earth. Since the apparent distance Δ and time τ are
        // not known in advance, they can be computed iteratively starting with τ = 0 (Δ = 0) until it converges to a
        // certain result within given accuracy.
        
        // lt0 is the value of light-time computed on the previous step and lt1 is the newly computed value of light
        // time on the current step; iterations continue until their difference in absolute value is not larger than
        // given precision; we start with lt0 = 0 (τ = 0) and lt1 being any number such that |lt0 - lt1| > p, where p
        // is the given precision (for instance, lt1 = 2p).
        lt0 = 0.0;
        lt1 = lt0 + 2 * pr;
        
        // computing correction due to light-time effect
        while (fabs(lt0 - lt1) > pr) {
            // computing heliocentric planetary position at an instant t - τ
            psp = heliocentric_planetary_position((t - lt0 / DAYS_IN_JULIAN_CENTURY) / 10.0, p);
            
            // computing rectangular coordinates of geocentric position of the planet from heliocentric position of the
            // planet and heliocentric position of the Earth 
            x = psp.distance * cos(psp.latitude) * cos(psp.longitude) - esp.distance * cos(esp.latitude) * cos(esp.longitude);
            y = psp.distance * cos(psp.latitude) * sin(psp.longitude) - esp.distance * cos(esp.latitude) * sin(esp.longitude);
            z = psp.distance * sin(psp.latitude) - esp.distance * sin(esp.latitude);
            
            // computing apparent distance from the planet to the Earth
            ds = sqrt(x * x + y * y + z * z);
            
            // saving previous light-time correction value
            lt1 = lt0;
            // computing new light-time correction value
            lt0 = 0.0057755183 * ds;
        }
        
        // converting geocentric position of the planet from rectangular coordinates to ecliptic
        ep.longitude = atan2(y, x);
        ep.latitude = atan2(z, sqrt(x * x + y * y));
        
        // since VSOP87 reference frame differs slightly from FK5, computing correction to ecliptic coordinates
        lp = ep.longitude - 1.397 * t * M_PI / 180.0 - 0.00031 * t * t * M_PI / 180.0;
        ep.longitude += (-0.09033 + 0.03916 * (cos(lp) + sin(lp) * tan(ep.latitude))) * M_PI / 64800.0;
        ep.latitude += (0.03916 * (cos(lp) - sin(lp))) * M_PI / 64800.0;
        
        // computing correction for aberration
        a = aberration(d, ep);
        
        // performing correction for aberration
        ep.longitude += a.longitude;
        ep.latitude += a.latitude;
        
        // performing correction for nutation in longitude
        ep.longitude += nutation_in_longitude(d);
        
        // shifting longitude to range [0, 2π)
        ep.longitude = fmod(ep.longitude, 2 * M_PI);
        if (ep.longitude < 0.0)
            ep.longitude += 2 * M_PI;
    }
    else {
        // if Earth is given as a planet to compute geocentric position of, then the only meaningful position may have
        // coordinates λ = 0, β = 0
        ep.longitude = 0.0;
        ep.latitude = 0.0;
    }
    
    return ep;
}

double planet_distance_to_sun(date d, Planet p)
{
    double t;                   // time in Julia millenia since beginning of the epoch J2000 til a given date
    spherical_point sp;         // true heliocentric position of a given planet
    
    // computing time in Julian millenia since the beginning of the epoch J2000 till a given date
    t = (julian_ephemeris_date(d) - J2000) / DAYS_IN_JULIAN_MILLENIUM;
    
    // computing true heliocentric position of the planet
    sp = heliocentric_planetary_position(t, p);
    
    return sp.distance;
}

double planet_distance_to_earth(date d, Planet p)
{
    double t;                   // time in Julia millenia since beginning of the epoch J2000 til a given date
    spherical_point psp;        // true heliocentric position of a given planet
    spherical_point esp;        // true heliocentric position of the Earth
    double x, y, z;             // rectangular coordinates of planet's geocentric true position
    double ds;                  // result distance from a given planet to the Earth
    
    if (p != EARTH){
        
        t = (julian_ephemeris_date(d) - J2000) / DAYS_IN_JULIAN_MILLENIUM;
        
        // computing heliocentric true position of Earth
        esp = heliocentric_planetary_position(t, EARTH);
        // computing heliocentric true position of a given planet
        psp = heliocentric_planetary_position(t, p);
        
        // computing rectangular coordinates of geocentric position of the planet from heliocentric position of the
        // planet and heliocentric position of the Earth
        x = psp.distance * cos(psp.latitude) * cos(psp.longitude) - esp.distance * cos(esp.latitude) * cos(esp.longitude);
        y = psp.distance * cos(psp.latitude) * sin(psp.longitude) - esp.distance * cos(esp.latitude) * sin(esp.longitude);
        z = psp.distance * sin(psp.latitude) - esp.distance * sin(esp.latitude);
        
        // computing distance from the planet to the Earth
        ds = sqrt(x * x + y * y + z * z);
    }
    else
        // if Earth is given as a planet to compute distance to, then only meaningful value can be 0
        ds = 0.0;
    
    return ds;
}

double planet_phase_angle(date d, Planet p)
{
    double pds;                 // distance from the Sun to a given planet
    double eds;                 // distance from the Earth to the Sun
    double pde;                 // distance from the Earth to a given planet
    double i;                   // result value of phase angle
    
    if (p != EARTH) {
        // computing distance from a given planet to the Sun
        pds = planet_distance_to_sun(d, p);
        // computing distance from the Earth to the Sun
        eds = sun_distance_to_earth(d);
        // computing distance from a given planet to th Earth
        pde = planet_distance_to_earth(d, p);
        
        // computing phase angle
        i = acos((pds * pds + pde * pde - eds * eds) / (2.0 * pds * pde));
    }
    else
        // if Earth is given as a planet to compute phase angle of, assigning error value
        i = -1;
    
    return i;
}

double planet_disk_illuminated_fraction(date d, Planet p)
{
    double i;                   // pahse angle of the planet
    double k;                   // result illuminated fraction of the planet disk
    
    if (p != EARTH) {
        // computing phase angle of a given planet
        i = planet_phase_angle(d, p);
        // computing illuminated fraction of the planet's disk
        k = (1.0 + cos(i)) / 2.0;
    }
    else
        // if Earth is given as a planet whoe illuminated fraction should be computed, assigning error value
        k = -1;
    
    return k;
}

/*
 * Computes saturnicentric position of the Earth referred to the plane of the ring at a givev date (d). This function
 * used later to determine apparent stellar magnitude of Saturn as seen from the Earth.
 *
 * Source: Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 301.
 */
static ecliptic_point saturnicentric_earth_position(date d)
{
    double t;                   // time in Julian centuries since the beginning of the epoch J2000 till a given date
    double i;                   // inclination of the plane of the ring(s) of Saturn (i)
    double o;                   // longitude of the ascending node referred to the mean equinox of date (☊)
    ecliptic_point gsp;         // geocentric apparent position of Saturn
    ecliptic_point sep;         // result satunicentric position of the Earth
    
    // computing time measured in Julian centuries since the beginning of the epoch J2000
    t = (julian_ephemeris_date(d) - J2000) / DAYS_IN_JULIAN_CENTURY;
    
    // computing inclination of the plane of the ring measured in degrees
    i = 28.075216 - 0.012998 * t + 0.000004 * t * t;
    // computing longitude of the ascending node measured in degrees
    o = 169.508470 + 1.394681 * t + 0.000412 * t * t;

    // converting inclination of the plane of the ring from degrees to radians (π = 180°)
    i *= M_PI / 180.0;
    // converting longitude of the ascending node from degrees to radians (π = 180°)
    o *= M_PI / 180.0;
    
    // computing geocentric apparent position of Saturn
    gsp = planet_apparent_position(d, SATURN);
    
    // computing saturnicentric position of the Earth referred to the plane of the ring from geocentric position of Saturn
    sep.longitude = atan2(sin(i) * sin(gsp.latitude) + cos(i) * cos(gsp.latitude) * sin(gsp.longitude - o),
                          cos(gsp.latitude) * cos(gsp.longitude - o));
    sep.latitude = asin(sin(i) * cos(gsp.latitude) * sin(gsp.longitude - o) - cos(i) * sin(gsp.latitude));
    
    // shifting longitude to interval [0, 2π)
    sep.longitude = fmod(sep.longitude, 2 * M_PI);
    if (sep.longitude < 0.0)
        sep.longitude += 2 * M_PI;
    
    return sep;
}

/*
 * Computes saturnicentric position of the Sun referred to the plane of the ring at a givev date (d). This function
 * used later to determine apparent stellar magnitude of Saturn as seen from the Earth.
 *
 * Source: Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 301.
 */
static ecliptic_point saturnicentric_sun_position(date d)
{
    double t;                   // time in Julian centuries since the beginning of the epoch J2000 till a given date
    double lt0, lt1;            // values of light-time correction during iterative computation of this effect
    spherical_point esp;        // heliocentric true position of the Earth
    spherical_point ssp;        // heliocentric true position of Saturn
    double x, y, z;             // rectangular coordinate of geocentric position of Saturn
    double ds;                  // apparent distance from Earth to Saturn (Δ)
    double i;                   // inclination of the plane of the ring(s) of Saturn (i)
    double o;                   // longitude of the ascending node referred to the plane of the ring(s) (☊)
    double n;                   // longitude of the ascending node of Saturn's orbit
    ecliptic_point ep;          // result saturnicentric position of the Sun
    
    // required order of precision when computing correction to the light-time effect
    const double pr = 10e-7;
    
    // computing time measured in Julian centuries since the beginning of the epoch J2000
    t = (julian_ephemeris_date(d) - J2000) / DAYS_IN_JULIAN_CENTURY;
    
    // VSOP87 calculations require time to be measured in Julian millenia instead of centuries, hence there is a
    // division by ten in the 
    
    // computing heliocentric position of Earth
    esp = heliocentric_planetary_position(t / 10.0, EARTH);
    
    // Light-correction is performed the way it is done when computing geocentric position of the planet. Refer to
    // the explanation in that function for more information on how light-time correction is performed
    
    // initializing light-time correction values
    lt0 = 0.0;
    lt1 = lt0 + 2 * pr;
    
    // computing correction due to light-time effect
    while (fabs(lt0 - lt1) > pr) {
        // computing heliocentric position of Saturn with current correction for light-time
        ssp = heliocentric_planetary_position((t - lt0 / DAYS_IN_JULIAN_CENTURY) / 10.0, SATURN);
        
        // computing rectangular coordinates of geocentric position of Saturn from heliocentric position of Saturn and
        // heliocentric position of the Earth 
        x = ssp.distance * cos(ssp.latitude) * cos(ssp.longitude) - esp.distance * cos(esp.latitude) * cos(esp.longitude);
        y = ssp.distance * cos(ssp.latitude) * sin(ssp.longitude) - esp.distance * cos(esp.latitude) * sin(esp.longitude);
        z = ssp.distance * sin(ssp.latitude) - esp.distance * sin(esp.latitude);
        
        // computing apparent distance from the planet to the Earth
        ds = sqrt(x * x + y * y + z * z);
        
        // saving previous light-time correction value
        lt1 = lt0;
        // computing new light-time correction value
        lt0 = 0.0057755183 * ds;
    }
    
    // computing inclination of the plane of the ring measured in degrees
    i = 28.075216 - 0.012998 * t + 0.000004 * t * t;
    // computing longitude of the ascending node measured in degrees
    o = 169.508470 + 1.394681 * t + 0.000412 * t * t;
    // computing longitude of the ascending node of Saturn's orbit measured in degrees
    n = 113.6655 + 0.8771 * t;
    
    // converting inclination of the plane of the ring from degrees to radians (π = 180°)
    i *= M_PI / 180.0;
    // converting longitude of the ascending node from degrees to radians (π = 180°)
    o *= M_PI / 180.0;
    // converting longitude of the ascending node of Saturn's orbit from degrees to radians (π = 180°)
    n *= M_PI / 180.0;
    
    // computing correction due to effect of aberration of the Sun as seen from Saturn
    ssp.longitude -= (0.01759 / ssp.distance) * M_PI / 180.0;
    ssp.latitude -= (0.000764 * cos(ssp.longitude - n) / ssp.distance) * M_PI / 180.0;
    
    // computing saturnicentric position of the Sun
    ep.longitude = atan2(sin(i) * sin(ssp.latitude) + cos(i) * cos(ssp.latitude) * sin(ssp.longitude - o),
                          cos(ssp.latitude) * cos(ssp.longitude - o));
    ep.latitude = asin(sin(i) * cos(ssp.latitude) * sin(ssp.longitude - o) - cos(i) * sin(ssp.latitude));
    
    // shifting longitude to interval [0, 2π)
    ep.longitude = fmod(ep.longitude, 2 * M_PI);
    if (ep.longitude < 0.0)
        ep.longitude += 2 * M_PI;
    
    return ep;
}

double planet_apparent_magnitude(date d, Planet p)
{
    double pds;                 // distance from the Sun to a given planet
    double pde;                 // distance from the Earth to a given planet
    double i;                   // illuminated fraction of the disk of a given planet
    ecliptic_point sep;         // saturnicentric apparent position of the Earth
    ecliptic_point ssp;         // saturnicentric apparent position of the Sun
    double b;                   // saturnicentric latitude of the Earth
    double du;                  // difference between saturnicentric longitudes of the Sun and the Earth
    double m;                   // result apparent stellar magnitude of a given planet
    
    if (p != EARTH) {
        // computing distance from the Sun to a given planet
        pds = planet_distance_to_sun(d, p);
        // computing distance from the Earth to a given planet
        pde = planet_distance_to_earth(d, p);
        // computing planet's phase angle
        i = planet_phase_angle(d, p);
        // converting value of the planet's phase angle from radians to degrees (π = 180°), since formulas for stellar
        // magnitude of a planet operates on the value expressed in degrees
        i *= 180.0 / M_PI;
        
        // computing apparent stellar magnitude of a given planet
        switch (p) {
            case MERCURY:
                m = -0.42 + 5 * log10(pds * pde) + 0.0380 * i - 0.000273 * i * i + 0.000002 * i * i * i;
                break;
            case VENUS:
                m = -4.40 + 5 * log10(pds * pde) + 0.0009 * i + 0.000239 * i * i - 0.00000065 * i * i * i;
                break;
            case MARS:
                m = -1.52 + 5 * log10(pds * pde) + 0.016 * i;
                break;
            case JUPITER:
                m = -9.40 + 5 * log10(pds * pde) + 0.005 * i;
                break;
            case SATURN:
                // in case of Saturn, apparent stellar magnitude depends also on the position of the ring(s), so these
                // extra computations are necessary
                
                // computing saturnicentric position of the Earth referred to the plane of the ring(s)
                sep = saturnicentric_earth_position(d);
                // computing saturnicentric position of the Sun referred to the plane of the ring(s)
                ssp = saturnicentric_sun_position(d);
                
                // computing absolute value of the saturnicentric latitude of the Earth
                b = fabs(sep.latitude);
                // computing absolute value of the difference between saturnicentric longitudes of the Earth and the Sun
                du = fabs(ssp.longitude - sep.longitude);
                
                m = -8.88 + 5 * log10(pds * pde) + 0.44 * du - 2.60 * sin(b) + 1.25 * sin(b) * sin(b);
                break;
            case URANUS:
                m = -7.19 + 5 * log10(pds * pde);
                break;
            case NEPTUNE:
                m = -6.87 + 5 * log10(pds * pde);
                break;
            default:
                break;
        }
    }
    
    // if Earth is given as a planet whose apparent stellar magnitude is to be computed, no meaningful value can be
    // given; in that case just return whatever value result variable was initialized with
    
    return m;
}