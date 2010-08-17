/*
 * geoobserver.c
 * Created by Serhii Tsyba (sertsy@gmail.com) on 19.07.10.
 *
 *
 * TODO: Consider a different method of determining rising and setting times.
 */

#include "geoobserver.h"
#include "earth.h"

#include <math.h>

double parallactic_angle(date d, geographic_point gp, equatorial_point ep)
{
    double gast;                // apparent siderial time at Greenwich (θ₀)
    double laha;                // local hour angle (H)
    double pa;                  // result parallactic angle (q)
    
    // calculating mean siderial time at Greenwich (θ₀)
    gast = greenwich_mean_siderial_time(d);
    // converting mean siderial time at Greenwich from hours to radians (π = 12ʰ)
    gast *= M_PI / 12.0;
    // calculatin apparent siderial time at Greenwich from mean
    gast += nutation_in_longitude(d) * cos(obliquity_of_ecliptic(d));
    
    // calculating local hour angle (H)
    laha = gast - gp.longitude - ep.right_ascension;
    
    // calculating parallactic angle (q)
    pa = sin(laha) / (tan(gp.latitude) * cos(ep.declination) - sin(ep.declination) * cos(laha));
    
    return pa;
}

static double transit(date d, geographic_point gp, equatorial_point ep, double sa)
{
    date d0;                    // date corresponding to 0ʰ of a given date
    double gast;                // apparent siderial time at Greenwich
    double m;                   // result time of rising
    
    // computing date corresponding to 0ʰ of a given date
    modf(d.day, &d0.day);
    d0.month = d.month;
    d0.year = d.year;
    
    // calculating apparent siderial time at Greenwich corresponding to 0ʰ of a given day
    // calculating mean siderial time at Greenwich
    gast = greenwich_mean_siderial_time(d0);
    // converting mean siderial time at Greenwich from hours to radians (π = 12ʰ)
    gast *= M_PI / 12.0;
    // converting mean siderial time at Greenwich to apparent
    gast += nutation_in_longitude(d) * cos(obliquity_of_ecliptic(d));
    
    // computing time of transit expressed in days (fraction of day)
    m = (ep.right_ascension + gp.longitude - gast) / (2 * M_PI);
    
    // shifting value to interval [0, 1] so that rise time is in interval [0ʰ, 24ʰ]
    m = fmod(m, 1.0);
    if (m < 0.0)
        m += 1.0;
    
    // converting days into hours (1ᵈ = 24ʰ)
    m *= 24.0;
    
    return m;
}

double rising(date d, geographic_point gp, equatorial_point ep, double sa)
{
    double cosh0;               // auxiliary computational variable
    double h0;                  // auxiliary computational variable
    double m;                   // result time of rising
    
    // performing auxiliary computations
    cosh0 = (sin(sa) - sin(gp.latitude) * sin(ep.declination)) / (cos(gp.latitude) * cos(ep.declination));
    
    // checking whether body rises at all; 
    if (fabs(cosh0) > 1)
        // body never rises;
        return -1;
    else {
        // body rises
        
        // performing auxiliary computations
        h0 = acos(cosh0);
        
        // computing rising time
        m = transit(d, gp, ep, sa) - 12 * h0 / M_PI;
        
        // fitting rising time into interval [0ʰ, 24ʰ]
        if (m < 0.0)
            m += 24.0;
        
        return m;
    }
}

double setting(date d, geographic_point gp, equatorial_point ep, double sa)
{
    double cosh0;               // auxiliary computational variable
    double h0;                  // auxiliary computational variable
    double m;                   // result time of setting
    
    // performing auxiliary computations
    cosh0 = (sin(sa) - sin(gp.latitude) * sin(ep.declination)) / (cos(gp.latitude) * cos(ep.declination));
    
    // checking whether body sets at all; 
    if (fabs(cosh0) > 1)
        // body never sets;
        return -1;
    else {
        // body sets
        
        // performing auxiliary computations
        h0 = acos(cosh0);
        
        // computing setting time
        m = transit(d, gp, ep, sa) + 12 * h0 / M_PI;
        
        // fitting setting time into interval [0ʰ, 24ʰ]
        if (m > 24.0)
            m -= 24.0;
        
        return m;
    }
}

double atmospheric_refraction(double a, double t, double p)
{
    double r;                   // result value of atmospheric refraction
    
    // computing value of atmospheric refraction
    r = 1.02 / (tan(a + 10.3 / (a + 5.11)));
    
    // performing correction of the value of atmospheric refraction to environmental conditions
    r *= (p / 101325.0) * (283.15 / t);
    
    return r;
}

equatorial_point diurnal_parallax(date d, geographic_point gp, double a, equatorial_point ep, double ehp)
{
    double u, s, c;             // auxiliary computational variables
    double da;                  // correction to right ascension due to diurnal parallax
    double gast;                // apparent siderial time at Greenwich
    double h;                   // geographic hour angle
    equatorial_point epp;       // result equatorial point corrected for diurnal parallax
    
    // performing auxiliary computations
    u = atan(EARTH_POLAR_RADIUS / EARTH_EQUATORIAL_RADIUS * tan(gp.latitude));
    s = EARTH_POLAR_RADIUS / EARTH_EQUATORIAL_RADIUS * sin(u) + a / EARTH_EQUATORIAL_RADIUS * sin(gp.latitude);
    c = cos(u) + a / EARTH_EQUATORIAL_RADIUS * cos(gp.latitude);
    
    // computing apparent siderial time at Greenwich
    // computing mean siderial time at Greenwich expressed in hours
    gast = greenwich_mean_siderial_time(d);
    // converting mean siderial time at Greenwich from hours to radians (π = 12ʰ)
    gast *= M_PI / 12.0;
    // computing apparent siderial time at Greenwich expressed in radians
    gast += nutation_in_longitude(d) * cos(obliquity_of_ecliptic(d) + nutation_in_obliquity(d));
    
    // computing geocentric hour angle
    h = gast - gp.longitude - ep.right_ascension;
    
    // computing correction to the right ascension of the given point due to diurnal parallax
    da = atan2(-c * sin(ehp) * sin(h), cos(ep.declination) - c * sin(ehp) * cos(h));
    // computing right ascension of the given point corrected for parallax
    epp.right_ascension = ep.right_ascension + da;
    // computing declination of the given point corrected for parallax
    epp.declination = atan2((sin(ep.declination) - s * sin(ehp)) * cos(da), cos(ep.declination) - c * sin(ehp) * cos(h));
    
    return epp;
}