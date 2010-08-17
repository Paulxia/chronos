/*
 * coordinates.c
 * Created by Serhii Tsyba (sertsy@gmail.com) on 28.06.10.
 */

#include "coordinates.h"
#include "calendar.h"
#include <math.h>

ecliptic_point equatorial_to_ecliptic(equatorial_point eqp, double e)
{
    ecliptic_point ecp;         // result point in ecliptical coordinate system
    
    // performing conversion to ecliptic coordinate system
    ecp.longitude = atan2(sin(eqp.right_ascension) * cos(e) + tan(eqp.declination) * sin(e), cos(eqp.right_ascension));
    ecp.latitude = asin(sin(eqp.right_ascension) * cos(e) - cos(eqp.declination) * sin(e) * sin(eqp.right_ascension));
    
    return ecp;
}

equatorial_point ecliptic_to_equatorial(ecliptic_point ecp, double e)
{
    equatorial_point eqp;       // result point in equatorial coordinate system
    
    // performing conversion to equatorial coordinate system
    eqp.right_ascension = atan2(sin(ecp.longitude) * cos(e) - tan(ecp.latitude) * sin(e), cos(ecp.longitude));
    eqp.declination = asin(sin(ecp.latitude) * cos(e) + cos(ecp.latitude) * sin(e) * sin(ecp.longitude));
    
    return eqp;
}

horizontal_point equatorial_to_horizontal(date d, geographic_point gp, equatorial_point ep)
{
    double lha;                 // local hour angle of equatorial coordinate system
    double gast;                // apparent siderial time at Greenwich for a given date
    horizontal_point hp;        // result point in horizontal coordinate system
    
    // calculating mean siderial time at Greenwich (θ₀)
    gast = greenwich_mean_siderial_time(d);
    // converting mean siderial time at Greenwich to degrees (1ʰ = 15°)
    gast *= 15.0;
    // converting mean siderial time at Greenwich to radians (π = 180°)
    gast *= M_PI / 180.0;
    
    // calculating local hour angle (H)
    lha = gast - gp.longitude - ep.right_ascension;
    
    // performing conversion to horizontal coordinate system
    hp.azimuth = atan2(sin(lha), cos(lha) * sin(gp.latitude) - tan(ep.declination) * cos(gp.latitude));
    hp.elevation = asin(sin(gp.latitude) * sin(ep.declination) + cos(gp.latitude) * cos(ep.declination) * cos(lha));
    
    return hp;
}

equatorial_point horizontal_to_equatorial(date d, geographic_point gp, horizontal_point hp)
{
    double lha;                 // local hour angle of equatorial coordinate system
    double gmst;                // mean siderial time at Greenwich for a given date
    equatorial_point ep;        // result point in equatorial coordinate system
    
    // computing Greenwich mean siderial time for a given time instant
    gmst = greenwich_mean_siderial_time(d);
    
    // performing conversion to equatorial coordinate system
    // computing local hour angle of equatorial coordinate system
    lha = atan2(sin(hp.azimuth), cos(hp.azimuth) * sin(gp.latitude) + tan(hp.elevation) * cos(gp.latitude));
    
    // computing right ascension from local hour angle
    ep.right_ascension = gmst - gp.longitude - lha;
    
    // computing declination
    ep.declination = asin(sin(gp.latitude) * sin(hp.elevation) - cos(gp.latitude) * cos(hp.elevation) * cos(hp.azimuth));
    
    return ep;
}