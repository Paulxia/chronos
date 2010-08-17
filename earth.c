/*
 * earth.c
 * Created by Serhii Tsyba (sertsy@gmail.com) on 24.06.10.
 */

#include "earth.h"
#include "sun.h"
#include "orbital.h"
#include <elp2000-82b/arguments.h>
#include <math.h>

/*
 * The following definition and two arrays define coefficients of the series for computation of the astronomical
 * nutation of the Earth's rotation axis adopted by IAU as 1980 IAU Theory of Nutation.
 * According to the Theory of Nutation, nutation in longitude (Δψ) and obliquity (Δε) are computed using following
 * expression
 *
 *                                      Δψ (Δε) = Σ (a + bt)cosφ
 *
 * where
 *
 *                                  φ = i₁l + i₂l' + i₃F + i₄D + i₅☊
 *
 * where l, l', F and D are Delaunay arguments and ☊ is the longitude of lunar ascending node referred to the mean
 * equinox of date; t is the time interval measured in Julian centuries since the beginning of the epoch J2000 till a
 * given date; multipliers iᵢ and coefficients a and b are constants provided in Theory of Nutation.
 *
 * Multipliers iᵢ are given in order specified above. Coefficients a and b are indexed at positions 2 and 3 (1 and 2
 * for 0-based indexing) for nutation in longitude and 4 and 5 (3 and 4) for nutation in obliquity. These coefficients
 * are measured in 10⁻⁵" (arcseconds). Coefficient at position 1 (0) is period measured in days, which is not used in
 * computations but is kept to preserve Theory of Nutation data structure.
 *
 * Source: P.K. Seidelman. 1980 IAU Theory of Nutation: The Final Report of the IAU Working Group on Nutation.
 *         U.S. Naval Observatory, Nautical Almanac Office, Washington, D.C. 20390, U.S.A., 1981.
 */

#define TOTAL_NUTATION_TERMS 106                                          // total amount of terms in series

static double nutation_multipliers[TOTAL_NUTATION_TERMS][5] = {
    {0, 0, 0, 0, 1},
    {0, 0, 2, -2, 2},
    {0, 0, 2, 0, 2},
    {0, 0, 0, 0, 2},
    {0, -1, 0, 0, 0},
    {1, 0, 0, 0, 0},
    {0, 1, 2, -2, 2},
    {0, 0, 2, 0, 1},
    {1, 0, 2, 0, 2},
    {0, -1, 2, -2, 2},
    {-1, 0, 0, 2, 0},
    {0, 0, 2, -2, 1},
    {-1, 0, 2, 0, 2},
    {1, 0, 0, 0, 1},
    {0, 0, 0, 2, 0},
    {-1, 0, 2, 2, 2},
    {-1, 0, 0, 0, 1},
    {1, 0, 2, 0, 1},
    {-2, 0, 0, 2, 0},
    {-2, 0, 2, 0, 1},
    {0, 0, 2, 2, 2},
    {2, 0, 2, 0, 2},
    {2, 0, 0, 0, 0},
    {1, 0, 2, -2, 2},
    {0, 0, 2, 0, 0},
    {0, 0, 2, -2, 0},
    {-1, 0, 2, 0, 1},
    {0, 2, 0, 0, 0},
    {0, 2, 2, -2, 2},
    {-1, 0, 0, 2, 1},
    {0, 1, 0, 0, 1},
    {1, 0, 0, -2, 1},
    {0, -1, 0, 0, 1},
    {2, 0, -2, 0, 0},
    {-1, 0, 2, 2, 1},
    {1, 0, 2, 2, 2},
    {0, -1, 2, 0, 2},
    {0, 0, 2, 2, 1},
    {1, 1, 0, -2, 0},
    {0, 1, 2, 0, 2},
    {-2, 0, 0, 2, 1},
    {0, 0, 0, 2, 1},
    {2, 0, 2, -2, 2},
    {1, 0, 0, 2, 0},
    {1, 0, 2, -2, 1},
    {0, 0, 0, -2, 1},
    {0, -1, 2, -2, 1},
    {2, 0, 2, 0, 1},
    {1, -1, 0, 0, 0},
    {1, 0, 0, -1, 0},
    {0, 0, 0, 1, 0},
    {0, 1, 0, -2, 0},
    {1, 0, -2, 0, 0},
    {2, 0, 0, -2, 1},
    {0, 1, 2, -2, 1},
    {1, 1, 0, 0, 0},
    {1, -1, 0, -1, 0},
    {-1, -1, 2, 2, 2},
    {0, -1, 2, 2, 2},
    {1, -1, 2, 0, 2},
    {3, 0, 2, 0, 2},
    {-2, 0, 2, 0, 2},
    {1, 0, 2, 0, 0},
    {-1, 0, 2, 4, 2},
    {1, 0, 0, 0, 2},
    {-1, 0, 2, -2, 1},
    {0, -2, 2, -2, 1},
    {-2, 0, 0, 0, 1},
    {2, 0, 0, 0, 1},
    {3, 0, 0, 0, 0},
    {1, 1, 2, 0, 2},
    {0, 0, 2, 1, 2},
    {1, 0, 0, 2, 1},
    {1, 0, 2, 2, 1},
    {1, 1, 0, -2, 1},
    {0, 1, 0, 2, 0},
    {0, 1, 2, -2, 0},
    {0, 1, -2, 2, 0},
    {1, 0, -2, 2, 0},
    {1, 0, -2, -2, 0},
    {1, 0, 2, -2, 0},
    {1, 0, 0, -4, 0},
    {2, 0, 0, -4, 0},
    {0, 0, 2, 4, 2},
    {0, 0, 2, -1, 2},
    {-2, 0, 2, 4, 2},
    {2, 0, 2, 2, 2},
    {0, -1, 2, 0, 1},
    {0, 0, -2, 0, 1},
    {0, 0, 4, -2, 2},
    {0, 1, 0, 0, 2},
    {1, 1, 2, -2, 2},
    {3, 0, 2, -2, 2},
    {-2, 0, 2, 2, 2},
    {-1, 0, 0, 0, 2},
    {0, 0, -2, 2, 1},
    {0, 1, 2, 0, 1},
    {-1, 0, 4, 0, 2},
    {2, 1, 0, -2, 0},
    {2, 0, 0, 2, 0},
    {2, 0, 2, -2, 1},
    {2, 0, -2, 0, 1},
    {1, -1, 0, -2, 0},
    {-1, 0, 0, 1, 1},
    {-1, -1, 0, 2, 1},
    {0, 1, 0, 1, 0}
};      // mutlipliers of arguments
static double nutation_coefficients[TOTAL_NUTATION_TERMS][5] = {
    {-6798.4, -171996.0, -174.2, 92025.0, 8.9},
    {182.6, -13187.0, -1.6, 5736.0, -3.1},
    {13.7, -2274.0, -0.2, 977.0, -0.5},
    {-3399.2, 2062.0, 0.2, -895.0, 0.5},
    {-365.3, -1426.0, 3.4, 54.0, -0.1},
    {27.6, 712.0, 0.1, -7.0, 0.0},
    {121.7, -517.0, 1.2, 224.0, -0.6},
    {13.6, -386.0, -0.4, 200.0, 0.0},
    {9.1, -301.0, 0.0, 129.0, -0.1},
    {365.2, 217.0, -0.5, -95.0, 0.3},
    {31.8, 158.0, 0.0, -1.0, 0.0},
    {177.8, 129.0, 0.1, -70.0, 0.0},
    {27.1, 123.0, 0.0, -53.0, 0.0},
    {27.7, 63.0, 0.1, -33.0, 0.0},
    {14.8, 63.0, 0.0, -2.0, 0.0},
    {9.6, -59.0, 0.0, 26.0, 0.0},
    {-27.4, -58.0, -0.1, 32.0, 0.0},
    {9.1, -51.0, 0.0, 27.0, 0.0},
    {-205.9, -48.0, 0.0, 1.0, 0.0},
    {1305.5, 46.0, 0.0, -24.0, 0.0},
    {7.1, -38.0, 0.0, 16.0, 0.0},
    {6.9, -31.0, 0.0, 13.0, 0.0},
    {13.8, 29.0, 0.0, -1.0, 0.0},
    {23.9, 29.0, 0.0, -12.0, 0.0},
    {13.6, 26.0, 0.0, -1.0, 0.0},
    {173.3, -22.0, 0.0, 0.0, 0.0},
    {27.0, 21.0, 0.0, -10.0, 0.0},
    {182.6, 17.0, -0.1, 0.0, 0.0},
    {91.3, -16.0, 0.1, 7.0, 0.0},
    {32.0, 16.0, 0.0, -8.0, 0.0},
    {386.0, -15.0, 0.0, 9.0, 0.0},
    {-31.7, -13.0, 0.0, 7.0, 0.0},
    {-346.6, -12.0, 0.0, 6.0, 0.0},
    {-1095.2, 11.0, 0.0, 0.0, 0.0},
    {9.5, -10.0, 0.0, 5.0, 0.0},
    {5.6, -8.0, 0.0, 3.0, 0.0},
    {14.2, -7.0, 0.0, 3.0, 0.0},
    {7.1, -7.0, 0.0, 3.0, 0.0},
    {-34.8, -7.0, 0.0, 0.0, 0.0},
    {13.2, 7.0, 0.0, -3.0, 0.0},
    {-199.8, -6.0, 0.0, 3.0, 0.0},
    {14.8, -6.0, 0.0, 3.0, 0.0},
    {12.8, 6.0, 0.0, -3.0, 0.0},
    {9.6, 6.0, 0.0, 0.0, 0.0},
    {23.9, 6.0, 0.0, -3.0, 0.0},
    {-14.7, -5.0, 0.0, 3.0, 0.0},
    {346.6, -5.0, 0.0, 3.0, 0.0},
    {6.9, -5.0, 0.0, 3.0, 0.0},
    {29.8, 5.0, 0.0, 0.0, 0.0},
    {411.8, -4.0, 0.0, 0.0, 0.0},
    {29.5, -4.0, 0.0, 0.0, 0.0},
    {-15.4, -4.0, 0.0, 0.0, 0.0},
    {-26.9, 4.0, 0.0, 0.0, 0.0},
    {212.3, 4.0, 0.0, -2.0, 0.0},
    {119.6, 4.0, 0.0, -2.0, 0.0},
    {25.6, -3.0, 0.0, 0.0, 0.0},
    {-3232.9, -3.0, 0.0, 0.0, 0.0},
    {9.8, -3.0, 0.0, 1.0, 0.0},
    {7.2, -3.0, 0.0, 1.0, 0.0},
    {9.4, -3.0, 0.0, 1.0, 0.0},
    {5.5, -3.0, 0.0, 1.0, 0.0},
    {1615.7, -3.0, 0.0, 1.0, 0.0},
    {9.1, 3.0, 0.0, 0.0, 0.0},
    {5.8, -2.0, 0.0, 1.0, 0.0},
    {27.8, -2.0, 0.0, 1.0, 0.0},
    {-32.6, -2.0, 0.0, 1.0, 0.0},
    {6786.3, -2.0, 0.0, 1.0, 0.0},
    {-13.7, -2.0, 0.0, 1.0, 0.0},
    {13.8, 2.0, 0.0, -1.0, 0.0},
    {9.2, 2.0, 0.0, 0.0, 0.0},
    {8.9, 2.0, 0.0, -1.0, 0.0},
    {9.3, 2.0, 0.0, -1.0, 0.0},
    {9.6, -1.0, 0.0, 0.0, 0.0},
    {5.6, -1.0, 0.0, 1.0, 0.0},
    {-34.7, -1.0, 0.0, 0.0, 0.0},
    {14.2, -1.0, 0.0, 0.0, 0.0},
    {117.5, -1.0, 0.0, 0.0, 0.0},
    {-329.8, -1.0, 0.0, 0.0, 0.0},
    {23.8, -1.0, 0.0, 0.0, 0.0},
    {-9.5, -1.0, 0.0, 0.0, 0.0},
    {32.8, -1.0, 0.0, 0.0, 0.0},
    {-10.1, -1.0, 0.0, 0.0, 0.0},
    {-15.9, -1.0, 0.0, 0.0, 0.0},
    {4.8, -1.0, 0.0, 0.0, 0.0},
    {25.4, -1.0, 0.0, 0.0, 0.0},
    {7.3, -1.0, 0.0, 1.0, 0.0},
    {4.7, -1.0, 0.0, 0.0, 0.0},
    {14.2, -1.0, 0.0, 0.0, 0.0},
    {-13.6, -1.0, 0.0, 0.0, 0.0},
    {12.7, 1.0, 0.0, 0.0, 0.0},
    {409.2, 1.0, 0.0, 0.0, 0.0},
    {22.5, 1.0, 0.0, -1.0, 0.0},
    {8.7, 1.0, 0.0, 0.0, 0.0},
    {14.6, 1.0, 0.0, -1.0, 0.0},
    {-27.3, 1.0, 0.0, -1.0, 0.0},
    {-169.0, 1.0, 0.0, 0.0, 0.0},
    {13.1, 1.0, 0.0, 0.0, 0.0},
    {9.1, 1.0, 0.0, 0.0, 0.0},
    {131.7, 1.0, 0.0, 0.0, 0.0},
    {7.1, 1.0, 0.0, 0.0, 0.0},
    {12.8, 1.0, 0.0, -1.0, 0.0},
    {-943.2, 1.0, 0.0, 0.0, 0.0},
    {-29.3, 1.0, 0.0, 0.0, 0.0},
    {-388.3, 1.0, 0.0, 0.0, 0.0},
    {35.0, 1.0, 0.0, 0.0, 0.0},
    {27.3, 1.0, 0.0, 0.0, 0.0}
};     // coefficients of terms


double geodesic_distance(geographic_point gp1, geographic_point gp2)
{
    double f, g, l;             // auxiliary computational variables
    double s, c;                // auxiliary computational variables
    double o;                   // auxiliary computational variable
    double r, d, h1, h2;        // auxiliary computational variables
    double gd;                  // result geodesic distance
    
    // performing auxiliary computations
    f = (gp1.latitude + gp2.latitude) / 2.0;
    g = (gp1.latitude - gp2.latitude) / 2.0;
    l = (gp1.longitude - gp2.longitude) / 2.0;
    
    s = sin(g)*sin(g) * cos(l)*cos(l) + cos(f)*cos(f) * sin(l)*sin(l);
    c = cos(g)*cos(g) * cos(l)*cos(l) + sin(f)*sin(f) * sin(l)*sin(l);
    
    o = atan2(sqrt(s), sqrt(c));
    
    r = sqrt(s * c) / o;
    
    d = 2 * o * EARTH_EQUATORIAL_RADIUS;
    h1 = (3 * r - 1) / 2 * c;
    h2 = (3 * r + 1) / 2 * s;
    
    // computing geodesic distance
    gd = d * (EARTH_FLATTERING * h1 * sin(f) * sin(f) * cos(g) * cos(g) + 1 - EARTH_FLATTERING * h2 * cos(f) * cos(f) * sin(g) * sin(g));
    
    return gd;
}

ecliptic_point precession(ecliptic_point ep, double jd0, double jd)
{
    double t0;                  // time interval measured in Julian centuries between starting epoch and J2000
    double t1;                  // time interval measured in Julian centuries between two given epochs
    double etha, pi, p;         // auxiliary computational variables
    double a, b, c;             // auxiliary computational variables
    ecliptic_point rep;         // result eqcliptic point corrected for precession
    
    // computing time intervals from given epochs
    t0 = (jd0 - J2000) / DAYS_IN_JULIAN_CENTURY;
    t1 = (jd - jd0) / DAYS_IN_JULIAN_CENTURY;
    
    // performing auxiliary computations; constants are measured in arcseconds
    etha = (47.0029 - 0.06603 * t0 + 0.000598 * t0 * t0) * t1 + (-0.03302 + 0.000598 * t0) * t1 * t1 + 0.000060 * t1 * t1 * t1;
    pi = 629554.9824 + 3289.4789 * t0 + 0.60622 * t0 * t0 - (869.8089 + 0.50491 * t0) * t1 + 0.03536 * t1 * t1;
    p = (5029.0966 + 2.22226 * t0 - 0.000042 * t0 * t0) * t1 + (1.11113 - 0.000042 * t0) * t1 * t1 - 0.000006 * t1 * t1 * t1;
    
    // converting computed values from arcesonds to radians (π = 648000")
    etha *= M_PI / (648000.0);
    pi *= M_PI / (648000.0);
    p *= M_PI / (648000.0);
    
    // performing auxiliary computations
    a = cos(etha) * cos(ep.latitude) * sin(pi - ep.longitude) - sin(etha) * sin(ep.latitude);
    b = cos(ep.latitude) * cos(pi - ep.longitude);
    c = cos(etha) * sin(ep.latitude) + sin(etha) * cos(ep.latitude) * sin(pi - ep.longitude);
    
    // computing ecliptic point reduced to the new epoch 
    rep.longitude = -atan2(a, b) + p + pi;
    rep.latitude = asin(c);
    
    return rep;
}

double nutation_in_longitude(date d)
{
    double t;                   // time measured in Julian centuries between given date and the beggingig of the epoch J2000 
    double da[TOTAL_DELAUNAY_ARGUMENTS];      // array holding Delaunay arguments
    double lan;                 // longitude of lunar ascending node (☊) referred to the mean equinox of date
    int i;                      // loop index variable
    double a;                   // argument of each term of the serie
    double n;                   // result nutation in longitude value
    
    // computing time since the beginning of the epoch J2000 to the given measured in Julian centuries
    t = (julian_date(d) - J2000) / DAYS_IN_JULIAN_CENTURY;
    
    // computing Delaunay arguments as given by semi-analytic lunar theory ELP
    compute_delaunay_arguments(t, FULL_SERIES_TOTAL_TERMS, da);
    
    // computing longitude of the lunar ascending node (☊)
    // Source: P.K. Seidelman. 1980 IAU Theory of Nutation: The Final Report of the IAU Working Group on Nutation,
    //         Celestial Mechanics, vol. 27, May 1982, p. 20.
    lan = 450160.28 - 6962890.539 * t + 7.455 * t * t + 0.008 * t * t * t;
    
    // computing nutation serie (see top comment for more information)
    for (i = 0, n = 0.0; i < TOTAL_NUTATION_TERMS; i++){
        // computing argument of the term of the serie
        a = nutation_multipliers[i][0] * da[L];
        a += nutation_multipliers[i][1] * da[LP];
        a += nutation_multipliers[i][2] * da[F];
        a += nutation_multipliers[i][3] * da[D];
        a += nutation_multipliers[i][4] * lan;
        
        // converting argument from arcseconds to radians (π = 648000")
        a *= M_PI / 648000.0;
        
        // accumulation nutation values
        n += (nutation_coefficients[i][1] + nutation_coefficients[i][2] * t) * sin(a);
    }
    
    // coefficients are given in 10⁻⁵", converting to degrees (1° = 36000000×10⁻⁵")
    n /= 36000000.0;
    
    // converting nutation values from arcseconds to radians (π = 180°)
    n *= M_PI / 180.0;
    
    return n;
}

double nutation_in_obliquity(date d)
{
    double t;                   // time measured in Julian centuries between given date and the beggingig of the epoch J2000 
    double da[TOTAL_DELAUNAY_ARGUMENTS];      // array holding Delaunay arguments
    double lan;                 // longitude of lunar ascending node (☊) referred to the mean equinox of date
    int i;                      // loop index variable
    double a;                   // argument of each term of the serie
    double n;                   // result nutation in obliquity value
    
    // computing time since the beginning of the epoch J2000 to the given measured in Julian centuries
    t = (julian_date(d) - J2000) / DAYS_IN_JULIAN_CENTURY;
    
    // computing Delaunay arguments as given by semi-analytic lunar theory ELP
    compute_delaunay_arguments(t, FULL_SERIES_TOTAL_TERMS, da);
    
    // computing longitude of the lunar ascending node (☊)
    // Source: P.K. Seidelman. 1980 IAU Theory of Nutation: The Final Report of the IAU Working Group on Nutation,
    //         Celestial Mechanics, vol. 27, May 1982, p. 20.
    lan = 450160.28 - 6962890.539 * t + 7.455 * t * t + 0.008 * t * t * t;
    
    // computing nutation serie (see top comment for more information)
    for (i = 0, n = 0.0; i < TOTAL_NUTATION_TERMS; i++){
        // computing argument of the term of the serie
        a = nutation_multipliers[i][0] * da[L];
        a += nutation_multipliers[i][1] * da[LP];
        a += nutation_multipliers[i][2] * da[F];
        a += nutation_multipliers[i][3] * da[D];
        a += nutation_multipliers[i][4] * lan;
        
        // converting argument from arcseconds to radians (π = 648000")
        a *= M_PI / 648000.0;
        
        // accumulation nutation values
        n += (nutation_coefficients[i][3] + nutation_coefficients[i][4] * t) * cos(a);
    }
    
    // coefficients are given in 10⁻⁵", converting to degrees (1° = 36000000×10⁻⁵")
    n /= 36000000.0;
    
    // converting nutation values from arcseconds to radians (π = 180°)
    n *= M_PI / 180.0;
    
    return n;    
}

ecliptic_point aberration(date d, ecliptic_point ep)
{
    ecliptic_point stp;         // true position of the Sun in ecliptic coordinates
    double eoe[TOTAL_ORBITAL_ELEMENTS];     // mean orbital elements of the Earth
    double e;                   // mean eccentricity of Earth's orbit
    double p;                   // mean logitude of perhelion of Earth's orbit
    double k;                   // aberration constant
    ecliptic_point dep;         // result change of longitude and latitude due to aberration
    
    // computing Sun's true position
    stp = sun_true_position(d);
    // computing mean orbital elements of the Earth
    compute_orbital_elements_of_date(d, EARTH, eoe);
    // taking values of eccentricity and longitude of perhelion for shorter notation in aberration formula
    e = eoe[ECCENTRICITY_OF_ORBIT];
    p = eoe[PERHELION_LONGITUDE];
    
    // converting value of aberration constant from arcseconds to radians (π = 648000")
    k *= M_PI / 648000.0;
    
    // computing change of ecliptic position due to effect of aberration
    dep.longitude = -k * (cos(stp.longitude - ep.longitude) - e * k * cos(p - ep.longitude)) / cos(ep.latitude);
    dep.latitude = -k * sin(ep.latitude) * (sin(stp.longitude - ep.longitude) - e * sin(p - ep.longitude));
    
    return dep;
}

double obliquity_of_ecliptic(date d)
{
    double t;                   // Julian centuries from a given to the beginning of the epoch J2000
    double e;                   // result obliquity of the ecliptic (ε)
    
    // computing time from the beginning of the epoch J2000 to a given date measured in Julian centuries
    t = (julian_date(d) - J2000) / DAYS_IN_JULIAN_CENTURY;
    
    // computing obliquity of the ecliptic (ε) using formula adopted by IAU; constants are measured in arcseconds
    // Source: J.H. Lieske, T. Lederle, W. Fricke and B. Morando. Expressions for the precession Quantities Based upon
    //         the IAU (1976) System of Astronomical Constants, Astronomy ans Astrophysics, vol. 58, 1977, p. 15.
    e = 84381.448 - 46.8150 * t - 0.00059 * t * t + 0.001813 * t * t * t;
    
    // converting arcseconds to radians (π = 648000")
    e *= M_PI / 648000.0;
    
    return e;
}