/*
 * calendar.c
 * Created by Serhii Tsyba (sertsy@gmail.com) on 23.06.10.
 */

#include "calendar.h"
#include <math.h>

/*
 * The following definitions and two arrays define values of difference between Dynamical Time and Universal time
 * commonly denoted as ΔT:
 *
 *                                  ΔT = Dynamical Time - Universal Time
 *
 * Data is split into two tables: for telescope era (1700 A.D. till present) and pre telescope era (1000 B.C. till 1700
 * A.D.). These values were computed by L.V. Morrison and F.R. Stephenson from records of lunar eclipses. Thus, ΔT
 * values given in the first table are more precise and later one are rather interpolations.
 * Each record consists of three values: year, at which value of ΔT is measured, the value of ΔT, measured in seconds
 * and value of uncertainty for given ΔT, also measured in seconds.
 * To find value of ΔT in intermediade years within tables limits, a linear interpolation technique is recommended.
 *
 * Source: L.V. Morrison, F.R. Stephenson. Historical values of the Earth's clock error ΔT and the calculation of
 *         eclipses. Journal for the History of Astronomy, vol. 35, 2004, pp. 327-336.
 *         L.V. Morrison, F.R. Stephenson. Addendum. Historical values of Earth's clock error. Journal for the History
 *         of Astronomy, vol. 36, 2005, p. 339.
 *         http://maia.usno.navy.mil/ser7/deltat.preds
 */
#define DELTA_T_TABLE_START_YEAR -1000          // lower limit of ΔT tables
#define DELTA_T_TABLE_END_YEAR 2020             // upper limit of ΔT tables

#define PRE_TELESCOPE_ERA_START_YEAR -1000      // start year of so called pre telescop era 1001 B.C.
#define PRE_TELESCOPE_ERA_YEAR_INTERVAL 100     // year interval at which values are given for pre telescope era
#define PRE_TELESCOPE_ERA_TOTAL_TERMS 28        // total amount of values given in table of pre telescope era

#define TELESCOPE_ERA_START_YEAR 1700           // start year of so called telescope era (or modern) 1700 A.D.
#define TELESCOPE_ERA_YEAR_INTERVAL 10          // year interval at which values are given for modern era
#define TELESCOPE_ERA_TOTAL_TERMS 33            // total amount of values given in table of modern era

static int deltat_pre_telescope_era[PRE_TELESCOPE_ERA_TOTAL_TERMS][3] = {
    {-1000, 25400, 640},
    {-900, 23700, 590},
    {-800, 22000, 550},
    {-700, 20400, 500},
    {-600, 18800, 460},
    {-500, 17190, 430},
    {-400, 15530, 390},
    {-300, 14080, 360},
    {-200, 12790, 330},
    {-100, 11640, 290},
    {0, 10580, 260},
    {100, 9600, 240},
    {200, 8640, 210},
    {300, 7680, 180},
    {400, 6700, 160},
    {500, 5710, 140},
    {600, 4740, 120},
    {700, 3810, 100},
    {800, 2960, 80},
    {900, 2200, 70},
    {1000, 1570, 55},
    {1100, 1090, 40},
    {1200, 740, 30},
    {1300, 490, 20},
    {1400, 320, 20},
    {1500, 200, 20},
    {1600, 120, 20},
    {1700, 9, 5}
};  // values of ΔT for pre telescope era
static int deltat_telescope_era[TELESCOPE_ERA_TOTAL_TERMS][3] = {
    {1700, 9, 5},
    {1710, 10, 3},
    {1720, 11, 3},
    {1730, 11, 3},
    {1740, 12, 2},
    {1750, 13, 2},
    {1760, 15, 2},
    {1770, 16, 2},
    {1780, 17, 1},
    {1790, 17, 1},
    {1800, 14, 1},
    {1810, 13, 1},
    {1820, 12, 1},
    {1830, 8, 1},
    {1840, 6, 0},
    {1850, 7, 0},
    {1860, 8, 0},
    {1870, 2, 0},
    {1880, -5, 0},
    {1890, -6, 0},
    {1900, -3, 0},
    {1910, 10, 0},
    {1920, 21, 0},
    {1930, 24, 0},
    {1940, 24, 0},
    {1950, 29, 0},
    {1960, 33, 0},
    {1970, 40, 0},
    {1980, 51, 0},
    {1990, 57, 0},
    {2000, 65, 0},
    {2010, 66, 0},
    {2020, 71, 4}
};          // values of ΔT for modern era

// contains two arrays of number of days in each month; first element being 0 given to coinside with Unknown_Month
// value of enumeration Month
static const int days_in_month[2][13] = {
    {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},    // for common year
    {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}     // for leap year
};

// first date of the Julian proleptic calendar (Julian date 0.0)
static const date JULIAN_START_DATE = {1.5, JANUARY, -4712};
// last date of the Julian calendar considered by Gregorian reform
static const date JULIAN_END_DATE = {4.0, OCTOBER, 1582};
// first date of the Gregorian calendar
static const date GREGORIAN_START_DATE = {15.0, OCTOBER, 1582};

/*
 * Compares two given dates (d1, d2) to see which one occurs later. Function will return 1 if dates occur in given
 * order, namely, d1 precedes d2, returns -1 if dates occur in opposite order, i.e. d2 occurs before d1 and will
 * return 0 if two dates are equal.
 */
static int compare_dates(date d1, date d2)
{
    if (d1.year == d2.year){
        if (d1.month == d2.month){
            if (d1.day == d2.day)
                return 0;
            else
                return d1.day < d2.day ? 1 : -1;
        }
        else
            return d1.month < d2.month ? 1 : -1;
    }
    else
        return d1.year < d2.year ? 1 : -1;
}

/*
 * Performs linear interpolation in two dimensional space given two points (x0, y0) and (x1, y1) and then computes
 * value of a y = f(x), given x, where f(x) is linear interpolant on [x0, x1].
 */
static double linear_interpolate(double x, double x0, double x1, double y0, double y1)
{
    double y;                   // result interpolated value for a given argument
    
    // computing value of interpolant
    y = y0 + (x - x0) * (y1 - y0) / (x1 - x0);
    
    return y;
}

int is_date_valid(date d)
{
    // checking if month number is in interval [1, 12]
    if (d.month < JANUARY || d.month > DECEMBER)
        return 0;
    
    // checking if day number is in interval [1, 30), [1, 31) or [1, 32), depending on the month number and whether
    // a given year is a common year or a leap year
    if (!(d.day >= 0 || d.day < days_in_month[is_leap_year(d.year)][d.month]) + 1)
        return 0;
    
    // checking whether a given date would result in a positive Julian date; although dates erlier than beginning of
    // Julian calendar are valid, some routines in Chronos library are invalid with negative Julian dates, namely,
    // calendar dates before 1.5 January 4713 B.C.
    if (compare_dates(d, JULIAN_START_DATE) > 0)
        return 0;
    
    // checking whether a given calendar date does not occur on days removed by Gregorian reform: 5 October 1582 till
    // 14 October 1582 inclusive
    if (compare_dates(JULIAN_END_DATE, d) > 0 &&
        compare_dates(GREGORIAN_START_DATE, d) < 0)
        return 0;
    
    return 1;
}

int is_leap_year(int y)
{
    // for Gregorian calendar a leap year is one divisible by four, but not by 100, except for 400
    if (y >= GREGORIAN_START_DATE.year)
        return (y % 4 == 0 && y % 100 != 0 || y % 400 == 0) ? 1 : 0;
    
    // for Julian calendar a leap year is one divisible by four
    else
        return (y % 4 == 0) ? 1 : 0;
}

Weekday day_of_week(date d)
{
    double jd;                  // Julian date of a given date
    Weekday wd;                 // result day of the week
    
    // computing Julian date from a given date
    jd = julian_date(d);
    
    // since Julian date 0 is Monday, then week day number of a given date is its Julian date modulo 7
    wd = (int)fmod(jd + 0.5, 7.0);
    
    // week day enumeration of this library assumes Monday being indexed at 1 and Sunday at 7
    wd += 1;
    
    return wd;
}

int day_of_year(date d)
{
    int l;                      // variable holding 1 for a leap year and 2 for a common year
    int yd;                     // result day of the year
    
    // initializing l with 2 for a leap year and 1 for a common year
    l = (is_leap_year(d.year)) ? 1 : 2;
    
    // computing day of the year
    yd = (int)(275.0 * (double)d.month / 9.0) - l * (int)(((double)d.month + 9.0) / 12.0) + (int)d.day - 30;
    
    return yd;
}

double julian_date(date d)
{
    int a, b;                   // auxiliary computational variables
    double jd;                  // result Julian date
    
    // if month is either January or February, assuming month to be December of the previous year
    if (d.month == JANUARY || d.month == FEBRUARY) {
        d.year--;
        d.month += 12;
    }
    
    // performing auxiliary computations
    // for Gregorian calendar date
    if (compare_dates(GREGORIAN_START_DATE, d) >= 0){
        a = d.year / 100;
        b = 2 - a + a / 4;
    }
    // for Julian calendar date
    else
        b = 0;
    
    // computing Julian date
    jd = (int)(365.25 * (double)(d.year + 4716)) + (int)(30.6001 * (double)(d.month + 1)) + d.day + b - 1524.5;
    
    return jd;
}

date calendar_date(double jdn)
{
    double a, b, c, d, e;       // auxiliary computational variables
    double ap;                  // auxiliary computational variable
    double z, f;                // integer and fractional parts of the Julian date
    date cd;                    // result calendar date
    
    // splitting Julian date into Julian day number and decimal hour part
    f = modf(jdn + 0.5, &z);
    
    // performing auxiliary computatations
    if (z < 2299161.0)
        a = z;
    else {
        modf((z - 1867216.25) / 36524.25, &ap);
        a = z + 1.0 + ap - (int)(ap / 4.0);
    }
    
    b = a + 1524.0;
    c = (int)((b - 122.1) / 365.25);
    d = (int)(365.25 * c);
    e = (int)((b - d) / 30.6001);
    
    // computing calendar date
    cd.day = b - d - (int)(30.6001 * e) + f;
    cd.month = e < 14.0 ? (int)e - 1 : (int)e - 13;
    cd.year = cd.month > 2 ? (int)c - 4716 : (int)c - 4715;
    
    return cd;
}

double julian_ephemeris_date(date d)
{
    double dt;                  // value of difference between Dynamical Time (DT) and Universal Time (UT), denoted as ΔT
    double jde;                 // result Julian ephemeris date (JDE)
    
    // computing difference between Dynamical Time and Universal Time
    dt = dynamical_time_difference(d);
    // converting ΔT value from seconds to Julian days (1ᵈ = 86400ˢ)
    dt /= 86400.0;
    
    // Julian Ephemeris date is Julian date in Dynamical Time
    jde = julian_date(d) + dt;
    
    return jde;
}

date date_of_easter(int y)
{
    int a, b, c, d, e, f;       // auxiliary computational variables
    date ed;                    // result Easter date
    
    // performing auxiliary computations
    
    // for Gregorian calendar;
    // in condition > (but not >=) is used because Gregorian calendar starts in October and and hence Easter of 1582 is
    // still computed by Julian calendar so the year must be at least 1583 to use these formulas
    if (y > GREGORIAN_START_DATE.year){
        a = y / 100;
        b = (a - ((a + 8) / 25) + 1) / 3;
        c = (19 * (y % 19) + a - (a / 4) - b + 15) % 30;
        d = (32 + 2 * (a % 4) + 2 * ((y % 100) / 4) - c - ((y % 100) % 4)) % 7;
        e = ((y % 19) + 11 * c + 22 * d) / 451;
        f = (c + d - 7 * e + 114);
    }
    // for Julian calendar
    else {
        a = (19 * (y % 19) + 15) % 30;
        b = (2 * (y % 4) + 4 * (y % 7) - a + 34) % 7;
        f = a + b + 114;
    }
    
    // computing date of Easter
    ed.day = f % 31 + 1;
    ed.month = f / 31;
    ed.year = y;
    
    return ed;
}

double dynamical_time_difference(date d)
{
    int i;                      // loop index variable
    double dy;                  // given date expressed in decimal years
    double dt;                  // result value of difference between Dynamical Time and Universal Time (ΔT)
    
    if (d.year < DELTA_T_TABLE_START_YEAR || d.year > DELTA_T_TABLE_END_YEAR)
        // if a given is beyond limits of the table of emperically computed values of ΔT then using proposed formula to
        // compute extrapolated approximate value
        dt = -20.0 + 32.0 * (d.year - 1820.0) * (d.year - 1820.0) / 10000.0;
    else {
        // if a value can be picked from a given table then converting given date to decimal years to interpolate
        // given year more precise
        dy = d.year + day_of_year(d) / (365.0 + is_leap_year(d.year));
        
        if (d.year < TELESCOPE_ERA_START_YEAR){
            // using table of years from pre telescope era: 1000 B.C. to 1700 A.D.
            
            // computing index of the year in the table that a given year follows
            i = (int)floor((d.year - PRE_TELESCOPE_ERA_START_YEAR) / PRE_TELESCOPE_ERA_YEAR_INTERVAL);
            
            // computing linearly interpolated value of ΔT
            dt = linear_interpolate(dy, deltat_pre_telescope_era[i][0], deltat_pre_telescope_era[i + 1][0],
                                    deltat_pre_telescope_era[i][1], deltat_pre_telescope_era[i + 1][1]);
        }
        else{
            // using table of years from modern era: 1700 A.D. to 2020 A.D.
            
            // computing index of the year in the table that a given year follows
            i = (int)floor((d.year - TELESCOPE_ERA_START_YEAR) / TELESCOPE_ERA_YEAR_INTERVAL);
            
            // computing linearly interpolated value of ΔT
            dt = linear_interpolate(dy, deltat_telescope_era[i][0], deltat_telescope_era[i + 1][0],
                                    deltat_telescope_era[i][1], deltat_telescope_era[i + 1][1]);
        }
    }
    
    return dt;
}

double greenwich_mean_siderial_time(date d)
{
    double jd;                  // Julian date of a given calendar date
    double t;                   // time interval measured in Julian centuries since the beginning of the epoch J2000
                                // till a given date
    double gmst;                // result Greenwich mean siderial time
    
    // computing Julian date from a given calendar date
    jd = julian_date(d);
    
    // computing amount of Julian centuries from the beginning of the epoch J2000 till a given calendar date
    t = (jd - J2000) / DAYS_IN_JULIAN_CENTURY;
    
    // computing Greenwich mean siderial time expressed in degrees
    gmst = 280.46061837 + 360.98564736629 * (jd - J2000) + 0.000387933 * t * t - t * t * t / 38710000;
    
    // shifting the value of mean siderial time at Greenwich to interval [0, 360);
    gmst = fmod(gmst, 360.0);
    if (gmst < 0.0)
        gmst += 360.0;
    
    // converting degrees to hours (1ʰ = 15°)
    gmst /= 15.0;
    
    return gmst;
}