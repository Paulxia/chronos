/*
 * calendar.h
 * Created by Serhii Tsyba (sertsy@gmail.com) on 23.06.10.
 *
 * This files contains data types and various routines related to time measurement and calendar phenomena.
 *
 * File defines date data type which is used throughout the library as the main type to represent time instance. Day
 * number of this data type is a decimal: integer part represents day number of the month and fractional part - hour
 * part of the day. Since there are 24 hours in a day, multiplying fractional part by 24 will result in obtaining hour
 * value of the day expressed in hours. For instance, date
 *
 *                                      12.55, December, 1900
 *
 * represents 12 December 1900 A.D. at 13ʰ12ᵐ00ˢ, since 0.55ᵈ × 24ʰ = 13,2ʰ = 13ʰ12ᵐ (0.2ʰ × 60ᵐ = 12ᵐ).
 *
 * Months are numbered from 1 to 12, 1 corresponding to January, 2 to February and so on. Values of provided Month
 * enumeration are advised to be used for convinience. Enumerated value of unknown month indexed at 0 is given to
 * indicate a potential error of a certain operation whose result is a month number.
 *
 * Year number is given as an integer. Positve years represent years A.D. and negative years represent years B.C.
 * However, since civil calendar contains no year 0 (1 B.C. is followed by 1 A.D.) value 0 corresponds to year 1 B.C.
 * Thus B.C. years are shifted one unit. For instance, date
 *
 *                                      12.55, December, -1900
 *
 * represents civil date 12 December 1901 B.C. at 13ʰ12ᵐ00ˢ.
 *
 * Convinient enumeration for week days is also given ranging from 1 to 7, where 1 corresponds to Monday, 2 to Tuesday
 * and so on. Week day value 0 is given for same purpose as in the case with months.
 *
 * Time instants must always be given in Universal Time, since library will adapt this value to Terrestrial Time when
 * needed. All output values are also given in Universal Time, unless stated otherwise.
 *
 * Each date used as an argument in any of the functions of this library must first be checked for validity within the
 * library. A provided function is_date_valid is given for this purpose. Even though this function checks that all
 * fields of date instance are in proper ranges, most importantly, it checks whether given date may be used with this
 * library. Some routines are not defined for negative Julian dates, i.e. calendar dates before 1 January 4713 B.C. at
 * noon. This is a limitation that should be kept in mind while using this library for hstoric times.
 */

#ifndef CALENDAR_H
#define CALENDAR_H

#define DAYS_IN_JULIAN_CENTURY 36525.0          // Julian century consists of exactly 36525 equal days
#define DAYS_IN_JULIAN_MILLENIUM 365250.0       // Julian millenium consists of exactly 365250 equal days

#define J2000 2451545.0                         // Julian date of the beginning of the standard epoch J2000

/*
 * An enumeration of week days provided for convinience. Indexing of week days starts at Monday, being at index 1 and
 * ends at Sunday being at index 7. Index 0 is given for potential undefined/error value in week day related routines
 * and is labeled unknown weekday.
 */
typedef enum {
    UNKNOWN_WEEKDAY = 0,
    MONDAY = 1,
    TUESDAY,
    WEDNESDAY,
    THURSDAY,
    FRIDAY,
    SATURDAY,
    SUNDAY = 7
} Weekday;

/*
 * An enumeration of months provided for convinience. Indexing starts at January, being at index 1 and ends at
 * December, being at index 12. Index 0 is given for potential undefined/error value in month related routines and is
 * labeled unknown month.
 */
typedef enum {
    UNKNOWN_MONTH = 0,
    JANUARY = 1,
    FEBRUARY,
    MARCH,
    APRIL,
    MAY,
    JUNE,
    JULY,
    AUGUST,
    SEPTEMBER,
    OCTOBER,
    NOVEMBER,
    DECEMBER = 12
} Month;

/*
 * A structure defining a calendar date consisting of day number, month and year. Day number is given as a decimal
 * value with its fractional part representing hour part of the day. Thus fractional part of the day number multiplied
 * by 24 will results in decimal hour value.
 * Month number is an integer in range [1, 12]. Using provided Month enumeration is advised for convinience.
 * Years BC are represented as negative numbers with 0 being 1 BC, -1 being 2 BC etc. since there is no year 0.
 */
typedef struct {
    double day;
    Month month;
    int year;
} date;

/*
 * Checks whether given date (d) is valid date for this library. A valid of this library is also a valid date in
 * general, however the opposite is not true (see below). Performs following checks:
 *  - day number is in range [1, 28/29/30/31], depending on the month and if a given year is common or leap;
 *  - month number is in range [1, 12];
 *  - date does not fall before Julian date 0 (1.5 January 4713 BC); even though such date is a valid calendar date,
 *    some routines in this library can not be used with negative Julian dates, and hence such restriction should be
 *    made;
 *  - date is not one of the dates removed by Gregorian reform (4 October 1582 - 15 October 1582);
 */
int is_date_valid(date d);

/*
 * Checks whether a given year (y) is common or leap. Works for both Gregorian and Julian calendar dates. Function
 * returns 0 if a given year is common and 1 if the year is a leap one.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 62.
 */
int is_leap_year(int y);

/*
 * Computes which day of the week is a given date (d).
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 65.
 */
Weekday day_of_week(date d);

/*
 * Computes day number of the year of a given date (d). Function returns value in range [1, 365] or [1, 366] depending
 * on the year of a given date being leap or common.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 65.
 */
int day_of_year(date d);

/*
 * Computes Julian date out of a given calendar date (d).
 * Julian date is the amount of time measured in days since 1.5 January 4713 BC. Julian day number is the integer part
 * of a Julian date.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 60.
 */
double julian_date(date d);

/*
 * Computes calendar date given out of a given Julian date (jd).
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 63.
 */
date calendar_date(double jd);

/*
 * Computes Julian Ephemeris date out of a given calendar date (d).
 * Julian date is the amount of time measured in days since 1.5 January 4713 BC but corrcted for Earth's clock error
 * (ΔT) at a given date. Julian day number is the integer part of a Julian Ephemeris date.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 177.
 */
double julian_ephemeris_date(date d);

/*
 * Computes date of Easter for a given year (y). Works for both Julian and Gregorian calendars. Day number of the
 * output date has no fractional part, since only day number is relevant.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 67.
 */
date date_of_easter(int y);

/*
 * Computes the Earth's clock error (ΔT), i.e. difference of Universal Time (UT) and Terresrial Time (TT) for a given
 * date (d).
 *
 * Source: L.V. Morrison, F.R. Stephenson. Historical values of the Earth's clock error ΔT and the calculation of
 *         eclipses. Journal for the History of Astronomy, vol. 35, 2004, pp. 327-336.
 *
 *         L.V. Morrison, F.R. Stephenson. Addendum. Historical values of Earth's clock error. Journal for the History
 *         of Astronomy, vol. 36, 2005, p. 339.
 *
 *         http://maia.usno.navy.mil/ser7/deltat.preds
 */
double dynamical_time_difference(date d);

/*
 * Computes mean siderial time at Greenwich meridian (GST) on a given date (d), i.e. the Greenwich hour angle of the
 * mean vernal point (the intersection  of the ecliptic of the date with the mean equator of the date). Output value is
 * expressed in decimal hours.
 *
 * Source: J. Meuss. Astronomical Algorothms. William-Bell, 1991, p. 84.
 */
double greenwich_mean_siderial_time(date d);

#endif // CALENDAR_H