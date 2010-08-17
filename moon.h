/*
 *  moon.h
 *  Created by Serhii Tsyba on 03.08.10.
 */

#include "calendar.h"
#include "coordinates.h"

ecliptic_point moon_true_position(date d);
ecliptic_point moon_apparent_position(date d);
double moon_distance_to_earth(date d);
double moon_phase_angle(date d);
double moon_disk_illuminated_fraction(date d);
double moon_bright_limb_position_angle(date d);