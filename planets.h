/*
 *  planets.h
 *  Created by Serhii Tsyba on 30.07.10.
 */

#ifndef PLANETS_H
#define PLANETS_H

#include "calendar.h"
#include "coordinates.h"
#include <vsop87d.h>

ecliptic_point planet_true_position(date d, Planet p);
ecliptic_point planet_apparent_position(date d, Planet p);
double planet_distance_to_sun(date d, Planet p);
double planet_distance_to_earth(date d, Planet p);
double planet_phase_angle(date d, Planet p);
double planet_disk_illuminated_fraction(date d, Planet p);
double planet_apparent_magnitude(date d, Planet p);

#endif // PLANETS_H