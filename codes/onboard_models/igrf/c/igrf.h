/*
  Implementation of IGRF Model.

  References:
    [1] Davis - Mathematical Modeling of Earth's Magnetic Field (2004)
    [2] Yang - Spacecraft Modeling Attitude Determination and Control (2019)

  Optimizations [a] and [b] are due to Alar Leibak from Estonia.

  Rishav (2021-12-26)
*/

#ifndef _IGRF_H_
#define _IGRF_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "igrf13.h"

#include <math.h>
#include <inttypes.h>

typedef struct date_time_t
{
  uint16_t year;
  uint8_t month;
  uint8_t day;
  uint8_t hour;
  uint8_t minute;
  uint8_t second;
}date_time;

uint8_t igrf(const date_time dt, const float x_sph[3], float b_ned[3]);
float igrf_get_inclination(const float b_ned[3]);
float igrf_get_declination(const float b_ned[3]);
float igrf_get_norm(const float b_ned[3]);

#ifdef __cplusplus
}
#endif
#endif // igrf.h
