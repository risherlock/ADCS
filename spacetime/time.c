#include <math.h>
#include "time.h"

double fix(double x)
{
  if (x > 0)
  {
    return floor(x);
  }

  return ceil(x);
}

// UTC to Julian day.
double time_julian_day(const utc_t t)
{
  double j0 = 367 * t.year - fix(7.0 * (t.year + fix((t.month + 9.0) / 12.0)) / 4.0) + fix(275.0 * t.month / 9.0) + t.day + 1721013.5;
  double ut = t.hour + t.minute / 60.0 + t.second / 3600.0;

  return j0 + ut / 24.0;
}

double normalize_zero_to_360(double angle)
{
  while (angle < 0)
  {
    angle += 360.0;
  }

  while (angle >= 360.0)
  {
    angle -= 360.0;
  }

  return angle;
}

// UTC to Greenwich sidereal time.
double time_greenwich_sidereal(const utc_t t)
{
  const double j0 = 367 * t.year - fix(7.0 * (t.year + fix((t.month + 9.0) / 12.0)) / 4.0) + fix(275.0 * t.month / 9.0) + t.day + 1721013.5;
  const double j = (j0 - 2451545.0) / 36525.0;
  const double j_sq = j * j;
  const double g0 = normalize_zero_to_360(100.4606184 + 36000.77004 * j + 0.000387933 * j_sq - 2.583e-8 * j_sq * j);
  const double ut = t.hour + t.minute / 60.0 + t.second / 3600.0;

  return g0 + 360.98564724 * ut / 24.0;
}

/**
 * @brief UTC to local sidereal time [deg].
 * @param t Time in UTC
 * @param el East longitude [deg]
 * @return Local sidereal time [deg]
 */
double time_local_sidereal_deg(const utc_t t, const double elon)
{
  const double lst = time_greenwich_sidereal(t) + elon;
  return lst - 360.0 * fix(lst / 360.0);
}

// UTC to local sidereal time [hr].
double time_local_sidereal_hr(const utc_t t, const double elon)
{
  return time_local_sidereal_deg(t, elon) / 15.0;
}
