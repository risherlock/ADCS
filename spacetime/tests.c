// Numerical validations
// 2024-12-25

#include <stdio.h>
#include "time.h"

int main(void)
{
  // Example 5.04. Ref [1].
  utc_t utc;
  utc.year = 2004;
  utc.month = 5;
  utc.day = 12;
  utc.hour = 14;
  utc.minute = 45;
  utc.second = 30;

  printf("Julian day number: %f\n", time_julian_day(utc)); // 2453138.115

  // Example 5.06. Ref [1].
  utc.year = 2004;
  utc.month = 3;
  utc.day = 3;
  utc.hour = 4;
  utc.minute = 30;
  utc.second = 0;

  const double el = 139.0 + 47.0 / 60.0 + 0.0 / 3600.0;

  printf("Local sidereal time [deg]: %f\n", time_local_sidereal_deg(utc, el)); // 8.57688
  printf("Local sidereal time [hr]: %f\n", time_local_sidereal_hr(utc, el)); // 0.571792

  return 0;
}
