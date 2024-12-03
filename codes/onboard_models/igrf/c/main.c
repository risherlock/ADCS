// Basic implementation of 'igrf.h'
// Rishav (2021-12-27)

#include <stdio.h>
#include "igrf.h"

#define R2D 57.2957795131
#define D2R 0.01745329251

int main()
{
  // Time
  date_time dt;
  dt.year = 2023;
  dt.month = 12;
  dt.day = 17;
  dt.hour = 0;
  dt.minute = 0;
  dt.second = 0;

  // 221B Baker Street
  const float latitude = 51.5238; // deg
  const float longitude = -0.1586; // deg
  const float height = 1000.0; // km
  const float x_sph[3] = {latitude, longitude, height};

  // Magnetic field in NED frame
  float b_ned[3] = {0.0};

  // Compute and print
  if(igrf(dt, x_sph, b_ned))
  {
    printf("Inputs:\n");
    printf("  Date: %d-%d-%d, %d:%d:%d\n", dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second);
    printf("  Latitude: %f deg\n", latitude);
    printf("  Longitude: %f deg\n", longitude);
    printf("  Height: %f km\n", height);
    printf("\nOutputs:\n");
    printf("  Bn: %f nT\n", b_ned[0]);
    printf("  Be: %f nT\n", b_ned[1]);
    printf("  Bd: %f nT\n", b_ned[2]);
    printf("  Magnitude: %f nT\n", igrf_get_norm(b_ned));
    printf("  Inclination: %f deg\n", igrf_get_inclination(b_ned) * R2D);
    printf("  Declination: %f deg\n", igrf_get_declination(b_ned) * R2D);
  }
  else
  {
    printf("Date error!\n");
  }

  return 0;
}
