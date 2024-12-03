#include "igrf.h"

#define PI 3.14159265358979323846
#define PI_2 1.57079632679489661923
#define R2D 57.2957795131
#define D2R 0.01745329251

// WGS84 params
static const double wgs84_f = 1 / 298.257223563;
static const double wgs84_a = 6378.137;
static const double wgs84_b = wgs84_a * (1 - wgs84_f);

// Constant parameters
static const double E2 = 1.0f - (wgs84_b / wgs84_a) * (wgs84_b / wgs84_a);
static const double E4 = E2 * E2;
static const double E6 = E4 * E2;
static const double E8 = E4 * E4;
static const double A21 = (512.0f * E2 + 128.0f * E4 + 60.0f * E6 + 35.0f * E8) / 1024.0f;
static const double A22 = (E6 + E8) / 32.0f;
static const double A23 = -3.0f * (4.0f * E6 + 3.0f * E8) / 256.0f;
static const double A41 = -(64.0f * E4 + 48.0f * E6 + 35.0f * E8) / 1024.0f;
static const double A42 = (4.0f * E4 + 2.0f * E6 + E8) / 16.0f;
static const double A43 = 15.0f * E8 / 256.0f;
static const double A44 = -E8 / 16.0f;
static const double A61 = 3.0f * (4.0f * E6 + 5.0f * E8) / 1024.0f;
static const double A62 = -3.0f * (E6 + E8) / 32.0f;
static const double A63 = 35.0f * (4.0f * E6 + 3.0f * E8) / 768.0f;
static const double A81 = -5.0f * E8 / 2048.0f;
static const double A82 = 64.0f * E8 / 2048.0f;
static const double A83 = -252.0f * E8 / 2048.0f;
static const double A84 = 320.0f * E8 / 2048.0f;

float get_years(const date_time dt)
{
  if ((dt.year < IGRF_START_YEAR) || (dt.year > IGRF_END_YEAR) ||
      (dt.month > 12) || (dt.day > 31) || (dt.day == 0) ||
      (dt.hour > 23) || (dt.minute > 60) || (dt.second > 60))
  {
    return -1;
  }

  // Days since IGRF_START_YEAR
  int years = dt.year - IGRF_START_YEAR;
  int days_arr[] = {0, 31, 59, 90, 120, 151, 182, 212, 243, 273, 304, 334};
  int is_leap = (((dt.year % 4) == 0) && (((dt.year % 100) != 0) || ((dt.year % 400) == 0)));
  double days = (days_arr[dt.month - 1] + dt.day + (dt.month > 2 ? is_leap : 0));
  double hours = dt.hour + (dt.minute / 60) + (dt.second / 3600);
  int total_days = is_leap ? 366 : 365;

  // Decimal years
  return years + days / total_days + hours / (24.0f * days);
}

/*
  Computes magnetic field strength [nT] in NED coordinates.

Inputs:
  dt = YYYY-MM-DD HH:MM:SS
  x_sph[3] = {geodetic latitude [deg], longitude [deg], surface height [km]}

Output:
  b_ned[3] = {bn, be, bd}, nT

Outline:
  1. Recursively calculate Legendre polynostatic const double wgs84_e2 = 6.69437999014e-3;
mials and its derivatives.
  2. Compute magnetic field strength using terms computed in 1.
*/
#include <stdio.h>
uint8_t igrf(const date_time dt, const float x_sph[3], float b_ned[3])
{
  const float a = 6371.2; // Radius of Earth, km
  const float theta = PI_2 - x_sph[0] * D2R;
  const float phi = x_sph[1] * D2R;
  const float h = x_sph[2];

  double ct, st;
  ct = cos(theta);
  st = sin(theta);

  // Geodetic to geocentric coordinates
  const double rho = hypot(wgs84_a * st, wgs84_b * ct);
  const double r = sqrt(h * h + 2 * h * rho + (pow(wgs84_a, 4) * st * st + pow(wgs84_b, 4) * ct * ct) / (rho * rho));
  const double cd = (h + rho) / r;
  const double sd = (wgs84_a * wgs84_a - wgs84_b * wgs84_b) / rho * ct * st / r;

  float temp_ct = ct;
  ct = ct * cd - st * sd;
  st = st * cd + temp_ct * sd;

  // Avoid singularity
  float epsilon = 1e-8;
  if (st < epsilon && st > -epsilon)
  {
    st = epsilon;
  }

  float years = get_years(dt);
  if (years < 0)
  {
    return 0;
  }

  // [a] Re-occurring power factors
  float ar_pow[IGRF_DEGREE];
  const float ar = a / r;
  ar_pow[0] = ar * ar * ar;
  for (uint8_t i = 1; i <= IGRF_DEGREE; i++)
  {
    ar_pow[i] = ar_pow[i - 1] * ar;
  }

  // [b] Re-occurring sines and cosines
  float sines[IGRF_DEGREE], cosines[IGRF_DEGREE];
  sines[0] = 0;
  cosines[0] = 1;
  sines[1] = sin(phi);
  cosines[1] = cos(phi);

  for (uint8_t i = 2; i < IGRF_DEGREE; i++)
  {
    if (i & 1)
    {
      sines[i] = sines[i - 1] * cosines[1] + cosines[i - 1] * sines[1];
      cosines[i] = cosines[i - 1] * cosines[1] - sines[i - 1] * sines[1];
    }
    else // even
    {
      sines[i] = 2 * sines[i >> 1] * cosines[i >> 1];
      cosines[i] = 2 * cosines[i >> 1] * cosines[i >> 1] - 1;
    }
  }

  float pnm, dpnm, p11, dp11, p10, dp10, p20, dp20;
  float br = 0.0f, bt = 0.0f, bp = 0.0f; // radial, theta, and phi components

  // n = m = 0
  p11 = 1;
  p10 = 1;
  dp11 = 0;
  dp10 = 0;

  int k = 0;
  for (uint8_t m = 0; m <= IGRF_DEGREE; m++)
  {
    for (uint8_t n = 1; n <= IGRF_DEGREE; n++)
    {
      if (m <= n)
      {
        if (n == m)
        {
          pnm = st * p11;
          dpnm = st * dp11 + ct * p11;

          p11 = pnm;
          dp11 = dpnm;
        }
        else
        {
          float Knm = 0;

          if (n != 1)
          {
            Knm = (float)((n - 1) * (n - 1) - m * m) / ((2.0f * n - 1) * (2.0f * n - 3));
          }

          pnm = ct * p10 - Knm * p20;
          dpnm = ct * dp10 - st * p10 - Knm * dp20;
        }

        p20 = p10;
        p10 = pnm;
        dp20 = dp10;
        dp10 = dpnm;

        // Linear interpolation of g and h
        k = (int)((0.5 * n * (n + 1) + m) - 1);
        const float g = g_val[k] + g_sv[k] * years;
        const float h = h_val[k] + h_sv[k] * years;

        float hsin = h * sines[m];
        float hcos = h * cosines[m];
        float gsin = g * sines[m];
        float gcos = g * cosines[m];

        br += ar_pow[n - 1] * (n + 1) * ((gcos + hsin) * pnm);
        bt += ar_pow[n - 1] * ((gcos + hsin) * dpnm);
        bp += ar_pow[n - 1] * (m * (-gsin + hcos) * pnm);
      }
    }
  }

  // Geocentric NED
  b_ned[0] = bt;
  b_ned[1] = -bp / st;
  b_ned[2] = -br;

  // Geocentric to geodetic coordinates
  const double GCLAT = (PI_2 - theta);
  const double SCL = sin(GCLAT);
  const double RI = a / r;
  const double A2 = RI * (A21 + RI * (A22 + RI * A23));
  const double A4 = RI * (A41 + RI * (A42 + RI * (A43 + RI * A44)));
  const double A6 = RI * (A61 + RI * (A62 + RI * A63));
  const double A8 = RI * (A81 + RI * (A82 + RI * (A83 + RI * A84)));
  const double CCL = sqrt(1 - SCL * SCL);
  const double S2CL = 2.0f * SCL * CCL;
  const double C2CL = 2.0f * CCL * CCL - 1.;
  const double S4CL = 2.0f * S2CL * C2CL;
  const double C4CL = 2.0f * C2CL * C2CL - 1.;
  const double S8CL = 2.0f * S4CL * C4CL;
  const double S6CL = S2CL * C4CL + C2CL * S4CL;
  const double DLTCL = S2CL * A2 + S4CL * A4 + S6CL * A6 + S8CL * A8;
  const double theta_gd = DLTCL + GCLAT;
  const double psi = sin(theta_gd) * sin(theta) - cos(theta_gd) * cos(theta);

  // Geodetic NED
  b_ned[0] = cos(psi) * bt - sin(psi) * br;
  b_ned[2] = -(sin(psi) * bt + cos(psi) * br);

  return 1;
}

// Magnetic declination in radians
float igrf_get_inclination(const float b_ned[3])
{
  const float norm = sqrt(b_ned[0] * b_ned[0] + b_ned[1] * b_ned[1]);
  return atan(b_ned[2] / norm);
}

// Magnetic declination in radians
float igrf_get_declination(const float b_ned[3])
{
  return atan2(b_ned[1], b_ned[0]);
}

// Magnitude of magnetic field
float igrf_get_norm(const float b_ned[3])
{
  return sqrt(b_ned[0] * b_ned[0] + b_ned[1] * b_ned[1] + b_ned[2] * b_ned[2]);
}
