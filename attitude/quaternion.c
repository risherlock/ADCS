#include <math.h>
#include "quaternion.h"

#define QUAT_COMPARISON_EPSILON 1e-5
#define QUAT_DIVISION_EPSILON 1e-9

void quat_unit(double q[4])
{
  q[0] = 1.0;
  q[1] = 0.0;
  q[2] = 0.0;
  q[3] = 0.0;
}

void quat_normalize(double q[4])
{
  double norm = quat_norm(q);

  if (norm > QUAT_DIVISION_EPSILON)
  {
    q[0] *= q[0] / norm;
    q[1] *= q[1] / norm;
    q[2] *= q[2] / norm;
    q[3] *= q[3] / norm;
  }
}

double quat_norm(const double q[4])
{
  return sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
}

void quat_conj(const double q[4], double q_conj[4])
{
  q_conj[0] = q[0];
  q_conj[1] = -q[1];
  q_conj[2] = -q[2];
  q_conj[3] = -q[3];
}

void quat_copy(const double q_src[4], double q_dest[4])
{
  q_dest[0] = q_src[0];
  q_dest[1] = q_src[1];
  q_dest[2] = q_src[2];
  q_dest[3] = q_src[3];
}

bool quat_is_equal(const double q1[4], const double q2[4])
{
  bool q0_eq = fabs(q1[0] - q2[0]) <= QUAT_COMPARISON_EPSILON;
  bool q1_eq = fabs(q1[1] - q2[1]) <= QUAT_COMPARISON_EPSILON;
  bool q2_eq = fabs(q1[2] - q2[2]) <= QUAT_COMPARISON_EPSILON;
  bool q3_eq = fabs(q1[3] - q2[3]) <= QUAT_COMPARISON_EPSILON;

  return q0_eq && q1_eq && q2_eq && q3_eq;
}

double quat_dot_prod(const double q1[4], const double q2[2])
{
  return q1[0] * q2[0] + q1[1] * q2[1] + q1[2] * q2[2] + q1[3] * q2[3];
}

void quat_prod(const double q1[4], const double q2[4], double q[4])
{
  q[0] = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3];
  q[1] = q1[1] * q2[0] + q1[0] * q2[1] + q1[2] * q2[3] - q1[3] * q2[2];
  q[2] = q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1];
  q[3] = q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0];
}

void quat_err(const double q1[4], const double q2[4], double qe[4])
{
  quat_conj(q2, qe);
  quat_prod(q1, qe, qe);
}

bool quat_to_axis_angle(const double q[4], double *psi, double v[3])
{
  bool status = true;
  *psi = 2.0 * acos(q[0]);
  double divisor = sqrt(1.0 - q[0] * q[0]);

  if (divisor > QUAT_DIVISION_EPSILON)
  {
    v[0] = q[1] / divisor;
    v[1] = q[2] / divisor;
    v[2] = q[3] / divisor;
  }
  else // Infinite solutions
  {
    // Arbitrary axis
    v[0] = 1;
    v[1] = 0;
    v[2] = 0;

    status = false;
  }

  return status;
}

void quat_axis_angle(const double psi, const double v[3], double q[4])
{
  double half_psi = psi / 2.0;
  double sin_half_psi = sin(half_psi);

  q[0] = cos(half_psi);
  q[1] = sin_half_psi * v[0];
  q[2] = sin_half_psi * v[1];
  q[3] = sin_half_psi * v[2];
}

void quat_rotate(const double q[4], const double v[3], double v_out[3])
{
  double q0_sq = q[0] * q[0];
  double q1_sq = q[1] * q[1];
  double q2_sq = q[2] * q[2];
  double q3_sq = q[3] * q[3];
  double q0q1x2 = 2.0 * q[0] * q[1];
  double q0q2x2 = 2.0 * q[0] * q[2];
  double q0q3x2 = 2.0 * q[0] * q[3];
  double q1q2x2 = 2.0 * q[1] * q[2];
  double q1q3x2 = 2.0 * q[1] * q[3];
  double q2q3x2 = 2.0 * q[2] * q[3];

  v_out[0] = q0_sq * v[0] + q0q2x2 * v[2] +
             q1_sq * v[0] - q0q3x2 * v[1] -
             q2_sq * v[0] + q1q2x2 * v[1] -
             q3_sq * v[0] + q1q3x2 * v[2];

  v_out[1] = q0_sq * v[1] - q0q1x2 * v[2] -
             q1_sq * v[1] + q0q3x2 * v[0] +
             q2_sq * v[1] + q1q2x2 * v[0] -
             q3_sq * v[1] + q2q3x2 * v[2];

  v_out[2] = q0_sq * v[2] + q0q1x2 * v[1] -
             q1_sq * v[2] - q0q2x2 * v[0] -
             q2_sq * v[2] + q1q3x2 * v[0] +
             q3_sq * v[2] + q2q3x2 * v[1];
}
