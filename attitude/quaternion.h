// Basic quaternion manipulations.
// 2024-12-22

#ifndef _ATTITUDE_QUAT_H_
#define _ATTITUDE_QUAT_H_

#include <stdbool.h>

void quat_unit(double q[4]);
void quat_normalize(double q[4]);
double quat_norm(const double q[4]);
void quat_conj(const double q[4], double q_conj[4]);
void quat_copy(const double q_src[4], double q_dest[4]);
bool quat_is_equal(const double q1[4], const double q2[4]);
double quat_dot_prod(const double q1[4], const double q2[2]);
void quat_prod(const double q1[4], const double q2[4], double q[4]);
void quat_err(const double q1[4], const double q2[4], double qe[4]);
bool quat_to_axis_angle(const double q[4], double *psi, double v[3]);
void quat_axis_angle(const double psi, const double v[3], double q[4]);
void quat_rotate(const double q[4], const double v[3], double v_out[3]);

#endif // quat.h
