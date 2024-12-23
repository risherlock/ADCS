// Maps one attitude parameterization to the other.
// 2024-12-23

#ifndef _ATTITUDE_ROTATION_H_
#define _ATTITUDE_ROTATION_H_

#include "euler.h"

void quat_to_dcm(const double q[4], double r[3][3]);
void dcm_to_quat(const double r[3][3], double q[4]);

void quat_to_euler(const double q[4], double e[3], const euler_seq_t es);
void euler_to_quat(const double e[3], const euler_seq_t es, double q[3]);

void euler_to_dcm(const double e[3], const euler_seq_t es, double r[3][3]);
void dcm_to_euler(const double r[3][3], double e[3], const euler_seq_t es);

#endif // rotation.h
