// Basic Direction Cosine Matrix manipulations.
// 2024-12-23

#ifndef _ATTITUDE_DCM_H_
#define _ATTITUDE_DCM_H_

void dcm_unit(double r[3][3]);
void dcm_prod(const double r1[3][3], const double r2[3][3]);
void dcm_rotate(const double r[3][3], const double v[3], double v_out[3]);

#endif // dcm.h
