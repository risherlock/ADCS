// Basic Euler angles manipulations.
// 2024-12-22

#ifndef _ATTITUDE_EULER_H_
#define _ATTITUDE_EULER_H_

typedef enum
{
  EULER_123,
  EULER_121,
  EULER_131,
  EULER_132,
  EULER_231,
  EULER_232,
  EULER_212,
  EULER_213,
  EULER_312,
  EULER_313,
  EULER_321,
  EULER_323
} euler_seq_t;

void euler_rotate(const double e[3], const euler_seq_t es, const double v[3], double v_out[3]);

#endif // euler.h
