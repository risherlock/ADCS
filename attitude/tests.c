#include <stdio.h>
#include "attitude.h"

void quaternion_tests(void);

int main(void)
{
  return 0;
}

void quaternion_tests(void)
{
  double q1[4];
  quat_unit(q1);

  double q2[4];
  const double psi = 1.27;
  const double v[3] = {1.0, 0.0, 0.0};
  quat_axis_angle(psi, v, q2);

  double v_rot[3];
  quat_rotate(q2, v, v_rot);
}
