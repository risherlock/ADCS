#include <iostream>
#include "time.h"

int main()
{
    time curr_time;
    curr_time.utc.year = 2021;
    curr_time.utc.month = 12;
    curr_time.utc.day = 31;
    curr_time.utc.hour = 0;
    curr_time.utc.min = 0;
    curr_time.utc.sec = 0;

    float decimal_days = curr_time.compute_decimal_days(2020, 1, 1);
    std::cout << decimal_days << std::endl;
}
