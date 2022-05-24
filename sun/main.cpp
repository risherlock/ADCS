// Sun vector model
//
// Source:
//      [1] Curtis - Orbital Mechanics For Engineering Students (2020)
//      [2] Vallado - Fundamentals of Astrodynamics and Applications (2013)
//
// Rishav (2021/08/26)

#include <iostream>
#include <iomanip>
#include <math.h>

int main()
{
    // Input
    uint16_t year = 2021;
    uint8_t month = 8;
    uint8_t day = 26;
    uint8_t hour = 0;
    uint16_t minute = 0;
    uint8_t second = 0;

    // Computes Julian day number at 0 UT for any year between 1900 and 2100
    double julian_days = 367 * year -
                         (int)(7 * (year + (int)((month + 9) / 12)) / 4) +
                         (int)(275 * month / 9) +
                         day + 1721013.5 +
                         hour / 24 + minute / 1440 + second / 86400;

    double T = (julian_days - 2451545.0) / 36525;
    double lambda_M = 280.460 + 36000.771 * T;
    double M = 357.5291092 + 35999.05034 * T;
    double lambda_ecliptic = lambda_M + 1.914666471 * sin(M) +
                             0.019994643 * sin(2 * M);
    double psi_ecliptic = 0;
    double epsilon = 23.439291 - 0.0130042 * T;

    double r = 1.000140612 - 0.016708617 * cos(M) - 0.000139589 * cos(2 * M);
    double r_x, r_y, r_z;

    r_x = r * cos(lambda_ecliptic);
    r_y = r * cos(epsilon) * sin(lambda_ecliptic);
    r_z = r * sin(epsilon) * sin(lambda_ecliptic);

    std::cout << std::endl;
    std::cout << "Input" << std::endl;
    std::cout << "=====" << std::endl;
    std::cout << "Year = " << +year << std::endl;
    std::cout << "Month = " << +month << std::endl;
    std::cout << "Day = " << +day << std::endl;
    std::cout << "Hour = " << +hour << std::endl;
    std::cout << "Minute = " << +minute << std::endl;
    std::cout << "Second = " << +second << std::endl;
    std::cout << std::endl;

    std::cout << std::fixed;
    std::cout << std::setprecision(5);
    std::cout << std::endl;
    std::cout << "Output" << std::endl;
    std::cout << "=======" << std::endl;
    std::cout << "norm(r) = " << r << std::endl;
    std::cout << "r_x = " << r_x << std::endl;
    std::cout << "r_y = " << r_y << std::endl;
    std::cout << "r_z = " << r_z << std::endl;
    std::cout << std::endl;

    return 0;
}
