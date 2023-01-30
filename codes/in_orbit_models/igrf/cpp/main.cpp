// Basic implementation of 'igrf.h'
// Rishav (2021/12/27)

#include <iostream>
#include "igrf.h"

float latitude = 28.3949;
float longitude = 84.1240;
float radius = 6371.2 + 1000;

uint16_t year = 2021;
uint8_t month = 12;
uint8_t day = 29;
uint8_t hour = 0;
uint8_t min = 0;
uint8_t sec = 0;

igrf compass;

int main()
{

    compass.set_date_time(year, month, day, hour, min, sec);
    compass.update(latitude * M_PI / 180.0f, longitude * M_PI / 180.0f, radius);

    // Results
    // Validation: http://www.geomag.bgs.ac.uk/data_service/models_compass/igrf_calc.html
    std::cout << "\r\n~~~ Inputs ~~~" << std::endl;
    std::cout << "Latitude: " << latitude << " degrees" << std::endl;
    std::cout << "Longitude: " << longitude << " degrees" << std::endl;
    std::cout << "Radius: " << radius << " km\r\n";
    std::cout << "YYYYMMDD: " << unsigned(year) << "/" << unsigned(month) << "/" << unsigned(day) << std::endl;
    std::cout << "HHMMSS:" << unsigned(hour) << ":" << unsigned(min) << ":" << unsigned(sec) << "\r\n"
              << std::endl;

    std::cout << "~~~ Outputs ~~~" << std::endl;
    std::cout << "North: " << compass.B_ned[0] << " nT\r\t" << std::endl;
    std::cout << "East: " << compass.B_ned[1] << " nT" << std::endl;
    std::cout << "Down: " << compass.B_ned[2] << " nT\r\n"
              << std::endl;

    std::cout << "Magnitude: " << compass.get_norm() << " nT" << std::endl;
    std::cout << "Horizontal intensity: " << compass.get_horizontal_intensity() << " nT" << std::endl;
    std::cout << "Declination: " << compass.get_declination() * 180.0 / M_PI << " degrees" << std::endl;
    std::cout << "Inclination: " << compass.get_inclination() * 180.0 / M_PI << " degrees\r\n"
              << std::endl;

    return 0;
}
