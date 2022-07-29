// Basic implementation of 'igrf.h'
// Rishav (2021/12/27)

#include <stdio.h>
#include <inttypes.h>
#include "igrf.h"

float latitude = 28.3949;
float longitude = 84.1240;
float radius = 6371.2 + 1000;

uint16_t year = 2021;
uint8_t month = 12;
uint8_t day = 29;
uint8_t hour = 10;
uint8_t min = 14;
uint8_t sec = 45;

int main()
{

    igrf_set_date_time(year, month, day, hour, min, sec);
    igrf_update(latitude * PI / 180.0f, longitude * PI / 180.0f, radius, 1);

    // Results
    // Validation: http: // www.geomag.bgs.ac.uk/data_service/models_compass/igrf_calc.html
    printf("\r\n~~~ Inputs ~~~\r\n");
    printf("Latitude: %f degrees\r\n", latitude);
    printf("Longitude: %f degrees\r\n", longitude);
    printf("Altitude: %f km\r\n", radius);
    printf("YYMMDD: %u/%u/%u\r\n", year, month, day);
    printf("HHMMSS: %u/%u/%u\r\n\r\n", hour, min, sec);

    printf("~~~ Outputs ~~~\r\n", latitude);
    printf("North (X): %f nT\r\n", B_ned[0]);
    printf("East (Y): %f  nT\r\n", B_ned[1]);
    printf("Down (Z): %f  nT\r\n\r\n", B_ned[2]);

    printf("Magnitude(F): %f nT\r\n", igrf_get_norm());
    printf("Horizontal intensity(H): %f  nT\r\n", igrf_get_horizontal_intensity());
    printf("Declination(D): %f  degrees\r\n", igrf_get_declination() * 180.0 / PI);
    printf("Inclination (I): %f  degrees\r\n\r\n", igrf_get_inclination() * 180.0 / PI);

    return 0;
}

