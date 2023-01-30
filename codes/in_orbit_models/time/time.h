// Time conversions for ADCS
// Rishav (2021/12/28)

class time
{
private:
public:
    // January 1, 1970, 0h UTC
    struct sutc
    {
        uint16_t year = 0;
        uint8_t month = 0;
        uint8_t day = 0;
        uint8_t hour = 0;
        uint8_t min = 0;
        uint8_t sec = 0;
    } utc;

    struct stle_appoch
    {
        uint16_t year = 0;
        float day = 0;
    } tle_appoch;

    time();

    void compute_utc_from_tle_apoch();
    void compute_tle_epoch_from_utc();
    float compute_local_sidereal_time();

    float compute_decimal_days(uint16_t year, uint8_t month, uint8_t day);
};

time::time()
{
}

float time::compute_decimal_days(uint16_t year, uint8_t month, uint8_t day)
{
    float decimal_days = 0;

    // If starting date is not Jan 1, count days to end of year
    int days[] = {0, 31, 59, 90, 120, 151, 182, 212, 243, 273, 304, 334};
    if ((month != 1 || day != 1) && (year < utc.year))
    {
        int isleap = (((year % 4) == 0) && (((year % 100) != 0) || ((year % 400) == 0)));
        int days_in_year = isleap ? 366 : 365;
        decimal_days = days_in_year - (days[month - 1] + day + (month > 2 ? isleap : 0));

        if (year < utc.year)
        {
            year = year + 1;
        }
    }

    // Count the number of days of complete years
    for (uint16_t i = year; i < utc.year; i++)
    {
        int isleap = (((i % 4) == 0) && (((i % 100) != 0) || ((i % 400) == 0)));
        int ith_year_days = isleap ? 366 : 365;
        decimal_days += ith_year_days;
    }

    // Count the number of days in last year
    int isleap = (((utc.year % 4) == 0) && (((utc.year % 100) != 0) || ((utc.year % 400) == 0)));
    double days_past_in_year = (days[utc.month - 1] + utc.day + (utc.month > 2 ? isleap : 0));
    double decimal_hours = utc.hour + (utc.min / 60) + (utc.sec / 3600);
    return decimal_days + days_past_in_year + decimal_hours / 24;
}
