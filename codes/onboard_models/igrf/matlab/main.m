% Rishav (2021/12/27)
clc
clear
close all

latitude = 28.3949*0;
longitude = 84.1240*0;
radius =  6371.2;
days = 0;

[Bn, Be, Bd] = igrf(latitude, longitude, radius, days);