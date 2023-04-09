function [Bn,Be,Bd] = igrf(theta,phi,r,days)
% Calculates magnetic field strength in local spherical coordinates
%
% Inputs:
% theta = Latitude measured in degrees positive from equator
% phi   = Longitude measured in degrees positive east from Greenwich
% r     = Geocentric radius, km
% days  = Decimal days since January 1, 2020
%
% Output:
%   Br, Bt, Bp = B in radial, theta and phi direction
%   Bn, Be, Bd = B in North, East and Down direction
%
% References:
% [1] Davis - Mathematical Modeling of Earth's Magnetic Field (2004)
% [2] Yaguang - Spacecraft Modeling Attitude Determination and Control:
%           Quaternion Based Approach (2019)
% 
% Note:
%   Variable names used in this code follows reference [2].
%
% Test: http://www.geomag.bgs.ac.uk/data_service/models_compass/igrf_calc.html
% Latest IGRF model: https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
%
% Rishav (2020/5/25)

% Check pole to avoid singularities
if (theta>-0.00000001 && theta<0.00000001)
    theta = 0.00000001;
elseif(theta<180.00000001 && theta>179.99999999)
    theta = 179.99999999;
end

deg2rad = pi/180;
theta = (90-theta)*deg2rad;
phi   = phi*deg2rad;
a     = 6371.2; % Radius of Earth, Km

% Read in the g and h Schmidt quasi-normalized coefficients
norm_g = fopen('./igrf_coeffs/norm_g_2020.txt');
norm_h = fopen('./igrf_coeffs/norm_h_2020.txt');
read_h = textscan(norm_h,'%f %f %f %f');
read_g = textscan(norm_g,'%f %f %f %f');

% Unpack the file values to variables
hn = read_h{1}; hm = read_h{2}; h_val = read_h{3}; h_svi = read_h{4};
gn = read_g{1}; gm = read_g{2}; g_val = read_g{3}; g_svi = read_g{4};

% Close files
fclose(norm_g);
fclose(norm_h);

% Gauss coefficients
N = max(gn);
g = zeros(N,N+1);
h = zeros(N,N+1);
for x=1:length(gn)
    g(gn(x),gm(x)+1) = g_val(x) + g_svi(x)*days/365;
    h(hn(x),hm(x)+1) = h_val(x) + h_svi(x)*days/365;
end

% 1. Calculate Legendre polynomials and its derivatives recursively
% 2. Compute magnetic field strength using terms computed in 1
Br  = 0; Bt = 0; Bp = 0;
Pnpnp = 1; Pnpm = 1;
dPnpnp = 0; dPnpm = 0;
for m = 0:N
    for n = 1:N
        if m <= n
            if (n == m)
                % Eqn(5.24b) and (5.25b) ...[2]
                Pnm = sin(theta)*Pnpnp;
                dPnm = sin(theta)*dPnpnp + cos(theta)*Pnpnp;
                
                Pnpnp  = Pnm; Pnpm = Pnpnp; Pnppm = 0;
                dPnpnp = dPnm; dPnpm = dPnpnp; dPnppm = 0;
                
            elseif (n == 1)
                % K = 0
                % Eqn(5.24c) and (5.25c) ...[2]
                Pnm = cos(theta)*Pnpm;
                dPnm = cos(theta)*dPnpm - sin(theta)*Pnpm;
                
                Pnppm = Pnpm; Pnpm = Pnm;
                dPnppm = dPnpm; dPnpm = dPnm;
                
            else % (n != m) && (n > 1)
                % Eqn(5.24e), (5.24c) and (5.25c) ...[2]
                Knm = ((n-1)^2-m^2)/((2*n-1)*(2*n-3));
                Pnm = cos(theta)*Pnpm - Knm*Pnppm;
                dPnm = cos(theta)*dPnpm - sin(theta)*Pnpm - Knm*dPnppm;
                
                Pnppm = Pnpm; Pnpm = Pnm;
                dPnppm = dPnpm; dPnpm = dPnm;
            end
            
            % Eqn(5.14a), (5.14b) and (5.14c) ...[2]
            Br = Br + (a/r)^(n+2)*(n+1)*...
                ((g(n,m+1)*cos(m*phi) + h(n,m+1)*sin(m*phi))*Pnm);
            
            Bt = Bt + (a/r)^(n+2)*...
                ((g(n,m+1)*cos(m*phi) + h(n,m+1)*sin(m*phi))*dPnm);
            
            Bp = Bp + (a/r)^(n+2)*...
                (m*(-g(n,m+1)*sin(m*phi) + h(n,m+1)*cos(m*phi))* Pnm);
        end
    end
end
Bt = -Bt;
Bp = -Bp/sin(theta);

Bn = -Bt;
Be = Bp;
Bd = -Br;
end
