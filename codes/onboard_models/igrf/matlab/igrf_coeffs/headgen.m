% Based on Davis - Mathematical Modeling of Earth’s Magnetic Field(2004)
% 2024-12-03

clc
clear
close all

igrf_start_year = 2025;

fid = fopen(strcat('igrf', num2str(igrf_start_year), '.csv'), 'r');
data = textscan(fid, '%s %f %f %f %f', 'Delimiter', ',', 'HeaderLines', 0);
fclose(fid);

gh = data{1};
n = data{2};
m = data{3};
val = data{4};
sv = data{5};

N = 13;
N = max(n);
g = zeros(N, N + 1);
h = zeros(N, N + 1);
hsv = zeros(N, N + 1);
gsv = zeros(N, N + 1);

for x = 1 : length(gh)
  if strcmp (gh(x), 'g')
    g(n(x), m(x) + 1) = val(x);
    gsv(n(x), m(x) + 1) = sv(x);
  else h(n(x), m(x) + 1) = val(x);
    hsv(n(x), m(x) + 1) = sv(x);
  end
end

count = 1;
S = zeros(N, N + 1);

for n = 1:N
  for m = 0:n
    if m > 1
      S(n, m+1) = S(n,m) * ((n-m+1)/(n+m))^0.5;
    elseif m > 0
      S(n, m + 1) = S(n, m) * (2 * (n - m + 1) / (n + m)) ^ 0.5;
    elseif n == 1
      S(n, 1) = 1;
    else
      S(n, 1) = S(n - 1, 1) * (2 * n - 1) / (n);
    end

    g_out(count) = g(n, m + 1) * S(n, m + 1);
    gsv_out(count) = gsv(n, m + 1) * S(n, m + 1);
    h_out(count) = h(n, m + 1) * S(n, m + 1);
    hsv_out(count) = hsv(n, m + 1) * S(n, m + 1);
    count = count + 1;
  end
end

count = count - 1;

header_file = strcat('igrf', num2str(igrf_start_year), '.h');
fid = fopen(header_file, 'w');

fprintf(fid, '// This file is generated using headgen.m\n\n');
fprintf(fid, strcat('#ifndef _IGRF', num2str(igrf_start_year), '_H_\n'));
fprintf(fid, strcat('#define _IGRF', num2str(igrf_start_year), '_H_\n\n'));
fprintf(fid, '#define IGRF_DEGREE %d\n', N);
fprintf(fid, '#define IGRF_END_YEAR %d\n', igrf_start_year + 5);
fprintf(fid, '#define IGRF_START_YEAR %d\n\n', igrf_start_year);
fprintf(fid, '// Schmidt quasi-normalized coefficients\n');
fprintf(fid, 'static const float g_val[] = {');
fprintf(fid, '%f, ', g_out(1 : end - 1));
fprintf(fid, '%f};\n', g_out(end));
fprintf(fid, 'static const float h_val[] = {');
fprintf(fid, '%f, ', h_out(1 : end - 1));
fprintf(fid, '%f};\n\n', h_out(end));
fprintf(fid, '// Secular variations\n');
fprintf(fid, 'static const float g_sv[] = {');
fprintf(fid, '%f, ', gsv_out(1 : end - 1));
fprintf(fid, '%f};\n', gsv_out(end, 4));
fprintf(fid, 'static const float h_sv[] = {');
fprintf(fid, '%f, ', hsv_out(1 : end - 1));
fprintf(fid, '%f};\n\n', hsv_out(end));
fprintf(fid, '#endif\n');
fclose(fid);

fprintf('%s generated successfully!\n', header_file);
