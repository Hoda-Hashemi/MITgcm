% Generate initial free-surface anomaly for Hashemi_primitive_ref
% Output: eta_init.bin (big-endian real*4)

nx = 1440;
ny = 720;

amp = 0.20;       % meters
sigLat = 10.0;    % deg
sigLon = 20.0;    % deg
lat0 = 0.0;       % deg
lon0 = 180.0;     % deg

lon = (0.5:nx)/nx*360.0;
lat = (0.5:ny)/ny*180.0 - 90.0;
[LON, LAT] = ndgrid(lon, lat);

% periodic longitude distance in degrees
 dlon = abs(LON-lon0);
 dlon = min(dlon, 360.0-dlon);

eta = amp*exp(-( (LAT-lat0).^2/(2*sigLat^2) + dlon.^2/(2*sigLon^2) ));
eta = eta - mean(eta(:));

fid = fopen('eta_init.bin','w','b');
fwrite(fid, eta, 'real*4');
fclose(fid);
