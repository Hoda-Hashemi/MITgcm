-----OUTDATED----
% Generate initial free-surface anomaly for vortexSphere_mitgcm_referenceCase
% Output: eta_init.bin (big-endian real*4)

nx = 1440;
ny = 720;

amp = 10.0;           % meters (same max anomaly scale as adjustment.cs-32x32x1)
lat0 = 0.0;           % deg
lon0 = 180.0;         % deg
cosineRadius = 0.04;  % nondimensional radius, as in adjustment.cs-32x32x1

lon = (0.5:nx)/nx*360.0;
lat = (0.5:ny)/ny*180.0 - 90.0;
[LON, LAT] = ndgrid(lon, lat);

% periodic longitude distance from anomaly center (degrees)
dlon = abs(LON-lon0);
dlon = min(dlon, 360.0-dlon);
dlat = LAT-lat0;

% Same compact cosine profile as adjustment.cs-32x32x1 gendata.m
x = 0.25*(dlon/360.0);
y = 0.25*(dlat/180.0);
R = sqrt(x.^2 + y.^2);
eta = amp*(0.5 + 0.5*cos(pi*min(R,0*R+cosineRadius)/cosineRadius));

% Depth-aware: initialize only where bathymetry indicates ocean.
fid = fopen('bathymetry.bin','r','b');
bathy = fread(fid,[nx ny],'real*4');
fclose(fid);
if min(bathy(:)) < 0
    ocean = bathy < 0;
else
    ocean = bathy > 0;
end
eta(~ocean) = 0.0;

fid = fopen('eta_init.bin','w','b');
fwrite(fid, eta, 'real*4');
fclose(fid);

