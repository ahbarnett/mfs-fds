% UK in shapefile format, from: http://www.diva-gis.org/gdata
% Needs: m_map from https://www.eoas.ubc.ca/~rich/map.html

addpath ~/matlab/m_map

s = m_shaperead('diva-gis/GBR_adm0');  % adm0 = coastal bdries
% adm2 = county bdries
ll = s.ncst{1};     % lat and long, as n*2 array
m_proj('stereographic','lat',55,'lon',0)  % center around UK
[x,y] = m_ll2xy(ll(:,1),ll(:,2));   % convert to xy via the projection
z = x+1i*y;
figure; plot(z); axis equal   % entire UK coastline, as 4e5 horiz & vert lines.
