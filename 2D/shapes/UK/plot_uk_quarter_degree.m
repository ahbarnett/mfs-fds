% get crude UK coastline at 1/4-degree from m_map, plot, and save
% Barnett 4/19/19

% Needs: m_map from https://www.eoas.ubc.ca/~rich/map.html
addpath ~/matlab/m_map

m_proj('stereographic','lat',55,'lon',0)  % center around UK
h = m_coast('line');           % load it in from m_map's shipped data
z = h.XData + 1i*h.YData;
figure; plot(z); axis equal     % world coastline, as around 5e3 pts

% break into "polygons"... (note NaNs separate lines in this data)
nns = find(isnan(z));  % inds of nans
lens = diff(nns)-2;    % length of all potential polys
ii = find(lens>0);     % inds within poly list of which viable
np = numel(ii);        % # polys
poly = cell(1,np);     % build cell array of poly coords in complx plane
for p=1:np
  i = ii(p);           % ind within viable list
  l = lens(i);         % len of this poly
  zind = nns(i);       % start index within z coord array
  poly{p} = z(zind+1:zind+l);
end
figure; for p=1:np, plot(poly{p},'.-'); hold on; end, axis equal; % all, colors

% extract one poly for UK
zin = 0 - 1i*0.06;
inp = nan(1,np);
for p=1:np
  z=poly{p}; inp(p) = inpolygon(real(zin),imag(zin),real(z),imag(z));
end
puk = find(inp);   % poly index of the UK
zuk = poly{puk};
zuk = zuk-mean(zuk);
zuk = zuk/max(abs(zuk));   % rescale to fit to unit disc
zuk = zuk(end:-1:1);       % make CCW orientation
figure; plot([zuk,zuk(1)],'.-'); axis equal  % plot closed polygon, index labels
hold on; for i=1:numel(zuk), text(real(zuk(i)),imag(zuk(i)),sprintf('%d',i));end
save uk_1poly_88.mat zuk

corangs = diff(angle(diff([zuk,zuk(1:2)])));  % changes in segment angles
corangs = mod(pi-corangs,2*pi);   % make interior corner angles in (0,2pi) range
[maxc,cmax] = max(corangs); [minc,cmin] = min(corangs);  % bad corners
cmax = cmax+1; cmin = cmin+1; % inds of bad corners
fprintf('range of corner angles = [%.3g,%.3g] rad = [%.3g,%.3g] pi\n',minc,maxc,minc/pi,maxc/pi)
figure; plot([zuk,zuk(1)],'.-'); axis equal; hold on;
plot(zuk(cmax),'*'); plot(zuk(cmin),'d'); title('the worst corners')

% could now get the 5 polys making up the UK...
