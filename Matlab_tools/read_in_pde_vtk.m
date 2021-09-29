function [X,Y,Z,M] = read_in_pde_vtk(varargin)

% Example use: [X,Y,Z,M] = read_in_pde_vtk('filename','filename.vtk')

for vac = 1:2:length(varargin)
    eval([varargin{vac},' = varargin{vac+1};']);
end

fid = fopen(filename);

str = fgetl(fid);
while ~contains(str,'POINTS')
    str = fgetl(fid);
end

npts = sscanf(str,'POINTS %d');
XY = textscan(fid,' %f %f %f',npts);
X = XY{1};
Y = XY{2};

while ~contains(str,'CELLS')
    str = fgetl(fid);
end

ncells = sscanf(str,'CELLS %d');
M = textscan(fid,' %f %f %f %f',ncells);
M = [M{2},M{3},M{4}];
M = M+1;

while ~contains(str,'POINT_DATA')
    str = fgetl(fid);
end

str = fgetl(fid);
str = fgetl(fid);

Z = textscan(fid,' %f %f %f %f %f %f %f %f %f %f',npts);
Z = [Z{1},Z{2},Z{3},Z{4},Z{5},Z{6},Z{7},Z{8},Z{9},Z{10}];
Z = Z';
Z = Z(:);
Z = Z(1:npts);

fclose(fid);

end