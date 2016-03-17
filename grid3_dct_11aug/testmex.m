% Author: Nick Zwart
% Date: 2011 aug 18
% Rev: 2011 aug 21
% A simple test of the grid3 mex compilation.

'Test Grid3:'
'    read crds'
% nVec: 3
% dim1: 1812
% dim2: 3
% dim3: 79
fid = fopen('spi_crds.raw');
    tmp = squeeze(fread(fid,inf, 'float32'));
    fclose(fid);
    size(tmp)
    tmp = reshape(tmp,[3,1812,3,79]);
    size(tmp) 
    crds = tmp;

'   read pre-weights'
% nVec: 1
% dim1: 1812
% dim2: 3
% dim3: 79
fid = fopen('soln_DCF.raw');
    tmp = squeeze(fread(fid,inf, 'float32'));
    fclose(fid);
    size(tmp)
    tmp = reshape(tmp,[1812,3,79]);
    size(tmp) 
    DCF = tmp;

'   read k-space data'
% nVec: 2
% dim1: 1812
% dim2: 3
% dim3: 79
fid = fopen('sim_kspc.raw');
    tmp = squeeze(fread(fid,inf, 'float32'));
    fclose(fid);
    size(tmp)
    tmp = reshape(tmp,[2,1812,3,79]);
    size(tmp) 
    data = tmp;

'   read solution image for comparison'
% nVec: 1
% dim1: 50
% dim2: 50
% dim3: 50
fid = fopen('soln_image.raw');
    tmp = squeeze(fread(fid,inf, 'float32'));
    fclose(fid);
    size(tmp)
    tmp = reshape(tmp,[50,50,50]);
    size(tmp) 
    soln_image = tmp;

'   Grid3 params:'
effMtx    = 75  % 1.5x the intended supported matrix
numThread = 8   % or 1 for no pthread lib exec

'   start Grid3 calc'
gdata = grid3_MAT(data,crds,DCF,effMtx,numThread);
size(gdata)

'   make a rolloff kernel'
delta = [1.0, 0.0]
k_not = [0.0, 0.0, 0.0]
DCF   = [1.0]
numThread = 1 % only have 1 data point
rokern = grid3_MAT(delta',k_not',DCF,effMtx,numThread);

'   fft into image space'
% DATA
% change to complex, fft, then shift
gdata = squeeze(gdata(1,:,:,:) + 1j*gdata(2,:,:,:));
gdata = fftn(gdata);
gdata = fftshift(gdata,1);
gdata = fftshift(gdata,2);
gdata = fftshift(gdata,3);

% ROLLOFF
% change to complex, shift, then fft
rokern = squeeze(rokern(1,:,:,:) + 1j*rokern(2,:,:,:));
rokern = fftn(rokern);
rokern = fftshift(rokern,1);
rokern = fftshift(rokern,2);
rokern = fftshift(rokern,3);
rokern = abs(rokern);

'   apply rolloff and crop'
gdata(rokern > 0) = gdata(rokern > 0) ./ rokern(rokern > 0);
xs = floor(effMtx/2 - effMtx/1.5/2)+1;
xe = floor(effMtx/2 + effMtx/1.5/2);
gdata = gdata(xs:xe,xs:xe,xs:xe);
gdata = single(abs(gdata)); % magnitude, float32

'   find difference between calc and soln'
diff = (gdata-soln_image).^2;
sum(sum(sum(diff)))

'   write output file'
tmp = gdata;
    size(tmp)
    fid = fopen('final_image.raw','w');
    fwrite(fid,tmp,'float32');
    fclose(fid);

