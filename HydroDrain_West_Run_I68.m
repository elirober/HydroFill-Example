%HydroDrain_West_Run_I68.m
%**************************************************************************
% MIT License
% 
% Copyright (c) [2025] [Robert N. Eli]
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
%**************************************************************************
clear;
clc;
tic;
%Set grid cell size = 1 meter:
dx=1.0;
%**************************************************************************
load('HydroFill_West_Run_I68.mat','nr','nc','IPMult', ...
        'idem','nbflag','ncul','maxncul','CIO','CIOflt','R','tifinfo');
%**************************************************************************
% CIO(1:ncul,1)=culvert modified # (culvert # +100 = 101,102,103,etc.)
% CIO(1:ncul,2)=culvert inlet n
% CIO(1:ncul,3)=culvert inlet m
% CIO(1:ncul,4)=culvert outlet n
% CIO(1:ncul,5)=culvert outlet m
% CIO(1:ncul,6)=culvert type: 1=circular cross-section, 2=box culvert
% CIO(1:ncul,7-8)=left open for future use.
%**************************************************************************
% CIOflt(1:ncul,1)=CulLength, m
% CIOflt(1:ncul,2)=CulSlope, m/m
% CIOflt(1:ncul,3)=culvert diameter, or width, m
% CIOflt(1:ncul,4)=ratio of width to height, m/m (circular=1.0)
% CIOflt(1:ncul,5)=culvert Mannings n
% CIOflt(1:ncul,6-8)=left open for future use
%**************************************************************************
%Convert idem (int32) to fdem (double). Note: f="floating-point"
fdem=(double(idem))/IPMult;
%Create the neighborhood box index matrix ma(6,4)
ma=int32([-1,1,-1,1;1,-1,1,-1;1,-1,1,-1;-1,1,1,-1;1,-1,-1,1;1,-1,-1,1]);
%Initialize "ma" selection counter:
ir=int32(0);
%Initialize Flow Length Index:
FLindx=int32(0);
%Initialize Culvert Length:
CulLength=double(0);
%Initialize Culvert Slope:
CulSlope=double(0);
%Set minimum culvert slope:
Culminslp=double(0.001);
%Initialize row and column offsets for downslope cell:
smaxnn=int32(0);
smaxmm=int32(0);
%Initialize storage for the row & column location of the downslope cell,
%and storage for maximum slope and ground distance.
dcell=zeros(nr,nc,2,'int32'); % 1=n,2=m (downslope cell)
dcellFS=zeros(nr,nc,'single'); % Friction Slope to downslope cell.
dcellGD=zeros(nr,nc,'single'); % Ground Distance to downslope cell.
bn=zeros(nr,nc,'int32'); %create basin number array and assign zero values.
basinct=int32(1); %basin number incrementer for boundary cell assignment.
n=1; %assign bn values in a clockwise direction around the box (2=1st val)
for m=1:nc
    basinct=basinct+1;
    bn(n,m)=basinct;
end
m=nc;
for n=2:nr
    basinct=basinct+1;
    bn(n,m)=basinct;
end
n=nr;
for m=nc-1:-1:1
    basinct=basinct+1;
    bn(n,m)=basinct;
end
m=1;
for n=nr-1:-1:2
    basinct=basinct+1;
    bn(n,m)=basinct;
end
%Compute downslope cells and corresponding distances and slopes:
for n=2:nr-1
    for m=2:nc-1
        if nbflag(n,m)==0   %no culvert entrance
            smax=0.0; % zero the map slope (ms) maximum
            smaxnn=0;
            smaxmm=0;
            smaxFS=0.0;
            smaxMS=0.0;
            if ir==4
                ir=1;
            else
                ir=ir+1;
            end
            for nn=ma(1,ir):ma(2,ir):ma(3,ir)
                for mm=ma(4,ir):ma(5,ir):ma(6,ir)
                    if nn~=0 || mm~=0
                        FLindx=nn*mm;
                        Ediff=fdem(n,m)-fdem(n+nn,m+mm);
                        if Ediff>0.0
                            if FLindx==0
                                mapdist=dx; %map distance to downslope cell
                                gd=sqrt(Ediff^2+mapdist^2); %=ground dist.
                            else
                                mapdist=1.4142*dx; %diagonal map dist.
                                gd=sqrt(Ediff^2+mapdist^2); %=ground dist.
                            end
                            fs=Ediff/gd; %Friction Slope
                            ms=Ediff/mapdist; %map slope
                            if ms>smax
                                smax=ms;
                                smaxnn=nn;
                                smaxmm=mm;
                                smaxMS=smax; %Map Slope to smax cell
                                smaxFS=fs; %Friction Slope to smax cell
                                smaxGD=gd; %Ground Distance to smax cell
                            end
                        end
                    end
                end
            end
            if smaxnn==0 && smaxmm==0
                disp('Error: Pit or Flat');
                disp(n);
                disp(m);
                return;
            end
            if smaxFS==0 || smaxMS==0
                disp('Error: zero slope smaxFS or zero slope smaxMS');
                disp(n);
                disp(m);
                return;
            end
            dcell(n,m,1)=n+smaxnn;
            dcell(n,m,2)=m+smaxmm;
            dcellFS(n,m)=smaxFS;
            dcellGD(n,m)=smaxGD;
        elseif nbflag(n,m)>100  %culvert entrance, set dcell at outlet:
            dcell(n,m,1)=CIO((nbflag(n,m)-100),4);
            dcell(n,m,2)=CIO((nbflag(n,m)-100),5);
            mapdist=sqrt(double((n-dcell(n,m,1))^2+(m-dcell(n,m,2))^2));
            Ediff=fdem(n,m)-fdem(dcell(n,m,1),dcell(n,m,2));
            if Ediff<=0.0
                CulLength=mapdist;
                CulSlope=Culminslp;
            else
                CulLength=sqrt(Ediff^2+mapdist^2);
                CulSlope=Ediff/mapdist;
            end
            dcellFS(n,m)=CulSlope;
            dcellGD(n,m)=CulLength;
            CIOflt((nbflag(n,m)-100),1)=CulLength;
            CIOflt((nbflag(n,m)-100),2)=CulSlope;
        else
            disp('Error: invalid nbflag(n,m) value')
            disp(n);
            disp(m);
            return;
        end
    end
end
%Begin Drainage Accumulation using the cell index method (cell area =1)
%total index area accumulates to the downstream cell along the max slope
%until reaching the grid boundary.
di=ones(nr,nc,'int32'); %drainage index matrix: initialize with value 1.
fthread=zeros((nr+nc)*10,2,'int32'); %flow thread storage stack.
ftct=int32(0); %zero flow thread count.
for n=2:nr-1
    disp(n);
    for m=2:nc-1
        ftct=ftct+1;
        fthread(ftct,1)=n;
        fthread(ftct,2)=m;
        nn=dcell(n,m,1);
        mm=dcell(n,m,2);
        di(nn,mm)=di(nn,mm)+1;
        ftct=ftct+1;
        fthread(ftct,1)=nn;
        fthread(ftct,2)=mm;
        while nn~=1 && nn~=nr && mm~=1 && mm~=nc
            nd=dcell(nn,mm,1);
            md=dcell(nn,mm,2);
            di(nd,md)=di(nd,md)+1;
            nn=nd;
            mm=md;
            ftct=ftct+1;
            fthread(ftct,1)=nn;
            fthread(ftct,2)=mm;
        end
        basin=bn(nn,mm);
        for k=1:ftct-1
            if bn(fthread(k,1),fthread(k,2))==basin
                break;
            else
                bn(fthread(k,1),fthread(k,2))=basin;
            end
        end
        ftct=0;
    end
end
%sum the drainage di to check against basin area in cell units:
diSum=int32(0);
basinarea=nr*nc;
n=1;
for m=1:nc
    diSum=diSum+di(n,m);
end
n=nr;
for m=1:nc
    diSum=diSum+di(n,m);
end
m=1;
for n=2:nr-1
    diSum=diSum+di(n,m);
end
m=nc;
for n=2:nr-1
    diSum=diSum+di(n,m);
end
%**************************************************************************
% Output drainage index (di) and basin number (bn) raster .tif files:
%**************************************************************************
save('HydroDrain_West_Run_I68.mat', ...
    'dx','nr','nc','di','bn','fdem','dcell','dcellFS','dcellGD','ncul', ...
    'R','tifinfo');
filename='HydroDrain_West_Run_I68_di.tif';
geotiffwrite(filename,di,R,'GeoKeyDirectoryTag', ...
    tifinfo.GeoTIFFTags.GeoKeyDirectoryTag);
filename='HydroDrain_West_Run_I68_bn.tif';
geotiffwrite(filename,bn,R,'GeoKeyDirectoryTag', ...
    tifinfo.GeoTIFFTags.GeoKeyDirectoryTag);
toc;
