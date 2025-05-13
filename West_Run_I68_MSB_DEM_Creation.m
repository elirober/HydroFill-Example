%West_Run_I68_MSB_DEM_Creation.m
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
clear
clc
%**************************************************************************
%import the DEM .tif file (float32 DEM, 2179x2015, 1m pixels):
[dem,R]=readgeoraster('West_Run_I68_Fill_corrected.tif');
tifinfo= geotiffinfo('West_Run_I68_Fill_corrected.tif');
%**************************************************************************
demsize=size(dem); %extract array size from the dem.
nr=int32(demsize(1,1)); %number of rows in the dem array.
nc=int32(demsize(1,2)); %number of columns in the dem array.
%**************************************************************************
%read in the Microsoft Buildings raster file for the West Run I68 fill.
% (float32 with 0.0 for no data and building number for structures: 1.0
% through 302.0, maximum). Some very small structures are not burned-in, or
% they are too small, so the actual structures to be included must be
% counted. Additionally, the float32 is converted to int32:
msb=int32(imread('WV_MSB_WestRun_I68_clip.tif'));
%**************************************************************************
%scan msb to compute the building areas and gather elevation statistics:
%**************************************************************************
maxid=int32(302); %load the building structure maximum id number.
msbstats=zeros(maxid,2,'int32');
%msbstats(1:maxid,1)=building structure id number (= 0 if missing or small)
%msbstats(1:maxid,2)=building area grid cell count (1x1 m grid cells)
%preload the building structure id numbers:
for k=1:maxid
    msbstats(k,1)=k;
end
msbelev=zeros(maxid,4,'single');
%msbelev(1:maxid,1)=foundation minimum elevation
%msbelev(1:maxid,2)=foundation maximum elevation
%msbelev(1:maxid,3)=foundation elevation difference
%msbelev(1:maxid,4)=roof elevation for the particular structure
%preload min and max elevations:
msbelev(1:maxid,1)=99999.9;
msbelev(1:maxid,2)=0.0;
%accumulate the building structure areas and update foundation elevations:
for n=1:nr
    for m=1:nc
        if msb(n,m)>0
            msbstats(msb(n,m),2)=msbstats(msb(n,m),2)+1;
            if dem(n,m)<msbelev(msb(n,m),1)
                msbelev(msb(n,m),1)=dem(n,m);
            elseif dem(n,m)>msbelev(msb(n,m),2)
                msbelev(msb(n,m),2)=dem(n,m);
            end
        end
    end
end
%filter out missing building structures and those < 10 sq.m.:
for k=1:maxid
    if msbstats(k,2)<10
        msbstats(k,1)=0;
    end
end
%compute building structure foundation elevation difference:
for k=1:maxid
    if msbstats(k,1)>0
        msbelev(k,3)=msbelev(k,2)-msbelev(k,1);
    else
        msbelev(k,1)=0.0;
        msbelev(k,2)=0.0;
        msbelev(k,3)=0.0;
    end
end
%compute the roof elevation:
for k=1:maxid
    if msbstats(k,1)>0
        msbelev(k,4)=msbelev(k,2)+(2+log10(single(msbstats(k,2))));
    else
        msbelev(k,4)=0.0;
    end
end
%generate the new dem with building structure roof elevations added:
for n=1:nr
    for m=1:nc
        if msb(n,m)>0
            if msbstats(msb(n,m),1)>0
                dem(n,m)=msbelev(msb(n,m),4);
            end
        end
    end
end
%Convert dem to idem for output within the .mat file:
%Load the Integer Precision Multiplier (IPMult) for 1 mm unit depth.
IPMult=1000.0;
idem=int32(dem*IPMult);
%**************************************************************************
%output dem as a .tif file, along with its geospatial data:
filename='West_Run_I68_MSB_DEM.tif';
geotiffwrite(filename,dem,R,'GeoKeyDirectoryTag', ...
    tifinfo.GeoTIFFTags.GeoKeyDirectoryTag);
%**************************************************************************
save('West_Run_I68_MSB_DEM_Creation.mat','nr','nc','maxid','msb', ...
    'msbstats','msbelev','dem','idem','IPMult','R','tifinfo');
