%West_Run_I68_MSB_Culvert_Import.m
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
load('West_Run_I68_MSB_DEM_Creation.mat','nr','nc','R','tifinfo');
%**************************************************************************
%initialize number of culverts to be stored, ncul:
ncul=int32(0);
%set culflag = 0 for no culverts, or =1 if culverts are to be processed:
culflag=int32(1);
%**************************************************************************
%create the cell status nbflag array for export as a geotiff format file.
nbflag=zeros(nr,nc,'int32');
%import the Culvert_IO file with matching positive and negative int32 #'s
% and the rest = 0, if culflag>0.
if culflag>0
    culvtIO=imread('West_Run_I68_Culvert_IO.tif');
    %count the number of culverts in culvtIO:
    for n=1:nr
        for m=1:nc
            if culvtIO(n,m)>0
                ncul=ncul+1;
            end
        end
    end
else
    culvtIO=zeros(nr,nc,'int32');
end
%**************************************************************************
%create the culvert Inlet-Outlet location array with #'s + 100
maxncul=int32(99); %Set maximum CIO array length (maximum # of culverts)
CIO=zeros(maxncul,8,'int32'); % Culvert data table with length maxncul
% CIO(1:ncul,1)=culvert modified # (culvert # +100 = 101,102,103,etc.)
% CIO(1:ncul,2)=culvert inlet n
% CIO(1:ncul,3)=culvert inlet m
% CIO(1:ncul,4)=culvert outlet n
% CIO(1:ncul,5)=culvert outlet m
% CIO(1:ncul,6)=culvert type: 1=circular cross-section, 2=box culvert
% CIO(1:ncul,7-8)=left open for future use.
%**************************************************************************
%create a companion culvert floating point table for storage of culvert
%descriptive information that must be stored as decimal.
CIOflt=zeros(maxncul,8,'double');
% CIOflt(1:ncul,1)=CulLength, m
% CIOflt(1:ncul,2)=CulSlope, m/m
% CIOflt(1:ncul,3)=culvert diameter, or width, m
% CIOflt(1:ncul,4)=ratio of width to height, m/m (circular=1.0)
% CIOflt(1:ncul,5)=culvert Mannings n
% CIOflt(1:ncul,6-8)=left open for future use
%**************************************************************************
% load the CIO table with culvert data:
for n=1:nr
    for m=1:nc
        if culvtIO(n,m)>0   %if inlet: store inlet # +100 & row=n,col=m
            CIO(culvtIO(n,m),1)=culvtIO(n,m)+100;
            CIO(culvtIO(n,m),2)=n;
            CIO(culvtIO(n,m),3)=m;
            CIO(culvtIO(n,m),6)=1; %assign culvert type as circular
        elseif culvtIO(n,m)<0   %if culvert exit: store row=n,col=m
            CIO(abs(culvtIO(n,m)),4)=n;
            CIO(abs(culvtIO(n,m)),5)=m;
        end
    end 
end
%load the nbflag array with the culvert inlet (+) and outlet (-), 101,102..
for nn=1:maxncul
    if CIO(nn,1)~=0
        nbflag(CIO(nn,2),CIO(nn,3))=CIO(nn,1);
        nbflag(CIO(nn,4),CIO(nn,5))=-CIO(nn,1);
    end
end
%**************************************************************************
%Output nbflag as a .tif file in order to export nbflag along with it's
% geospatial data. Zeros except for Culvert I-O: +101,-101,+102,-102, etc.
filename='West_Run_I68_nbflag.tif';
geotiffwrite(filename,nbflag,R,'GeoKeyDirectoryTag', ...
    tifinfo.GeoTIFFTags.GeoKeyDirectoryTag);
%**************************************************************************
save('West_Run_I68_MSB_Culvert_Import.mat','nbflag','nr','nc', ...
    'ncul','maxncul','CIO','CIOflt','R','tifinfo');
