%HydroFill_Ponding_West_Run_I68.m
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
%start timer:
tic;
%Load the input workspace:
load('West_Run_I68_MSB_DEM_Creation.mat','nr','nc','idem','R','tifinfo');
%Set cell size, dx, in meters:
dx=1.0; %for informational perposes (not used).
%set minimum idem delta = 1 in units used (mm or 1/10 mm)
mindelt=int32(1);
%**************************************************************************
nit=int32(2500); %upper limit on number of iterations
%**************************************************************************
%Create water depth matrix:
wdepth=zeros(nr,nc,'int32'); %*********************************************
%Compute stack size estimate equal to 1% of DEM size
stksize=int32(nr*nc*0.01);
%create abstack to store row and column numbers (parallel stacks a & b), &
%fstack to store current fill depths: 
abstack=zeros(2,2,stksize,'int32');
fstack=zeros(stksize,3,'int32');
%Initialize the stack order number:
a=int32(1);
b=int32(2);
%Initialize the Stack counters:
na=int32(0);
nb=int32(0);
nf=int32(0);
%create the cell status flag array:
nbflag=zeros(nr,nc,'uint16');
%nbflag=0 cell is not stored in the nb pstack
%nbflag=1 cell is stored in the nb pstack
%nbflag=2 cell is a new pit or a flat, and stored in the nb pstack
PorF=zeros(1,3,'int32'); %PorF = Pit or Flat criteria.
%PorF(1)= # of neighbor cells delt>0 (>0 is higher than focus cell)
%PorF(2)= # of neighbor cells delt=0
%PorF(3)= # of neighbor cells delt<0 (<0 is lower than focus cell)
LoPlus=int32(0);  %LoPlus = Lowest positive Delta
output=zeros(nit,6,'int32'); %create output matrix with the variables:
%output(nit,1)=k=iteration number
%output(nit,2)=npits=# of pits remaining
%output(nit,3)=nflats=# of cells with a neighbor of equal elev. remaining
%output(nit,4)=na=current size of the "a" stack
%output(nit,5)=nb=current size of the "b" stack
%output(nit,6)=nf=current size of the "f" fillstack
%zero the number of pits and flats:
npits=int32(0);
nflats=int32(0);
%cycle through the entire dem and locate pits and flats
% and store in in the "a" stack.
for n=2:nr-1
    for m=2:nc-1
        LoPlus=99999999;
        PorF(:)=0;
        for nn=-1:1:1
            for mm=-1:1:1
                if nn~=0 || mm~=0
                    delt=idem(n+nn,m+mm)-idem(n,m);
                    if delt>0
                        PorF(1)=PorF(1)+1;
                        if delt<LoPlus
                            LoPlus=delt;
                        end
                    elseif delt==0
                        PorF(2)=PorF(2)+1;
                    elseif delt<0
                        PorF(3)=PorF(3)+1;
                    end
                end
            end
        end
        if PorF(1)==8 %Pit discovered, store cell in the "a"stack.
            na=na+1;
            abstack(a,1,na)=n;
            abstack(a,2,na)=m;
        elseif PorF(3)==0 && PorF(2)>0 %Flat discovered, store in "a" stack
            na=na+1;
            abstack(a,1,na)=n;
            abstack(a,2,na)=m;
        end
    end 
end
%initial loading of the "a" stack is complete, begin iterative filling.
for k=1:nit
    disp(k);
    npits=0;
    nflats=0;
    for stkinc=1:na %stack increment, stkinc.
        n=abstack(a,1,stkinc);
        m=abstack(a,2,stkinc);
        LoPlus=99999999;
        PorF(:)=0;
        for nn=-1:1:1
            for mm=-1:1:1
                if nn~=0 || mm~=0
                    delt=idem(n+nn,m+mm)-idem(n,m);
                    if delt>0
                        PorF(1)=PorF(1)+1;
                        if delt<LoPlus
                            LoPlus=delt;
                        end
                    elseif delt==0
                        PorF(2)=PorF(2)+1;
                    elseif delt<0
                        PorF(3)=PorF(3)+1;
                    end
                end
            end
        end
        if PorF(1)==8 %Pit discovered, raise elevation to LoPlus value:
            nf=nf+1;
            fstack(nf,1)=n;
            fstack(nf,2)=m;
            fstack(nf,3)=LoPlus; %store fill depth for addition to dem.
            wdepth(n,m)=wdepth(n,m)+LoPlus; %add to the fill depth array.
            if nbflag(n,m)==0
                nb=nb+1;
                abstack(b,1,nb)=n;
                abstack(b,2,nb)=m;
            end
            nbflag(n,m)=2;
            npits=npits+1;
        elseif PorF(3)==0 && PorF(2)>0 %Flat discovered,raise by mindelt:
            nf=nf+1;
            fstack(nf,1)=n;
            fstack(nf,2)=m;
            fstack(nf,3)=mindelt;  %store minimum fill depth.
            wdepth(n,m)=wdepth(n,m)+mindelt;
            if nbflag(n,m)==0
                nb=nb+1;
                abstack(b,1,nb)=n;
                abstack(b,2,nb)=m;
            end
            nbflag(n,m)=2;
            nflats=nflats+1;
        end
        if nbflag(n,m)==2 %new pit or flat: check surrounding cells
            for nn=-1:1:1
                for mm=-1:1:1
                    if nn~=0 || mm~=0
                        if n+nn~=1 && n+nn~=nr && m+mm~=1 && m+mm~=nc
                            if nbflag(n+nn,m+mm)==0
                                delt=idem(n+nn,m+mm)- ...
                                (idem(n,m)+fstack(nf,3));
                                if delt==0 || delt==mindelt
                                    nb=nb+1;
                                    abstack(b,1,nb)=n+nn;
                                    abstack(b,2,nb)=m+mm;
                                    nbflag(n+nn,m+mm)=1;
                                end
                            end
                        end
                    end
                end
            end
            nbflag(n,m)=1;
        end
    end
    %add fill depth increases to dem array:
    for nfinc=1:nf
        idem(fstack(nfinc,1),fstack(nfinc,2))= ...
           idem(fstack(nfinc,1),fstack(nfinc,2))+fstack(nfinc,3);
    end
    output(k,1)=k;
    output(k,2)=npits;
    output(k,3)=nflats;
    output(k,4)=na;
    output(k,5)=nb;
    output(k,6)=nf;
    %If nb counter is zero, iterations are complete - terminate the job:
    if nb==0
        break     %breaks continuation of the (for K=1:nit) loop.
    end
    %Switch stack order, including counters. Zero 2nd stack counter:
    if k<nit
    aa=a;
    a=b;
    b=aa;
    na=nb;
    nb=0;
    nf=0;
    nbflag(:,:)=0;
    end
end
save('HydroFill_Ponding_West_Run_I68.mat','nr','nc', ...
    'idem','wdepth','nb','nbflag','nit','R','tifinfo');
filename='HydroFill_Ponding_West_Run_I68_wdepth.tif';
geotiffwrite(filename,wdepth,R,'GeoKeyDirectoryTag', ...
    tifinfo.GeoTIFFTags.GeoKeyDirectoryTag);
toc;
