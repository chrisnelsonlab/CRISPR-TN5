%--------------------------------------------------------------------------
%Plot profile of AAV integration
%Code author: Chris Nelson
%Public Version 1.0
%Updated 2018-12-30
%--------------------------------------------------------------------------
function AAV_Plot
clc
clearvars
disp('starting');
d=dir('output\heartAAV\*aav.dat'); %Change to the correct directory
AAVsize1=4465;  %Change to size of AAV genome
AAVsize2=3665;  %Change to size of  AAV genome
%SET UP FOUR DATA STRUCTURES
AAVarray1 = zeros(AAVsize1,1);
AAVarray2 = zeros(AAVsize1,1);
AAVarray3 = zeros(AAVsize2,1);
AAVarray4 = zeros(AAVsize2,1);

%Loop iterates through files to tabulate AAV genome integration profile.
%The file format is based on the output from CRISPR_TN5.m
for j=1:length(d)
    fid = fopen(strcat('output\',d(j).name));  %This line will need to change if the output folder has a different name
    disp(strcat('Opening:',d(j).name));
    while true
        thisline = fgetl(fid);
        if ~ischar(thisline)
           break;
        end
        
        C =strsplit(thisline);
        if(length(C) == 5)
            AAVGenome = string(C(1,3));
            AAVstart = cellfun(@str2num,C(1,4));
            AAVend = cellfun(@str2num,C(1,5));
           
            if(AAVGenome == "AAVGENOME1")
                AAVarray1 = incrementAAV(AAVarray1,AAVstart, AAVend);
            elseif(AAVGenome == "AAVGENOME1R")
                AAVarray2 = incrementAAV(AAVarray2,AAVstart, AAVend);
            elseif(AAVGenome == "AAVGENOME2")
                AAVarray3 = incrementAAV(AAVarray3,AAVstart, AAVend);
            elseif(AAVGenome == "AAVGENOME2R")
                AAVarray4 = incrementAAV(AAVarray4,AAVstart, AAVend);
            else
            end
        end
 
    end
     
end

AAVarray2flip = flip(AAVarray2); %Array 2 will be reverse complement
AAVarray4flip = flip(AAVarray4); %Array 4 will be reverse complement
AAVarray1Total = AAVarray1 + AAVarray2flip; %Add two arrays together
AAVarray2Total = AAVarray3 + AAVarray4flip; %Add two arrays together

%Plot the AAV integration data
subplot(3,4,[1,2]);
plot(AAVarray1)
ylabel('AAVGenome1 Sense')


subplot(3,4,[5,6]);
plot(AAVarray2)
ylabel('AAVGenome1 Antisense')

subplot(3,4,[3,4]);
plot(AAVarray3)
ylabel('AAVGenome2 Sense')

subplot(3,4,[7,8]);
plot(AAVarray4)
ylabel('AAVGenome2 Antisense')

subplot(3,4,[9,10]);
plot(AAVarray1Total)
ylabel('AAVGenome1 Total')

subplot(3,4,[11,12]);
plot(AAVarray2Total)
ylabel('AAVGenome2 Total')

fclose('all');
end

%Function to tabulate the AAV integration into one array.
function AAVarray = incrementAAV(AAVarray, AAVstart, AAVend)
    for m = AAVstart:AAVend
        AAVarray(m) = AAVarray(m)+1;
    end
end