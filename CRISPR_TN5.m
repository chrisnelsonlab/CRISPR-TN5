%--------------------------------------------------------------------------
%Analysis of multiplex-CRISPR samples enriched with TN5 (NEXTERA) and PCR
%Code author: Chris Nelson
%Public Version 1.0
%Updated 2018-12-27
%--------------------------------------------------------------------------
function main
clc
clearvars

%Create a summary file to append data for analysis. 
Summaryoutputfile = 'READCOUNTS'; %Name of summary output file. Data appends to this file
fileID = fopen(strcat(Summaryoutputfile,'_quantiative.dat'),'a'); %Open summary output file
fprintf(fileID,strcat('Started: ',string(datetime('now')),'\n')); %Write date and time
fprintf(fileID,'Filename\t Unedited\t Indel\t NHEJ\t Substitution\t Delteion\t Deletion with indel\t Inversion\t Inversion with indel\t AAV\t Mystery\n'); %Write data headers
fclose(fileID); %Close file

%-------------------EDIT HERE----------------------------------------------
%TN5sort function has the following input (READ1, READ2, Output file name, Summaryoutput file (from above), parameter file)
%Example commands. Copy and paste as many commands as needed to process data
TN5sort('EX_READ1_truncated.fastq','EX_READ2_truncated.fastq','EXAMPLEg2',Summaryoutputfile,'pargRNA1.dat'); %Command for example data from github
TN5sort('EX_READ1_truncated.fastq','EX_READ2_truncated.fastq','EXAMPLEg2',Summaryoutputfile,'pargRNA2.dat'); %Command for example data from github
fprintf(fileID,strcat('Started: ',string(datetime('now')),'\n')); %Write date and time
%Copy and paste three lines above as needed...
%--------------------------------------------------------------------------

fprintf(fileID,strcat('Finsihed: ',string(datetime('now')),'\n')); %Write date and time
flcose(fileID);
end

function TN5sort(fwdread,revread,outfile1,outfile2,parfile)
disp(strcat('Starting: ', outfile1)); 

[T, T2] = readzippedfastq(fwdread,revread); %READ FASTQ - Check if zipped or not. Unzip if needed - T1=Read1 and T2=Reverse complement of Read2.
disp('Done reading fastq');

%Read in experiment specific parameters-------------------------------------------------------------------
fid = fopen(parfile);                   %Open parameter file
Intact = fgetl(fid);                    %Line  1 - Amplicon sequence if no editing takes place
Deletion = fgetl(fid);                  %Line  2 - Amplicon sequence if a targeted deletion occurs
Inversion = fgetl(fid);                 %Line  3 - Amplicon sequence if a targeted inversion occurs
DSBbp = str2num(fgetl(fid));            %Line  4 - DSB location in bp
WINDOWbp = str2num(fgetl(fid));         %Line  5 - BP window to check for indels on each side of DSB. Suggest 5.
MINSCORE = str2num(fgetl(fid));         %Line  6 - Minimum score for aligning to gene products   %%MOVE TO OPTIONAL AND CALCULATE
SHORTMINSCORE = str2num(fgetl(fid));    %Line  7 - Minimum score for aligning short portion      %%MOVE TO OPTIONAL AND CALCULATE[Roughly = 1.24*DSBbp]
AAVSCORE = str2num(fgetl(fid));     	%Line  8 - Minimum alignment score for the AAV genome    %%MOVE TO OPTIONAL AND CALCULATE
dysinv =  fgetl(fid);                   %Line  9 - Large portion of target gene with inversion for unexpected modifications (e.g. chewback)
dysrev = fgetl(fid);                    %Line 10 - Large portion of target gene for unexpected modifications (e.g. chewback)
minsize = str2num(fgetl(fid));          %Line 11 - minimum amplicon size. Reject shorter amplicons
AAVGENOME1 = fgetl(fid);                %Line 12 - AAV Genome 1
AAVGENOME2 = fgetl(fid);                %Line 13 - AAV genome 2
DYSGapOpenValue = str2num(fgetl(fid));  %Line 14 - Value for gap opening when aligning to the dystrophin gene (Default 28)
AAVGAPOPENVALUE = str2num(fgetl(fid));  %Line 15 - Value for gap opening when aligning to the AAV genome (default 8)
DYSEXTENDGAP = str2num(fgetl(fid));     %Line 16 - Value for gap extending gaps when aligning (Default 0 - no penalty)
fclose(fid);
%END FILE READ----------------------------------------------------------

%Define remaining variables-----------------------------------------------
llim = 1;  %Start bp of amplicon is always 1
winleft = DSBbp-WINDOWbp; %Window around DSB
winright = DSBbp+WINDOWbp; %Window around DSB
newStr = extractBetween(Intact,winleft+1,winright); %temp
intexact = newStr{1}; %Exact region around the DSB
newStr = extractBetween(Deletion,winleft+1,winright);
delexact = newStr{1}; 
newStr = extractBetween(Inversion,winleft+1,winright);
invexact = newStr{1};
newStr = extractBetween(Intact,1,DSBbp);
LeftAlign = newStr{1}; %Get DNA sequence from the start of the primer to the location of the DSB
AAVGENOME1R = seqrcomplement(AAVGENOME1); %Reverse complement of AAV genome 1
AAVGENOME2R = seqrcomplement(AAVGENOME2); %Reverse compleent of AAV genome 2
%--------------------------------------------------------------------------

%initialize Variables------------------------------------------------------
delcount = 0;       %Number of deletions reads
invcount = 0;       %Number of inversions reads
intcount = 0;       %Number of intact reads
aavcount = 0;       %Number of AAV integration events
mysterycount = 0;   %Reads that did not align
intactsubcount = 0; %Count the intact reads that have a substituion within the window
intactNHEJcount=0;  %count the intact reads that have an insertion/deletion/substitution
delindelcount = 0;  %Count the deletion reads that have an insertion or deletion
invindelcount = 0;  %Count the inverion reads that have an insertion or deletion
intactindelcount=0; %Count indels in intact reads

intactmatrix = cell(50000,1);    %output table for intact reads
deletionmatrix = cell(10000,1);  %output table for deletion reads
inversionmatrix = cell(10000,1); %Output table for inverted reads
aavintmatrix = cell(2000,1);    %Output table for aav integration reads
mysterymatrix = cell(10000,1);   %Output table for unalinged reads
indelmatrix = cell(10000,1);    %output table for indels
delindelmatrix = cell(10000,1);  %output table for deletion-indels
invindelmatrix = cell(10000,1);  %output table for inversion-indels
%END initialize variables

%Align each read with expected products
disp('Aligning');
for j = 1:height(T)
    temp = T{j,1};  %Forward read
    S = char(temp);
    temp2 = T2{j,1}; %Reverse read
    S2 = char(temp2);
       
    if(cellfun('length',temp) > minsize) %Skip short amplicons
        [s0, Alignment] = swalign(S(llim:DSBbp), LeftAlign,'Alphabet', 'NT'); 
        
        if(s0 > SHORTMINSCORE)  %Excludes reads that don't align to the target locus
            [sf, af] = nwalign(dysrev,S(llim:length(S)), 'GapOpen', DYSGapOpenValue, 'Glocal', true,'ExtendGap', DYSEXTENDGAP,'Alphabet', 'NT');
            [sr, ar] = nwalign(dysrev,S2(llim:length(S2)), 'GapOpen', DYSGapOpenValue, 'Glocal', true,'ExtendGap', DYSEXTENDGAP,'Alphabet', 'NT');
            [si, ai] = nwalign(dysinv,S2(llim:length(S2)), 'GapOpen', DYSGapOpenValue, 'Glocal', true,'ExtendGap', DYSEXTENDGAP,'Alphabet', 'NT');
            [tempstart, tempend] = findminmax(af);
            [tempstart2, tempend2] = findminmax(ar);
            [tempstart3, tempend3] = findminmax(ai);

            %align to reference and discard poor scores
            [s1, AlignInt] = swalign(S, Intact,'Alphabet', 'NT');
            [s2, AlignDel] = swalign(S, Deletion,'Alphabet', 'NT');
		    [s3, AlignInv] = swalign(S, Inversion,'Alphabet', 'NT');

           
           %CONDITION 1 - Read is an intact read-------------------------------------------------------------
           if (s1>s2)&&(s1>s3)&&(s1>MINSCORE)
                if(strfind(S,intexact)< 151)  %PERFECT INTACT
                     intcount=intcount+1;
                     s=sprintf('%s\t%s\t%s\t%d\t%d\t%d\t%d\t', S,S2,'intact',tempstart,tempend,tempstart2,tempend2);
                     intactmatrix(intcount) = cellstr(s);
                else  %INDELS AND SUBSTITUTION
                    intactNHEJcount = intactNHEJcount + 1; %COUNT ALL SUBSTITIONS AND INDELS AS NHEJ
                    DashFIND = strfind(AlignInt(1,:),'-');
                    DashFIND2 = strfind(AlignInt(3,:),'-');
                    DashFINDlength = length(DashFIND);
                    DashFINDlength2 = length(DashFIND2);
                    %SUBSITUTION
                    if(DashFINDlength==0)&&(DashFINDlength2==0)
                        intactsubcount = intactsubcount + 1;
                        s=sprintf('%s\t%s\t%s\t%d\t%d\t%d\t%d\t',S,S2,'SUBSTITUTION',tempstart,tempend,tempstart2,tempend2);
                        indelmatrix(intactNHEJcount) = cellstr(s);
                    %INSERTION OR DELETION (INDEL)
                    else %record indel in table and count
                        if ((min([DashFIND 151])<winright)&&(max([DashFIND 1])>winleft)) || ((min([DashFIND2 151])<winright)&&(max([DashFIND2 1])>winleft))
                            s=sprintf('%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t',S,S2,'INDEL',tempstart,tempend,tempstart2,tempend2, DashFINDlength, DashFINDlength2, min(DashFIND), max(DashFIND), min(DashFIND2), max(DashFIND2));
                            %s=strcat(S,':',S2,':Arealindel:',num2str(tempstart),':',num2str(tempend),':',num2str(tempstart2),':',num2str(tempend2),':',num2str(DashFINDlength),':',num2str(DashFINDlength2),':',num2str(min(DashFIND)),':',num2str(max(DashFIND)),':',num2str(min(DashFIND2)),':',num2str(max(DashFIND2)));
                            intactindelcount = intactindelcount + 1;
                        else %out of window
                           s=sprintf('%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t',S,S2,'INDEL_OUTSIDE_WINDOW',tempstart,tempend,tempstart2,tempend2, DashFINDlength, DashFINDlength2, min(DashFIND), max(DashFIND), min(DashFIND2), max(DashFIND2));
                           %s=strcat(S,':',S2,':outsideindel:',num2str(tempstart),':',num2str(tempend),':',num2str(tempstart2),':',num2str(tempend2),':',num2str(DashFINDlength),':',num2str(DashFINDlength2),':',num2str(min(DashFIND)),':',num2str(max(DashFIND)),':',num2str(min(DashFIND2)),':',num2str(max(DashFIND2)));
                        end    
                        indelmatrix(intactNHEJcount) = cellstr(s);
                    end
                end
            %END INTACT----------------------------------------------------
                    
            %CONDITION 2 - Read is a Deletion------------------------------
            elseif(s2>s3)&&(s2>MINSCORE)
                delcount=delcount+1;
                s=sprintf('%s\t%s\t%d\t%d\t%d\t%d\t', S,S2,tempstart,tempend,tempstart2,tempend2);
                deletionmatrix(delcount) = cellstr(s);
                disp('deletion');
                if(strfind(S,delexact)< 151) 
                      %do nothing
                else
                    %record indel in table and count
                    delindelcount = delindelcount + 1;
                    delindelmatrix(delindelcount) = cellstr(S);
                    
                end
            %END DELETION -------------------------------------------------
            
            %CONDITION 3 - Read is an Inversion----------------------------
            elseif(s3>MINSCORE)
                invcount=invcount+1;
                s=sprintf('%s\t%s\t%d\t%d\t%d\t%d\t', S,S2,tempstart,tempend,tempstart3,tempend3);
                inversionmatrix(invcount) = cellstr(s);
                disp('inversion');
                if(strfind(S,invexact)< 151) 
                      %do nothing
                else
                    invindelcount = invindelcount + 1;
                    invindelmatrix(invindelcount) = cellstr(S);
                end
            %END INVERSION-------------------------------------------------
            
            %CONDITION 4 - includes AAV integration, large deletion, etc 
            else
                %CHECK FOR VECTOR GENOME-----------------------------------
                [s4, a4] = nwalign(AAVGENOME1,S(DSBbp:length(S)), 'GapOpen', AAVGAPOPENVALUE,'Glocal', true,'Alphabet', 'NT');
                [s5, a5] = nwalign(AAVGENOME1R,S(DSBbp:length(S)), 'GapOpen', AAVGAPOPENVALUE,'Glocal', true,'Alphabet', 'NT');
                [s6, a6] = nwalign(AAVGENOME2,S(DSBbp:length(S)), 'GapOpen', AAVGAPOPENVALUE,'Glocal', true,'Alphabet', 'NT');
                [s7, a7] = nwalign(AAVGENOME2R,S(DSBbp:length(S)), 'GapOpen', AAVGAPOPENVALUE,'Glocal', true,'Alphabet', 'NT');
                
               %MYSTERY SEQUENCE ------------------------------------------
               %TO DO: Add code to check genome wide for translocations or
               %large deletions. Right now this is done manually
               if(max([s4,s5,s6,s7])<AAVSCORE)
                        mysterycount = mysterycount+1;
                        s=sprintf('%s\t%s\t%d\t%d\t%d\t%d\t',S,S2,tempstart,tempend,tempstart2,tempend2);
                        mysterymatrix(mysterycount) = cellstr(s); 
               %END MYSTERY SEQUENCE--------------------------------------
               
               %Found vector integration----------------------------------
               else
                   %AAV function (stringwhich is aav) (temploc),
                   [RETURNEDSTRING,RETURNEDSTART,RETURNEDEND]= SCOREAAVGENOME(s4,s5,s6,s7,a4,a5,a6,a7);
                   aavcount=aavcount+1;
                   s=sprintf('%s\t%s\t%s\t%d\t%d\t', S,S2,RETURNEDSTRING,RETURNEDSTART,RETURNEDEND);
                   aavintmatrix(aavcount) = cellstr(s);
                   disp(RETURNEDSTRING);
               end
               %END VECTOR GENOME-----------------------------------------
            end
        end
    end
end

%Output all numerical data and aligned reads into separate files
disp(table(intcount,intactindelcount,delcount,invcount,aavcount,mysterycount));
fileID = fopen(strcat(outfile2,'_quantiative.dat'),'a'); 
fprintf(fileID,'%s\t %6d\t %6d\t %6d\t %6d\t %6d\t %6d\t %6d\t %6d\t %6d\t %6d\n',outfile1,[intcount intactindelcount intactNHEJcount intactsubcount delcount delindelcount invcount invindelcount aavcount mysterycount]);
fclose(fileID);

Tabdata1 = cell2table(mysterymatrix);
writetable(Tabdata1,strcat(outfile1,'_mystery.dat'))
Tabdata2 = cell2table(aavintmatrix);
writetable(Tabdata2,strcat(outfile1,'_aav.dat'))
Tabdata3 = cell2table(deletionmatrix);
writetable(Tabdata3,strcat(outfile1,'_deletion.dat'))
Tabdata4 = cell2table(indelmatrix);
writetable(Tabdata4,strcat(outfile1,'_indel.dat'))
Tabdata5 = cell2table(inversionmatrix);
writetable(Tabdata5,strcat(outfile1,'_inversion.dat'))
Tabdata6 = cell2table(intactmatrix);
writetable(Tabdata6,strcat(outfile1,'_intact.dat'))
end

%Checks if zipped, unzips, and assigns fastq, strips out sequence data, and places into matlab table format
function [T, T2] = readzippedfastq(fwdread,revread)
    if(endsWith(fwdread,'.gz'))
        fwdread = char(gunzip(fwdread));
    end
    if(endsWith(revread, '.gz'))
        revread = char(gunzip(revread));
    end
    forwardreads = fastqread(fwdread);  %Read in forward reads
    reversereads = fastqread(revread);  %Read in reverse reads
    structsize = size(forwardreads);
    revstructsize = size(reversereads);
    for i = 1:structsize(1,2)
        forwardreadS(i).Sequence = forwardreads(i).Sequence;
    end
    for rccount = 1:revstructsize(1,2)
        reversereadS(rccount).Sequence = seqrcomplement(reversereads(rccount).Sequence);
    end
    T = struct2table(forwardreadS);
    T2 = struct2table(reversereadS);
end

%Finds beginning and end of alignment
function [seqstart, seqend] = findminmax(alignment)
    Arr = strfind(alignment(3,:),'A');
    Trr = strfind(alignment(3,:),'T');
    Crr = strfind(alignment(3,:),'C');
    Grr = strfind(alignment(3,:),'G');
    seqstart = min([Arr,Trr,Crr,Grr]);
    seqend = max([Arr,Trr,Crr,Grr]);
end

%Selects best aligning AAV genome. In cases where the score is equal, it randomly selects from the best scoring alignments.
function [AAVGENOME, AAVstartleft, AAVstartright] = SCOREAAVGENOME(s4,s5,s6,s7,a4,a5,a6,a7)
    if (s4>max([s5,s6,s7]))%S4 is the highest
        AAVselect = 4;
    elseif (s5>max([s4,s6,s7])) %S5 is the highest
        AAVselect = 5;
    elseif (s6>max([s4,s5,s7])) %S6 is the highest
        AAVselect = 6;
    elseif (s7>max([s4,s5,s6])) %S7 is the highest
        AAVselect = 7;
    elseif (s4==s5)&&(s5==s6)&&(s6==s7) %All equal
        AAVselect = randi([4 7],1,1);
    elseif (s4==s5)&&(s5==s6)%3 are equal
        AAVselect = randi([4 6],1,1);
    elseif (s5==s6)&&(s6==s7)
        AAVselect = randi([5 7],1,1);
    elseif (s4==s6)&&(s6==s7)
        temprand = randi([1 3],1,1);
        if(temprand == 1)
            AAVselect = 4;
        elseif(temprand == 2)
            AAVselect = 6;
        else
            AAVselect = 7;
        end
    elseif (s4==s5)&&(s5==s7)
        temprand = randi([1 3],1,1);
        if(temprand == 1)
            AAVselect = 4;
        elseif(temprand == 2)
            AAVselect = 5;
        else
            AAVselect = 7;
        end
    elseif (s4==s5)&&(s4>max([s6,s7]))
        AAVselect = randi([4,5],1,1);
    elseif (s4==s6)&&(s4>max([s5,s7]))
        temprand = randi([1 2],1,1);
        if(temprand == 1)
            AAVselect = 4;
        else
            AAVselect = 6;
        end
    elseif (s4==s7)&&(s4>max([s5,s6]))
        temprand = randi([1 2],1,1);
        if(temprand == 1)
            AAVselect = 4;
        else
            AAVselect = 7;
        end
    elseif (s5==s6)&&(s5>max([s4,s7]))
        AAVselect = randi([5,6],1,1);
    elseif (s5==s7)&&(s5>max([s4,s6]))
        temprand = randi([1 2],1,1);
        if(temprand == 1)
            AAVselect = 5;
        else
            AAVselect = 7;
        end
    elseif (s6==s7)&&(s6>max([s4,s5]))
        AAVselect = randi([6,7],1,1);
    else
        disp('f');
    end
    
    %ASSIGN VARIABLES
    if(AAVselect == 4)
        AAVGENOME = 'AAVGENOME1';
        [AAVstartleft, AAVstartright] = findminmax(a4);
    elseif(AAVselect == 5)
        AAVGENOME = 'AAVGENOME1R';
        [AAVstartleft, AAVstartright] = findminmax(a5);
    elseif(AAVselect == 6)
        AAVGENOME = 'AAVGENOME2';
        [AAVstartleft, AAVstartright] = findminmax(a6);
    elseif(AAVselect == 7)
        AAVGENOME = 'AAVGENOME2R';
        [AAVstartleft, AAVstartright] = findminmax(a7);
    else
        AAVstartleft = 0;
        AAVstartright =0;
    end
end