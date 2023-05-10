function CheckingPairingMyron(data, stats, graph, cycle, title1, title2)
%HOW TO USE IT
%Function takes a cell array type input with Prefixes (movies) in each cell; for
%example CheckingPairingMyron({[Prefix1], [Prefix2]}) where each of the
%prefixes finishes the path to the folder with the ParticlesPlus data. Can also separate
%by condition for example CheckingStickiness({[Prefix1];[Prefix2]}) where conditions
%are separated by the semi-colon and an infinite number of movies can be grouped in each condition

%IMPORTANT: Conditions have to match in length so if they don't put in [NaN] values as
%placeholders until they do (the function will skip through these. Really
%only supports two conditions - not really multi-dimensional. 

%WHAT IT DOES
%This function will end up outputing the mean distance between
%homologus loci for the rest of the movie after the loci are considered pairing, which is currently
%defined as when in any given 4 minutes or 7 33s frames in a movie, the 3D distance 
%between homologous chromosome loci remains less than 0.65um. It will also
%output the percent of nuclei that have pairing chromosomes which is the ratio of pairing traces
%over total number of nuclei scored. This is for every movie for each condition. 
%Conditions will be separated and compared. Each nucleus that has pairing chromosomes 
%in any given movie will be recorded onlu once. 

%EXTRA VARIABLES
%Furthermore by typing 'StatMe' as an input variable for stats when there
%are two conditions present the function will perform a fisher's exact test between the
%percent of pairing traces for the two conditions and output the p value.
%IMPORTANT: Fisher's test is one-tailed so gypsy has to be the first
%condition and no gypsy must be the second condition - it will tell you if
%the first condition is significantly bigger than the second condition

%Finally, by typing 'GraphMe' as an input variable for graph a histogram of 
%the percent pairing traces +/- error. This will appear for the two conditions with a p value and 
%each condition can be named whatever by inputing a string for the title1 and
%title2 variable ex: 'Gypsy' and 'No Gypsy'.

%Cycle defines whether you use pairing parameters for cycle 13 (i.e
%'Cycle13') or cycle 14 (i.e. 'Cycle 14)

MasterArray = [];
ForStats = [];
%Creates two arrays, a MasterArray to store data and one to perform stats

ogdirectory = 'P:/BatemanUG1/LivemRNA'
%Stores the all time current working directory before all this code as a variable
%IMPORTANT: Have to define cycle

if exist('cycle', 'var') && strcmp(cycle, 'Cycle13')
    pairing_param = 5;
elseif exist('cycle', 'var') && strcmp(cycle, 'Cycle14')
    pairing_param = 6;
else
    cycle = 'null';
    pairing_param = 6;
    disp('Cycle not specified pairing paramter used is > 6 frames under 0.65um')
end
%Defines cycle parameters to use for pairing definition

for r = 1:length(data(:,1))
    nuclei_count = 0;
    %Loops through each condition specified in the input data
    Sticky = [];
    %Makes an empty array
    disp(strcat('Data for condition number', num2str(r)))
    %Tells you what condition the data (Sticky array) is that you're
    %looking at
    framecount = 0;
    MinorStats = [];
    for i = 1:length(data(1,:))
        MinorSticky = [];
        %Loops through the cell array
        if isa(data{r, i}, 'char')
            [~,~,DropboxFolder,~,~]=...
            DetermineLocalFolders(data{r,i});
            DataFolder=[DropboxFolder,filesep,data{r,i}];
            load([DataFolder,filesep,'ParticlesPlus.mat'])
            MovieIndex = i;
        else
            continue;
        end
        a = 1;
        amax = length(ParticlesPlus{1});
        %Sets parameters for the while loop below where amax is equal to
        %the total number of nuclei in the ParticlesPlus Data of Interest
        while a < amax + 1
        %Loops through each nucleus
             if getfield(ParticlesPlus{1}, {a}, 'Approved') == 1
             %If the nucleus was approved in the tracking then this function 
             %includes the nucleus in the analysis
             nuclei_count = nuclei_count + 1;
             GaussDists = getfield(ParticlesPlus{1}, {a}, 'Gauss3DDistWithZ');
             %Grabs the 3DGaussian Distance between chromosomes with Z and
             %stores it as a single row vector variable GaussDists
             framecount = framecount + length(GaussDists(~isnan(GaussDists)));
                 for b = 4:(length(GaussDists)-4)
                 %loops through the 3DGaussian Distance data within an
                 %index where we can look 3 frames before and 3 after for a
                 %total of 7 frames or ~4 minutes
                 Pairing_Check = GaussDists(b-3: b+3);
                 count_paired = 0;
                 for distance = 1:length(Pairing_Check)
                     if Pairing_Check(distance) < 0.65
                         count_paired = count_paired + 1;
                     else
                         continue
                     end
                 end
                     if count_paired > pairing_param
                     %If the distance between loci drops below 0.65um 
                     %for ~4 minutes we're calling this a pairing trace
                         Stickyelem = [];
                         %Creates an empty array to store data of which
                         %movie, in what nucleus, and what frame the
                         %pairing trace occured and the average distance between loci 
                         %after pairing occured (and how many frames before
                         %or after pairing was established are missing)
                         Stickyelem = [Stickyelem MovieIndex];
                         %Adds the index corresponding to the movie at that index in the data (cell array) that is inputed into  
                         %the function to denote in which movie the contact occured 
                         Stickyelem = [Stickyelem getfield(ParticlesPlus{1}, {a}, 'Nucleus')];
                         %Adds the nucleus number where pairing was determined to occur to the array
                         Stickyelem = [Stickyelem b];
                         %Then adds the frame number where the pairing was established, in that nucleus, to the same row

                         GaussDistsTrunc = GaussDists(b+1:length(GaussDists));
                         GaussDistswoNaN = (GaussDistsTrunc(~isnan(GaussDistsTrunc)));
                         %Creates an array of the distances after pairing
                         %is determined and removes any NaN values 
                         GaussDistsb4Contact = GaussDists(1:b);
                         GaussDistsb4ContactwoNaN = (GaussDistsb4Contact(~isnan(GaussDistsb4Contact)));
                         %Creates an array of the distances before pairing
                         %and removes any NaN values

                         frame_num = length(GaussDistsTrunc);
                         frame_num_woNaN = length(GaussDistswoNaN);
                         %Creates a variable of how many frames there are after pairing and how
                         %many have data points

                         frame_num_before = length(GaussDistsb4Contact);
                         frame_num_beforewoNaN = length(GaussDistsb4ContactwoNaN);
                         %Creates a variable of how many frames there are before pairing and how
                         %many have data points
                         
                         Stickyelem = [Stickyelem mean(GaussDistswoNaN)];
                         Stickyelem = [Stickyelem frame_num];
                         Stickyelem = [Stickyelem frame_num_woNaN];
                         Stickyelem = [Stickyelem frame_num_before];
                         Stickyelem = [Stickyelem frame_num_beforewoNaN];
                         %Adds all these values to the pairing trace array
                         
                         MinorSticky = [MinorSticky; Stickyelem];
                         %Adds this individual contact data to a master array
                     else
                         continue;
                     end
                 end
             end
             a = a+1;
             % Moves onto the next nucleus and repeats the analysis
        end
        if length(MinorSticky) == 0
            continue
        else
            [~, uidx] = unique(MinorSticky(:, 2), 'stable');
            MinorSticky = MinorSticky(uidx, :); 
            Sticky = [Sticky; MinorSticky];
            %removes all repeat pairing traces for nuclei that may satisfy
            %the condition more than once and adds the set of pairing traces for that movie
            %to a master array for the whole condition you're looking at 
        end
    end

    %At this point the program has looped through all the ParticlePlus datasets
    %of each Prefix you've specificied and put each set of pairing traces for each movie and the
    %accompanying data into this master Sticky Array

    cd(ogdirectory)
    %Returns us to the original working directory, which for us is
    %'P:\BatemanUG1\LivemRNA'
    
    if length(Sticky) == 0
        disp('There were no pairing traces in this set of movies')
        Sticky = zeros([1, 8]);
        %Assures Sticky has data but just a row of 0s if there was no
        %pairing traces to be used later
        PostContactData = [];
        %Creates an empty array
        for row = 1:length(Sticky(:, 4))
        %Loops through each pairing occurence
            PostContactData = [PostContactData; Sticky(row,4)];
            %Takes the mean of the distance between loci post pairing being established 
            %and adds the means for each pairing occurence to an array
        end
    else
        PostContactData = [];
        %Creates an empty array
        for row = 1:length(Sticky(:, 4))
        %Loops through each pairing trace
            PostContactData = [PostContactData; Sticky(row,4)];
            %Takes the mean of the frames post contact and adds the
            %means for each contact to an array
        end
    end
    meanPC = mean(PostContactData);
    %Takes the mean distance after pairing for all pairing traces
    %from the data for each condition listed in the initial function
    stderror = std(PostContactData) / sqrt(length(PostContactData));
    %Takes the SEM (Standard error of the mean) 
    
    StickyTable = array2table(Sticky);
    StickyTable.Properties.VariableNames = {'MovieIndex', 'Nucleus#', 'Frame#', 'Dist After Pairing', '# Frames After Pairing', '# After woNaN', '# Frames Cefore', '# Cefore wo NaN'};
    %Creates a Table in which the noted variables are displayed for an
    %entire condition
    
    disp(StickyTable)
    disp('The total number of nuclei for this condition is')
    disp(nuclei_count)
    disp('The number of frames analyzed is')
    disp(framecount)
    if Sticky(1,1) == 0
        disp('The number of pairing traces is')
        disp(0)
        disp('The proportion of nuclei having pairing is')
        disp(0)
    else
        disp('The number of pairing traces is')
        disp(length(PostContactData))
        disp('The proportion of nuclei having pairing is')
        disp(length(PostContactData)/nuclei_count)
    end
    disp('The mean distance in loci after pairing is established is ')
    disp(meanPC)
    disp(stderror)
    disp('Units are um')
    %Displays the initial array with movie# as inputted in fxn, nucleus, frame of pairing, and distances between
    %chromosomes after pairing and finally the mean and stderror of this mean distance dataset
    
    MinorArray = [meanPC, stderror];
    MasterArray = [MasterArray; MinorArray];
    %Append an array of the mean and stderror for the condition you're
    %looking at onto a master array holding data for all conditions
    if Sticky(1,1) == 0
        p_master = 0;
        MinorStats = [p_master nuclei_count p_master sqrt((p_master * (1-p_master))/nuclei_count)];
        ForStats = [ForStats; MinorStats];
    else
        p_master = length(Sticky(:, 2))/nuclei_count;
        MinorStats = [length(Sticky(:,2)) nuclei_count-length(Sticky(:,2)) p_master sqrt((p_master * (1-p_master))/nuclei_count)];
        ForStats = [ForStats; MinorStats];
    end
    %Append the necessary information to perform the fisher's statistics
    %test to an array MinorStats which holds data for a specific condition
    %that then adds to a master array ForStats which holds data for all
    %conditions
    %These are the number of nuclei that have paired homologs, the
    %total number of nuclei, the proportion of nuclei that have pairing to
    %all nuclei and the standard deviation for each proportion 
        
end
%This end causes the function to loop back and goes to the next condition and does everything
%over again

PairingArray = ForStats(:, [3,4]) * 100;
FisherTable = ForStats(:, [1,2]);
%Creates a PairingArray with the proportion and std to be used in graphing
%and a printed table as percentages
%Creates a FisherTable to be used for Fisher's statistics test 

if exist('title1') && exist('title2')
    MasterTable = array2table(MasterArray);
    MasterTable.Properties.RowNames = {title1, title2};
    MasterTable.Properties.VariableNames = {'Distance', 'Error'};
    disp(' ')
    disp('Data for both conditions:')
    disp(MasterTable)
    PairingTable = array2table(PairingArray);
    PairingTable.Properties.RowNames = {title1, title2};
    PairingTable.Properties.VariableNames = {'Percent Paired', 'Error'};
    disp(' ')
    disp(PairingTable)
else
    title1 = 'null';
    title2 = 'null';
    MasterTable = array2table(MasterArray);
    MasterTable.Properties.VariableNames = {'Distance', 'Error'};
    disp(' ')
    disp('Data for both conditions:')
    disp(MasterTable)
    PairingTable = array2table(PairingArray);
    PairingTable.Properties.VariableNames = {'Percent Paired', 'Error'};
    disp(' ')
    disp(PairingTable)
end
%Displays mean distance post pairing and stderror as well as percentage of nuclei having pairing
%and standard deviation for two conditions and will title them if titles are given

if exist('stats', 'var') && strcmp(stats, 'StatMe')
    [~, p_val, ~] = fishertest(FisherTable, 'Tail', 'right', 'Alpha', 0.05);
    fprintf('Fishers exact test')
    disp(' ')
    fprintf('p-value = %f\n', p_val)
    if p_val >0.05
        disp('Not significant')
    else
        disp('Significant')
    end
else
    stats = 'null';
end
%Performs a fisher's test on the pairing data between two conditions if there's 
%two conditions present and you ask for it through 'StatMe'

if exist('graph', 'var') && strcmp(graph, 'GraphMe')
    barvalues = PairingArray(:,1);
    barerrors = PairingArray(:,2);
    bar(barvalues)
    hold on
    errorbar(barvalues, barerrors)
    set(gca, 'XTickLabel', {title1, title2})
    ylabel('Percent Paired')
    if exist('p_val', 'var')
        legend(['p_val = ' num2str(p_val)])
    end
else
    graph = 'null';
end
%Graphs the percent of nuclei having pairing and standard deviation for each condition as a histogram. If the p
%value was created via StatMe (only works if you have two conditions) it
%will display this too. Will also set the titles for two conditions if you
%named them.
    
end