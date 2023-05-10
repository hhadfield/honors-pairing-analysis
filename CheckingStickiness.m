function CheckingStickiness(data, stats, graph, title1, title2)
%HOW TO USE IT
%Function takes a cell array type input with Prefixes (movies) in each cell; for
%example CheckingStickiness({[Prefix1], [Prefix2]}) where each of the
%prefixes finishes the path to the folder with the ParticlesPlus data. Can also separate
%by condition for example CheckingStickiness({[Prefix1];[Prefix2]}) where conditions
%are separated by the semi-colon and an inifinite number of movies can be grouped in each condition
%Conditions have to match in length so if they don't put in [NaN] values as
%placeholders until they do (the function will skip through these.

%WHAT IT DOES
%This function will end up outputing the percent of nuclei which have
%homologs that have experienced a pairing chance which is defined as at
%least one instance of interhomolog distance less than 0.65um. 

%EXTRA VARIABLES
%Furthermore by typing 'StatMe' as an input variable for stats when there
%are two conditions present the function will perform a one-sided Fisher's exact test between the
%mean distance for the two conditions and output the p value. Finally, by typing 'GraphMe' as an
%input variable for graph a histogram of the rates of pairing chances will appear
%for the two conditions with p value and each condition can be named
%whatever by inputing a string for title1 and title2 variable

MasterArray = [];
ForStats = [];
%Creates a MasterArray to store data and a ForStats array to perform stats

ogdirectory = 'P:/BatemanUG1/LivemRNA'
%Stores the all time current working directory before all this code as a variable
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
            %Finds the folder in which the movie/prefix you're looking at
            %is
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
                 for b = 2:(length(GaussDists)-3)
                 %loops through the 3DGaussian Distance data within an
                 %index where we can look one frame before and two after
                     if (GaussDists(b) <= 0.65) 
                     %If the distance between loci drops below or equal to
                     %0.65 microns adds the loci/nucleus to a data array 
                         Stickyelem = [];
                         %Creates an empty array to store data of which movie, in what nucleus, and what frame the contact
                         %occured and the distances between loci (as well as erases the previous contact data)
                         Stickyelem = [Stickyelem MovieIndex];
                         %Adds the index corresponding to the movie at that index in the data (cell array) that is inputed into  
                         %the CheckingStickiness function to denote in which movie the contact occured 
                         Stickyelem = [Stickyelem getfield(ParticlesPlus{1}, {a}, 'Nucleus')];
                         %Adds the nucleus number where the contact occured to the array
                         Stickyelem = [Stickyelem b];
                         %Then adds the frame number where the contact occured, in that nucleus, to the same row
                         Stickyelem = [Stickyelem GaussDists(b-1)];
                         %Then adds the distance 30s pre-contact to the array
                         Stickyelem = [Stickyelem GaussDists(b)];
                         %Then adds the contact distance (must be less than
                         %0.6um)
                         Stickyelem = [Stickyelem GaussDists(b+1)];
                         %Adds the distance 30s after the initial contact 
                         Stickyelem = [Stickyelem GaussDists(b+2)];
                         %Adds the distance 1 minute after initial contact
                         MinorSticky = [MinorSticky; Stickyelem];
                     else
                         continue;
                     end
                 end
             end
             a = a+1;
             % Moves onto the next nucleus and repeats the analysis
        end
         if length(MinorSticky)==0
             continue
         else
             [~, uidx] = unique(MinorSticky(:, 2), 'stable');
             MinorSticky = MinorSticky(uidx, :);
             Sticky = [Sticky; MinorSticky];
             %Erases any repeat nuclei which had loci scored more than once as
             %having a chance to pair because they had multiple instances
             %of interhomolog distance less than 0.65um
         end
    end

    %At this point the program has looped through all the ParticlePlus datasets
    %of each Prefix you've specificied and put each individual contact and the
    %accompanying data into this master Sticky Array

    cd(ogdirectory)
    %Returns us to the original working directory, which for us is
    %'P:\BatemanUG1\LivemRNA'

    if length(Sticky) == 1
        disp('There were no contacts in this set of movies')
    else
        PostContactData = [];
        %Creates an empty array
        for row = 1:length(Sticky(:, 4))
        %Loops through each contact occurence
            PostContactArray = [Sticky(row, 6);Sticky(row, 7)];
            MeanPC = mean(PostContactArray);
            PostContactData = [PostContactData; MeanPC];
            %Takes the mean of the two frames post contact and adds the
            %means for each contact to an array
        end
    end
    PostContactData = PostContactData(~isnan(PostContactData));
    %Removes all NaN values (did it at this point and not earlier to not
    %corrupt the data)
    meanPC = mean(PostContactData);
    %Takes the mean of the change in distance after contact for all contacts
    %from all the data listed in the initial function
    stderror = std(PostContactData) / sqrt(length(PostContactData));
    %Takes the SEM (Standard error of the mean) 
    
    StickyTable = array2table(Sticky);
    StickyTable.Properties.VariableNames = {'MovieIndex', 'Nucleus#', 'Frame#', 'Dist 30s Pre', 'Dist @ Contact', 'Dist 30s After', 'Dist 1min After'};
    %Creates a Table in which the noted variables are displayed for an
    %entire condition
    
    disp('The total number of nuclei for this condition is')
    disp(nuclei_count)
    disp(StickyTable)
    disp('The number of recordable contacts between loci is')
    disp(length(Sticky(:, 1)))
    disp('The number of recordable contacts analyzed is')
    disp(length(PostContactData))
    disp('The number of frames analyzed is')
    disp(framecount)
    disp('The number of contacts per nuclei is')
    disp(length(PostContactData)/nuclei_count)
    disp('The mean distance in loci for the minute post contact is ')
    disp(meanPC)
    disp(stderror)
    disp('Units are um')
    disp(framecount/nuclei_count * (33/60))
    %Displays the initial array with movie# as inputted in fnctn, nucleus, frame, and distances between
    %chromosomes; the mean distances between chromosomes over the minute post contact
    %and finally the mean and stderror of this mean distance dataset
    
    MinorArray = [meanPC, stderror];
    MasterArray = [MasterArray; MinorArray];
    %Append an array of the mean and stderror for the condition you're
    %looking at onto a master array
    if Sticky(1,1) == 0
        p_master = 0;
        MinorStats = [p_master nuclei_count p_master sqrt((p_master * (1-p_master))/nuclei_count)];
        ForStats = [ForStats; MinorStats];
    else
        p_master = length(Sticky(:, 2))/nuclei_count;
        MinorStats = [length(Sticky(:,2)) nuclei_count-length(Sticky(:,2)) p_master sqrt((p_master * (1-p_master))/nuclei_count)];
        ForStats = [ForStats; MinorStats];
    end
    %Makes proportions for fisher's test with total nuclei with a pairing chance, nuclei without, the percent of nuclei that had
    %loci with a pairing chance of total nuclei, and finally stderror of a proportion
end
%This end loops back and goes to the next condition and does everything
%over again

if exist('title1') && exist('title2')
    MasterTable = array2table(MasterArray);
    MasterTable.Properties.RowNames = {title1, title2};
    MasterTable.Properties.VariableNames = {'Distance', 'Error'};
    disp(' ')
    disp('Data for both conditions:')
    disp(MasterTable)
else
    title1 = 'null';
    title2 = 'null';
    MasterTable = array2table(MasterArray);
    MasterTable.Properties.VariableNames = {'Distance', 'Error'};
    disp(' ')
    disp('Data for both conditions:')
    disp(MasterTable)
end
%Displays mean distance minute post contact and stderror for two conditions
%and will title them if titles are given

PairingArray = ForStats(:, [3,4]) * 100;
FisherTable = ForStats(:, [1,2]);
%Creates a PairingArray with the proportion and std to be used in graphing
%and a printed table as percentages
%Creates a FisherTable to be used for Fisher's statistics test 

if exist('title1') && exist('title2')
    disp('Data for both conditions:')
    PairingTable = array2table(PairingArray);
    PairingTable.Properties.RowNames = {title1, title2};
    PairingTable.Properties.VariableNames = {'Percent Nuclei with Pairing Chance', 'Error'};
    disp(' ')
    disp(PairingTable)
else
    title1 = 'null';
    title2 = 'null';
    disp('Data for both conditions:')
    PairingTable = array2table(PairingArray);
    PairingTable.Properties.VariableNames = {'Percent Nuclei with Pairing Chance', 'Error'};
    disp(' ')
    disp(PairingTable)
end
%Displays a table with percent of nuclei with loci with a pairing chance
%and the associated error 

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
%Performs a fisher's test if there's two conditions present and you ask for it
%through 'StatMe'

if exist('graph', 'var') && strcmp(graph, 'GraphMe')
    barvalues = PairingArray(:,1);
    barerrors = PairingArray(:,2);
    bar(barvalues)
    hold on
    errorbar(barvalues, barerrors)
    set(gca, 'XTickLabel', {title1, title2})
    ylabel('Percent of Nuclei with Pairing Chance')
    if exist('p_val', 'var')
        legend(['p_val = ' num2str(p_val)])
    end
else
    graph = 'null';
end
%Graphs the percent of nuclei with a pairing chance and stderror for each condition as a histogram. If the p
%value was created via StatMe (only works if you have two conditions) it
%will display this too. Will also set the titles for two conditions if you
%named them

end