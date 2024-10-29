
%% generate fake files suitable for annotation by HOMER annotatePeaks
makePerms=false; %only needed on first run

if makePerms
homeTable=readtable('/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/HomertDP43Binding.csv'); % read in your data (output from Homer)
outputDir='/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/HomerPermutationTest/'; % generate a directory where you want your permutations to go (code will produce this if you haven't)
permNumber=2000; %number of Permutations to run
homerPermGenerator(homeTable,outputDir,permNumber)
end

% outputting our actual ratio values
refTable=readtable('/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/TDP43_BindingSites_Annotated.xlsx');
annotation=string(refTable.Annotation(:)); chromLoc=string(refTable.Chr(:)); detAnnotation=string(refTable.DetailedAnnotation(:));
[intergenicRat,nonCodeRat,promoteRat, lineRat, aluRat, intronLineRat, intronALURat,alphaRat,repRat, intronRepRat,intronalphaRat]...
    = annoRatioMaker(annotation,detAnnotation);

%%
alphaList=[]; repList=[];
for tt=1:height(refTable)
     chrNum= str2num(refTable.Chr{tt}( (find(refTable.Chr{tt}=='r')+1):end));

    if isempty(chrNum) && contains(refTable.Chr{tt},'X' )
        chrNum=23;
    elseif isempty(chrNum) && contains(refTable.Chr{tt},'Y' )
            chrNum=24;
    end

    if contains(refTable.DetailedAnnotation{tt}, 'ALR/Alpha|Satellite|centr')
alphaList=[alphaList,chrNum ];
    elseif contains(refTable.DetailedAnnotation{tt}, '(')
repList=[repList, chrNum];
    end

end

figure
subplot(1,2,1)
histogram(alphaList, 'BinWidth',1)
title('Alpha Satellite Centers')
xlabel('Chromosome Number')
ylabel('Peak Count')
subplot(1,2,2)
histogram(repList, 'BinWidth',1, 'FaceColor','m')
title('Simple Repetitive Elements')
xlabel('Chromosome Number')
ylabel('Peak Count')
sgtitle('Chromsomal Distribution of Repetive elements')


%% Hits are a bit strange for just the HEK293 output... I just screen for enhancers on the enhancer atlas site now
% %% checkin for enhancers just for fun
% enhancerLocations= readtable("/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/Consensus_enhancers_HEK293.txt");
% % this comes from http://www.enhanceratlas.org/data/download/enhancer/hs/HEK293.bed
% annoCatch=[];
% 
% for tt=1: height(refTable)
%     chr=refTable.Chr{tt}; startr= refTable.Start(tt); endr= refTable.End(tt); 
%     %sigh definitely could implement this more elegantly
%     for dd=1:height(enhancerLocations)
%         if strcmp(chr, enhancerLocations.Var1{dd}) && startr>= enhancerLocations.Var2(dd) && endr<= enhancerLocations.Var3(dd)
%         annoCatch=[annoCatch,dd];
%         fprintf(['enhancer at ' chr,])
%         end
%     end
% end

%% generating permutation comparators

 permList=dir('/Users/pennerc/HomerPermOutput'); permList={permList.name}; permList(1:5)=[];

ammendr= cell(1,6); %not sure exactly why this happens but the last six are replaced by various meta data... irritating 
startPt=1994;
for dd=1:6
ammendr{dd}= ['tdp43bindingoutput_TDP43BindingPerm_', num2str(startPt+dd),'.txt'];
end
permList=[permList,ammendr];

intergenicRatPerms=nan(1, length(permList));
nonCodeRatPerms=nan(1, length(permList));
promoteRatPerms= nan(1, length(permList));
lineRatPerms= nan(1, length(permList));
aluRatPerms= nan(1, length(permList));
intronLineRatPerms= nan(1, length(permList));
intronALURatPerms=nan(1, length(permList));
alphaRatPerms=nan(1, length(permList));
repRatPerms=nan(1, length(permList));
intronRepRatPerms=nan(1, length(permList));
intronalphaRatPerms=nan(1, length(permList));
for zz=1:length(permList)
tableAtPlay=readtable(['/Users/pennerc/HomerPermOutput/',permList{zz}]);
if ~isempty(tableAtPlay)
annotation=string(tableAtPlay.Annotation(:)); chromLoc=string(tableAtPlay.Chr(:)); detAnnotation=string(tableAtPlay.DetailedAnnotation(:));
[intergenicRatPerms(zz),nonCodeRatPerms(zz),promoteRatPerms(zz), lineRatPerms(zz), aluRatPerms(zz), intronLineRatPerms(zz), intronALURatPerms(zz),alphaRatPerms(zz), repRatPerms(zz),intronRepRatPerms(zz), intronalphaRatPerms(zz)] = ...
    annoRatioMaker(annotation,detAnnotation);
end
end

%% plotting resultant histograms

homerPermPlotter(intergenicRat,intergenicRatPerms, 'Intergenic')

homerPermPlotter(nonCodeRat,nonCodeRatPerms, 'Intronic/Non-Coding')

homerPermPlotter(promoteRat,promoteRatPerms, 'Promoter')

homerPermPlotter(lineRat,lineRatPerms, 'LINE regions')

homerPermPlotter(aluRat,aluRatPerms, 'ALU regions')

homerPermPlotter(intronALURat,intronALURatPerms, 'alu regions within Introns')

homerPermPlotter(intronLineRat,intronLineRatPerms, 'Line regions within Intronic')


homerPermPlotter(alphaRat,alphaRatPerms, 'Alpha Satelites')
*homerPermPlotter(intronRepRat,intronRepRatPerms, 'repettive Sequences  within Introns')
homerPermPlotter(intronalphaRat,intronalphaRatPerms, 'alpha satelites within Introns')



