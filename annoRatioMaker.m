function [intergenicRat,nonCodeRat,promoteRat, lineRat, aluRat, intronLineRat, intronALURat,alphaRat,repRat, intronRepRat,intronalphaRat] = annoRatioMaker(annotation,detAnnotation)
%% cleaning up annotations
warning('off','all')

try

for dd=1:length(annotation   )
    annotationAtPlay=convertStringsToChars(annotation(dd)); %convert to char 
    detAtPlay=convertStringsToChars(detAnnotation(dd));
    if ~strcmp(annotationAtPlay, 'Intergenic')
        annotation(dd)= convertCharsToStrings(annotationAtPlay(1: (find(isspace(annotationAtPlay),1)-1)   )) ; %remove space
    end
    if ~strcmp(detAtPlay, 'Intergenic')
        if detAtPlay(1:2)=='L1' | detAtPlay(1:2)=='L2' | detAtPlay(1:2)=='L3'
            detAnnotation(dd)="LINE"; %change to LINE 
        elseif (detAtPlay(end-2:end)=='Alu' )
            detAnnotation(dd)="Alu";
        elseif (detAtPlay(1))=='('
            detAnnotation(dd)="repeat";
        elseif sum(detAtPlay=='|') >=1
            detAnnotation(dd)= convertCharsToStrings(detAtPlay(1: (find(detAtPlay=='|',1)-1)  )); 
        elseif sum(detAtPlay==' ')
        detAnnotation(dd)=convertCharsToStrings(detAtPlay(1: (find(isspace(detAtPlay),1)-1)   ));
        end
    end
end


catch
b=1;
end
%% outputting metrics 
intergenicRat=sum(annotation=="Intergenic")/ length(annotation);
nonCodeRat=sum(annotation=="intron" | annotation=="non-coding")/length(annotation);
promoteRat= sum(annotation=="promoter-TSS")/length(annotation);
lineRat= sum(detAnnotation=='LINE')/ length(detAnnotation);
aluRat=sum(detAnnotation=='Alu')/length(detAnnotation);
intronLineRat= sum(annotation=="intron" & detAnnotation=="LINE")/ sum(annotation=='intron');
intronALURat=sum(annotation=="intron" & detAnnotation=="Alu")/ sum(annotation=='intron');
alphaRat=sum(detAnnotation=='ALR/Alpha')/length(detAnnotation);
repRat=sum(detAnnotation=='repeat')/length(detAnnotation);
intronalphaRat= sum(annotation=="intron" & detAnnotation=="'ALR/Alpha'")/ sum(annotation=='intron');
intronRepRat=sum(annotation=="intron" & detAnnotation=='repeat')/ sum(annotation=='intron');
    

end
