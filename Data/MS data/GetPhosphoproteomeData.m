function [ data ] = GetPhosphoproteomeData()
load('Humphrey_2013.mat', 'matchedData')

usefulColumns=[1 3 4 5 8 10 11 50:52 33 20:22];

LY=matchedData(:,25);
MK=matchedData(:,30);

LY.LYDirection = sign(matchedData.LYExpsInsulinStarvedMedianlog2);
MK.MKDirection = sign(matchedData.MKExpsInsulinStarvedMedianlog2);

LY{:,1}=2.^LY{:,1};
MK{:,1}=2.^MK{:,1};

matchedData=[matchedData(:,usefulColumns) LY MK];
matchedData.LYMean=nanmean(LY{:,1},2);
matchedData.LYSEM=nanSEM(LY{:,1},2);
matchedData.MKMean=nanmean(MK{:,1},2);
matchedData.MKSEM=nanSEM(MK{:,1},2);

matchedData.meanValues=[ones(size(matchedData,1),1) matchedData.meanValues];
matchedData.SEMValues=[zeros(size(matchedData,1),1) matchedData.SEMValues];

matchedData = sortrows(matchedData,'IPI368PositioninProtein','ascend');
matchedData = sortrows(matchedData,'Genenamesprimary','ascend');

matchedData.SEMValues(matchedData.SEMValues<1e-10)=inf;

data=matchedData;
data(~any(~isnan(data.meanValues(:,2:end)),2),:)=[]; %Remove data with only nans

end