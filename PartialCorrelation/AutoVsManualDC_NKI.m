function [DiceCoef1, DiceCoef2, AvgDice] = AutoVsManualDC_NKI(subject, ext)

subject = char(subject);

auto_results = importdata(['/scr/murg2/MachineLearning/partialcorr/20comps_results/NKI_ICA_indiv/' subject ext '.1D']);
auto44 = auto_results;
auto44(auto44==1)=0;
auto44(auto44==2)=1;
auto45 = auto_results;
auto45(auto45==2)=0;

manual_results = importdata(['/scr/murg2/NKI/' subject '/manual_results_' subject '.1D']);
manual44 = manual_results;
manual44(manual44==1)=0;
manual44(manual44==2)=1;
manual45 = manual_results;
manual45(manual45==2)=0;

label = manual44+manual45;

manual45(manual45>0)=2;
auto45(auto45>0)=3;
Overlap45 = auto45-manual45;
[r] = find(Overlap45==1);
OverlapSize45=size(r);
[r1] = find(manual45==2);
labelsize_manual45=size(r1);
[r2] = find(auto45==3);
labelsize_auto45=size(r2);
DiceCoef1 = 2*OverlapSize45/(labelsize_manual45+labelsize_auto45);
manual44(manual44>0)=2;
auto44(auto44>0)=3;
Overlap44 = auto44-manual44;
[r3] = find(Overlap44==1);
OverlapSize44=size(r3);
[r4] = find(manual44==2);
labelsize_manual44=size(r4);
[r5] = find(auto44==3);
labelsize_auto44=size(r5);
DiceCoef2 = 2*OverlapSize44/(labelsize_manual44+labelsize_auto44);

AvgDice = (DiceCoef1 + DiceCoef2)/2;

end