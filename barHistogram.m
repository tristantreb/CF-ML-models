function barHistogram(data, plottitle, xtitle)
% bar plot with interventions grouped by duration
% Input: vector that need groupby, xlabel, title
% Output: bar plot
arguments
    data (:,1) double
    plottitle (1,1) string = "Distribution of the data inputs"
    xtitle (1,1) string = "Data"
end

Intr = table(data,'VariableNames',{'data'});
Intr = groupcounts(Intr, 'data');
bar(Intr.data, Intr.GroupCount)
grid('on')
title(plottitle)
xlabel(xtitle)
ylabel('Frequency')
% automatic xticks increment selection
xspan=max(Intr.data)-min(Intr.data);
if xspan <= 10
    xticksincrement=1;
elseif xspan <= 20
    xticksincrement=2;
else
    xticksincrement=5;
end
xticks(min(Intr.data):xticksincrement:max(Intr.data))
% automatic yticks increment selection
if max(Intr.GroupCount) <= 10
    yticksincrement=1;
elseif max(Intr.GroupCount) <= 20
    yticksincrement=2;
elseif max(Intr.GroupCount) <= 50
    yticksincrement=5;
else
    yticksincrement=10;
end
yticks(0:yticksincrement:max(Intr.GroupCount)+yticksincrement)
end

