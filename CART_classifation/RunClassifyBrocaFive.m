subject_list = importdata(['/scr/murg2/MachineLearning/FiveNewSubjects.txt']);

cost = [];

for i=1:length(subject_list)
[Cost] = ClassifyBroca(subject_list(i));

cost = horzcat(cost, Cost);

fprintf('Subject # %u\n',i)

end

save('/scr/murg2/MachineLearning/ClassCost5subjects.mat', '-v7.3', 'cost');