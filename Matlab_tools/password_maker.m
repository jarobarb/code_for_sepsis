list1 = char(64+[1:26]);
list2 = char(96+[1:26]);
list3 = num2str([0:9]')';
list = [list1,list2,list3];

inds = ceil(numel(list)*rand(8,1));

list(inds)