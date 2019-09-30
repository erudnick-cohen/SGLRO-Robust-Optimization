totalC = 0;
for i = 2:length(countList)
    totalC = totalC+countList{i}.counter;
end
conTotal = totalC
objTotal = countList{1}.counter