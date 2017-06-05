function model = addMetcomp(model)
    comps = [];
    for i = 1:length(model.mets)
        curComp = strsplit(model.mets{i}, '[');
        curComp = strsplit(curComp{2}, ']');
        comps = [comps; curComp{1}];
    end
    model.comps = cellstr(unique(comps));
    model.metComps = zeros(length(comps),1);

    for i = 1:length(comps)
        compNr = find(ismember(model.comps, comps(i)));
        model.metComps(i) = compNr;
    end
end
