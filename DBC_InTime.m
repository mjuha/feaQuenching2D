function DBC_InTime(t)

global nodeSet lineNode DBCSet U

[m,~] = size(lineNode);
numericCell = nodeSet(:,1);
numericVector = cell2mat(numericCell);
% check nodeSet
for i=1:m
    phytag = lineNode(i,1);
    % search in nodeSet
    [row,~] = find(numericVector==phytag);
    if  isempty(row)
        continue
    end
    % find in DBCSet first
    name = nodeSet{row,2};
    % find in DBCSet first
    [row,~] = find(strcmp(DBCSet,name),3);
    for j=1:length(row)
        dof = DBCSet{row(j),2};   
        if strcmp(dof,'TFunction')
            value = DBCSet{row(j),3};
            if value(1) == 1 % linear function
                val = TLinear(t,value(2),value(3)); % time, To and m
            elseif value(1) == 2 % data table
                val = TTable(t);
            end
            U(1,lineNode(i,2:3)) = val;       
        end
    end
end

end