function T = TTable(t)

global TTableData

T = interp1(TTableData(:,1), TTableData(:,2), t);

end