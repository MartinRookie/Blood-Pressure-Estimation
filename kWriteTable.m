function kWriteTable(data,RowName,ColumnNames,savepath)
% 一行代码实现表格数据的写入的封装函数

% data：表格中的数据
% RowNames：行名称，cell类型，内容为字符串。需要与data的行数相同
% ColumnNames：列名称，cell类型，内容为字符串。需要与data的列数相同
% savepath：包含文件名的文件保存路径，可以为相对路径或绝对路径，比如"data.csv"或'F:\data.csv'

%  Copyright (c) 2023 Mr.看海 All rights reserved.
%  原文链接 https://zhuanlan.zhihu.com/p/603604365
%  代码地址：http://khsci.com/docs/index.php/2023/02/05/export/

%% 1.格式化输入数据方向
RowName = RowName(:);

%% 2.生成数据表格T
collen = length(ColumnNames);
for i = 1:collen
    eval([ColumnNames{i},'=data(:,i)']);
end

tableStr = 'T = table(';
for i = 1:collen
    tableStr = [tableStr,ColumnNames{i},','];
end
tableStr = [tableStr,'''RowNames'',RowName)'];
eval(tableStr);
%% 3.写入文件
writetable(T,savepath,'WriteRowNames',true)  %写入文件

end