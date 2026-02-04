clc
close all

full_data = [ts', X];
data_table = table(full_data);
writetable(data_table, './assets/output_data.csv');