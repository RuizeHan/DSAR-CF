function [max_response,PSR] = get_PSR(response_map)

[max_resp_row, max_rows] = max(response_map, [], 1);
[max_response, max_col] = max(max_resp_row, [], 2);
max_row=max_rows(max_col);
% size_r=size(response_map);
% maxregion_index = [max_row-5:max_row+5;max_col-5:max_col+5];


response_map(max(1,max_row-5):max_row+5,max(1,max_col-5):max_col+5)=2;
response_map((response_map==2))=[];
mean_sidelobe=mean(response_map);
std_response_map=std(response_map);
PSR= (max_response-mean_sidelobe)/std_response_map;

