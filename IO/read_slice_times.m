function time = read_slice_times()
    % READ_SLICE_TIMES  return the times associated with each slice file
    %                   temporal data is stored in slice_times.txt
    %
    % David Deepwell, 2019

    dat = importdata('slice_times.txt', ',', 1);
    time = dat.data(:,2);
end
