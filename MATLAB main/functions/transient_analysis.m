function transient_analysis(y,NR_periods)
    y_reshape = reshape(y,[],NR_periods);
    y_diff = y_reshape - y_reshape(:,end);
    y_diff_rms = rms(y_diff);
    figure
    plot(1:NR_periods-1,y_diff_rms(1:NR_periods-1),'.--','Markersize',15)
    xticks(1:NR_periods-1)
    grid on
    title("RMS of difference between p-th and last period")
    xlabel("Period p")
end

