

for Iter = 1: size(alldata,2)

    temp = alldata{1,Iter}(1,:);
    sqtemp = temp .* temp;
    sqtemp(isnan(sqtemp))=0.000000001;
    meansqtemp = mean(sqtemp(1,:));
    stdsqtemp = std(sqtemp(1,:));
    SM_Spike_Data = [];

    for i = 1:size(sqtemp,2)
        if sqtemp(1,i) > (meansqtemp + 0.2*stdsqtemp); % 0.00001;
            SM_Spike_Data(Iter,1,i) = 1;
        else
            SM_Spike_Data(Iter,1,i) = 0;
        end
    end

    % Downsample to 200 Hz for plotting with downsampled newdata
    a = []; aa = [];
    a = squeeze(SM_Spike_Data(Iter,1,:));
        
    for i = 1: 50: floor(size(SM_Spike_Data,3))-50
        aa(end+1) = mean(a(i:i+50));
    end
    
    figure('Name',['Iter_' num2str(Iter)]);
    plot(aa);
    hold on;
    plot(newdata(Iter,:));

end

autoArrangeFigures(3,3)

save('SM_Spike_Data','SM_Spike_Data');

