function msc = analysis10x_mean(name,dt,TH)
%cd ../stacks
%clear all
name = [name,'/'];
cfiles=dir([name,'img*CFP-CFP_00*.tif']);yfiles=dir([name,'img_*CFP-YFP_00*.tif']);nn=length(cfiles);
h=fspecial('average',[10 10]);
n=1;t=1;
TH = TH -100;
%get data from stack-files (fret-concentration and locations of cells)

    for t=1:nn
        imc=double(imread([name,cfiles(t).name]))-100;
        imy=double(imread([name,yfiles(t).name]))-100;
        %imc = imc(36:46,206:218);
        %imy = imy(36:46,206:218);
        mask=imy>TH;
        %mask=imfilter(imy>TH,h);
        msc(t)=mean(mean(imc(mask)./imy(mask)));
        msc0(t,1)=mean(mean(imc./imy));
        msc0(t,2)=mean(mean(imc));
        msc0(t,3)=mean(mean(imy));
        msc0(t,4)=mean(imc(mask)./imy(mask));
        if mod(t,10)==0,disp(t),end
    end
save([name,'data.mat'], 'msc*','dt'); %, msc* dt
%----------------------
% make compound figure


    %clear all
    load([name,'data']);
    figure;
    plot(dt*[1:size(msc,2)],(msc(:)));
    xlabel('Time(min)');ylabel('FRET-Intensity');
    saveas(gcf,[name,'fret.fig']);
    
    figure;
    plot(dt*[1:size(msc,2)],(msc0(:,2)/mean(msc0(:,2))));
    hold all;
    plot(dt*[1:size(msc,2)],(msc0(:,3)/mean(msc0(:,3))));
    xlabel('Time(min)');ylabel('CFP and YFP Traces');
    saveas(gcf,[name,'cfp-yfp.fig']);
% 
%     init_thresh = (mean(msc) + max(msc))/2;
%     [peaks2,time2, intervals, time3] = getIntervals(msc, 3, dt, init_thresh);
%     figure;
%     plot(dt*[1:length(msc)],msc,'.-',time2,peaks2,'ro');
%     title('Peaks of Signal');
%     xlabel('Time(min)');
%     saveas(gcf, [name,'peaks.fig']);
%     
%     figure;
%     plot(time3,intervals,'.-');
%     title('Intervals/Period between Peaks');
%     xlabel('Time after Acquisition (min)');
%     ylabel('Period (min)');
%     saveas(gcf, [name,'intervals.fig']);
%     
%     density(1) = 0;
%     for i=1:length(intervals)
%         density(i) = getDensity(name, dt, time3(i), TH);
%     end
%     
%     figure;
%     plot(time3,density,'.-');
%     title('Density vs. Time');
%     xlabel('Time after Acquisition (min)');
%     ylabel('Density ~ cells/area');
%     saveas(gcf, [name,'density.fig']);
%     
%     figure;
%     plot(intervals,density,'o');
%     corrmat = corr([density' , intervals']);
%     title(['Density versus Intervals (Correlation = ', num2str(corrmat(1,2)), ')']);
%     xlabel('Density (Cells/Area)');
%     ylabel('Period between Peaks (min)');
%     saveas(gcf, [name,'densityvsintervals.fig']);
end