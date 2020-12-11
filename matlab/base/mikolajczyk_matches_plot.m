%function [] = mikolajczyk_matches_plot(imgsA_, real_matches, correct_matches, varargin)
%
%conf = struct('debugLevel', 0, ...
%    'detectors', {'sift'}, ...
%    'descriptors', {'sift'}, ...
%    'dataset', '', ...
%    'tempSaveDir', '.', ...
%    'bestMatchMax', 1,...
%    'coverages', 1,...
%    'saveEPS', 0, ... %AGB
%    'saveImg', 1,... % AGB
%    'imgFormat', 'png',...% AGB
%    'showResults', true); % AGB
%conf = mvpr_getargs(conf,varargin);
function [ r,c,classes ] = mikolajczyk_matches_plot(classes,classnames,methods,real_matches,correct_matches,varargin)

conf = struct('debugLevel', 0, ...
              'dataset', '', ...
              'legends', [], ...
              'imgSaveDir', './', ...
              'imgSaveFile', 'repeatabilityplot_IMG_SAVE', ...
              'bestMatchMax',1,...
              'saveAlsoFig', 1,...
              'coverages',10,...
              'imgFormat', 'png',...
              'results',[],...
              'configures',[]);
conf = mvpr_getargs(conf,varargin);

if isempty(conf.legends)
  conf.legends = methods;
end;

% Add all -row for results over all classes
classnames{numel(classnames)+1,1} = 'all';

n_methods = numel(methods);
n_classes = numel(classnames);
%n_bestMatch = size(real_matches{1},1);
n_bestMatch = conf.bestMatchMax;
means=zeros(n_methods, n_classes);
medians=zeros(n_methods, n_classes);
% covperimg=zeros(n_methods, imgTot);

for method=1:n_methods
    for c=1:n_classes-1
        %medians(method,c) = median(real_matches{method}(conf.bestMatchMax,classes==c));
        %means(method,c) = mean(real_matches{method}(conf.bestMatchMax,classes==c));
        medians(method,c) = median(real_matches{method}(n_bestMatch,classes==c));
        means(method,c) = mean(real_matches{method}(n_bestMatch,classes==c));
        
        for m=1:conf.coverages
          %coverages(method,c,m) = sum(real_matches{method}(conf.bestMatchMax,classes==c)>=m);
          coverages(method,c,m) = sum(real_matches{method}(n_bestMatch,classes==c)>=m);
        end
        
        if conf.debugLevel
            classnames{c}
            medians(method,c)
            means(method,c)
        end
    end
    
    %medians(method,n_classes) = median(real_matches{method}(conf.bestMatchMax,:));
    %means(method,n_classes) = mean(real_matches{method}(conf.bestMatchMax,:));
    medians(method,n_classes) = median(real_matches{method}(n_bestMatch,:));
    means(method,n_classes) = mean(real_matches{method}(n_bestMatch,:));
    
    for m=1:conf.coverages
        %coverages(method,n_classes,m) = sum(real_matches{method}(conf.bestMatchMax,:)>=m)/(n_classes-1);
        coverages(method,n_classes,m) = sum(real_matches{method}(n_bestMatch,:)>=m)/(n_classes-1);
    end
    
    if conf.debugLevel
        ['all']
        medians(method,n_classes)
        means(method,n_classes)
    end
        descMedian(method) = medians(method,n_classes);
        descMeans(method) = means(method,n_classes);
    %     cpi = real_matches{method}(conf.bestMatchMax,:);
    %     % Coverage per image over 10 will be displayed as 10
    %     cpi(cpi>10)=10;
    %     covperimg(method,:) = cpi; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fprintf('\tSaving %d coverage files:', n_coverages);

%
% Upper limit of coverage: number of image pairs for each class
%

xmaxCover = numel(classes) / (n_classes-1);

for c = 1:conf.coverages
    % Print N
    fprintf(' %d', c);
    
    % Output file base name
    coverageFile = fullfile(conf.imgSaveDir,...
             ['result-figure-' conf.dataset '-coverage-' num2str(c) '-besmatches-' ...
              num2str(n_bestMatch)]); 
    %coverageFile = ip_get_fullname('result-figure', ...
    %                               ['coverage-' num2str(c) '-b' num2str(conf.bestMatchMax)], ...
    %conf.dataset, '', ...
    %    'tempSaveDir', conf.tempSaveDir);
    
    % Coverages
    if conf.debugLevel > 0
      figure('visible', 'on', 'color', 'white');
    else
      figure('visible', 'off', 'color', 'white');
    end;      
    hbar = barh(coverages(:,:,c)');
    axis([0 1.1*xmaxCover 0 n_classes+1]);
    set(gca, 'ytick', 1:n_classes);
    set(gca,'yticklabel', classnames);
    set(gca, 'FontSize', 14);
    %xlabel(['Coverage-' num2str(c) ' (K=' num2str(conf.bestMatchMax) ')']);
    xlabel(['Coverage-' num2str(c) ' (K=' num2str(n_bestMatch) ')']);

    if conf.debugLevel > 0
      input('DEBUG[1]: Coverage plot <RETURN>');
    end;
      
    % Export
    %if conf.saveEPS, saveas(gcf, [coverageFile '.eps']); end
    export_fig([coverageFile '.' conf.imgFormat]);
    
    %% Only save legend once, i.e. after processing last coverage
    %if c==n_coverages
       % Legend file
       coverageLegendFile = ip_get_fullname('result-figure', ...
          ['coverage-b' num2str(conf.bestMatchMax) '-legend'], ...
           conf.dataset, '', ...
           'tempSaveDir', conf.tempSaveDir);
       % Add legend to last figure
       legend(conf.legends, 'location', 'bestoutside');
       % Hide bars and axes
       set([hbar gca], 'visible', 'off');
       % Export
       if conf.saveEPS, saveas(gcf, [coverageLegendFile '.eps']); end
       if conf.saveImg, export_fig(coverageLegendFile, imgFormat); end
       % Restore
       % TODO: Not needed since we don't present coverages in the end
       % set([hbarcover gca], 'visible', 'on');        
    %end
    
    % Close down figure
    close(gcf);
end
fprintf(' done!\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\tSaving mean/median files...\n');

%
% AGB: Heuristically set a common x-axis limit
%
%xmaxMeanMedian = 1;

% Labels
%xlabelmean = ['Mean #(correct matches) (K=' num2str(conf.bestMatchMax) ')'];
%xlabelmedian = ['Median #(correct matches) (K=' num2str(conf.bestMatchMax) ')'];
xlabelmean = ['Mean #(correct matches) (K=' num2str(n_bestMatch) ')'];
xlabelmedian = ['Median #(correct matches) (K=' num2str(n_bestMatch) ')'];

% Files
%meanFile = ip_get_fullname('result-figure', ['mean-all-b' num2str(conf.bestMatchMax)], ...
%    conf.dataset, '', ...
%    'tempSaveDir', conf.tempSaveDir);
meanFile = fullfile(conf.imgSaveDir,...
                    ['result-figure-' conf.dataset '-mean-all-besmatches-' ...
                    num2str(n_bestMatch)]); 
%medianFile = ip_get_fullname('result-figure', ['median-all-b' num2str(conf.bestMatchMax)], ...
%    conf.dataset, '', ...
%    'tempSaveDir', conf.tempSaveDir);
medianFile = fullfile(conf.imgSaveDir,...
                      ['result-figure-' conf.dataset '-median-all-besmatches-' ...
                    num2str(n_bestMatch)]); 
%{
% Number of matches average
if conf.debugLevel > 0
  hmean = figure('visible', 'on', 'color', 'white');
else
  hmean = figure('visible', 'on', 'color', 'white');
end;
hbarmean = barh(means');
set(gca,'yticklabel', classnames);
set(gca, 'FontSize', 14);
xlabel(xlabelmean);
%axis([0 1.1*xmaxMeanMedian 0 n_classes+1])
%if conf.saveEPS, saveas(gcf, [meanFile '.eps']); end
legend(conf.legends, 'Location', 'Best','interpreter','none');
export_fig([meanFile '.' conf.imgFormat]);
% Add legend
legend(conf.legends, 'Location', 'Best','interpreter','none');
% Hide bars and axes
set([hbarmean gca], 'visible', 'off');
% Export
%export_fig([meanFile '-legend'], imgFormat);
export_fig([meanFile '-legend.' conf.imgFormat]);
% Restore
set([hbarmean gca], 'visible', 'on');
if conf.debugLevel > 0
  input('DEBUG[1]: Average matching plot <RETURN>');
end;

% Number of matches median
hmedian = figure('visible', 'off', 'color', 'white');
hbarmedian = barh(medians');
set(gca,'yticklabel', classnames);
set(gca, 'FontSize', 14);
%axis([0 1.1*xmaxMeanMedian 0 n_classes+1])
%title('median of matches per class');
xlabel(xlabelmedian);
%ylabel('dataset');
% legend(method_strs, 'location', 'bestoutside');
%if conf.saveEPS, saveas(gcf, [medianFile '.eps']); end
legend(conf.legends, 'Location', 'Best','interpreter','none');
export_fig([medianFile '.' conf.imgFormat]);
% Add legend
legend(conf.legends, 'Location', 'Best','interpreter','none');
% Hide bars and axes
set([hbarmedian gca], 'visible', 'off');
% Export
export_fig([medianFile '-legend' '.' conf.imgFormat]);
% Restore
set([hbarmedian gca], 'visible', 'on');

%
% Save table
%
h = figure('name', 'Overall results table');
uitable('RowName', conf.legends,...
    'ColumnName', {'Mean', 'Median'},...
    'Data', [means(:,end) medians(:,end)],...
        'Position',[20 20 400 100]);
tableFile = fullfile(conf.imgSaveDir,...
                     ['result-figure-' conf.dataset '-mean-median-all-besmatches-' ...
                    num2str(n_bestMatch) '-table']); 
%tableFile = ip_get_fullname('result-figure', ['mean-median-all-b' num2str(conf.bestMatchMax) '-table'], ...
%    conf.dataset, '', ...
%    'tempSaveDir', conf.tempSaveDir);
%if conf.saveEPS, saveas(h, [tableFile '.eps']); end
%export_fig([tableFile '.' conf.imgFormat]);
saveas(h, [tableFile '.png']);
close(h);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AGB: After 2012.11.15, this code suddenly crashed for n_methods==4 !?
%
% fprintf('\tSaving number of features...\n');

% % AGB: This is not medians!?
% medianFile = ip_get_fullname('result-figure', ['num-features-b' num2str(conf.bestMatchMax) '.png'], ...
%                              conf.dataset, '', ...
%                              'tempSaveDir', conf.tempSaveDir);

%keyboard
for d = 1:n_methods % Loop over descriptor types
     	for m = 1:size(correct_matches{d},2) % Loop over descriptors in A
     		n_features(d,m) = size(correct_matches{d}{m},2);
     	end
    x = correct_matches{d}'; % [features_A bestMatchMax] logical matrix
         x2 = correct_matches{d}; % [bestMatchMax features_A] logical matrix
    % AGB: If number of best matches > 1 here, the above statement does not extract
    % the right information, since x will contain a cell of a many
    % [bmm M] logicals (one per image), and later we expect [1 M] logicals
    for idx=1:numel(x) % Loop over all images
        % x{idx} = x{idx}(conf.bestMatchMax,:); % AGB: Take out only the right row of logical array
        x{idx} = x{idx}(n_bestMatch,:); % AGB: Take out only the right row of logical array
    end
    
    % n_features_per_class:
    %   (unused)
    % n_features_per_image:
    %   [n_methods n_images] matrix of #(correspondences) per image
    for c=1:n_classes-1 % Loop over all classes, except 'all'
         		n_features_per_class(d,c+1) = size(cell2mat(x2(classes==c)),2) / sum(classes==c);
        y = x(classes==c); %
        idxs = find((classes==c)==1);
        %        for i = 1:size(y,2)
        for i = 1:size(y,1) % over all images of this class
            % 			n_features_per_image(d,idxs(i)) = size(y{i},2);
            % AGB: Corrected this to not just take the size, but only the
            % true's
            % NO: That is equal to the number of matches
            n_features_per_image(d,idxs(i)) = size(y{i},2);
        end
        
    end
    %sum all features from different classes and get the average
    numOfRegions(d) = (sum(n_features_per_class(d,:)))/(n_classes-1);
    
    display(['NumOfRegions ',num2str(numOfRegions(d))]);
    display(['Median ',num2str(descMedian(d))]);
    display(['Mean ',num2str(descMeans(d))]);

      %SAVE corresp_avg and repeat_avg
    dirPath = '/home/antti/Desktop/Work/descriptor_vocbenchmark/matlab/DetectorTesting/';
    filePath = [dirPath,'DescriptorResults_',conf.legends{d},'.mat'];
    
    
  
  
    value = struct('numOfRegions', numOfRegions(d),...
                        'medians', descMedian(d), ...
                        'means', descMeans(d));%,...
    %                    'configures', conf.configures(d));
    if 0~= exist(filePath)   
        load(filePath,'-mat');
        temp = cell2mat(values);
        %check if the value already exists
        if ~any([temp.means] == value.means) || ... 
                ~any([temp.medians] == value.medians)
        values{end+1} = value;
        save(filePath,'values');
        end
    
    else
      values={};
      values{1} = value;
%      save(filePath,'values');
    end
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\tSaving averaged means...\n');

%
% AGB: Heuristically set a common x-axis limit
%
xmaxMeanAllAvg = 0.3;

% Labels
%xlabelmeanall = ['Mean #(correct matches)/#(features) (K=' num2str(conf.bestMatchMax) ')'];
xlabelmeanall = ['Mean #(correct matches)/#(features) (K=' num2str(n_bestMatch) ')'];

% Files
%meanFile = ip_get_fullname('result-figure', ['mean-all-averaged-b' num2str(conf.bestMatchMax)], ...
%meanFile = ip_get_fullname('result-figure', ['mean-all-averaged-b' num2str(n_bestMatch)], ...
%                           conf.dataset, '', ...
%                           'tempSaveDir', conf.tempSaveDir);
meanFile = fullfile(conf.imgSaveDir,...
                    ['result-figure-' conf.dataset '-mean-all-averaged-besmatches-' ...
                    num2str(n_bestMatch)]); 

for method=1:n_methods
%    matches = ( real_matches{method}(conf.bestMatchMax,:) ./ n_features_per_image(method,:) );
    matches = ( real_matches{method}(n_bestMatch,:) ./ n_features_per_image(method,:) );
    for c=1:n_classes-1
        m_means(method,c) = mean(matches(classes==c));
    end
end
% Number of matches average
%keyboard
% figure('Visible', 'off');
hmeanall = figure('visible', 'off', 'color', 'white');
%barh((means ./ n_features_per_class)');
hbarmeanall = barh(m_means');
set(gca,'yticklabel', classnames);
set(gca, 'FontSize', 14);
xlabel(xlabelmeanall);
% AGB: Set max probability < 1
axis([0 1.1*xmaxMeanAllAvg 0 n_classes+1])
% legend(method_strs, 'location', 'bestoutside');

%if conf.saveEPS, saveas(gcf, [meanFile '.eps']); end
export_fig([meanFile '.' conf.imgFormat]);
% Add legend
legend(conf.legends, 'Location', 'Best','interpreter','none');
% Hide bars and axes
set([hbarmeanall gca], 'visible', 'off');
% Export
export_fig([meanFile '-legend'], conf.imgFormat);
% Restore
set([hbarmeanall gca], 'visible', 'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% AGB: Present relevant results
%if conf.showResults
    %fprintf('Showing descriptor results for %d nearest neighbor(s)\n', conf.bestMatchMax);
    fprintf('Showing descriptor results for %d nearest neighbor(s)\n', n_bestMatch);
    % Generate a figure for showing all results
%    hfig=figure('name', ['Results for ' num2str(conf.bestMatchMax) ' nearest neighbor(s)'], 'color', 'white');
    hfig=figure('name', ['Results for ' num2str(n_bestMatch) ' nearest neighbor(s)'], 'color', 'white');
    % Show all bar charts in subfigures
    h131=subplot(131);
    h132=subplot(132);
    h133=subplot(133);
    copyobj(hbarmean, h131);
    copyobj(hbarmedian, h132);
    copyobj(hbarmeanall, h133);
    % Add horizontal labels
    xlabel(h131, xlabelmean);
    xlabel(h132, xlabelmean);
    xlabel(h133, xlabelmeanall);
    % Add vertical labels
%    axis(h131,[0 1.1*xmaxMeanMedian 0 n_classes+1]);
%    axis(h132,[0 1.1*xmaxMeanMedian 0 n_classes+1]);
%    axis(h133, [0 1.1*xmaxMeanAllAvg 0 n_classes+1])
    set(h131,'ytick',1:n_classes);
    set(h132,'ytick',1:n_classes);
    set(h133,'ytick',1:n_classes);
    set(h131,'yticklabel', classnames);
    set(h132,'yticklabel', cell(1,n_classes+1));
    set(h133,'yticklabel', cell(1,n_classes+1));
    % Add plot titles
    title(h131, get(hmean,'Name'));
    title(h132, get(hmedian,'Name'));
    title(h133, get(hmeanall,'Name'));
    % Add legend
    hleg=legend(hbarmeanall, conf.legends,'interpreter','none');%, 'location', 'bestoutside');
    hleg=copyobj(hleg,hfig);
    set(hleg, 'units','normalized','position',[0.85,0.825,.05,.05]);%[.78,.85,.06,.06]);
    % Set font
    set(h131, 'FontSize', 14);
    set(h132, 'FontSize', 14);
    set(h133, 'FontSize', 14);
    % Close original figures
    close([hmean hmedian hmeanall]);
%end
%}
end % function
