%MIKOLAJCZYK_REPEATABILITY_PLOT Plots repeatability results
%
% See also MIKOLAJCZYK_REPEATABILITY.M .

function [ r,c,classes ] = mikolajczyk_repeatability_plot(classes,classnames,methods,results,varargin)

conf = struct('debugLevel', 0, ...
              'dataset', '', ...
              'legends', [], ...
              'imgSaveDir', './', ...
              'imgSaveFile', 'repeatabilityplot_IMG_SAVE', ...
              'saveAlsoFig', 1,...
              'settings', []);
conf = mvpr_getargs(conf,varargin);

if isempty(conf.legends)
  conf.legends = methods;
end;

% Add all -row for results over all classes
classnames{numel(classnames)+1,1} = 'all';

% There may be some problems with some detectors, so we convert 
% possible NaNs to 0
for i=1:numel(results)
    results(i).repeatability(isnan(results(i).repeatability)) = 0;
end

repeat_avg = zeros(numel(methods), size(classnames,1));
repeat_med = zeros(numel(methods), size(classnames,1));

corresp_med = zeros(numel(methods), size(classnames,1));
corresp_avg = zeros(numel(methods), size(classnames,1));

% numel(results) == number of methods that results are available
for i=1:numel(results)
  rep = results(i).repeatability;
  cor = results(i).n_correspondences;

  % For each class
  for c=1:numel(classnames)-1
    repeat_med(i, c) = median(rep(classes==c));
    repeat_avg(i, c) = mean(rep(classes==c));
    corresp_med(i, c) = median(cor(classes==c));
    corresp_avg(i, c) = mean(cor(classes==c));
  end

  % All
  repeat_med(i, numel(classnames)) = median(rep);
  repeat_avg(i, numel(classnames)) = mean(rep);
  corresp_med(i, numel(classnames)) = median(cor);
  corresp_avg(i, numel(classnames)) = mean(cor);

  fprintf([conf.legends{i} ': avg corresp=%.2f, avg repeat=%.2f\n'],...
          corresp_avg(i,end),repeat_avg(i,end));
  corresp_avg(i,end)/(repeat_avg(i,end)/100)  
  
  
  %SAVE corresp_avg and repeat_avg
  dirPath = '/home/antti/Desktop/Work/descriptor_vocbenchmark/matlab/DetectorTesting/';
  filePath = [dirPath,'DetectorResults_',conf.legends{i},'.mat'];
  
  value = struct('corresp_avg', corresp_avg(i,end),...
                        'repeat_avg', repeat_avg(i,end));%, ...
                        %'configures', conf.settings(i));
  if 0~= exist(filePath)   
    load(filePath,'-mat');
    temp = cell2mat(values);
    %chech if the value already exists
    if ~any([temp.corresp_avg] == value.corresp_avg) || ... 
            ~any([temp.repeat_avg] == value.repeat_avg)
    values{end+1} = value;
    save(filePath,'values');
    end
    
  else
      values ={};
      values{1} = value;
      %save(filePath,'values');
  end

end

% Plot average repeatability for all the groups and methods
hrepavg = figure('Name', 'Average repeatability', 'Visible', 'off');
set(gcf, 'color', 'white');
hbarrepmean = barh(repeat_avg');
set(gca, 'YTickLabel', classnames); 
set(gca, 'FontSize', 14);
% title(['Average repeatability rates per image for each dataset']);
xlabelrep = 'Repeatability [%]: #(correspondences)/#(regions)';
xlabel(xlabelrep);
% ylabel('Dataset');
hold on;
legend(conf.legends, 'Location', 'Best','interpreter','none');
%fileOut = ip_get_fullname('result-figure', ...
%                          'detrep-avg.png', ...
%                           conf.dataset, '', ...
%                          'tempSaveDir', conf.tempSaveDir);
fileOut = fullfile(conf.imgSaveDir,[conf.imgSaveFile '-detrep-avg.png']);
if conf.saveAlsoFig
    saveas(gcf,strrep(fileOut,'.png', '.fig'));
end
export_fig(fileOut);
% close all;
% close(h);

% Plot median repeatability for all the groups and methods
hrepmed=figure('Name', 'Median repeatability', 'Visible', 'off');
set(gcf, 'color', 'white');
hbarrepmed = barh(repeat_med');
set(gca, 'YTickLabel', classnames); 
set(gca, 'FontSize', 14);
% title(['Median repeatability rates per image for each dataset']);
xlabel(xlabelrep);
% ylabel('Dataset');
hold on;
legend(conf.legends, 'Location', 'Best','interpreter','none');
%fileOut = ip_get_fullname('result-figure', ...
%                          'detrep-median.png', ...
%                           conf.dataset, '', ...
%                          'tempSaveDir', conf.tempSaveDir);
%
fileOut = fullfile(conf.imgSaveDir,[conf.imgSaveFile '-detrep-median.png']);
if conf.saveAlsoFig
    saveas(gcf,strrep(fileOut,'.png', '.fig'));
end
export_fig(fileOut);
% close;
% close(h);



% Plot average correspondences
% figure('Visible', 'off');
hcoravg = figure('Name', 'Average number of correspondences', 'Visible', 'off');
set(gcf, 'color', 'white');
hbarcormean = barh(corresp_avg');
set(gca, 'YTickLabel', classnames);
set(gca, 'FontSize', 14);
xlabelcor = '#(correspondences)';
xlabel(xlabelcor);
% title(['Average number of correspondences']);
% ylabel('Dataset');
hold on;
legend(conf.legends, 'Location', 'Best','interpreter','none');
%fileOut = ip_get_fullname('result-figure', ...
%                          'detcor-avg.png', ...
%                           conf.dataset, '', ...
%                          'tempSaveDir', conf.tempSaveDir);
%
fileOut = fullfile(conf.imgSaveDir,[conf.imgSaveFile '-detcor-avg.png']);
if conf.saveAlsoFig
    saveas(gcf,strrep(fileOut,'.png', '.fig'));
end
export_fig(fileOut);
% close;


% Plot median correspondences
hcormed=figure('Visible', 'off', 'Name', 'Median number of correspondences');
% figure('Name', 'Median number of correspondences per image for each dataset');
set(gcf, 'color', 'white');
hbarcormed = barh(corresp_med');
set(gca, 'YTickLabel', classnames);
set(gca, 'FontSize', 14);
xlabel(xlabelcor);
% title(['Median number of correspondences per image for each dataset']);
% ylabel('Dataset');
hold on;
legend(conf.legends, 'Location', 'Best','interpreter','none');
%fileOut = ip_get_fullname('result-figure', ...
%                          'detcor-median.png', ...
%                           conf.dataset, '', ...
%                          'tempSaveDir', conf.tempSaveDir);
%
fileOut = fullfile(conf.imgSaveDir,[conf.imgSaveFile '-detcor-median.png']);
if conf.saveAlsoFig
    saveas(gcf,strrep(fileOut,'.png', '.fig'));
end
export_fig(fileOut);
% close(h);


% Build matrices which have all the repeatabilities and correspondences
repeat = [];
corresp = [];
for i=1:numel(methods)
    repeat = [ repeat results(i).repeatability ];
    corresp = [ corresp results(i).n_correspondences ];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AGB: Present relevant results

fprintf('Showing detector results\n');
% Generate a figure for showing all results
hfig=figure('Name', 'Detector results', 'Color', 'white');
% Show all bar charts in subfigures
h221=subplot(221);
h222=subplot(222);
h223=subplot(223);
h224=subplot(224);
copyobj(hbarrepmean, h221);
copyobj(hbarrepmed, h222);
copyobj(hbarcormean, h223);
copyobj(hbarcormed, h224);
% Add horizontal labels
xlabel(h221, xlabelrep);
xlabel(h222, xlabelrep);
xlabel(h223, xlabelcor);
xlabel(h224, xlabelcor);
% Add vertical labels
% ylim(h221,[0 n_classes+1]);
% axis(h132,[0 100 0 n_classes+1]);
n_classes = numel(classnames);
set(h221,'ytick',1:n_classes);
set(h222,'ytick',1:n_classes);
set(h223,'ytick',1:n_classes);
set(h224,'ytick',1:n_classes);
set(h221,'yticklabel', classnames);
set(h222,'yticklabel', cell(1,n_classes));
set(h223,'yticklabel', classnames);
set(h224,'yticklabel', cell(1,n_classes));
% Add plot titles
title(h221, get(hrepavg,'Name'));
title(h222, get(hrepmed,'Name'));
title(h223, get(hcoravg,'Name'));
title(h224, get(hcormed,'Name'));
% Add legend
hleg=legend(hbarcormed, conf.legends,'interpreter','none');%, 'location', 'bestoutside');
hleg=copyobj(hleg,hfig);
set(hleg, 'units','normalized','position',[0.85,0.85,.05,.05]);%[.78,.85,.06,.06]);
% Set font
set(h221, 'FontSize', 14);
set(h222, 'FontSize', 14);
set(h223, 'FontSize', 14);
set(h224, 'FontSize', 14);
% Close original figures

close([hrepavg hrepmed hcoravg hcormed]);
end % function