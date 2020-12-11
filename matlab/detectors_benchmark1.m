% DETECTORS_BENCHMARK1 VOC benchmark detectors
%
% Runs the detector benchmark described in ref. [1] & [2].
% Evaluation is based on the
% configuration file settings in DETECTOR_BENCHMARK1_CONF.M .
%
% References
% [1] Lankinen, J.; Kangas, V.; Kamarainen, J.-K.: A comparison of
%     local feature detectors and descriptors for visual object
%     categorization by intra-class repeatability and matching, In 21th
%     Int. Conf. on Pattern Recognition (ICPR2012), Tsukuba Science
%     City, Japan, 2012.
% [2] Lankinen, J.; Buch, A.G.; Kangas, V.; Kamarainen,
%     J.-K.;Krueger N.: A Comparison of Feature Detectors and
%     Descriptors for Visual Object Class Matching,
%     Submitted to IEEE Trans. on Image Processing.
%
% See also DETECTOR_BENCHMARK1_CONF.M and DESCRIPTORS_BENCHMARK1.M .

% Clear display and start timing
clearvars -except confParameter timeFunction
close all;
tic

% 0. Config & Initialisation
% Run the configuration function
conf = detectors_benchmark1_conf();

imgTotA = mvpr_lcountentries(fullfile(conf.filelistRoot,conf.imgsA));
imgTotB = mvpr_lcountentries(fullfile(conf.filelistRoot,conf.imgsB));
if (imgTotA ~= imgTotB)
  error(['The files ''' conf.imgsA ''' and ''' conf.imgsB ''' contain different numbers of images']);
else
  imgTot = imgTotA;
end;
imgTot = 1; %FOR REDUCING THE TIME TO DEGUG
% 1.
% Compute and store the true transformations using the ground truth
% landmark locations (used to call ip_compare_write_transforms())
fhA = mvpr_lopen(fullfile(conf.filelistRoot,conf.landmarksA),'read');
fhB = mvpr_lopen(fullfile(conf.filelistRoot,conf.landmarksB),'read');
for imgInd = 1:imgTot
    
  if 0 == mod(imgInd,30)
     fprintf(['\r1. Computing ground truth transformations ' ...
           'between image pairs: %4d/%4d'],imgInd,imgTot);
  end
  fileEntryA = mvpr_lread(fhA);
  fileEntryB = mvpr_lread(fhB);

  if strcmp(fileEntryA, fileEntryB)
    warning(['Using the same file in A list and B list makes no sense:'...
             fileEntryA]);
  end

  [lmDir lmFile lmExt] = fileparts(fileEntryA{1});
  [lmDirB lmFileB lmExtB] = fileparts(fileEntryB{1});
  lmSavePath = fullfile(conf.tempSaveDir,conf.landmarksRoot,lmDir);
  if (~exist(lmSavePath,'dir'))
    mkdir(lmSavePath)
  end;
  lmSaveFile = [lmFile '_AtoB_' lmFileB '_H.txt'];

  if conf.enforceTransformations == true || ~ ...
        exist(fullfile(lmSavePath,lmSaveFile), 'file')
    % Load landmarks (groundtruth)
    dA = load(fullfile(conf.landmarksRoot,fileEntryA{1}));
    dB = load(fullfile(conf.landmarksRoot,fileEntryB{1}));
    % Compute transformation
    [H] = mvpr_h2d_corresp (dA', dB', 'hType', conf.hType);
    % Save transformation
    save(fullfile(lmSavePath,lmSaveFile), 'H', '-ascii');
  end;
end;
fprintf(' ...Done!\n');
mvpr_lclose(fhA);
mvpr_lclose(fhB);

% 2.
% Extract locate features from each image and store the results
for detdescNo = 1:length(conf.detdesc)
  %fprintf('\r2.%1d Run detector and descriptor: %s\n', detdescNo,conf.detdesc{detdescNo});
  % Open image list file and create sift features
  timeTotalElapsed = 0;
  timeElapsed = 0;
  fhA = mvpr_lopen(fullfile(conf.filelistRoot,conf.imgsA),'read');
  fhB = mvpr_lopen(fullfile(conf.filelistRoot,conf.imgsB),'read');
  for imgInd = 1:imgTot
    timeETA = (timeTotalElapsed / (imgInd-1)) * (imgTot-imgInd);
    if 0 == mod(imgInd,30)
    fprintf('\r  Extract image descriptors: %4d/%4d (ETA: %.1fs)', ...
            imgInd,imgTot, timeETA);
    end
    tic;
      fileEntryA = mvpr_lread(fhA);
      fileEntryB = mvpr_lread(fhB);
      
      [imgDirA imgFileA imgExtA] = fileparts(fileEntryA{1});
      [imgDirB imgFileB imgExtB] = fileparts(fileEntryB{1});
      
      fSavePath = fullfile(conf.tempSaveDir,conf.imageRoot,imgDirA);
      if (~exist(fSavePath,'dir'))
        mkdir(fSavePath)
      end;
      fSaveFileA = [imgFileA '_' conf.detdesc{detdescNo} '_features.mat'];
      fSaveFileB = [imgFileB '_' conf.detdesc{detdescNo} '_features.mat'];
      
      % Extract local feature (if not exist or enforced)
      if conf.enforceFeatures == true || ~exist(fullfile(fSavePath,fSaveFileA), 'file') ...
            || ~exist(fullfile(fSavePath,fSaveFileB), 'file')

        imgA = mvpr_imread(fullfile(conf.imageRoot, fileEntryA{1}),'range', [0 1]);
        imgB = mvpr_imread(fullfile(conf.imageRoot, fileEntryB{1}),'range', [0 1]);
        
        % Extract features (with the default parameters)
        [framesA descriptorsA configures] = ...
            mvpr_feature_extract_new(imgA, ...
                                     'method', conf.detdesc{detdescNo},...
                                     'debugLevel',conf.debugLevel);
        [framesB descriptorsB configures] = ...
            mvpr_feature_extract_new(imgB, ...
                                     'method', conf.detdesc{detdescNo},...
                                     'debugLevel',conf.debugLevel);                
        if imgInd == 1
            configures
        end
                                 
        % Save frames (detected regions) and their descriptors
        frames = framesA; descriptors = descriptorsA;
        save(fullfile(fSavePath,fSaveFileA),'frames', 'descriptors');
        frames = framesB; descriptors = descriptorsB;
        save(fullfile(fSavePath,fSaveFileB),'frames', 'descriptors');
      end;
      timeElapsed = toc;
      timeTotalElapsed = timeTotalElapsed + timeElapsed;
  end
%    confSet(detdescNo) = configures;
    fprintf(' ...Done!\n');
    mvpr_lclose(fhA);
    mvpr_lclose(fhB);
    
    %fprintf('\rDone! Total time elapsed: %.2fs, %.2fs per file.  \n', ...
    %        timeTotalElapsed, timeTotalElapsed / imgTot);
end;

% 3.
% Mask detected features using object foreground masks (if exist)
if (~isempty(conf.masksA) && ~isempty(conf.masksB))
  for detdescNo = 1:length(conf.detdesc)
    fprintf('\r3.%1d (%s) Applying foreground masks to detected regions\n', detdescNo,conf.detdesc{detdescNo});
    timeTotalElapsed = 0;
    timeElapsed = 0;
    fhA_im = mvpr_lopen(fullfile(conf.filelistRoot,conf.imgsA),'read');
    fhB_im = mvpr_lopen(fullfile(conf.filelistRoot,conf.imgsB),'read');
    fhA_ma = mvpr_lopen(fullfile(conf.filelistRoot,conf.masksA),'read');
    fhB_ma = mvpr_lopen(fullfile(conf.filelistRoot,conf.masksB),'read');
    for imgInd = 1:imgTot
      timeETA = (timeTotalElapsed / (imgInd-1)) * (imgTot-imgInd);
      if 0 == mod(imgInd,30)
      fprintf('\r  Image pair: %4d/%4d (ETA: %.1fs)', ...
              imgInd,imgTot, timeETA);
      end
      tic;
        fileEntryA_im = mvpr_lread(fhA_im);
        fileEntryB_im = mvpr_lread(fhB_im);
        fileEntryA_ma = mvpr_lread(fhA_ma);
        fileEntryB_ma = mvpr_lread(fhB_ma);
        
        [imgDirA imgFileA imgExtA] = fileparts(fileEntryA_im{1});
        [imgDirB imgFileB imgExtB] = fileparts(fileEntryB_im{1});
        
        % Names of feature files (from 2.)
        fSavePath = fullfile(conf.tempSaveDir,conf.imageRoot,imgDirA);
        fSaveFileA = [imgFileA '_' conf.detdesc{detdescNo} '_features.mat'];
        fSaveFileB = [imgFileB '_' conf.detdesc{detdescNo} '_features.mat'];

        % Masked feature files
        fSaveFileA_new = [imgFileA '_' conf.detdesc{detdescNo} '_features_masked.mat'];
        fSaveFileB_new = [imgFileB '_' conf.detdesc{detdescNo} '_features_masked.mat'];

        if conf.enforceFeatures == true || ~exist(fullfile(fSavePath,fSaveFileA_new), 'file') ...
            || ~exist(fullfile(fSavePath,fSaveFileB_new), 'file')

          % Load detected regions
          load(fullfile(fSavePath,fSaveFileA),'frames', 'descriptors');
          framesA = frames; descriptorsA = descriptors;
          load(fullfile(fSavePath,fSaveFileB),'frames', 'descriptors');
          framesB = frames; descriptorsB = descriptors;
          
          % Load the object foreground masks
          maskA = load(fullfile(conf.masksRoot, fileEntryA_ma{1}));
          maskB = load(fullfile(conf.masksRoot, fileEntryB_ma{1}));
          
          switch conf.maskType
           case 'caltech101',
            % Image A
            b = maskA.box_coord;
            c = maskA.obj_contour;
            c(1,:) = c(1,:) + b(3);
            c(2,:) = c(2,:) + b(1);
            in = inpolygon (framesA(1,:), framesA(2,:), c(1,:), c(2,:));
            framesA = framesA(:,in);
            descriptorsA = descriptorsA(:,in);
            
            % Image B
            b = maskB.box_coord;
            c = maskB.obj_contour;
            c(1,:) = c(1,:) + b(3);
            c(2,:) = c(2,:) + b(1);
            in = inpolygon (framesB(1,:), framesB(2,:), c(1,:), c(2,:));
            framesB = framesB(:,in);
            descriptorsB = descriptorsB(:,in);

           case 'r-caltech101',
            imgAinfo = imfinfo(fullfile(conf.imageRoot,fileEntryA_im{1}));      
            imgBinfo = imfinfo(fullfile(conf.imageRoot,fileEntryB_im{1}));      
            imgASz = [imgAinfo.Width imgAinfo.Height];
            imgBSz = [imgBinfo.Width imgBinfo.Height];
            % Must be [1,imgSize] both dims
            boxA = max(maskA.box_coord,ones(size(maskA.box_coord)));
            boxA = min(boxA, [imgASz(1)*ones(size(maskA.box_coord,1),1) ...
                              imgASz(2)*ones(size(maskA.box_coord,1),1)]); 
            boxB = max(maskB.box_coord,ones(size(maskB.box_coord)));
            boxB = min(boxB, [imgBSz(1)*ones(size(maskB.box_coord,1),1) ...
                              imgBSz(2)*ones(size(maskB.box_coord,1),1)]); 
            in = inpolygon (framesA(1,:), framesA(2,:), boxA(:,1), boxA(:,2));
            framesA = framesA(:,in);
            descriptorsA = descriptorsA(:,in);
            in = inpolygon (framesB(1,:), framesB(2,:), boxB(:,1), boxB(:,2));
            framesB = framesB(:,in);
            descriptorsB = descriptorsB(:,in);
            
           otherwise,
            error('Undefined mask type, don''t know how to process?');
          end;
          % Save frames (detected regions) and their descriptors
          frames = framesA; descriptors = descriptorsA;
          save(fullfile(fSavePath,fSaveFileA_new),'frames', 'descriptors');
          frames = framesB; descriptors = descriptorsB;
          save(fullfile(fSavePath,fSaveFileB_new),'frames', 'descriptors');
        end;
        timeElapsed = toc;
        timeTotalElapsed = timeTotalElapsed + timeElapsed;
      end
      fprintf(' ...Done!\n');
      mvpr_lclose(fhA_im);
      mvpr_lclose(fhB_im);
      mvpr_lclose(fhA_ma);
      mvpr_lclose(fhB_ma);
      fprintf('Done! Total time elapsed: %.2fs, %.2fs per file.  \n', ...
              timeTotalElapsed, timeTotalElapsed / imgTot);
  end;
else
  error(['The case that no masks available not implemented (only ' ...
         'Caltech-101 used so far)']);
end;
  
% 4.
% Read transformations (1.) and detecter areas (2.) and compute
% their matching scores
% configures
for detdescNo = 1:length(conf.detdesc)
  fprintf('\r4.%1d (%s) Matching detected regions\n', detdescNo,conf.detdesc{detdescNo});
  % Open image list file and create sift features
  timeTotalElapsed = 0;
  timeElapsed = 0;
  fhA_im = mvpr_lopen(fullfile(conf.filelistRoot,conf.imgsA),'read');
  fhB_im = mvpr_lopen(fullfile(conf.filelistRoot,conf.imgsB),'read');
  fhA_lm = mvpr_lopen(fullfile(conf.filelistRoot,conf.landmarksA),'read');
  fhB_lm = mvpr_lopen(fullfile(conf.filelistRoot,conf.landmarksB),'read');
  fhA_ma = mvpr_lopen(fullfile(conf.filelistRoot,conf.masksA),'read');
  fhB_ma = mvpr_lopen(fullfile(conf.filelistRoot,conf.masksB),'read');
  conf_curr = conf;
  % If save exists, then read, check if config is the same how
  % many files already computed
  matchSaveDir = [conf.tempSaveDir];
  matchSaveFile = ['SAVE_detectors_benchmark1_' ...
                   conf.conf_sets{conf.conf_setno}...
                   '_' conf.detdesc{detdescNo} '.mat'];
  if exist(fullfile(matchSaveDir,matchSaveFile),'file') && ~conf.enforceMatching
    load(fullfile(matchSaveDir,matchSaveFile),'res','conf','imgInd');
    startInd = imgInd;
    fprintf(['\n   ===> Using existing matching results in ' ...
             fullfile(matchSaveDir,matchSaveFile) ' \n']);
    conf = conf_curr; % not checked anyway so return old
  else
    startInd = 1;
  end;
  for imgInd = 1:imgTot
    timeETA = (timeTotalElapsed / (imgInd-1)) * (imgTot-imgInd);
    
    if 0 == mod(imgInd,30)
        fprintf('\r  Matching detectors of image pair: %4d/%4d (ETA: %.1fs)', ...
               imgInd,imgTot, timeETA);
    end
    tic;
      
    fileEntryA_im = mvpr_lread(fhA_im);
    fileEntryB_im = mvpr_lread(fhB_im);
    fileEntryA_lm = mvpr_lread(fhA_lm);
    fileEntryB_lm = mvpr_lread(fhB_lm);
    fileEntryA_ma = mvpr_lread(fhA_ma);
    fileEntryB_ma = mvpr_lread(fhB_ma);
    
    if startInd > imgInd
      continue;
    else
      [imgDirA imgFileA imgExtA] = fileparts(fileEntryA_im{1});
      [imgDirB imgFileB imgExtB] = fileparts(fileEntryB_im{1});
      
      % Name of the transformation matrix file
      [lmDir lmFile lmExt] = fileparts(fileEntryA_lm{1});
      [lmDirB lmFileB lmExtB] = fileparts(fileEntryB_lm{1});
      lmSavePath = fullfile(conf.tempSaveDir,conf.landmarksRoot,lmDir);
      lmSaveFile = [lmFile '_AtoB_' lmFileB '_H.txt'];
      
      % Names of feature files (from 2.)
      fSavePath = fullfile(conf.tempSaveDir,conf.imageRoot,imgDirA);
      fSaveFileA = [imgFileA '_' conf.detdesc{detdescNo} '_features_masked.mat'];
      fSaveFileB = [imgFileB '_' conf.detdesc{detdescNo} '_features_masked.mat'];
      
      % Load the groundtruth transformation
      H = load(fullfile(lmSavePath,lmSaveFile));
      
      % Load detected regions
      load(fullfile(fSavePath,fSaveFileA),'frames', 'descriptors');
      framesA = frames; descriptorsA = descriptors;
      load(fullfile(fSavePath,fSaveFileB),'frames', 'descriptors');
      framesB = frames; descriptorsB = descriptors;
      
      % Request image sizes
      imgAinfo = imfinfo(fullfile(conf.imageRoot,fileEntryA_im{1}));      
      imgBinfo = imfinfo(fullfile(conf.imageRoot,fileEntryB_im{1}));      
      imgASz = [imgAinfo.Width imgAinfo.Height];
      imgBSz = [imgBinfo.Width imgBinfo.Height];
      
      % THis function wraps functionality copied from Mikolajczyk
      [erro,repeat,corresp,match_score,matches, twi] = ...
          mikolajczyk_repeatability(framesA,framesB,descriptorsA,descriptorsB,H,imgASz,imgBSz,1);
      
      %if mvpr_lcountentries(fileA) <= 3 || mvpr_lcountentries(fileB) <= 3
      %  res.erro(imgInd)                = 0;
      %  res.repeatability(imgInd)       = 0;
      %  res.n_matches(imgInd)           = 0;
      %  res.n_correspondences(imgInd)   = 0;
      %  res.twi{imgInd}                 = 0;
      %else
      %conf.detectorOverlap = 50
      switch conf.detectorOverlap
       case 40,
        idx = 4;
       case 50,
        idx = 5;
       otherwise,
        error(['Unknown overlap criterion (bounded by those supported ' ...
               'in Mikolajczyk code']);
      end
      res.erro(imgInd)                = erro(idx);
      res.repeatability(imgInd)       = repeat(idx);
      res.n_correspondences(imgInd)   = corresp(idx);
      res.match_score(imgInd)         = match_score;
      res.n_matches(imgInd)           = matches;
      res.twi{imgInd}                 = twi; % AGB: I guess this is not needed for detector tests?
      
      %% AGB: Save
      %reserro = res.erro(imgInd);
      %resrepeatability = res.repeatability(imgInd);
      %resn_correspondences = res.n_correspondences(imgInd);
      %resmatch_score = res.match_score(imgInd);
      %resn_matches = res.n_matches(imgInd);
      %restwi = res.twi{imgInd};
      %save(dataFile,...
      %     'reserro', 'resrepeatability',...
      %     'resn_correspondences', 'resmatch_score',...
      %     'resn_matches', 'restwi');
      
      save(fullfile(matchSaveDir,matchSaveFile),'res','conf','imgInd');
      
      timeElapsed = toc;
      timeTotalElapsed = timeTotalElapsed + timeElapsed;
    end
  end;
  results(detdescNo) = res;
  fprintf(' ...Done!\n');
  mvpr_lclose(fhA_im);
  mvpr_lclose(fhB_im);
  mvpr_lclose(fhA_lm);
  mvpr_lclose(fhB_lm);
  mvpr_lclose(fhA_ma);
  mvpr_lclose(fhB_ma);
  fprintf('\rDone! Total time elapsed: %.2fs, %.2fs per file.  \n', ...
          timeTotalElapsed, timeTotalElapsed / imgTot);
end;

% 5.
% Plot the matching results
fprintf('\r5. All done, plotting and saving the results\n');
imgPairClass = {};
fhA_im = mvpr_lopen(fullfile(conf.filelistRoot,conf.imgsA),'read');
fhB_im = mvpr_lopen(fullfile(conf.filelistRoot,conf.imgsB),'read');
classnames = {};
classes = zeros(imgTot,1);
for imgInd = 1:imgTot
    fileEntryA_im = mvpr_lread(fhA_im);
    fileEntryB_im = mvpr_lread(fhB_im);
    classA = fileEntryA_im{2};
    classB = fileEntryB_im{2};
    if ~strcmp(classA,classB)
      warning(['Class mismatch for line ' num2str(imgInd) 'in '...
               fullfile(conf.filelistRoot,conf.imgsA)]);
    end;
    class_n = find(strcmp(classnames,classA));
    if isempty(class_n) % create new
        class_n = numel(classnames)+1;
        classnames{class_n,1} = classA;
    end
    classes(imgInd) = class_n;
end;

mikolajczyk_repeatability_plot(classes,classnames, conf.detdesc, results, ...
                               'dataset', conf.conf_sets{conf.conf_setno},...
                               'imgSaveDir', conf.tempSaveDir,...
                               'imgSaveFile', ['PLOT_detectors_benchmark1' conf.conf_sets{conf.conf_setno}],...
                               'saveAlsoFig', conf.saveAlsoFig,...
                               'legends', conf.detdescLabels,...
                               'debugLevel',conf.debugLevel);%,...
                               %'settings', confSet);
%end

disp 'All done!'
toc
