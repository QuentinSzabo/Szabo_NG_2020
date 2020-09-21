%% FISH - DAPI 2D analysis

% - Diameters of chromatin nanodomains within FISH probes or within DAPI staining

% - Input files: folders containing FISH or DAPI channel-separated images
% (3D ROIs surrounding a single FISH locus or individual z-slices of a single nucleus) named '*C1.tif'

%% Paths and variables to define
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder = cd;
inputfolder = [folder '\input\']; % Contains input files
outputfolder = [folder '\output\'];
if ~exist(outputfolder, 'dir')
  mkdir(outputfolder);
end
close all
%%%%%%%%%%%%%%%%%%%%%%
f = 1; % 1 for FISH analysis
species = 1; % 1 for mouse; 2 for drosophila
%%%%%%%%%%%%%%%%%%%%%%

channel = 'C1'; % Channel to be analyzed
xypixelsize = 0.04; % xy pixel size in µm
zsize = 0.125; % z size in µm

if f == 1
    if species == 1
        minvolume = 200; % Minimum FISH volume (in voxels) in mouse
    elseif species == 2
        minvolume = 40; % Minimum FISH volume (in voxels) in drosophila
    end
else
    % Graphical scalebar parameters for DAPI
    if species == 1
        LnW = 4; scalebar = [21 145; 21 21];
    elseif species == 2
        LnW = 8; scalebar = [18 142; 18 18]; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output summary files summing data from all input folders

allwatershed = table; allwatershedname = 'all_watershed.csv'; % Watershed segmentation
allfull = table; allfullname = 'all_full.csv'; % Full segmentation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
list = dir(inputfolder);
for k = 1:length(list) - 2
%% Input/output folders and files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inputsubfolder = [inputfolder list(k+2).name];
outputsubfolder = [outputfolder list(k+2).name];
if ~exist(outputsubfolder, 'dir')
  mkdir(outputsubfolder);
end

filepatternC1=fullfile(inputsubfolder, ['*' channel '.tif']);
filesC1=dir(filepatternC1);
nfiles=length(filesC1);

% Output summary files
summarywatershed = table; summarywatershedname = 'summary_watershed.csv';
summaryfull = table; summaryfullname = 'summary_full.csv' ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
for i=1:nfiles
%% Image information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     inputnameC1=filesC1(i).name;
     inputfoldernameC1=fullfile(inputsubfolder, inputnameC1);

     imdata = imfinfo(inputfoldernameC1);
     widthC = [imdata.Width]; nxpixels = widthC(1);
     height = [imdata.Height]; nypixels = height(1);
     nzslices = numel(imdata); 
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%% Reading, filtering and segmentation of images 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If FISH, 3D segmentation and extraction of a single z-slice
if f == 1

     im3dC1 = zeros(nypixels, nxpixels, nzslices);
     for zsliceC1 = 1 : nzslices
            imsliceC1 = imread(inputfoldernameC1, zsliceC1);
               im3dC1(:,:,zsliceC1) = imsliceC1;
     end
     
     imfC1 = imgaussfilt3(im3dC1, .5); imf16C1 = uint16(imfC1);
     projC1 = max(imf16C1, [], 3);
     
     threshC1 = graythresh(imf16C1); maskC1 = imbinarize(imf16C1,threshC1); 
     projmaskC1 = max(maskC1, [], 3);

     % Extraction of 3D segmented object
     connectedcC1 = bwconncomp(maskC1);
     prestatsC1 = regionprops3(connectedcC1, 'Volume');
     idxC1 = find([prestatsC1.Volume] > minvolume);    
     mask2C1 = ismember(labelmatrix(connectedcC1), idxC1);
     mask2C1 = imclearborder(mask2C1);
     projmask2C1 = max(mask2C1, [], 3);
     statsC1 = regionprops3(mask2C1, 'Volume', 'Centroid', 'ConvexHull');
     
     % Only keeps images with 1 segmented object per channel
     centroidC1 = statsC1.Centroid;
     nobjcentroidC1 = size(centroidC1,1);
     if nobjcentroidC1 ~= 1
        continue
     end   
    
     % Extraction of a random 2D slice
	 maxz = round(max(statsC1.ConvexHull{1,1}(:,3)));
	 minz = round(min(statsC1.ConvexHull{1,1}(:,3)));
	 zrange = maxz - minz;
	 zROI = randi([minz minz+zrange],1);
	 im2DC1 = im3dC1(:,:,zROI);
	 im2DC1 = uint16(im2DC1);
	 mask2DC1 = mask2C1(:,:,zROI);
	 im2DC1(~mask2DC1) = 0;
      
else
% If DAPI

     im2DC1 = imread(inputfoldernameC1);    
     im2DC1 =  uint16(im2DC1);
     imf2DC1 = imgaussfilt(im2DC1, .5); imf162DC1 = uint16(imf2DC1);
     thresh2DC1 = graythresh(imf162DC1);
     mask2DC1 = imbinarize(imf162DC1,thresh2DC1); 
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Extraction of segmented objects - Full segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     connectedc2DC1 = bwconncomp(mask2DC1);
     prestats2DC1 = regionprops(connectedc2DC1, 'Area');
     idx2DC1 = find([prestats2DC1.Area] > 9);    
     mask22DC1 = ismember(labelmatrix(connectedc2DC1), idx2DC1);
     stats2DC1 = regionprops('table',mask22DC1, 'Area', 'Centroid','EquivDiameter');
     
     pixarea = xypixelsize^2;

     areafull = stats2DC1{:,1}.* pixarea;
     diameterfull = stats2DC1{:,3}.*xypixelsize;
     nsubobjfull = size(stats2DC1,1);
     id = repmat(i, nsubobjfull,1);

     statstablefull = table(id, areafull, diameterfull);
     summaryfull = [summaryfull; statstablefull];    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Watershed segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     im = im2DC1; im(~mask22DC1) = 0;
     I = mat2gray(im);
     imA = -I; imA(~mask22DC1) = Inf;
     A = watershed(imA); A(~mask22DC1) = 0;
     
     connectedcA = bwconncomp(A);
     prestatsA = regionprops(connectedcA, 'Area');
     idxA = find([prestatsA.Area] > 9); 
     maskA = ismember(labelmatrix(connectedcA), idxA);
     maskA = imclearborder(maskA);
     statsA = regionprops('table', maskA, 'Area', 'Centroid', 'EquivDiameter');
     
     areawatershed = statsA{:,1}.* pixarea;
     diameterwatershed = statsA{:,3}.*xypixelsize;
     nsubobjwatershed = size(statsA,1);
     id = repmat(i, nsubobjwatershed,1);

     statstable = table(id, areawatershed, diameterwatershed);
     summarywatershed = [summarywatershed; statstable];
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FISH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if f == 1
%%  Image max projection

     adjprojC1  = imadjust(mat2gray(projC1),[0.05; 1]);
     bound1 = bwboundaries(projmask2C1);
    
     close all           
     fig = figure('pos',[600 200 400 400]);
     set(gcf,'color','w'); set(gca,'xtick',[], 'ytick', []); ax = gca; outerpos = ax.OuterPosition; ti = ax.TightInset; 
     left = outerpos(1) + ti(1); bottom = outerpos(2) + ti(2); ax_width = outerpos(3) - ti(1) - ti(3);
     ax_height = outerpos(4) - ti(2) - ti(4); ax.Position = [left bottom ax_width ax_height]; fig.InvertHardcopy = 'off';
       
     imshow(adjprojC1);
     hold on
     for p = 1:length(bound1)
           boundaryC1 = bound1{p};
           plot(boundaryC1(:,2), boundaryC1(:,1), 'Color', 'y', 'LineWidth', 2.25)
     end 

     imname = [strrep(inputnameC1, [channel '.tif'], [channel '_maxproj_segmented']), '.tif']; 
     outputimname = fullfile(outputsubfolder, imname); 
     print(outputimname, '-dtiff', '-r100')
    
%% Image 2D full segmentation  

     adjC1 = im2DC1;
     bound1 = bwboundaries(mask22DC1);
     
     close all    
     fig = figure('pos',[600 200 400 400]);
     set(gcf,'color','w'); set(gca,'xtick',[], 'ytick', []); ax = gca; outerpos = ax.OuterPosition; ti = ax.TightInset; 
     left = outerpos(1) + ti(1); bottom = outerpos(2) + ti(2); ax_width = outerpos(3) - ti(1) - ti(3);
     ax_height = outerpos(4) - ti(2) - ti(4); ax.Position = [left bottom ax_width ax_height]; fig.InvertHardcopy = 'off';
      
     imshow(adjC1, []);
     hold on
     for p = 1:length(bound1)
           boundaryC1 = bound1{p};
           plot(boundaryC1(:,2), boundaryC1(:,1), 'Color', 'y', 'LineWidth', 2.25)
     end 

     imname = [strrep(inputnameC1, [channel '.tif'], [channel '_slice_segmented']), '.tif']; 
     outputimname = fullfile(outputsubfolder, imname); 
     print(outputimname, '-dtiff', '-r100')

%% Image 2D watershed segmentation

     close all
     fig = figure('pos',[600 200 400 400]);
     set(gcf,'color','w'); set(gca,'xtick',[], 'ytick', []); ax = gca; outerpos = ax.OuterPosition; ti = ax.TightInset; 
     left = outerpos(1) + ti(1); bottom = outerpos(2) + ti(2); ax_width = outerpos(3) - ti(1) - ti(3);
     ax_height = outerpos(4) - ti(2) - ti(4); ax.Position = [left bottom ax_width ax_height]; fig.InvertHardcopy = 'off';
      
     imshow(adjC1, []);
     hold on
     for j = 1:size(statsA,1)
             connectedcA1 = bwconncomp(maskA);
             idxAj = find([statsA.Area] == statsA{j,1}); 
             maskAj = ismember(labelmatrix(connectedcA1), idxAj);
             boundj = bwboundaries(maskAj);
         for p = 1:length(boundj)
                 boundaryC1j = boundj{p};
                 plot(boundaryC1j(:,2), boundaryC1j(:,1), 'Color', 'y', 'LineWidth', 2.25)
         end        
     end

     imname = [strrep(inputnameC1, [channel '.tif'], 'slice_watershed'), '.tif']; 
     outputimname = fullfile(outputsubfolder, imname); 
     print(outputimname, '-dtiff', '-r100')
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
%% DAPI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Image
                  
     if species == 1
     adjC1  = imadjust(mat2gray(imf162DC1),[.1; .9]);
     elseif species == 2
     adjC1  = imadjust(mat2gray(imf162DC1),[.1; .85]);  
     end
     close all
     fig = figure('pos',[600 200 400 400]);
     set(gcf,'color','w'); set(gca,'xtick',[], 'ytick', []); ax = gca; outerpos = ax.OuterPosition; ti = ax.TightInset; 
     left = outerpos(1) + ti(1); bottom = outerpos(2) + ti(2); ax_width = outerpos(3) - ti(1) - ti(3);
     ax_height = outerpos(4) - ti(2) - ti(4); ax.Position = [left bottom ax_width ax_height]; fig.InvertHardcopy = 'off';

     imshow(adjC1);
     hold on
     line(scalebar(1,:), scalebar(2,:), 'Color', 'w', 'LineWidth', LnW)
    
     imname = [strrep(inputnameC1, '.tif','_slice'), '.tif']; 
     outputimname = fullfile(outputsubfolder, imname); 
     print(outputimname, '-dtiff', '-r600')
    
 %% Image full segmentation   
          
     bound1 = bwboundaries(mask22DC1);
    
     close all
     fig = figure('pos',[600 200 400 400]);
     set(gcf,'color','w'); set(gca,'xtick',[], 'ytick', []); ax = gca; outerpos = ax.OuterPosition; ti = ax.TightInset; 
     left = outerpos(1) + ti(1); bottom = outerpos(2) + ti(2); ax_width = outerpos(3) - ti(1) - ti(3);
     ax_height = outerpos(4) - ti(2) - ti(4); ax.Position = [left bottom ax_width ax_height]; fig.InvertHardcopy = 'off';
    
     imshow(adjC1);
     hold on
         for p = 1:length(bound1)
                 boundaryC1 = bound1{p};
                 plot(boundaryC1(:,2), boundaryC1(:,1), 'Color', 'y', 'LineWidth', .5)
         end 

     imname = [strrep(inputnameC1, '.tif','_slice_segmented'), '.tif']; 
     outputimname = fullfile(outputsubfolder, imname); 
     print(outputimname, '-dtiff', '-r300')
     hold off

    
%% Image 2D watershed segmentation

     close all
     fig = figure('pos',[600 200 400 400]);
     set(gcf,'color','w'); set(gca,'xtick',[], 'ytick', []); ax = gca; outerpos = ax.OuterPosition; ti = ax.TightInset; 
     left = outerpos(1) + ti(1); bottom = outerpos(2) + ti(2);
     ax_width = outerpos(3) - ti(1) - ti(3); ax_height = outerpos(4) - ti(2) - ti(4);
     ax.Position = [left bottom ax_width ax_height]; fig.InvertHardcopy = 'off';
      
     imshow(adjC1);
     hold on
     for j = 1:size(statsA,1)
            connectedcA1 = bwconncomp(maskA);
            idxAj = find([statsA.Area] == statsA{j,1}); 
            maskAj = ismember(labelmatrix(connectedcA1), idxAj);
            boundj = bwboundaries(maskAj);

            for p = 1:length(boundj)
                 boundaryC1j = boundj{p};
                 plot(boundaryC1j(:,2), boundaryC1j(:,1), 'Color', 'y', 'LineWidth', .5)
            end 
     end

     imname = [strrep(inputnameC1, '.tif','_slice_watershed'), '.tif']; 
     outputimname = fullfile(outputsubfolder, imname); 
     print(outputimname, '-dtiff', '-r600')
     hold off

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     close all;

end

%% Summary files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    writetable(summaryfull, fullfile(outputsubfolder, summaryfullname));
    writetable(summarywatershed, fullfile(outputsubfolder, summarywatershedname));
       
    allwatershed = [allwatershed; summarywatershed];
    allfull = [allfull; summaryfull];
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%% Figure with all images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figallimages = figure;
    imagespattern=fullfile(outputsubfolder, '*.tif');
    imagesfiles=dir(imagespattern);
    fileNames = ({imagesfiles.name});
    inputfileNames=fullfile(outputsubfolder, fileNames);
    allimages = montage(inputfileNames);
    
    allimname = 'all_images'; 
    outputallimname = fullfile(outputsubfolder, allimname); 
    print(outputallimname, '-dtiff', '-r800')
    print(outputallimname, '-dpdf', '-fillpage')
    
    close all
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%% Files summing data from all input folders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    writetable(allwatershed, fullfile(outputfolder, allwatershedname));
    writetable(allfull, fullfile(outputfolder, allfullname));

