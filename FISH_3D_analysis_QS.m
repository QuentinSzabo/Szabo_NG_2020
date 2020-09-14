%% FISH 3D analysis QS

% - Intermingling between FISH probes if two channels are analyzed
% - FISH probe structure if one channel is analyzed

% - Input files: folders containing FISH channel-separated images (3D ROIs surrounding FISH loci)
% named '*C1.tif' and/or '*C2.tif'

%% Paths and variables to define
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder = cd;
inputfolder = [folder '\input_1C_C2\']; % Contains input files
outputfolder = [folder '\output_1C_C2\'];
if ~exist(outputfolder, 'dir')
  mkdir(outputfolder);
end
close all
%%%%%%%%%%%%%%%%%%%%%%
typeM = 1; % Type of microscopy: 1 for 3D-SIM; 0 for conv WF
nC = 1; % Number of channels to be analyzed: 2 for probes intermingling; 1 for probe structure
if nC == 1
    channel = 'C2'; % Channel to be analyzed if one color analysis
%%%%%%%%%%%%%%%%%%%%%% 

    minsubdomainvolume = 36; % Minimum FISH watershed segmented volume in voxels
else
    channel = 'C1';
end
if typeM == 1
    xypixelsize = 0.04; % xy pixel size in µm for 3D-SIM
    xysize = 55; % xy image size for averaged image for 3D-SIM
else
    xypixelsize = 0.08; % xy pixel size in µm for conv WF
    xysize = 28; % xy image size for averaged image for conv WF
end
zsize = 0.125; % z size in µm
minvolume = 200; % Minimum FISH segmented volume in voxels

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
if nC == 2
filepatternC2=fullfile(inputsubfolder, '*C2.tif');
filesC2=dir(filepatternC2);
end
nfiles=length(filesC1);

% Output summary files
if nC == 1
volumesummary = table; volumesummaryname = 'summary_volumes.csv' ;
principalaxislengthsummary = table; principalaxislengthsummaryname = 'summary_axes.csv';
sphericitysummary = table; sphericitysummaryname = 'summary_sphericity_.csv';
subdomainvolumesummary = table; subdomainvolumesummaryname = 'summary_subdomain_volume.csv';
nsubdomainssummary = table; nsubdomainssummaryname = 'summary_n_subdomains.csv';
end

if nC == 2
overlapfractionsummary = table; overlapfractionsummaryname = 'summary_overlap_fraction.csv';
distancecentroidsummary = table;distancecentroidsummaryname = 'summary_distance_centroids.csv';
mergeimC1 = zeros(xysize, xysize); mergeimC2 = zeros(xysize, xysize);
alldistmatrixoverlap = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
for i=1:nfiles
%% Image information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     inputnameC1=filesC1(i).name;
     inputfoldernameC1=fullfile(inputsubfolder, inputnameC1);
     if nC == 2
     inputnameC2=filesC2(i).name;
     inputfoldernameC2=fullfile(inputsubfolder, inputnameC2);
     end
   
     imdata = imfinfo(inputfoldernameC1);
     widthC = [imdata.Width]; nxpixels = widthC(1);
     height = [imdata.Height]; nypixels = height(1);
     nzslices = numel(imdata); 
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%% Reading, filtering and segmentation of images 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     im3dC1 = zeros(nypixels, nxpixels, nzslices);
     for zsliceC1 = 1 : nzslices
            imsliceC1 = imread(inputfoldernameC1, zsliceC1);
               im3dC1(:,:,zsliceC1) = imsliceC1;
     end
     if nC == 2
     im3dC2 = zeros(nypixels, nxpixels, nzslices);
     for zsliceC2 = 1 : nzslices
            imsliceC2 = imread(inputfoldernameC2, zsliceC2);
               im3dC2(:,:,zsliceC2) = imsliceC2;
     end
     end
     
     imfC1 = imgaussfilt3(im3dC1, .5); imf16C1 = uint16(imfC1);
     projC1 = max(imf16C1, [], 3);
     if nC == 2
     imfC2 = imgaussfilt3(im3dC2, .5); imf16C2 = uint16(imfC2); 
     projC2 = max(imf16C2, [], 3);
     end
          
     threshC1 = graythresh(imf16C1); maskC1 = imbinarize(imf16C1,threshC1); 
     projmaskC1 = max(maskC1, [], 3);
     if nC == 2
     threshC2 = graythresh(imf16C2); maskC2 = imbinarize(imf16C2,threshC2); 
     projmaskC2 = max(maskC2, [], 3);
     end
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Extraction of segmented objects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     connectedcC1 = bwconncomp(maskC1);
     prestatsC1 = regionprops3(connectedcC1, 'Volume');
     idxC1 = find([prestatsC1.Volume] > minvolume);    
     mask2C1 = ismember(labelmatrix(connectedcC1), idxC1);
     mask2C1 = imclearborder(mask2C1);
     projmask2C1 = max(mask2C1, [], 3);
     statsC1 = regionprops3(mask2C1, 'Volume', 'Centroid','PrincipalAxisLength','SurfaceArea');
     
     if nC == 2
     connectedcC2 = bwconncomp(maskC2);
     prestatsC2 = regionprops3(connectedcC2, 'Volume');
     idxC2 = find([prestatsC2.Volume] > minvolume);    
     mask2C2 = ismember(labelmatrix(connectedcC2), idxC2);
     mask2C2 = imclearborder(mask2C2);
     projmask2C2 = max(mask2C2, [], 3);
     statsC2 = regionprops3(mask2C2, 'Centroid');
     end
     
     % Only keeps images with 1 segmented object per channel
     if size(statsC1.Centroid,1) ~= 1
        continue
     end
     if nC == 2
        if size(statsC2.Centroid,1) ~= 1
            continue
        end
     end
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Intermingling between C1 and C2 if 2 channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     if nC == 2
     maskmerge = (mask2C1 + mask2C2);
     maskoverlap = maskmerge > 1;
     projmaskoverlap = max(maskoverlap, [], 3);
     
     % Overlap fraction
     overlapfraction = jaccard(mask2C1, mask2C2);
     overlapfractionsummary = [overlapfractionsummary; table(overlapfraction)];
     
     % 3D distance between centroids in µm
     xdist = (statsC2.Centroid(1) - statsC1.Centroid(1)) * xypixelsize;
     ydist = (statsC2.Centroid(2) - statsC1.Centroid(2)) * xypixelsize;
     zdist = (statsC2.Centroid(3) - statsC1.Centroid(3)) * zsize;
     distance = sqrt(xdist^2 + ydist^2 + zdist^2);
     distancecentroidsummary = [distancecentroidsummary; table(distance)];
     
     % Matrix overlap
     statsprojmask2C1 = regionprops('table', projmask2C1, 'Area', 'Centroid'); 
     xshiftC1 = round(nxpixels/2 - statsprojmask2C1.Centroid(1));
     yshiftC1 = round(nypixels/2 - statsprojmask2C1.Centroid(2));
     shiftmatrixC1 = circshift(projmask2C1,[yshiftC1 xshiftC1]);
     mergeimC1 = imadd(mergeimC1, double(shiftmatrixC1));
     
     statsprojmask2C2 = regionprops('table', projmask2C2, 'Area', 'Centroid'); 
     xshiftC2 = round(nxpixels/2 - statsprojmask2C2.Centroid(1));
     yshiftC2 = round(nypixels/2 - statsprojmask2C2.Centroid(2));
     shiftmatrixC2 = circshift(projmask2C2,[yshiftC2 xshiftC2]);
     mergeimC2 = imadd(mergeimC2, double(shiftmatrixC2));
     
     dist2D = sqrt((statsprojmask2C1.Centroid(1)-statsprojmask2C2.Centroid(1))^2 + (statsprojmask2C1.Centroid(2)-statsprojmask2C2.Centroid(2))^2);
     alldistmatrixoverlap = [alldistmatrixoverlap; dist2D];
     end
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%% Probe structure analysis if 1 channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     if nC == 1

     voxelvol = xypixelsize * xypixelsize * zsize;
     
     % Volume in µm3
     volume = statsC1.Volume(1) * voxelvol;
     volumesummary = [volumesummary; table(volume)];
     
     % PrincipalAxisLength in µm
     xdistC1 = statsC1.PrincipalAxisLength(1) * xypixelsize;
     ydistC1 = statsC1.PrincipalAxisLength(2) * xypixelsize;
     zdistC1 = statsC1.PrincipalAxisLength(3) * zsize;
     principalaxislength = sqrt(xdistC1^2 + ydistC1^2 + zdistC1^2);
     principalaxislengthsummary = [principalaxislengthsummary; table(principalaxislength)];

     % Sphericity
     sphericity = ((pi^(1/3))*((6*statsC1.Volume)^(2/3)))/(statsC1.SurfaceArea);
     sphericitysummary = [sphericitysummary; table(sphericity)];
  
  
%% Watershed segmentation

     im = im3dC1; im(~mask2C1) = 0; 
     I = mat2gray(im); 
     imA = -I; imA(~mask2C1) = Inf;
     A = watershed(imA); A(~mask2C1) = 0;
     
     connectedcA = bwconncomp(A);
     prestatsA = regionprops3(connectedcA, 'Volume');
     idxA = find([prestatsA.Volume] > minsubdomainvolume); 
     maskA = ismember(labelmatrix(connectedcA), idxA);
     statsA = regionprops3(maskA, 'Volume', 'Centroid','ConvexHull');
     subdomaincentroid = statsA.Centroid;
     
     % Subdomain volumes in µm3
     subdomainvolume = statsA{:,1}.* voxelvol;
     subdomainvolumesummary = [subdomainvolumesummary; table(subdomainvolume)];
     
     % Number of of subdomains
     nsubdomains = size(statsA,1);
     nsubdomainssummary = [nsubdomainssummary; table(nsubdomains)];

     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Image with segmentation borders

       adjprojC1  = imadjust(mat2gray(projC1),[.05; 1]);
       if nC == 2
       adjprojC2  = imadjust(mat2gray(projC2),[.05 ; 1]);
       overlayC1C2 = imfuse(adjprojC2, adjprojC1);
       end
   
       bound1 = bwboundaries(projmask2C1);
       if nC == 2       
       bound2 = bwboundaries(projmask2C2);
       bound3 = bwboundaries(projmaskoverlap);
       of = ['OF =' ' ' num2str(overlapfraction,3)];
       end
       
      close all      
      fig = figure('pos',[600 200 400 400]);
      set(gcf,'color','w'); set(gca,'xtick',[], 'ytick', []); ax = gca; outerpos = ax.OuterPosition; ti = ax.TightInset; 
      left = outerpos(1) + ti(1); bottom = outerpos(2) + ti(2); ax_width = outerpos(3) - ti(1) - ti(3);
      ax_height = outerpos(4) - ti(2) - ti(4); ax.Position = [left bottom ax_width ax_height]; fig.InvertHardcopy = 'off';
       
       if nC == 1
       imshow(adjprojC1);
       colb = 'y';
       else
       imshow(overlayC1C2);
       colb = 'magenta';
       end
         hold on
         for p = 1:length(bound1)
                 boundaryC1 = bound1{p};
                 plot(boundaryC1(:,2), boundaryC1(:,1), 'Color', colb, 'LineWidth', 2.25)
         end 
         if nC == 2
         for q = 1:length(bound2)
                 boundaryC2 = bound2{q};
                 plot(boundaryC2(:,2), boundaryC2(:,1), 'Color', 'green', 'LineWidth', 2.25)
         end   
         for r = 1:length(bound3)
                 boundaryC3 = bound3{r};
                 plot(boundaryC3(:,2), boundaryC3(:,1), 'Color', 'white', 'LineWidth', 3.25)
         end  
         text(5, 9, of, 'Color', 'white', 'FontSize', 24)
         end
         
       if nC == 2
       imname = [strrep(inputnameC1, ['_' channel '.tif'], '_C1_C2_segmented'), '.tif']; 
       else
       imname = [strrep(inputnameC1, ['_' channel '.tif'], ['_' channel '_segmented']), '.tif']; 
       end
       outputimname = fullfile(outputsubfolder, imname); 
       print(outputimname, '-dtiff', '-r300')
    
%% Image 3D watershed segmentation if 1 channel

    if nC == 1
      close all
      fig = figure('pos',[600 200 400 400]);
      set(gcf,'color','w'); set(gca,'xtick',[], 'ytick', []); ax = gca; outerpos = ax.OuterPosition; ti = ax.TightInset; 
      left = outerpos(1) + ti(1); bottom = outerpos(2) + ti(2); ax_width = outerpos(3) - ti(1) - ti(3);
      ax_height = outerpos(4) - ti(2) - ti(4); ax.Position = [left bottom ax_width ax_height]; fig.InvertHardcopy = 'off';
      set(gca,'Ydir','reverse');
      hold on

      i1 = patch(isosurface(maskA));
      daspect([3.125 3.125 1])
      xlim([0 nxpixels]); ylim([0 nypixels]); zlim([0 nzslices])
      i1.FaceColor = [0.5 0.75 1];
      i1.EdgeColor = 'none';
      i1.FaceAlpha = .25;
      camlight 
      lighting gouraud
    
        for j = 1:nsubdomains
            scatter3(subdomaincentroid(j,1),subdomaincentroid(j,2),subdomaincentroid(j,3), 100,'filled')
        end

      imname = [strrep(inputnameC1, [channel '.tif'], [channel '_subdomains_3D_segmented']), '.tif']; 
      outputimname = fullfile(outputsubfolder, imname); 
      print(outputimname, '-dtiff', '-r300')
      
    end      
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
      close all

%% Summary files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     if nC == 1
     writetable(volumesummary, fullfile(outputsubfolder, volumesummaryname));
     writetable(principalaxislengthsummary, fullfile(outputsubfolder, principalaxislengthsummaryname));
     writetable(sphericitysummary, fullfile(outputsubfolder, sphericitysummaryname));
     writetable(subdomainvolumesummary, fullfile(outputsubfolder, subdomainvolumesummaryname));
     writetable(nsubdomainssummary, fullfile(outputsubfolder, nsubdomainssummaryname));
     end  
     
     if nC ==2
     writetable(overlapfractionsummary, fullfile(outputsubfolder, overlapfractionsummaryname));
     writetable(distancecentroidsummary, fullfile(outputsubfolder, distancecentroidsummaryname)); 
     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%% Averaged image if 2 channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if nC ==2
     shiftoverlap = round(mean(alldistmatrixoverlap)/2);

     matrixC1 = circshift(mergeimC1,shiftoverlap);    
     matrixC2 = circshift(mergeimC2,-shiftoverlap);
     
     overlay = imfuse(matrixC2, matrixC1);
     imagesc(overlay)
     overlayname = [('averaged_overlay'), '.tif']; 
     outputoverlayname = fullfile(outputsubfolder, overlayname); 
     imwrite(overlay, outputoverlayname);
     close all
     end
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%% Figure with all images     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     fig = figure;
     set(gcf,'color','w');     

     imagespattern=fullfile(outputsubfolder, '*segmented.tif');

     imagesfiles=dir(imagespattern);
     fileNames = ({imagesfiles.name});
     inputfileNames = fullfile(outputsubfolder, fileNames);
     allimages = montage(inputfileNames);
 
     allimname = 'all_images'; 
     outputallimname = fullfile(outputsubfolder, allimname); 
     print(outputallimname, '-dtiff', '-r800')
     print(outputallimname, '-dpdf', '-fillpage')
    
     close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if  typeM == 1
      microscopy = "3D-SIM";
  else
      microscopy = "Conventional";
  end

  if nC == 2
      parameters = table(microscopy, nC, xypixelsize, zsize, minvolume);
      writetable(parameters, fullfile(outputsubfolder, 'parameters.txt'));
  else
     parameters = table(microscopy, nC, string(channel), xypixelsize, zsize, minvolume, minsubdomainvolume);
     parameters.Properties.VariableNames(3) = "channel";
     writetable(parameters, fullfile(outputsubfolder, 'parameters.txt'));
  end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

