function saliency_img_SF=get_saliency_wCtr(srcImg)

    addpath(genpath('Funcs'));

    %% 1. Parameter Settings
    doFrameRemoving = true;
    useSP = true;           %You can set useSP = false to use regular grid for speed consideration

        %% Pre-Processing: Remove Image Frames

        if doFrameRemoving
            [noFrameImg, frameRecord] = removeframe_wCtr(srcImg, 'sobel');
            [h, w, chn] = size(noFrameImg);
        else
            noFrameImg = srcImg;
            [h, w, chn] = size(noFrameImg);
            frameRecord = [h, w, 1, h, 1, w];
        end

        %% Segment input rgb image into patches (SP/Grid)
        pixNumInSP = 600;                           %pixels in each superpixel
        spnumber = round( h * w / pixNumInSP );     %super-pixel number for current image

        if useSP
            [idxImg, adjcMatrix, pixelList] = SLIC_Split(noFrameImg, spnumber);
        else
            [idxImg, adjcMatrix, pixelList] = Grid_Split(noFrameImg, spnumber);        
        end
        %% Get super-pixel properties
        spNum = size(adjcMatrix, 1);
        meanRgbCol = GetMeanColor(noFrameImg, pixelList);
        meanLabCol = colorspace('Lab<-', double(meanRgbCol)/255);
        meanPos = GetNormedMeanPos(pixelList, h, w);
        bdIds = GetBndPatchIds(idxImg);
        colDistM = GetDistanceMatrix(meanLabCol);
        posDistM = GetDistanceMatrix(meanPos);
        [clipVal, geoSigma, neiSigma] = EstimateDynamicParas(adjcMatrix, colDistM);

        %% Saliency Optimization
        [bgProb, bdCon, bgWeight] = EstimateBgProb(colDistM, adjcMatrix, bdIds, clipVal, geoSigma);
        wCtr = CalWeightedContrast(colDistM, posDistM, bgProb);
        optwCtr = SaliencyOptimization(adjcMatrix, bdIds, colDistM, neiSigma, bgWeight, wCtr);

        saliency_img=SaveSaliencyMap(optwCtr, pixelList,frameRecord,true);
        
       %% Saliency Filter
        [cmbVal, contrast, distribution] = SaliencyFilter(colDistM, posDistM, meanPos);
    
        %smapName=fullfile(savepath_SF, strcat(noSuffixName, '_SF.png'));
        saliency_img_SF = SaveSaliencyMap(cmbVal, pixelList, frameRecord, true);  
        
       %% Geodesic Saliency
        geoDist = GeodesicSaliency(adjcMatrix, bdIds, colDistM, posDistM, clipVal);
    
        %smapName=fullfile(savepath_GS, strcat(noSuffixName, '_GS.png'));
        saliency_img_GS =SaveSaliencyMap(geoDist, pixelList, frameRecord, true);
        
        
end