% folderName = 'C:\Users\CHENG-BANG\Dropbox\ImageData\Healthcare\BreastHistopathologyImages\IDCimage';
% folder_to_output = 'C:\Users\CHENG-BANG\Dropbox\ImageData\Healthcare\BreastHistopathologyImages\HRQA';
% pc = 0.5;
function ConvertAllpics(folderName, folder_to_output, imgType, l, v, Ub, Lb, pc)
    % find all the pictures in the current folder
    folder_to_search = folderName;
    filetype = fullfile(folder_to_search, imgType);
    filelist = dir(filetype);
    % Convert all pictures into SHRQA
    for Readinfile = 1 : length(filelist)
      File2Read = fullfile(folder_to_search, filelist(Readinfile).name);
      img_tmp = imread(File2Read); 
      if size(img_tmp,1) == size(img_tmp,2)
          ST_tmp = SpatialTraverse(img_tmp, l, v);
          % 
          Chk_ST_tmp = sum(ST_tmp, 2);
          if sum(Chk_ST_tmp(:)>=700) < pc*(2^(2*l));
            [s, ~] = SymbG(ST_tmp, Ub, Lb);
            [IdxM, HRR, HMean, HVar, HSkew, HKurt, HEnt, HGini] = HRQA(s, size(Ub,1),1);
            Tmp_HRQA = horzcat(IdxM, HRR, HMean, HVar, HSkew, HKurt, HEnt, HGini);
            filename = filelist(Readinfile).name;filename = filename(1:end-4);
            File2Write = fullfile(folder_to_output, strcat(filename,'.csv'));
            csvwrite(File2Write,Tmp_HRQA);
          end
      end
    end
    % find all the folders in the current folder
    Filelist = dir(folder_to_search);
    folderFlags = [Filelist.isdir];
    subFolders = Filelist(folderFlags);
    for fd = 1:length(subFolders)
        if subFolders(fd).name ~= '.' | subFolders(fd).name ~= '.'
            folderName1 = strcat(folderName,'\',subFolders(fd).name);
            ConvertAllpics(folderName1, folder_to_output, imgType, l, v, Ub, Lb, pc);
        end
    end
end