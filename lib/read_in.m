function [timeVec, trueStrain, expStress] = read_in(filename,cyclic)

if cyclic == 0
    fileID = fopen(filename);
    data = textscan(fileID,'%f %f %f %f','Delimiter',',','HeaderLines',1);
    timeVec = data{1};
    trueStrain = data{2};
    expStress = data{3};
else
    Area = 2.0 * 9.0;
    fileID = fopen(filename);
    data = textscan(fileID,'%f %f %f','Delimiter',',','HeaderLines',3);
    timeVec = data{1};
    forces = data{2};
    displ = data{3};
    
    expStress = forces./Area;
    trueStrain = displ
end

end

