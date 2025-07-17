function [Nodes, Elements] = Readmesh( fname )

fid = fopen(fname, 'rt');  
S = textscan(fid, '%s', 'Delimiter', '\n'); 
S = S{1};

idxS = strfind(S, 'Coordinates');             
idx1 = find(not(cellfun(@isempty, idxS)));  

idxS = strfind(S, '% Elements (triangles)');          
idx2 = find(not(cellfun(@isempty, idxS)));    

idxS = strfind(S, '% Elements (triangles)');      
idx3 = find(not(cellfun(@isempty, idxS)));    

idxS = strfind(S, 'END');
idx4 = find(not(cellfun(@isempty, idxS)));

pause(0.01)
Nodes = S(idx1(1)+1:idx2(1)-1);  
Nodes = cell2mat(cellfun(@str2num, Nodes, 'UniformOutput', false)); 
Elements = S(idx3(1)+1:idx4(1)-1);
Elements = cell2mat(cellfun(@str2num, Elements, 'UniformOutput', false));

end