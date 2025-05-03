function [sourceTile, detTile] = getTileNumbers(padname)

    switch padname
        case {'GA00274'}
        
            sourceTile{1} = [34,35,36];
            sourceTile{2} = [31,32,33];
            sourceTile{3} = [28,29,30];
            sourceTile{4} = [19,20,21];
            sourceTile{5} = [22,23,24];
            sourceTile{6} = [25,26,27];
            sourceTile{7} = [1,2,3];
            sourceTile{8} = [4,5,6];
            sourceTile{9} = [7,8,9];
            sourceTile{10} = [16,17,18];
            sourceTile{11} = [13,14,15];
            sourceTile{12} = [10,11,12];
            
            detTile{1} = [45,46,47,48];
            detTile{2} = [41,42,43,44];
            detTile{3} = [37,38,39,40];
            detTile{4} = [25,26,27,28];
            detTile{5} = [29,30,31,32];
            detTile{6} = [33,34,35,36];
            detTile{7} = [1,2,3,4];
            detTile{8} = [5,6,7,8];
            detTile{9} = [9,10,11,12];
            detTile{10} = [21,22,23,24];
            detTile{11} = [17,18,19,20];
            detTile{12} = [13,14,15,16];
        
        case {'GA00438_NF', 'GA00440_NF'}
        
            sourceTile{1} = [1,2,3];
            sourceTile{2} = [4,5,6];
            sourceTile{3} = [7,8,9];
            sourceTile{4} = [34,35,36];
            sourceTile{5} = [31,32,33];
            sourceTile{6} = [28,29,30];
            sourceTile{7} = [16,17,18];
            sourceTile{8} = [13,14,15];
            sourceTile{9} = [10,11,12];
            sourceTile{10} = [19,20,21];
            sourceTile{11} = [22,23,24];
            sourceTile{12} = [25,26,27];
            
            detTile{1} = [1,2,3,4];
            detTile{2} = [5,6,7,8];
            detTile{3} = [9,10,11,12];
            detTile{4} = [45,46,47,48];
            detTile{5} = [41,42,43,44];
            detTile{6} = [37,38,39,40];
            detTile{7} = [21,22,23,24];
            detTile{8} = [17,18,19,20];
            detTile{9} = [13,14,15,16];
            detTile{10} = [25,26,27,28];
            detTile{11} = [29,30,31,32];
            detTile{12} = [33,34,35,36];
    
        case {'GA00440', 'GA00438', 'GA00439'}
    
            sourceTile{1} = [1,2,3];
            sourceTile{2} = [4,5,6];
            sourceTile{3} = [7,8,9];
            sourceTile{4} = [43,44,45];
            sourceTile{5} = [40,41,42];
            sourceTile{6} = [37,38,39];
            sourceTile{7} = [22,23,24];
            sourceTile{8} = [19,20,21];
            sourceTile{9} = [16,17,18];
            sourceTile{10} = [25,26,27];
            sourceTile{11} = [28,29,30];
            sourceTile{12} = [31,32,33];
            sourceTile{13} = [10,11,12];
            sourceTile{14} = [34,35,36];
            sourceTile{15} = [13,14,15];
            
            detTile{1} = [1,2,3,4];
            detTile{2} = [5,6,7,8];
            detTile{3} = [9,10,11,12];
            detTile{4} = [57,58,59,60];
            detTile{5} = [53,54,55,56];
            detTile{6} = [49,50,51,52];
            detTile{7} = [29,30,31,32];
            detTile{8} = [25,26,27,28];
            detTile{9} = [21,22,23,24];
            detTile{10} = [33,34,35,36];
            detTile{11} = [37,38,39,40];
            detTile{12} = [41,42,43,44];
            detTile{13} = [13,14,15,16];
            detTile{14} = [45,46,47,48];
            detTile{15} = [17,18,19,20];
    
        case {'GA00370', 'GA00369', 'GA00351'}
    
            sourceTile{1} = [1,2,3];
            sourceTile{2} = [4,5,6];
            sourceTile{3} = [7,8,9];
            sourceTile{4} = [34,35,36];
            sourceTile{5} = [31,32,33];
            sourceTile{6} = [19,20,21];
            sourceTile{7} = [16,17,18];
            sourceTile{8} = [13,14,15];
            sourceTile{9} = [22,23,24];
            sourceTile{10} = [25,26,27];
            sourceTile{11} = [10,11,12];
            sourceTile{12} = [28,29,30];
            
            detTile{1} = [1,2,3,4];
            detTile{2} = [5,6,7,8];
            detTile{3} = [9,10,11,12];
            detTile{4} = [45,46,47,48];
            detTile{5} = [41,42,43,44];
            detTile{6} = [25,26,27,28];
            detTile{7} = [21,22,23,24];
            detTile{8} = [17,18,19,20];
            detTile{9} = [29,30,31,32];
            detTile{10} = [33,34,35,36];
            detTile{11} = [13,14,15,16];
            detTile{12} = [37,38,39,40];
    end
end
