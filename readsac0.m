
function [A,B,C,D]=readsac0(sacfl)    
    fid=fopen(sacfl,'r',  'ieee-le');
    A = fread( fid, [ 70 1 ], 'float32' );
    B = fread( fid, [ 40 1 ], 'int32' );
    C = char(fread ( fid, [ 1 192 ], 'char' ));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read SAC file data.
    %
    D = fread( fid, 'float32' );

    A(A==-12345.0) = NaN;
    B(B==-12345) = NaN;
    C = cellstr(reshape(C,8,24)');
    C(strmatch('-12345',C)) = {''};
    fclose(fid);
