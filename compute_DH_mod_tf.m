function [output] = compute_DH_mod_tf(ai_1,alphai_1,di,ttai)
%COMPUTE_TRANSFORM computes symbolic modified DH transform
output = [ cos(ttai)                      , -sin(ttai)                      , 0                   , ai_1;
                sin(ttai)*cos(alphai_1) , cos(ttai)*cos(alphai_1) , -sin(alphai_1) , -sin(alphai_1)*di;
                sin(ttai)*sin(alphai_1)  , cos(ttai)*sin(alphai_1) , cos(alphai_1) , cos(alphai_1)*di;
                0                                , 0                                , 0                   , 1                       ];

end

