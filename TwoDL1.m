%elements in the matrix
A=1
B=2
C=3


% null element


%opeartions
% ih iv dh dv
ih=1;
iv=2;
dh=3;
dv=4;

%cost vector 
costT=[1 1  1  1];

% two example matrices
MA=[
        C , A , A , A  
        A , A , A , A  
        A , B , A , A    
];

MB=[
        C , A , A , A  
        A , A , A , A  
        A , B , A , A    
];


%matching 1

vv=[
1 1;
1 1;
1 1;
1 1;
1 1;
];

tv=[
ih;
ih;
ih;
iv;
iv;
];

insv=[
5
6
7
8
9
];

MB1=transformation(MA,vv, tv,insv);

MB1
MB
%-----------------------

% 
% vv=[
%     1 1;
%     2 1;
%     3 1;
%     2 2;
% ]
% 
% 
% 
% 8
% 
% 
% iv=1


