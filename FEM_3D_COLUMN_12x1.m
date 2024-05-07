%% Matlab mesh
%% FEM_3D_COLUMN, Created by Gmsh
%% ASCII
clear msh;
msh.nbNod = 48;
msh.POS = [
0 0 300;
0 0 0;
0 25 300;
0 25 0;
25 0 300;
25 0 0;
25 25 300;
25 25 0;
0 0 27.27272727272727;
0 0 54.54545454545455;
0 0 81.81818181818181;
0 0 109.0909090909091;
0 0 136.3636363636364;
0 0 163.6363636363636;
0 0 190.9090909090909;
0 0 218.1818181818182;
0 0 245.4545454545455;
0 0 272.7272727272727;
0 25 27.27272727272727;
0 25 54.54545454545455;
0 25 81.81818181818181;
0 25 109.0909090909091;
0 25 136.3636363636364;
0 25 163.6363636363636;
0 25 190.9090909090909;
0 25 218.1818181818182;
0 25 245.4545454545455;
0 25 272.7272727272727;
25 0 27.27272727272727;
25 0 54.54545454545455;
25 0 81.81818181818181;
25 0 109.0909090909091;
25 0 136.3636363636364;
25 0 163.6363636363636;
25 0 190.9090909090909;
25 0 218.1818181818182;
25 0 245.4545454545455;
25 0 272.7272727272727;
25 25 27.27272727272727;
25 25 54.54545454545455;
25 25 81.81818181818181;
25 25 109.0909090909091;
25 25 136.3636363636364;
25 25 163.6363636363636;
25 25 190.9090909090909;
25 25 218.1818181818182;
25 25 245.4545454545455;
25 25 272.7272727272727;
];
msh.MAX = max(msh.POS);
msh.MIN = min(msh.POS);
msh.LINES =[
 2 9 0
 9 10 0
 10 11 0
 11 12 0
 12 13 0
 13 14 0
 14 15 0
 15 16 0
 16 17 0
 17 18 0
 18 1 0
 1 3 0
 4 19 0
 19 20 0
 20 21 0
 21 22 0
 22 23 0
 23 24 0
 24 25 0
 25 26 0
 26 27 0
 27 28 0
 28 3 0
 2 4 0
 6 29 0
 29 30 0
 30 31 0
 31 32 0
 32 33 0
 33 34 0
 34 35 0
 35 36 0
 36 37 0
 37 38 0
 38 5 0
 5 7 0
 8 39 0
 39 40 0
 40 41 0
 41 42 0
 42 43 0
 43 44 0
 44 45 0
 45 46 0
 46 47 0
 47 48 0
 48 7 0
 6 8 0
 2 6 0
 1 5 0
 4 8 0
 3 7 0
];
msh.QUADS =[
 2 9 19 4 0
 9 10 20 19 0
 10 11 21 20 0
 11 12 22 21 0
 12 13 23 22 0
 13 14 24 23 0
 14 15 25 24 0
 15 16 26 25 0
 16 17 27 26 0
 17 18 28 27 0
 18 1 3 28 0
 6 8 39 29 0
 29 39 40 30 0
 30 40 41 31 0
 31 41 42 32 0
 32 42 43 33 0
 33 43 44 34 0
 34 44 45 35 0
 35 45 46 36 0
 36 46 47 37 0
 37 47 48 38 0
 38 48 7 5 0
 2 6 29 9 0
 9 29 30 10 0
 10 30 31 11 0
 11 31 32 12 0
 12 32 33 13 0
 13 33 34 14 0
 14 34 35 15 0
 15 35 36 16 0
 16 36 37 17 0
 17 37 38 18 0
 18 38 5 1 0
 4 19 39 8 0
 19 20 40 39 0
 20 21 41 40 0
 21 22 42 41 0
 22 23 43 42 0
 23 24 44 43 0
 24 25 45 44 0
 25 26 46 45 0
 26 27 47 46 0
 27 28 48 47 0
 28 3 7 48 0
 2 4 8 6 0
 1 5 7 3 0
];
msh.HEXAS =[
 2 9 29 6 4 19 39 8 0
 9 10 30 29 19 20 40 39 0
 10 11 31 30 20 21 41 40 0
 11 12 32 31 21 22 42 41 0
 12 13 33 32 22 23 43 42 0
 13 14 34 33 23 24 44 43 0
 14 15 35 34 24 25 45 44 0
 15 16 36 35 25 26 46 45 0
 16 17 37 36 26 27 47 46 0
 17 18 38 37 27 28 48 47 0
 18 1 5 38 28 3 7 48 0
];
msh.PNT =[
 1 0
 2 0
 3 0
 4 0
 5 0
 6 0
 7 0
 8 0
];
