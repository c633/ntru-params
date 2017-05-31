read("ntru_params.gp");

test() = {
    my(maps, ps, psplit, fail);
    maps = [[401, 113, 27, 693, 0.095, -45.4, -0.6],
    [541, 49, 15, 800, 0.149, -26.9, -13.1],
    [659, 38, 13,902, 0.175, -21.9, -17.7],
    [449, 134, 35, 770, 0.100, -49.0, -0.3],
    [613, 55, 17, 905, 0.142, -31.5, -14.9],
    [761, 42, 15, 1026, 0.183, -23.1, -20.9],
    [677, 157, 45, 1129, 0.096, -67.4, -2.0],
    [887, 81, 27, 1294, 0.143, -43.9, -21.9],
    [1087, 63, 23, 1464, 0.175, -34.2, -31.9]
    ];
    for(i = 1, #maps,
        ps = calculatePs(maps[i][1], maps[i][2], maps[i][3], maps[i][4], maps[i][5]);
        psplit = calculatePSplit(maps[i][1], maps[i][2], maps[i][3], maps[i][4]);
        if(ps != maps[i][6] || psplit != maps[i][7],
            fail += 1;
            printf("Test fail %s, \n\tgot (%g, %g)\n",
                maps[i], ps, psplit);
        );
    );
    printf("FAIL: %d/%d", fail, #maps);
}