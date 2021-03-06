import numpy as np

substance_p= {'Q': 3,
              'WH': 0,
              'WV': 300,
              'fasta': 'RPKPQQFFGLM',
              'modifications': {'11': {'C_carbo': {'H': 1, 'N': 1, 'O': -1}}},
              'name': 'substanceP',
              'spectrum': (np.array([     191.932,
                                          271.183,
                                          415.8,
                                          425.948,
                                          444.232,
                                          444.359,
                                          444.568,
                                          445.236,
                                          449.284,
                                          449.909,
                                          450.247,
                                          450.576,
                                          450.912,
                                          451.073,
                                          451.23,
                                          451.766,
                                          451.892,
                                          452.795,
                                          453.08,
                                          453.197,
                                          453.301,
                                          453.803,
                                          454.11,
                                          454.848,
                                          454.901,
                                          455.24,
                                          455.571,
                                          455.833,
                                          455.915,
                                          456.016,
                                          456.069,
                                          456.2,
                                          456.318,
                                          456.399,
                                          456.443,
                                          456.88,
                                          457.062,
                                          457.892,
                                          457.979,
                                          458.062,
                                          458.259,
                                          458.383,
                                          459.07,
                                          459.885,
                                          459.972,
                                          475.927,
                                          496.342,
                                          497.326,
                                          498.322,
                                          501.138,
                                          539.586,
                                          557.435,
                                          581.378,
                                          596.367,
                                          597.292,
                                          598.302,
                                          601.319,
                                          624.384,
                                          625.393,
                                          626.407,
                                          626.474,
                                          635.076,
                                          656.359,
                                          666.305,
                                          697.133,
                                          707.855,
                                          708.433,
                                          725.358,
                                          726.348,
                                          727.335,
                                          731.066,
                                          741.324,
                                          746.009,
                                          747.33,
                                          748.33,
                                          751.438,
                                          752.442,
                                          753.442,
                                          754.462,
                                          771.36,
                                          789.906,
                                          790.264,
                                          792.298,
                                          801.086,
                                          805.496,
                                          810.221,
                                          813.224,
                                          815.206,
                                          816.789,
                                          818.37,
                                          835.176,
                                          839.268,
                                          843.906,
                                          844.998,
                                          849.017,
                                          849.099,
                                          852.21,
                                          852.448,
                                          852.945,
                                          853.416,
                                          854.437,
                                          855.416,
                                          855.527,
                                          856.112,
                                          856.422,
                                          856.514,
                                          860.391,
                                          862.153,
                                          867.943,
                                          869.413,
                                          871.409,
                                          871.827,
                                          873.756,
                                          874.337,
                                          875.383,
                                          877.395,
                                          883.424,
                                          884.484,
                                          890.205,
                                          896.435,
                                          897.193,
                                          898.188,
                                          898.49,
                                          898.568,
                                          899.513,
                                          900.513,
                                          900.611,
                                          901.209,
                                          901.522,
                                          903.475,
                                          904.368,
                                          905.264,
                                          905.522,
                                          906.361,
                                          907.38,
                                          908.375,
                                          910.359,
                                          913.492,
                                          914.024,
                                          914.104,
                                          915.107,
                                          916.149,
                                          918.505,
                                          922.498,
                                          924.483,
                                          927.183,
                                          936.384,
                                          978.048,
                                          983.471,
                                          1003.604,
                                          1004.6,
                                          1016.452,
                                          1019.84,
                                          1027.544,
                                          1030.571,
                                          1030.661,
                                          1037.525,
                                          1044.54,
                                          1045.567,
                                          1046.474,
                                          1046.602,
                                          1047.585,
                                          1048.597,
                                          1049.606,
                                          1050.588,
                                          1054.595,
                                          1055.565,
                                          1059.788,
                                          1060.606,
                                          1062.016,
                                          1062.579,
                                          1067.493,
                                          1068.564,
                                          1069.551,
                                          1070.594,
                                          1072.575,
                                          1076.041,
                                          1077.067,
                                          1078.563,
                                          1078.725,
                                          1079.565,
                                          1080.569,
                                          1081.535,
                                          1081.614,
                                          1082.551,
                                          1087.606,
                                          1088.602,
                                          1090.557,
                                          1091.664,
                                          1094.547,
                                          1095.577,
                                          1096.55,
                                          1097.787,
                                          1100.508,
                                          1101.575,
                                          1102.251,
                                          1102.609,
                                          1103.386,
                                          1103.611,
                                          1104.614,
                                          1105.614,
                                          1106.226,
                                          1106.618,
                                          1107.609,
                                          1107.932,
                                          1111.083,
                                          1111.543,
                                          1125.562,
                                          1126.57,
                                          1139.623,
                                          1140.595,
                                          1149.595,
                                          1173.691,
                                          1174.682,
                                          1175.688,
                                          1183.661,
                                          1189.731,
                                          1190.607,
                                          1191.639,
                                          1200.667,
                                          1201.667,
                                          1202.685,
                                          1216.468,
                                          1216.688,
                                          1217.704,
                                          1218.452,
                                          1218.703,
                                          1219.645,
                                          1219.735,
                                          1220.662,
                                          1220.751,
                                          1221.71,
                                          1222.625,
                                          1231.665,
                                          1234.075,
                                          1234.647,
                                          1237.703,
                                          1238.695,
                                          1239.64,
                                          1240.688,
                                          1245.64,
                                          1246.633,
                                          1247.648,
                                          1248.652,
                                          1257.658,
                                          1262.688,
                                          1262.815,
                                          1263.649,
                                          1271.693,
                                          1277.679,
                                          1278.719,
                                          1284.727,
                                          1285.761,
                                          1287.688,
                                          1288.683,
                                          1289.676,
                                          1290.71,
                                          1291.701,
                                          1292.726,
                                          1293.689,
                                          1296.738,
                                          1300.759,
                                          1301.763,
                                          1303.702,
                                          1304.616,
                                          1304.718,
                                          1305.731,
                                          1306.72,
                                          1307.691,
                                          1308.717,
                                          1310.389,
                                          1312.624,
                                          1312.685,
                                          1314.688,
                                          1315.687,
                                          1316.661,
                                          1317.647,
                                          1320.738,
                                          1321.77,
                                          1322.146,
                                          1322.368,
                                          1322.791,
                                          1330.152,
                                          1331.713,
                                          1332.11,
                                          1332.722,
                                          1332.943,
                                          1333.554,
                                          1333.733,
                                          1334.063,
                                          1334.412,
                                          1334.727,
                                          1335.734,
                                          1336.13,
                                          1336.724,
                                          1337.651,
                                          1339.748,
                                          1342.713,
                                          1344.678,
                                          1345.769,
                                          1345.918,
                                          1346.708,
                                          1347.735,
                                          1348.534,
                                          1348.736,
                                          1349.038,
                                          1349.318,
                                          1349.501,
                                          1349.744,
                                          1349.859,
                                          1350.743,
                                          1351.221,
                                          1351.315,
                                          1351.75,
                                          1352.433,
                                          1352.764,
                                          1353.623,
                                          1353.734,
                                          1354.661,
                                          1354.749,
                                          1355.331,
                                          1355.681,
                                          1356.562,
                                          1362.728,
                                          1363.725,
                                          1364.722,
                                          1365.737,
                                          1366.748,
                                          1367.749,
                                          1368.744,
                                          1369.802,
                                          1370.69,
                                          1370.783,
                                          1371.707,
                                          1371.853,
                                          1372.738,
                                          1373.701,
                                          1373.871,
                                          1374.683,
                                          1376.705,
                                          1380.695,
                                          1381.75,
                                          1397.64,
                                          1401.604,
                                          1402.682,
                                          1411.073,
                                          1416.633,
                                          1418.701,
                                          1419.77,
                                          1420.786,
                                          1444.694,
                                          1461.786,
                                          1464.65,
                                          1483.735,
                                          1499.602]),
                           np.array(     [17.36,
                                          98.33,
                                          17.23,
                                          15.21,
                                          208.4,
                                          6.41,
                                          117.6,
                                          44.26,
                                          19.72,
                                          553.6,
                                          474.8,
                                          151.3,
                                          221.8,
                                          28.98,
                                          87.63,
                                          9.758,
                                          48.64,
                                          33.83,
                                          35.56,
                                          29.08,
                                          20.22,
                                          22.9,
                                          16.5,
                                          22.39,
                                          21.45,
                                          164.9,
                                          172.5,
                                          13.84,
                                          10.05,
                                          51.68,
                                          38.28,
                                          18.27,
                                          21.73,
                                          25.66,
                                          20.05,
                                          26.55,
                                          47.55,
                                          11.86,
                                          30.0,
                                          6.661,
                                          14.77,
                                          12.98,
                                          10.42,
                                          8.798,
                                          62.61,
                                          11.28,
                                          65.6,
                                          21.18,
                                          35.55,
                                          8.599,
                                          4.456,
                                          16.09,
                                          21.41,
                                          23.39,
                                          155.4,
                                          84.63,
                                          54.76,
                                          273.3,
                                          92.99,
                                          32.37,
                                          0.9093,
                                          45.41,
                                          32.02,
                                          13.9,
                                          33.97,
                                          13.9,
                                          35.38,
                                          130.8,
                                          80.73,
                                          30.41,
                                          8.852,
                                          23.18,
                                          7.039,
                                          48.43,
                                          44.57,
                                          10.08,
                                          547.5,
                                          248.5,
                                          37.49,
                                          32.98,
                                          21.34,
                                          18.52,
                                          31.01,
                                          15.94,
                                          11.97,
                                          32.66,
                                          28.39,
                                          29.66,
                                          24.17,
                                          32.19,
                                          7.584,
                                          17.86,
                                          23.36,
                                          24.66,
                                          17.4,
                                          3.753,
                                          17.09,
                                          43.39,
                                          19.05,
                                          119.2,
                                          46.81,
                                          9.646,
                                          26.98,
                                          10.41,
                                          1.973,
                                          29.86,
                                          39.88,
                                          34.42,
                                          14.56,
                                          19.51,
                                          21.51,
                                          8.911,
                                          7.999,
                                          8.507,
                                          54.05,
                                          8.153,
                                          21.16,
                                          9.917,
                                          31.86,
                                          20.68,
                                          24.9,
                                          10.19,
                                          87.3,
                                          1.528,
                                          527.6,
                                          225.3,
                                          10.77,
                                          12.18,
                                          81.92,
                                          15.36,
                                          4.048,
                                          25.55,
                                          36.45,
                                          20.54,
                                          24.86,
                                          24.41,
                                          21.69,
                                          59.49,
                                          8.518,
                                          0.2358,
                                          16.32,
                                          5.675,
                                          14.82,
                                          20.53,
                                          11.95,
                                          8.285,
                                          36.41,
                                          31.1,
                                          22.29,
                                          28.73,
                                          10.66,
                                          30.55,
                                          23.9,
                                          16.82,
                                          41.17,
                                          12.27,
                                          21.23,
                                          55.41,
                                          159.3,
                                          21.85,
                                          411.5,
                                          329.7,
                                          188.6,
                                          5.514,
                                          14.11,
                                          28.84,
                                          4.459,
                                          7.404,
                                          48.16,
                                          6.394,
                                          26.25,
                                          5.535,
                                          31.04,
                                          21.39,
                                          21.23,
                                          42.16,
                                          5.675,
                                          4.626,
                                          982.0,
                                          19.78,
                                          329.6,
                                          128.0,
                                          9.565,
                                          25.02,
                                          39.75,
                                          64.58,
                                          61.91,
                                          6.055,
                                          13.12,
                                          73.9,
                                          51.92,
                                          62.57,
                                          3.831,
                                          11.42,
                                          26.0,
                                          5.25,
                                          39.94,
                                          9.986,
                                          1140.0,
                                          961.9,
                                          444.3,
                                          7.547,
                                          86.0,
                                          53.31,
                                          4.365,
                                          6.723,
                                          45.73,
                                          76.99,
                                          30.66,
                                          24.29,
                                          7.007,
                                          17.81,
                                          51.6,
                                          8.726,
                                          21.9,
                                          35.54,
                                          12.34,
                                          28.92,
                                          16.49,
                                          64.39,
                                          20.25,
                                          36.97,
                                          30.4,
                                          364.0,
                                          609.8,
                                          8.031,
                                          294.2,
                                          43.2,
                                          87.76,
                                          6.82,
                                          19.23,
                                          5.59,
                                          10.11,
                                          9.036,
                                          7.274,
                                          20.51,
                                          10.66,
                                          25.28,
                                          31.91,
                                          9.263,
                                          12.02,
                                          19.39,
                                          17.31,
                                          25.48,
                                          42.61,
                                          15.94,
                                          3.819,
                                          8.122,
                                          11.01,
                                          9.179,
                                          31.13,
                                          19.34,
                                          22.01,
                                          30.25,
                                          40.2,
                                          47.07,
                                          201.7,
                                          132.4,
                                          38.05,
                                          51.51,
                                          28.16,
                                          19.65,
                                          20.11,
                                          70.05,
                                          29.18,
                                          166.6,
                                          196.3,
                                          107.7,
                                          34.04,
                                          24.42,
                                          8.539,
                                          26.87,
                                          45.37,
                                          46.75,
                                          93.41,
                                          23.72,
                                          14.56,
                                          112.5,
                                          37.32,
                                          5.574,
                                          20.16,
                                          9.237,
                                          28.94,
                                          406.3,
                                          4.938,
                                          1397.0,
                                          3.838,
                                          3.057,
                                          999.9,
                                          13.53,
                                          8.78,
                                          395.1,
                                          85.22,
                                          4.626,
                                          44.7,
                                          6.261,
                                          26.28,
                                          11.6,
                                          17.45,
                                          4.364,
                                          3.831,
                                          21.36,
                                          719.3,
                                          9.618,
                                          1858.0,
                                          4.709,
                                          10.95,
                                          23.88,
                                          2407.0,
                                          84.74,
                                          1610.0,
                                          10.03,
                                          4.469,
                                          564.1,
                                          10.82,
                                          134.0,
                                          1.624,
                                          68.8,
                                          37.62,
                                          43.07,
                                          7.076,
                                          54.78,
                                          8.174,
                                          8.002,
                                          60.06,
                                          244.2,
                                          436.9,
                                          283.7,
                                          102.0,
                                          72.51,
                                          34.43,
                                          130.1,
                                          3.728,
                                          181.3,
                                          2.423,
                                          133.9,
                                          32.43,
                                          12.26,
                                          7.073,
                                          21.17,
                                          19.7,
                                          16.48,
                                          20.15,
                                          26.87,
                                          19.73,
                                          8.122,
                                          35.5,
                                          34.83,
                                          44.43,
                                          6.782,
                                          19.79,
                                          25.82,
                                          15.09,
                                          26.2,
                                          24.87]))}
