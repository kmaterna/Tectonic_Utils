"""functions to support the conversion between utm coordinate system (meters) and lat/lon (degrees)"""

import numpy as np


def deg2utm(lat, lon):
    """
    [x,y,utmzone] = deg2utm(Lat,Lon)
    Description: Function to convert lat/lon vectors into UTM coordinates (WGS84).
    Some code has been extracted from UTM.m function by Gabriel Ruiz Martinez.

    Inputs:
     Lat: Latitude vector (np array).   Degrees.  +ddd.ddddd  WGS84
     Lon: Longitude vector (np array).  Degrees.  +ddd.ddddd  WGS84
    Outputs:
     x, y , utmzone.   See example

    Example 1:
     Lat=[40.3154333; 46.283900; 37.577833; 28.645650; 38.855550; 25.061783];
     Lon=[-3.4857166; 7.8012333; -119.95525; -17.759533; -94.7990166; 121.640266];
     [x,y,utmzone] = deg2utm(Lat,Lon);

     fprintf('%7.0f ',x)
        458731  407653  239027  230253  343898  362850
     fprintf('%7.0f ',y)
        4462881 5126290 4163083 3171843 4302285 2772478
     utmzone =
        30 T
        32 T
        11 S
        28 R
        15 S
        51 R

    Example 2: If you have Lat/Lon coordinates in Degrees, Minutes and Seconds
     LatDMS=[40 18 55.56; 46 17 2.04];
     LonDMS=[-3 29  8.58;  7 48 4.44];
     Lat=dms2deg(mat2dms(LatDMS)); %convert into degrees
     Lon=dms2deg(mat2dms(LonDMS)); %convert into degrees
     [x,y,utmzone] = deg2utm(Lat,Lon)

    Author:
    Rafael Palacios
    Universidad Pontificia Comillas
    Madrid, Spain
    Version: Apr/06, Jun/06, Aug/06, Aug/06
    Aug/06: fixed a problem (found by Rodolphe Dewarrat) related to southern hemisphere coordinates.
    Aug/06: corrected m-Lint warnings
    """
    # Argument checking
    n1 = len(lat);
    n2 = len(lon);
    if n1 != n2:
        raise Exception('Lat and Lon vectors should have the same length');

    # Earth shape constants
    sa = 6378137.000000;
    sb = 6356752.314245;

    # Memory pre-allocation
    x = np.zeros(np.shape(lon));
    y = np.zeros(np.shape(lat));
    utmzone = ['60 X' for x in lon];

    # Main loop
    for i in range(len(lat)):
        la = lat[i];
        lo = lon[i];

        # %e = (((sa ^ 2) - (sb ^ 2)) ^ 0.5) / sa;
        e2 = np.sqrt((sa * sa) - (sb * sb)) / sb;
        e2cuadrada = e2 * e2;
        c = (sa * sa) / sb;
        # alpha = (sa - sb) / sa; % f \
        # ablandamiento = 1 / alpha; % 1 / f

        lat_rad = la * (np.pi / 180);
        lon_rad = lo * (np.pi / 180);

        Huso = np.fix((lo / 6) + 31);
        S = ((Huso * 6) - 183);
        deltaS = lon_rad - (S * (np.pi / 180));

        if la < -72:
            Letra = 'C';
        elif la < -64:
            Letra = 'D';
        elif la < -56:
            Letra = 'E';
        elif la < -48:
            Letra = 'F';
        elif la < -40:
            Letra = 'G';
        elif la < -32:
            Letra = 'H';
        elif la < -24:
            Letra = 'J'
        elif la < -16:
            Letra = 'K';
        elif la < -8:
            Letra = 'L';
        elif la < 0:
            Letra = 'M';
        elif la < 8:
            Letra = 'N';
        elif la < 16:
            Letra = 'P';
        elif la < 24:
            Letra = 'Q';
        elif la < 32:
            Letra = 'R';
        elif la < 40:
            Letra = 'S';
        elif la < 48:
            Letra = 'T';
        elif la < 56:
            Letra = 'U';
        elif la < 64:
            Letra = 'V';
        elif la < 72:
            Letra = 'W';
        else:
            Letra = 'X';

        a = np.cos(lat_rad) * np.sin(deltaS);
        epsilon = 0.5 * np.log((1 + a) / (1 - a));
        nu = np.arctan(np.tan(lat_rad) / np.cos(deltaS)) - lat_rad;
        denominator = np.power((1 + (e2cuadrada * np.power(np.cos(lat_rad), 2))), 0.5);
        v = c / denominator * 0.9996;
        ta = (e2cuadrada / 2) * np.power(epsilon, 2) * np.power((np.cos(lat_rad)), 2);
        a1 = np.sin(2 * lat_rad);
        a2 = a1 * (np.cos(lat_rad)) * (np.cos(lat_rad));
        j2 = lat_rad + (a1 / 2);
        j4 = ((3 * j2) + a2) / 4;
        j6 = ((5 * j4) + (a2 * np.power(np.cos(lat_rad), 2))) / 3;
        alfa = (3 / 4) * e2cuadrada;
        beta = (5 / 3) * alfa * alfa;
        gama = (35 / 27) * np.power(alfa, 3);
        Bm = 0.9996 * c * (lat_rad - alfa * j2 + beta * j4 - gama * j6);
        xx = epsilon * v * (1 + (ta / 3)) + 500000;
        yy = nu * v * (1 + ta) + Bm;

        if yy < 0:
            yy = 9999999 + yy;

        x[i] = xx;
        y[i] = yy;
        utmzone[i] = str(Huso)+' '+Letra;

    return [x, y, utmzone];


def utm2deg(xx, yy, utmzone):
    """
    [Lat, Lon] = utm2deg(x, y, utmzone)

    Description: Function to convert vectors of UTM coordinates into Lat/Lon vectors (WGS84).
    Some code has been extracted from UTMIP.m function by Gabriel Ruiz Martinez.

    Inputs:
        x, y , utmzone.

    Outputs:
        Lat: Latitude vector.   Degrees.  +ddd.ddddd  WGS84
        Lon: Longitude vector.  Degrees.  +ddd.ddddd  WGS84

    Example 1:
        x=[ 458731;  407653;  239027;  230253;  343898;  362850];
        y=[4462881; 5126290; 4163083; 3171843; 4302285; 2772478];
        utmzone=['30 T'; '32 T'; '11 S'; '28 R'; '15 S'; '51 R'];
        [Lat, Lon]=utm2deg(x,y,utmzone);

        fprintf('%11.6f ',lat)
            40.315430   46.283902   37.577834   28.645647   38.855552   25.061780
        fprintf('%11.6f ',lon)
            -3.485713    7.801235 -119.955246  -17.759537  -94.799019  121.640266

    Example 2: If you need Lat/Lon coordinates in Degrees, Minutes and Seconds
        [Lat, Lon]=utm2deg(x,y,utmzone);
        LatDMS=dms2mat(deg2dms(Lat))
        LatDMS =
        40.00         18.00         55.55
        46.00         17.00          2.01
        37.00         34.00         40.17
        28.00         38.00         44.33
        38.00         51.00         19.96
        25.00          3.00         42.41

        LonDMS=dms2mat(deg2dms(Lon))
        LonDMS =
        -3.00         29.00          8.61
        7.00          48.00          4.40
        -119.00        57.00         18.93
        -17.00         45.00         34.33
        -94.00         47.00         56.47
        121.00         38.00         24.96

    Author:
    Rafael Palacios
    Universidad Pontificia Comillas
    Madrid, Spain
    Version: Apr/06, Jun/06, Aug/06
    Aug/06: corrected m-Lint warnings
    """

    # Argument checking
    n1 = len(xx);
    n2 = len(yy);
    n3 = len(utmzone);
    if n1 != n2 or n3 != n2:
        raise Exception('x, y, and utmzone vectors should have the same number or rows');

    # Earth constants
    sa = 6378137.000000;
    sb = 6356752.314245;

    # Memory pre-allocation
    lat_return = np.zeros(len(xx),);
    lon_return = np.zeros(len(xx),);

    # Main Loop
    for i in range(len(xx)):
        if utmzone[i][3] > 'X' or utmzone[i][3] < 'C':
            print('utm2deg: Warning utmzone should be a vector of strings like "30 T", not "30 t"\n');
        if utmzone[i][3] > 'M':
            hemis = 'N';  # Northern hemisphere
        else:
            hemis = 'S';
        x = xx[i];
        y = yy[i];
        zone = float(utmzone[i][0:2]);

        # % e = (((sa ^ 2) - (sb ^ 2)) ^ 0.5) / sa;
        e2 = np.power( ( np.power(sa, 2) - np.power(sb, 2) ), 0.5) / sb;
        e2cuadrada = np.power(e2, 2);
        c = np.power(sa, 2) / sb;

        X = x - 500000;
        if hemis == 'S' or hemis == 's':
            Y = y - 10000000;
        else:
            Y = y;

        S = ((zone * 6) - 183);
        lat = Y / (6366197.724 * 0.9996);
        v = (c / np.power( 1 + (e2cuadrada * np.power(np.cos(lat), 2)), 0.5) ) * 0.9996;
        a = X / v;
        a1 = np.sin(2 * lat);
        a2 = a1 * np.power(np.cos(lat), 2);
        j2 = lat + (a1 / 2);
        j4 = ((3 * j2) + a2) / 4;
        j6 = ((5 * j4) + (a2 * np.power(np.cos(lat), 2))) / 3;
        alfa = (3 / 4) * e2cuadrada;
        beta = (5 / 3) * np.power(alfa, 2);
        gama = (35 / 27) * np.power(alfa, 3);
        Bm = 0.9996 * c * (lat - alfa * j2 + beta * j4 - gama * j6);
        b = (Y - Bm) / v;
        Epsi = ((e2cuadrada * np.power(a, 2) ) / 2) * np.power(np.cos(lat), 2);
        Eps = a * (1 - (Epsi / 3));
        nab = (b * (1 - Epsi)) + lat;
        senoheps = (np.exp(Eps) - np.exp(-Eps)) / 2;
        Delt = np.arctan(senoheps / (np.cos(nab)));
        TaO = np.arctan(np.cos(Delt) * np.tan(nab));
        longitude = (Delt * (180 / np.pi)) + S;
        latitude = (lat + (1 + e2cuadrada * np.power(np.cos(lat), 2) - (3 / 2) * e2cuadrada * np.sin(lat) *
                           np.cos(lat) * (TaO - lat)) * (TaO - lat)) * (180 / np.pi);

        lat_return[i] = latitude;
        lon_return[i] = longitude;

    return [lat_return, lon_return];
