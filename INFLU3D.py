import numpy as np




# Define constants
PI = np.pi
TWOPI = 2.0 * PI
EPS = 1e-06
PNDS = 1.0

# Define source and doublet strengths
DUB = 1.0 / (4.0 * PI)
SIG = 1.0 / (4.0 * PI)

# Define panel coordinates
X = np.array([-0.5, 0.5, 0.5, -0.5, -0.5])
Y = np.array([-0.5, -0.5, 0.5, 0.5, -0.5])
Z = np.array([0.0, 0.0, 0.0, 0.0, 0.0])

# Define point of interest
PXI, PYI, PZI = map(float, input('ENTER POINT OF INTEREST (X,Y,Z): ').split())
PXO, PYO, PZO = 0.0, 0.0, 0.0

while True:
    RJ31 = 0.0
    CJK1 = 0.0
    VXS = 0.0
    VYS = 0.0
    VZS = 0.0
    VXD = 0.0
    VYD = 0.0
    VZD = 0.0
    VX = 0.0
    VY = 0.0
    VZ = 0.0

    # Calculate relative coordinates
    PX = PXI - PXO
    PY = PYI - PYO
    PZ = PZI - PZO
    RDIST = np.sqrt(PX * PX + PY * PY + PZ * PZ)
    PNLX = np.sum(X[:4]) * 0.25
    PNLY = np.sum(Y[:4]) * 0.25
    PNLZ = np.sum(Z[:4]) * 0.25
    PNX = PX - PNLX
    PNY = PY - PNLY
    PNZ = PZ - PNLZ
    PNS = np.sqrt(PNX * PNX + PNY * PNY + PNZ * PNZ)
    D1X = X[2] - X[0]
    D1Y = Y[2] - Y[0]
    D1Z = Z[2] - Z[0]
    D2X = X[3] - X[1]
    D2Y = Y[3] - Y[1]
    D2Z = Z[3] - Z[1]
    CRX = D1Y * D2Z - D2Y * D1Z
    CRY = D2X * D1Z - D1X * D2Z
    CRZ = D1X * D2Y - D2X * D1Y
    CRSQ = np.sqrt(CRX * CRX + CRY * CRY + CRZ * CRZ)
    AREA = CRSQ / 2.0
    CNX = CRX / CRSQ
    CNY = CRY / CRSQ
    CNZ = CRZ / CRSQ
    PNN = CNX * PNX + CNY * PNY + CNZ * PNZ
    TCMX = (X[2] + X[3]) / 2.0 - PNLX
    TCMY = (Y[2] + Y[3]) / 2.0 - PNLY
    TCMZ = (Z[2] + Z[3]) / 2.0 - PNLZ
    TMS = np.sqrt(TCMX * TCMX + TCMY * TCMY + TCMZ * TCMZ)
    CMX = ((X[2] + X[3]) / 2.0 - PNLX) / TMS
    CMY = ((Y[2] + Y[3]) / 2.0 - PNLY) / TMS
    CMZ = ((Z[2] + Z[3]) / 2.0 - PNLZ) / TMS
    CLX = CMY * CNZ - CNY * CMZ
    CLY = CNX * CMZ - CMX * CNZ
    CLZ = CMX * CNY - CNX * CMY

    for J in range(4):
        K = J + 1
        AX = PX - X[J]
        AY = PY - Y[J]
        AZ = PZ - Z[J]
        BX = PX - X[K]
        BY = PY - Y[K]
        BZ = PZ - Z[K]
        SX = X[K] - X[J]
        SY = Y[K] - Y[J]
        SZ = Z[K] - Z[J]
        A = np.sqrt(AX * AX + AY * AY + AZ * AZ)
        B = np.sqrt(BX * BX + BY * BY + BZ * BZ)
        S = np.sqrt(SX * SX + SY * SY + SZ * SZ)
        SM = SX * CMX + SY * CMY + SZ * CMZ
        SL = SX * CLX + SY * CLY + SZ * CLZ
        AM = AX * CMX + AY * CMY + AZ * CMZ
        AL = AX * CLX + AY * CLY + AZ * CLZ
        BM = BX * CMX + BY * CMY + BZ * CMZ
        ALL = AM * SL - AL * SM

        if A + B - S > 0.0 and S > 0.0:
            RJ3 = np.log((A + B + S) / (A + B - S)) / S
        else:
            RJ3 = 0.0

        PA = PNZ * PNZ * SL + ALL * AM
        PB = PA - ALL * SM
        RNUM = SM * PNZ * (B * PA - A * PB)
        DNOM = PA * PB + PNZ * PNZ * A * B * SM * SM
        if abs(PNZ) < EPS:
            DE = 0.0
        else:
            if RNUM != 0:
                DE = np.arctan2(RNUM, DNOM)
            else:
                DE = 0.0
        RJ31 = RJ31 - SIG * ALL * RJ3
        CJK1 = CJK1 - DUB * DE
        VXS = VXS + SIG * (RJ3 * (SM * CLX - SL * CMX) + DE * CNX)
        VYS = VYS + SIG * (RJ3 * (SM * CLY - SL * CMY) + DE * CNY)
        VZS = VZS + SIG * (RJ3 * (SM * CLZ - SL * CMZ) + DE * CNZ)
        AVBX = AY * BZ - AZ * BY
        AVBY = AZ * BX - AX * BZ
        AVBZ = AX * BY - AY * BX
        ADB = AX * BX + AY * BY + AZ * BZ
        VMOD = (A + B) / (A * B * (A * B + ADB))
        VXD = VXD + DUB * VMOD * AVBX
        VYD = VYD + DUB * VMOD * AVBY
        VZD = VZD + DUB * VMOD * AVBZ

    DTT = TWOPI
    if RDIST > 0.0:
        PNDS = PNZ ** 2 / RDIST
    if PNDS < EPS and RDIST > EPS:
        DTT = PNZ * AREA / np.sqrt(RDIST) / RDIST
    if abs(DTT) < abs(CJK1):
        CJK1 = DTT
    if RDIST < EPS * EPS:
        CJK1 = -TWOPI

    CJK = CJK1
    BJK = RJ31 - PNZ * CJK1
    VX = VXD + VXS
    VY = VYD + VYS
    VZ = VZD + VZS
    TVS = np.sqrt(VXS * VXS + VYS * VYS + VZS * VZS)
    TVD = np.sqrt(VXD * VXD + VYD * VYD + VZD * VZD)
    TV = np.sqrt(VX * VX + VY * VY + VZ * VZ)

    print('AREA OF PANEL =', AREA)
    print('SOURCE (POTENTIAL) =', BJK)
    print('SOURCE (VELOCITY):')
    print('VXS =', VXS)
    print('VYS =', VYS)
    print('VZS =', VZS)
    print('TVS =', TVS)
    print('DOUBLET (POTENTIAL):', CJK)
    print('DOUBLET (VELOCITY):')
    print('VXD =', VXD)
    print('VYD =', VYD)
    print('VZD =', VZD)
    print('TVD =', TVD)
    print('TOTAL VELOCITY:')
    print('VX =', VX)
    print('VY =', VY)
    print('VZ =', VZ)
    
    another_try = input('DO YOU WANT ANOTHER TRY? 1:YES, 2:NO ')
    if another_try != '1':
        break
