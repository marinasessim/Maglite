def legendreP(n,x):
    if n == 0:
        lp = 1
    elif n == 1:
        lp = x
    elif n == 2:
        lp = 0.5 * (3 * pow(x, 2) - 1)
    elif n == 3:
        lp = 0.5 * (5 * pow(x, 3) - 3 * x)
    elif n == 4:
        lp = 0.125 * (35 * pow(x, 4) - 30 * pow(x, 2) + 3)
    elif n == 5:
        lp = 0.125 * (63 * pow(x, 5) - 70 * pow(x, 3) + 15 * x)
    elif n == 6:
        lp =  0.0625 * (231 * pow(x, 6) - 315 * pow(x, 4) + 105 * pow(x, 2) - 5)
    elif n == 7:
        lp =  0.0625 * (429 * pow(x, 7) - 693 * pow(x, 5) + 315 * pow(x, 3) - 35 * x)
    elif n == 8:
        lp =  0.0625 * (6435 * pow(x,8) - 12012 * pow(x,6) + 6930 * pow(x, 4) - 1260 * pow(x, 2) + 35)
    elif n == 9:
        lp =  0.0078125 * (12155 * pow(x, 9) - 25740 * pow(x, 7) + 18018 * pow(x,5) - 4620 * pow(x, 3) + 315 * x)
    elif n == 10:
        lp =  0.00390625 * (46189 * pow(x, 10) - 109395 * pow(x, 8) + 90090 * pow(x, 6) - 30030 * pow(x, 4) + 3465 * pow(x, 2) - 63)
    else:
        print("Legendre expansion goes up to l = 10")
    return lp

def shiftedLP(n,x):
    if n == 0:
        lp = 1
    elif n == 1:
        lp = 2 * x - 1
    elif n == 2:
        lp = 6 * pow(x, 2) - 6 * x + 1
    elif n == 3:
        lp = 20 * pow(x, 3) - 30 * pow(x, 2) + 12 * x - 1
    elif n == 4:
        lp = 70 * pow(x, 4) - 140 * pow(x, 3) + 90 * pow(x, 2) - 20 * x + 1
    elif n == 5:
        lp = -1 + 30 * x - 210 * pow(x, 2) + 560 * pow(x, 3) - 630 * pow(x, 4) + 252 * pow(x, 5)
    elif n == 6:
        lp = 1 - 42 * x + 420 * pow(x, 2) - 1680 * pow(x, 3) + 3150 * pow(x, 4) - 2772 * pow(x, 5) + 924 * pow(x, 6)
    elif n == 7:
        lp = -1 + 56 * x - 756 * pow(x, 2) + 4200 * pow(x, 3) - 11550 * pow(x, 4) + 16632 * pow(x, 5) - 12012 * pow(x, 6) + 3432 * pow(x, 7)
    elif n == 8:
        lp = 1 - 72 * x + 1260 * pow(x, 2) - 9240 * pow(x, 3) + 34650 * pow(x, 4) - 72072 * pow(x, 5) + 84084 * pow(x, 6) - 51480 * pow(x, 7) + 12870 * pow(x, 8)
    else:
        print("Legendre expansion goes up to n = 8")
    return lp
