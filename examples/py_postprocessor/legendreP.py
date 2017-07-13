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
        lp =  0.0625 * (6435 * pow(x, 8) - 12012 * pow(x, 6) + 6930 * pow(x, 4) - 1260 * pow(x, 2) + 35)
    elif n == 9:
        lp =  0.0078125 * (12155 * pow(x, 9) - 25740 * pow(x, 7) + 18018 * pow(x,5) - 4620 * pow(x, 3) + 315 * x)
    elif n == 10:
        lp =  0.00390625 * (46189 * pow(x, 10) - 109395 * pow(x, 8) + 90090 * pow(x, 6) - 30030 * pow(x, 4) + 3465 * pow(x, 2) - 63)
    else:
        print("Legendre expansion goes up to l = 10")
    return lp
