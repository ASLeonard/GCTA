/*
 * Implementations of the statistical functions
 *
 * 2010 by Jian Yang <jian.yang.qt@gmail.com>
 * Modernized 2026 with C++20/23/26 features
 *
 * This file is distributed under the GNU General Public
 * License, Version 3.  Please see the file LICENSE for more
 * details
 */

#include "StatFunc.h"
#include <numeric>
#include <random>
#include <ranges>
#include <numbers>
#include <algorithm>
#include "Logger.h"

////////// P-value Calculatiion Functions Start ////////////////

double StatFunc::t_prob(double df, double t_value, bool two_tail) {
    double p = StatFunc::betai(df * 0.5, 0.5, df / (df + (t_value * t_value)));

    if (two_tail) return p;

    return p * 0.5;
}

double StatFunc::F_prob(double df_1, double df_2, double F_value) {
    return StatFunc::betai(df_2 * 0.5, df_1 * 0.5, df_2 / (df_2 + df_1 * F_value));
}

double StatFunc::betai(double a, double b, double x) {
    double bt;

    if (x < 0.0 || x > 1.0) LOGGER.e(0, "bad x in routine betai!");

    if (x == 0.0 || x == 1.0) bt = 0.0;
    else bt = exp(gammln(a + b) - gammln(a) - gammln(b) + a * log(x) + b * log(1.0 - x));

    if (x < (a + 1.0) / (a + b + 2.0)) return bt * betacf(a, b, x) / a;
    else return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
}

double StatFunc::gammln(double xx) noexcept {
    static constexpr std::array<double, 6> cof = {
        76.18009172947146, -86.50532032941677,
        24.01409824083091, -1.231739572450155, 
        0.1208650973866179e-2, -0.5395239384953e-5
    };

    double y = xx;
    double x = xx;
    double tmp = x + 5.5;
    tmp -= (x + 0.5) * std::log(tmp);
    double ser = 1.000000000190015;
    
    for (const auto c : cof) {
        ser += c / ++y;
    }
    
    return -tmp + std::log(2.5066282746310005 * ser / x);
}

double StatFunc::betacf(double a, double b, double x) {
    constexpr int MAXIT = 300;
    constexpr double Eps = 1.0e-08;
    constexpr double FPMIN = std::numeric_limits<double>::min() / std::numeric_limits<double>::epsilon();

    const double qab = a + b;
    const double qap = a + 1.0;
    const double qam = a - 1.0;
    
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (CommFunc::Abs(d) < FPMIN) d = FPMIN;
    d = 1.0 / d;
    double h = d;
    
    for (int m = 1; m <= MAXIT; ++m) {
        const int m2 = 2 * m;
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (CommFunc::Abs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (CommFunc::Abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (CommFunc::Abs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (CommFunc::Abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        const double del = d * c;
        h *= del;
        if (CommFunc::Abs(del - 1.0) <= Eps) break;
    }
    return h;
}

double StatFunc::chi_prob(double df, double chi_sqr_val) {
    return 1 - StatFunc::gammp(df * 0.5, chi_sqr_val * 0.5);
}

double StatFunc::gammp(double a, double x) {
    if (x < 0.0 || a <= 0.0) {
        LOGGER.e(0, "invalid arguments in routine gammp");
    }

    double gamser, gammcf, gln;
    if (x < a + 1.0) {
        gser(gamser, a, x, gln);
        return gamser;
    } else {
        gcf(gammcf, a, x, gln);
        return 1.0 - gammcf;
    }
}

void StatFunc::gser(double &gamser, double a, double x, double &gln) {
    constexpr int ITMAX = 500;
    constexpr double Eps = 1.0e-08;

    gln = gammln(a);
    if (x <= 0.0) {
        if (x < 0.0) {
            LOGGER.e(0, "x is less than 0 in routine gser");
        }
        gamser = 0.0;
        return;
    }
    
    double ap = a;
    double sum = 1.0 / a;
    double del = sum;
    
    for (int n = 0; n < ITMAX; ++n) {
        ++ap;
        del *= x / ap;
        sum += del;
        if (CommFunc::Abs(del) < CommFunc::Abs(sum) * Eps) {
            gamser = sum * std::exp(-x + a * std::log(x) - gln);
            return;
        }
    }
    LOGGER.e(0, "a is too large, and ITMAX is too small in routine gser");
}

void StatFunc::gcf(double &gammcf, double a, double x, double &gln) {
    constexpr int ITMAX = 100;
    constexpr double EPS = std::numeric_limits<double>::epsilon();
    constexpr double FPMIN = std::numeric_limits<double>::min() / EPS;

    gln = gammln(a);
    double b = x + 1.0 - a;
    double c = 1.0 / FPMIN;
    double d = 1.0 / b;
    double h = d;
    
    for (int i = 1; i <= ITMAX; ++i) {
        const double an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if (CommFunc::Abs(d) < FPMIN) d = FPMIN;
        c = b + an / c;
        if (CommFunc::Abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        const double del = d * c;
        h *= del;
        if (CommFunc::Abs(del - 1.0) <= EPS) {
            gammcf = std::exp(-x + a * std::log(x) - gln) * h;
            return;
        }
    }
    LOGGER.e(0, "a is too large, and ITMAX is too small in gcf");
}
////////// P-value Calculatiion Functions End ////////////////

///////// Random Number Generation Functions Start ////////

void StatFunc::gasdev_seq(int &idum, std::vector<double> &vec, int size, double means, double var) {
    vec.clear();
    vec.resize(size);

    const double sd = std::sqrt(var);
    for (auto& v : vec) {
        v = means + sd * gasdev(idum);
    }
}

void StatFunc::gasdev_seq(int &idum, std::vector<double> &vec, int size, double var) {
    if (size < 2) {
        LOGGER.e(0, "invalid size. StatFunc::gasdev_seq");
    }
    if (CommFunc::FloatEqual(var, 0.0)) {
        vec.assign(size, 0.0);
        return;
    }

    const double n = static_cast<double>(size);
    double s_square = 0.0;

    while (s_square < 0.01 * var) {
        gasdev_seq(idum, vec, size, 0.0, var);
        
        // Calculate and subtract mean
        const double ave = std::accumulate(vec.begin(), vec.end(), 0.0) / n;
        std::ranges::for_each(vec, [ave](double& v) { v -= ave; });
        
        // Calculate variance
        s_square = std::transform_reduce(
            vec.begin(), vec.end(), 0.0, std::plus<>{},
            [](double v) { return v * v; }
        ) / (n - 1.0);
    }
    
    const double c = std::sqrt(var / s_square);
    std::ranges::for_each(vec, [c](double& v) { v *= c; });
}

double StatFunc::gasdev(int &idum) {
    static int iset = 0;
    static double gset;

    if (idum < 0) iset = 0;
    
    if (iset == 0) {
        double v1, v2, rsq;
        do {
            v1 = 2.0 * ran1(idum) - 1.0;
            v2 = 2.0 * ran1(idum) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        
        const double fac = std::sqrt(-2.0 * std::log(rsq) / rsq);
        gset = v1 * fac;
        iset = 1;
        return v2 * fac;
    } else {
        iset = 0;
        return gset;
    }
}

double StatFunc::UniformDev(double a, double b, int &idum) {
    if (a >= b) LOGGER.e(0, "b must be larger than a. StatFunc::UniformDev");
    if (idum > 0) idum *= -1;
    return a + (b - a) * ran1(idum);
}

static std::random_device ran_device;
static std::mt19937 gen(ran_device());
static std::uniform_real_distribution<double> u_dis(0, 1);

double StatFunc::ran1(int &idum) {
    return u_dis(gen);
}
/*
    const int IA = 16807, IM = 2147483647, IQ = 127773, IR = 2836, NTAB = 32;
    const int NDIV = (1 + (IM - 1) / NTAB);
    const double EPS = 3.0e-16, AM = 1.0 / IM, RNMX = (1.0 - EPS);
    static int iy = 0;
    static vector<int> iv(NTAB);
    int j, k;
    double temp;

    if (idum <= 0 || !iy) {
        if (-idum < 1) idum = 1;
        else idum = -idum;
        for (j = NTAB + 7; j >= 0; j--) {
            k = idum / IQ;
            idum = IA * (idum - k * IQ) - IR*k;
            if (idum < 0) idum += IM;
            if (j < NTAB) iv[j] = idum;
        }
        iy = iv[0];
    }
    k = idum / IQ;
    idum = IA * (idum - k * IQ) - IR*k;
    if (idum < 0) idum += IM;
    j = iy / NDIV;
    iy = iv[j];
    iv[j] = idum;
    if ((temp = AM * iy) > RNMX) return RNMX;
    else return temp;
}
*/

double StatFunc::chidev(int &idum, double df) {
    if (df > 2.0) {
        return 2.0 * cheng_gamdev(idum, df * 0.5);
    }
    LOGGER.e(0, "invalid degree of freedom. StatFunc::chidev");
    return 0.0;
}

double StatFunc::cheng_gamdev(int &idum, double alpha) {
    const double beta = alpha - 1.0;
    
    while (true) {
        const double u1 = UniformDev(0.0, 1.0, idum);
        const double u2 = UniformDev(0.0, 1.0, idum);
        
        const double nu = (alpha - 1.0 / (6.0 * alpha)) * u1 / (beta * u2);
        const double d_buf = 2.0 * (u2 - 1.0) / beta + nu + 1.0 / nu;
        
        if (d_buf < 2.0) {
            return beta * nu;
        }
        
        const double d_buf2 = 2.0 * std::log(u2) / beta - std::log(nu) + nu;
        if (d_buf2 < 1.0) {
            return beta * nu;
        }
    }
}

int StatFunc::RandAbs(int a, int b, int &seed) {
    int stk = b - a + 1;
    int stl = 2;
    while (stl < stk) {
        stl = stl + stl;
    }
    
    const int modul = 4 * stl;
    stk = seed;
    int sti = 1;
    int rand = 0;
    
    while (sti <= 1) {
        stk = 5 * stk;
        stk = stk % modul;
        stl = stk / 4 + a;
        if (stl <= b) {
            rand = stl;
            ++sti;
        }
    }
    seed = stk;
    return rand;
}

///////// Random Number Generation Functions End ////////

double StatFunc::chi_val(double df, double prob) {
    constexpr double eps = 1.0e-08;
    double walk = 100.0;
    double chi_val = walk;
    double way = 0.0;
    double preway = 0.0;
    double prob_buf = chi_prob(df, chi_val);

    if (CommFunc::Abs(prob_buf - prob) < eps) {
        return chi_val;
    }

    if (prob_buf > prob) {
        preway = way = 1.0;
        chi_val += walk;
    } else {
        preway = way = -1.0;
        chi_val -= walk;
    }

    while (true) {
        prob_buf = chi_prob(df, chi_val);
        way = (prob_buf > prob) ? 1.0 : -1.0;

        if (CommFunc::Abs(preway - way) > eps) {
            walk *= 0.5;
        }
        chi_val += walk * way;
        preway = way;

        if (walk < eps) break;
    }

    return chi_val;
}

double StatFunc::t_val(double df, double prob) {
    constexpr double eps = 1.0e-08;
    double walk = 100.0;
    double t_val = walk;
    double way = 0.0;
    double preway = 0.0;
    double prob_buf = t_prob(df, t_val, false);

    if (CommFunc::Abs(prob_buf - prob) < eps) {
        return t_val;
    }

    if (prob_buf > prob) {
        preway = way = 1.0;
        t_val += walk;
    } else {
        preway = way = -1.0;
        t_val -= walk;
    }

    while (true) {
        prob_buf = t_prob(df, t_val, false);
        way = (prob_buf > prob) ? 1.0 : -1.0;

        if (CommFunc::Abs(preway - way) > eps) {
            walk *= 0.5;
        }
        t_val += walk * way;
        preway = way;

        if (walk < eps) break;
    }

    return t_val;
}

double StatFunc::F_val(double df_1, double df_2, double prob) {
    constexpr double eps = 1.0e-08;
    double walk = 100.0;
    double F_val = walk;
    double way = 0.0;
    double preway = 0.0;
    double prob_buf = F_prob(df_1, df_2, F_val);

    if (CommFunc::Abs(prob_buf - prob) < eps) {
        return F_val;
    }

    if (prob_buf > prob) {
        preway = way = 1.0;
        F_val += walk;
    } else {
        preway = way = -1.0;
        F_val -= walk;
    }

    while (true) {
        prob_buf = F_prob(df_1, df_2, F_val);
        way = (prob_buf > prob) ? 1.0 : -1.0;

        if (CommFunc::Abs(preway - way) > eps) {
            walk *= 0.5;
        }
        F_val += walk * way;
        preway = way;

        if (walk < eps) break;
    }

    return F_val;
}

double StatFunc::ControlFDR(std::span<const double> P_Value, double alpha, bool Restrict) {
    const int Size = static_cast<int>(P_Value.size());
    if (Size <= 1) {
        LOGGER.e(0, "invalid size. StatFunc::ControlFDR");
    }
    
    std::vector<double> P_ValueBuf(P_Value.begin(), P_Value.end());
    std::ranges::stable_sort(P_ValueBuf);
    
    double FDR_Threshold = 0.0;
    for (int i = 0; i < Size; ++i) {
        const double d_Buf = (static_cast<double>(i + 1) * alpha) / static_cast<double>(Size);
        
        if (i == Size - 1) {
            FDR_Threshold = (P_ValueBuf[i] <= d_Buf || Restrict) ? d_Buf : -1.0;
            break;
        }
        
        if (P_ValueBuf[i] <= d_Buf && P_ValueBuf[i + 1] > d_Buf) {
            FDR_Threshold = d_Buf;
            break;
        } else if (i == Size - 1) {
            FDR_Threshold = -1.0;
        }
    }
    return FDR_Threshold;
}

double StatFunc::ControlFDR_Zou(std::span<const double> GenePValue, double FDR) {
    const int Size = static_cast<int>(GenePValue.size());
    std::vector<double> Temp(GenePValue.begin(), GenePValue.end());
    std::ranges::sort(Temp);
    
    double PValue = 0.0;
    for (int i = 0; i < Size; ++i) {
        const double threshold = (i + 1) * FDR / Size;
        
        if (Temp[i] <= threshold && Temp[i + 1] > threshold) {
            PValue = threshold;
            break;
        } else if (Temp[i] <= threshold && i == Size - 1) {
            PValue = threshold;
            break;
        } else if (i == Size - 1) {
            PValue = FDR;
        }
    }
    return PValue;
}

double StatFunc::ControlFDR_Storey(std::vector<double> &P_Value, std::vector<double> &Q_Value, double CrtQ, double &FDR) {
    if (P_Value.empty()) {
        return CrtQ;
    }

    // Initialize Lambda
    std::vector<double> Lambda;
    for (double dBuf = 0.0; dBuf < 0.9001; dBuf += 0.05) {
        Lambda.push_back(dBuf);
    }

    // Calculate Pi0
    std::ranges::stable_sort(P_Value);
    const double Pi0 = CalcuPi0(P_Value, Lambda);

    // Calculate q-value
    const int m = static_cast<int>(P_Value.size());
    const double m0 = Pi0 * static_cast<double>(m);
    
    Q_Value.clear();
    Q_Value.resize(m);
    Q_Value[m - 1] = Pi0 * P_Value[m - 1];
    
    for (int i = m - 2; i >= 0; --i) {
        Q_Value[i] = (m0 * P_Value[i]) / static_cast<double>(i + 1);
        if (Q_Value[i] > Q_Value[i + 1]) {
            Q_Value[i] = Q_Value[i + 1];
        }
    }

    // Calculate FDR-adjusted critical p-value
    bool Flag = false;
    double CrtVal = 1.0;
    for (int i = 0; i < m - 1; ++i) {
        if (Q_Value[i] < CrtQ && Q_Value[i + 1] > CrtQ) {
            CrtVal = 0.5 * (P_Value[i] + P_Value[i + 1]);
            Flag = true;
            break;
        }
    }
    if (!Flag) {
        CrtVal = 0.5 * (P_Value[0] + P_Value[1]);
    }

    // Estimate FDR
    const int iBuf = std::ranges::count_if(P_Value, [CrtVal](double p) { return p < CrtVal; });
    FDR = (CrtVal * m0) / static_cast<double>(iBuf);
    if (FDR > 1.0) {
        FDR = 1.0;
    }

    return CrtVal;
}

double StatFunc::CalcuPi0(std::vector<double> &P_Value, std::vector<double> &Lambda) {
    const int P_Size = static_cast<int>(P_Value.size());
    const int LambdaSize = static_cast<int>(Lambda.size());
    
    std::vector<double> Pi(LambdaSize);
    std::vector<double> NewPi(LambdaSize);
    
    for (int i = 0; i < LambdaSize; ++i) {
        const int Count = std::ranges::count_if(P_Value, [lambda = Lambda[i]](double p) {
            return p > lambda;
        });
        Pi[i] = Count / (P_Size * (1.0 - Lambda[i]));
    }
    
    double PrePi = 0.0;
    spline(Lambda, Pi, 0.0, 0.0, NewPi);
    splint(Lambda, Pi, NewPi, Lambda[LambdaSize - 1], PrePi);
    return std::min(PrePi, 1.0);
}

void StatFunc::spline(std::vector<double> &x, std::vector<double> &y, double yp1, double ypn, std::vector<double> &y2) {
    const int n = static_cast<int>(y2.size());
    std::vector<double> u(n - 1);
    
    if (yp1 > 0.99e30) {
        y2[0] = u[0] = 0.0;
    } else {
        y2[0] = -0.5;
        u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
    }
    
    for (int i = 1; i < n - 1; ++i) {
        const double sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
        const double p = sig * y2[i - 1] + 2.0;
        y2[i] = (sig - 1.0) / p;
        u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
        u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
    }
    
    double qn, un;
    if (ypn > 0.99e30) {
        qn = un = 0.0;
    } else {
        qn = 0.5;
        un = (3.0 / (x[n - 1] - x[n - 2])) * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
    }
    
    y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);
    for (int k = n - 2; k >= 0; --k) {
        y2[k] = y2[k] * y2[k + 1] + u[k];
    }
}

void StatFunc::splint(std::vector<double> &xa, std::vector<double> &ya, std::vector<double> &y2a, double x, double &y) {
    const int n = static_cast<int>(xa.size());
    int klo = 0;
    int khi = n - 1;
    
    while (khi - klo > 1) {
        const int k = (khi + klo) >> 1;
        if (xa[k] > x) {
            khi = k;
        } else {
            klo = k;
        }
    }
    
    const double h = xa[khi] - xa[klo];
    if (h == 0.0) {
        LOGGER.e(0, "bad xa input to routine splint");
    }
    
    const double a = (xa[khi] - x) / h;
    const double b = (x - xa[klo]) / h;
    y = a * ya[klo] + b * ya[khi] + 
        ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h) / 6.0;
}

double StatFunc::erf(double x) noexcept {
    constexpr double a1 = 0.070523084;
    constexpr double a2 = 0.0422820123;
    constexpr double a3 = 0.0092705272;
    constexpr double a4 = 0.0001520143;
    constexpr double a5 = 0.0002765672;
    constexpr double a6 = 0.0000430638;
    
    const double tmp = 1.0 + a1 * x + a2 * std::pow(x, 2) + a3 * std::pow(x, 3) + 
                       a4 * std::pow(x, 4) + a5 * std::pow(x, 5) + a6 * std::pow(x, 6);
    return 1.0 - std::pow(tmp, -16);
}

/*
n <- length(p)
lp <- length(p)
i <- lp:1L
o <- order(p, decreasing = TRUE)
ro <- order(o) #This is the index of the initial p-values, highest to lowest
pmin(1, cummin(n/i * p[o]))[ro]
*/

std::vector<double> StatFunc::ControlFDR_BH(std::span<const double> p_value) {
    const int n = static_cast<int>(p_value.size());
    
    std::vector<std::pair<double, int>> pval_buf(n);
    for (int i = 0; i < n; ++i) {
        pval_buf[i] = {p_value[i], i};
    }
    
    std::ranges::stable_sort(pval_buf, [](const auto& a, const auto& b) {
        return a.first > b.first;
    });
    
    std::vector<double> fdr(n);
    double min_val = 1.0;
    
    for (int i = 0; i < n; ++i) {
        const double c = static_cast<double>(n) / static_cast<double>(n - i) * pval_buf[i].first;
        min_val = std::min(c, min_val);
        fdr[pval_buf[i].second] = std::min(1.0, min_val);
    }
    
    return fdr;
}

double StatFunc::pnorm(double x) noexcept {
    const double z = (x > 0) ? -x : x;
    constexpr double sqrt2pi = 2.50662827463;
    
    const double t0 = 1.0 / (1.0 + 0.2316419 * std::fabs(z));
    const double z1 = std::exp(-0.5 * z * z) / sqrt2pi;
    const double p0 = z1 * t0 * (0.31938153 + t0 * (-0.356563782 + 
                      t0 * (1.781477937 + t0 * (-1.821255978 + 1.330274429 * t0))));
    
    return (x >= 0) ? p0 : 1.0 - p0;
}

double StatFunc::dnorm(double x) noexcept {
    constexpr double inv_sqrt_2pi = 0.39894228;
    return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}

double StatFunc::qnorm_sub(double x, double y) {
    const double y2 = y * y;
    const double y3 = y2 * y;
    const double y4 = y3 * y;
    const double x2 = x * x;
    const double x3 = x2 * x;
    
    return y + 0.5 * x * y2 + (2 * x2 + 1) * y3 / 6.0 + 
           (6 * x3 + 7 * x) * y4 / 12.0;
}

double StatFunc::qnorm(double p, bool upper) {
    if (upper) {
        p = 1.0 - p;
    }
    
    double x = 0.0;
    for (int i = 0; i < 4; ++i) {
        x = x + qnorm_sub(x, (pnorm(x) - p) / dnorm(x));
    }
    return x;
}

double StatFunc::pchisq(double x, double df) {
    if (x < 0) return -9;

    double p, q;
    int st = 0; // error variable
    int w = 1; // function variable
    double bnd = 1; // boundary function

    // NCP is set to 0
    cdfchi(&w, &p, &q, &x, &df, &st, &bnd);

    // Check status
    if (st != 0) return -9;

    // Return p-value
    return q;
}

double StatFunc::qchisq(double q, double df) {
    if (q < 0) return -9;
    else if (q >= 1) return 0;

    double x;
    double p = 1 - q;
    int st = 0; // error variable
    int w = 2; // function variable
    double bnd = 1; // boundary function

    // NCP is set to 0
    cdfchi(&w, &p, &q, &x, &df, &st, &bnd);

    // Check status
    if (st != 0) return -9;

    // Return p-value
    return x;
}

//#################
// functions to calculate pchisqsum 

double StatFunc::pchisqsum(double x, Eigen::VectorXd lambda) {
    double pval = psadd(x, lambda);
    if (pval > 1.0) pval = psatt(x, lambda);
    return pval;
}

double StatFunc::psadd(double x, Eigen::VectorXd lambda) {
    double d = lambda.maxCoeff();
    if (d <= 0.0) return 2.0;
    lambda = lambda.array() / d;
    x = x / d;

    double lmin = 0.0;
    double m = lambda.minCoeff();
    if (m < 0.0) lmin = 0.499995 / m;
    else if (x > lambda.sum()) lmin = -0.01;
    else lmin = -0.5 * (double) lambda.size() / x;
    double lmax = 0.499995 / lambda.maxCoeff();

    double hatzeta = Brents_Kp_min_x(lambda, x, lmin, lmax, 1e-08);
    if(hatzeta > lmax + 9) return 2.0;
    double sign = (hatzeta < 0.0) ? -1.0 : 1.0;
    double w = sign * sqrt(2 * (hatzeta * x - K(hatzeta, lambda)));
    double v = hatzeta * sqrt(Kpp(hatzeta, lambda));

    // debug
    //LOGGER<<"hatzeta = "<<hatzeta<<endl;
    //LOGGER<<"w = "<<w<<endl;
    //LOGGER<<"v = "<<v<<endl;

    
    if (fabs(hatzeta) < 1e-04) return 2.0;
    else return pnorm(w + log(v / w) / w);
}

double StatFunc::psatt(double x, Eigen::VectorXd lambda) {
    double sum=lambda.sum();
    if(CommFunc::FloatEqual(sum,0.0)) return 2.0;
  
    double sq_sum=lambda.dot(lambda);
    double sum_sq=sum*sum;
    double a=sq_sum/sum;
    double b=sum_sq/sq_sum;
    
    return pchisq(x/a, b);
}

double StatFunc::K(double zeta, Eigen::VectorXd &lambda) {
    return -0.5*(1.0 - 2.0 * zeta * lambda.array()).log().sum();
}

double StatFunc::Kp(double zeta, Eigen::VectorXd &lambda) {
    return (lambda.array() / (1.0 - 2.0 * zeta * lambda.array())).sum();
}

double StatFunc::Kpp(double zeta, Eigen::VectorXd &lambda) {
    return 2.0 * (lambda.array().square() / (1.0 - 2.0 * zeta * lambda.array()).array().square()).sum();
}

double StatFunc::Kp_min_x(double zeta, Eigen::VectorXd &lambda, double x) {
    return Kp(zeta, lambda) - x;
}

double StatFunc::Brents_Kp_min_x(Eigen::VectorXd &lambda, double x, double lowerLimit, double upperLimit, double errorTol) {
    double a = lowerLimit;
    double b = upperLimit;
    double c = 0;
    double d = 1.7976931348623157E+308;

    double fa = Kp_min_x(a, lambda, x);
    double fb = Kp_min_x(b, lambda, x);

    double fc = 0;
    double s = 0;
    double fs = 0;
    
    // if f(a) f(b) >= 0 then error-exit
    if (fa * fb >= 0) {
        if (fa < fb)
            return a;
        else
            return b;
    }

    // if |f(a)| < |f(b)| then swap (a,b) end if
    if (fabs(fa) < fabs(fb)) {
        double tmp = a;
        a = b;
        b = tmp;
        tmp = fa;
        fa = fb;
        fb = tmp;
    }

    c = a;
    fc = fa;
    bool mflag = true;
    int i = 0;

    while (!(fb == 0) && (fabs(a - b) > errorTol)) {
        if ((fa != fc) && (fb != fc))
            // Inverse quadratic interpolation
            s = a * fb * fc / (fa - fb) / (fa - fc) + b * fa * fc / (fb - fa) / (fb - fc) + c * fa * fb / (fc - fa) / (fc - fb);
        else
            // Secant Rule
            s = b - fb * (b - a) / (fb - fa);

        double tmp2 = (3 * a + b) / 4;
        if ((!(((s > tmp2) && (s < b)) || ((s < tmp2) && (s > b)))) || (mflag && (fabs(s - b) >= (fabs(b - c) / 2))) || (!mflag && (fabs(s - b) >= (fabs(c - d) / 2)))) {
            s = (a + b) / 2;
            mflag = true;
        } else {
            if ((mflag && (fabs(b - c) < errorTol)) || (!mflag && (fabs(c - d) < errorTol))) {
                s = (a + b) / 2;
                mflag = true;
            } else
                mflag = false;
        }
        fs = Kp_min_x(s, lambda, x);
        d = c;
        c = b;
        fc = fb;
        if (fa * fs < 0) {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }

        // if |f(a)| < |f(b)| then swap (a,b) end if
        if (fabs(fa) < fabs(fb)) {
            double tmp = a;
            a = b;
            b = tmp;
            tmp = fa;
            fa = fb;
            fb = tmp;
        }
        i++;
        if (i > 1000) return upperLimit+10;
    }
    return b;
}


std::vector<size_t> StatFunc::sort_re_index(std::span<const double> x) {
    std::vector<size_t> x_index(x.size());
    std::iota(x_index.begin(), x_index.end(), 0);

    std::ranges::sort(x_index, [&x](size_t index1, size_t index2) {
        return x[index1] > x[index2];
    });

    return x_index;
}
