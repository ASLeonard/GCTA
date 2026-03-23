/*
 * Implementations of the statistical functions
 *
 * 2010 by Jian Yang <jian.yang.qt@gmail.com>
 *
 * This file is distributed under the GNU General Public
 * License, Version 3.  Please see the file LICENSE for more
 * details
 */

#include "StatFunc.h"
#include "numeric"
#include <random>
#include "Logger.h"

namespace {

// Returns ln(Q(a, x)) — the log of the upper regularized incomplete gamma —
// using Boost for moderate x and a log-space asymptotic expansion for extreme x
// where gamma_q(a, x) would underflow to zero.
static double log_gamma_q_upper(double a, double x) {
    double q = boost::math::gamma_q(a, x);
    if (q > 0.0) return std::log(q);
    // Asymptotic expansion coefficients: Q(a,x) ~ exp(-x + (a-1)*ln(x) - lgamma(a)) * S
    // where S = 1 + (a-1)/x + (a-1)(a-2)/x^2 + ... (truncate before diverging)
    double log_leading = -x + (a - 1.0) * std::log(x) - std::lgamma(a);
    double s = 1.0, term = 1.0, min_abs = 1.0;
    for (int k = 1; k <= 50; ++k) {
        term *= (a - k) / x;
        if (std::abs(term) >= min_abs) break;
        min_abs = std::abs(term);
        s += term;
    }
    return log_leading + std::log(std::abs(s));
}

} // anonymous namespace

////////// P-value Calculatiion Functions Start ////////////////

double StatFunc::t_prob(double df, double t_value, bool two_tail) {
    double p = boost::math::cdf(
        boost::math::complement(boost::math::students_t(df), std::abs(t_value)));
    return two_tail ? 2.0 * p : p;
}

double StatFunc::F_prob(double df_1, double df_2, double F_value) {
    return boost::math::cdf(
        boost::math::complement(boost::math::fisher_f(df_1, df_2), F_value));
}

double StatFunc::chi_prob(double df, double chi_sqr_val, bool log_p) {
    if (log_p) return log_gamma_q_upper(df * 0.5, chi_sqr_val * 0.5);
    return boost::math::cdf(
        boost::math::complement(boost::math::chi_squared(df), chi_sqr_val));
}

////////// P-value Calculatiion Functions End ////////////////

///////// Random Number Generation Functions Start ////////

void StatFunc::gasdev_seq(int &idum, vector<double> &vec, int size, double means, double var) {
    vec.clear();
    vec.resize(size);

    int i = 0;
    for (i = 0; i < size; i++) vec[i] = means + sqrt(var) * StatFunc::gasdev(idum);
}

void StatFunc::gasdev_seq(int &idum, vector<double> &vec, int size, double var) {
    if (size < 2) LOGGER.e(0, "invalid size. StatFunc::gasdev_seq");
    if (CommFunc::FloatEqual(var, 0.0)) {
        vec.clear();
        vec.resize(size);
        return;
    }

    int i = 0;
    double ave = 0.0, s_square = 0.0;
    double n = (double) size;

    while (s_square < 0.01 * var) {
        StatFunc::gasdev_seq(idum, vec, size, 0.0, var);
        for (i = 0; i < size; i++) ave += vec[i];
        ave /= n;
        for (i = 0; i < size; i++) vec[i] -= ave;
        for (i = 0, s_square = 0.0; i < size; i++) s_square += vec[i] * vec[i];
        s_square = s_square / (n - 1.0);
    }
    double c = sqrt(var / s_square);
    for (i = 0; i < size; i++) vec[i] *= c;
}

double StatFunc::gasdev(int &idum) {
    static int iset = 0;
    static double gset;
    double fac, rsq, v1, v2;

    if (idum < 0) iset = 0;
    if (iset == 0) {
        do {
            v1 = 2.0 * ran1(idum) - 1.0;
            v2 = 2.0 * ran1(idum) - 1.0;
            rsq = v1 * v1 + v2*v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1*fac;
        iset = 1;
        return v2*fac;
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

double StatFunc::chidev(int &idum, const double df) {
    if (df > 2.0) return 2.0 * cheng_gamdev(idum, df * 0.5);
    else LOGGER.e(0, "invalid degree of freedom. StatFunc::chidev");
    return 0;
}

double StatFunc::cheng_gamdev(int &idum, const double alpha) {
    double u1 = 0.0, u2 = 0.0, nu = 0.0, d_buf = 0.0;
    double beta = (alpha - 1.0);
    while (1) {
        u1 = StatFunc::UniformDev(0.0, 1.0, idum);
        u2 = StatFunc::UniformDev(0.0, 1.0, idum);
        ;
        nu = (alpha - 1.0 / (6.0 * alpha)) * u1;
        nu /= (beta * u2);
        d_buf = 2.0 * (u2 - 1.0) / beta + nu + 1.0 / nu;
        if (d_buf < 2.0) return beta * nu;
        else {
            d_buf = 2.0 * log(u2) / beta - log(nu) + nu;
            if (d_buf < 1.0) return beta * nu;
        }
    }
}

int StatFunc::RandAbs(int a, int b, int &seed) { //a,bΪ������Ҷ˵㣬seedΪ����
    int rand;
    int stk = b - a + 1;
    int stl = 2;
    while (stl < stk) stl = stl + stl;
    int modul = 4 * stl;
    stk = seed;
    int sti = 1;
    while (sti <= 1) {
        stk = 5 * stk;
        stk = stk % modul;
        stl = stk / 4 + a;
        if (stl <= b) {
            rand = stl;
            sti = sti + 1;
        }
    }
    seed = stk;
    return (rand);
}

///////// Random Number Generation Functions End ////////

double StatFunc::chi_val(double df, double prob) {
    return boost::math::quantile(
        boost::math::complement(boost::math::chi_squared(df), prob));
}

double StatFunc::t_val(double df, double prob) {
    return boost::math::quantile(
        boost::math::complement(boost::math::students_t(df), prob));
}

double StatFunc::F_val(double df_1, double df_2, double prob) {
    return boost::math::quantile(
        boost::math::complement(boost::math::fisher_f(df_1, df_2), prob));
}

double StatFunc::ControlFDR(const vector<double> &P_Value, double alpha, bool Restrict) {
    int i = 0, Size = P_Value.size();
    if (Size <= 1) LOGGER.e(0, "invalid size. StatFunc::ControlFDR");
    double FDR_Threshold = 0.0;
    vector<double> P_ValueBuf(Size);
    for (i = 0; i < Size; i++) P_ValueBuf[i] = P_Value[i];
    stable_sort(P_ValueBuf.begin(), P_ValueBuf.end());
    for (i = 0; i < Size; i++) {
        double d_Buf = ((double) (i + 1)) * alpha / ((double) Size);
        if (i == Size - 1) {
            if (P_ValueBuf[i] <= d_Buf || Restrict) FDR_Threshold = d_Buf;
            else FDR_Threshold = -1.0;
            break;
        }
        if (P_ValueBuf[i] <= d_Buf && P_ValueBuf[i + 1] > d_Buf) {
            FDR_Threshold = d_Buf;
            break;
        } else if (i == Size - 1) FDR_Threshold = -1.0;
    }
    return FDR_Threshold;
}

double StatFunc::ControlFDR_Zou(const vector<double> &GenePValue, double FDR) {
    int i = 0, Size = GenePValue.size();
    double PValue = 0.0;
    vector<double> Temp(Size);
    for (i = 0; i < Size; i++) Temp[i] = GenePValue[i];
    sort(Temp.begin(), Temp.end());
    for (i = 0; i < Size; i++) {
        if (Temp[i] <= (i + 1) * FDR / Size && Temp[i + 1]>(i + 1) * FDR / Size) {
            PValue = (i + 1) * FDR / Size;
            break;
        } else if (Temp[i] <= (i + 1) * FDR / Size && i == Size - 1) {
            PValue = (i + 1) * FDR / Size;
            break;
        } else if (i == Size - 1) {
            PValue = FDR;
        }
    }
    return PValue;
}

double StatFunc::ControlFDR_Storey(vector<double> &P_Value, vector<double> &Q_Value, double CrtQ, double &FDR) {
    if (P_Value.empty()) return CrtQ;

    // Initialize Lambda
    double dBuf = 0.0;
    vector<double> Lambda;
    for (dBuf = 0.0; dBuf < 0.9001; dBuf += 0.05) Lambda.push_back(dBuf);

    // Calculate Pi0
    stable_sort(P_Value.begin(), P_Value.end());
    double Pi0 = CalcuPi0(P_Value, Lambda);

    // Calculate q-value
    int i = 0, m = P_Value.size();
    double m0 = Pi0 * (double) m;
    Q_Value.clear();
    Q_Value.resize(m);
    Q_Value[m - 1] = Pi0 * P_Value[m - 1];
    for (i = m - 2; i >= 0; i--) {
        Q_Value[i] = (m0 * P_Value[i]) / (double) (i + 1);
        if (Q_Value[i] > Q_Value[i + 1]) Q_Value[i] = Q_Value[i + 1];
    }

    // Calculate FDR-adjusted critical p-value
    bool Flag = false;
    double CrtVal = 1.0;
    for (i = 0; i < m - 1; i++) {
        if (Q_Value[i] < CrtQ && Q_Value[i + 1] > CrtQ) {
            CrtVal = 0.5 * (P_Value[i] + P_Value[i + 1]);
            Flag = true;
            break;
        }
    }
    if (!Flag) CrtVal = 0.5 * (P_Value[0] + P_Value[1]);

    // Estimate FDR
    int iBuf = 0;
    for (i = 0; i < m; i++) {
        if (P_Value[i] < CrtVal) iBuf++;
    }
    FDR = (CrtVal * m0) / (double) iBuf;
    if (FDR > 1.0) FDR = 1.0;

    return CrtVal;
}

double StatFunc::CalcuPi0(vector<double> &P_Value, vector<double> &Lambda) {
    int i = 0, j = 0, P_Size = P_Value.size(), LambdaSize = Lambda.size();
    vector<double> Pi(LambdaSize), NewPi(LambdaSize);
    for (i = 0; i < LambdaSize; i++) {
        int Count = 0;
        for (j = 0; j < P_Size; j++) {
            if (P_Value[j] > Lambda[i]) Count++;
        }
        Pi[i] = Count / (P_Size * (1.0 - Lambda[i]));
    }
    double PrePi = 0.0;
    spline(Lambda, Pi, 0.0, 0.0, NewPi);
    splint(Lambda, Pi, NewPi, Lambda[LambdaSize - 1], PrePi);
    return PrePi < 1.0 ? PrePi : 1.0;
}

void StatFunc::spline(vector<double> &x, vector<double> &y, const double yp1, const double ypn, vector<double> &y2) {
    int i, k;
    double p, qn, sig, un;

    int n = y2.size();
    vector<double> u(n - 1);
    if (yp1 > 0.99e30)
        y2[0] = u[0] = 0.0;
    else {
        y2[0] = -0.5;
        u[0] = (3.0 / (x[1] - x[0]))*((y[1] - y[0]) / (x[1] - x[0]) - yp1);
    }
    for (i = 1; i < n - 1; i++) {
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
        p = sig * y2[i - 1] + 2.0;
        y2[i] = (sig - 1.0) / p;
        u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
        u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
    }
    if (ypn > 0.99e30)
        qn = un = 0.0;
    else {
        qn = 0.5;
        un = (3.0 / (x[n - 1] - x[n - 2]))*(ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
    }
    y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);
    for (k = n - 2; k >= 0; k--)
        y2[k] = y2[k] * y2[k + 1] + u[k];
}

void StatFunc::splint(vector<double> &xa, vector<double> &ya, vector<double> &y2a, const double x, double &y) {
    int k;
    double h, b, a;

    int n = xa.size();
    int klo = 0;
    int khi = n - 1;
    while (khi - klo > 1) {
        k = (khi + klo) >> 1;
        if (xa[k] > x) khi = k;
        else klo = k;
    }
    h = xa[khi] - xa[klo];
    if (h == 0.0) LOGGER.e(0, "bad xa input to routine splint");
    a = (xa[khi] - x) / h;
    b = (x - xa[klo]) / h;
    y = a * ya[klo] + b * ya[khi]+((a * a * a - a) * y2a[klo]
            +(b * b * b - b) * y2a[khi])*(h * h) / 6.0;
}

// Default: upper-tail P(Z > x)
double StatFunc::pnorm(double x) {
    return boost::math::cdf(
        boost::math::complement(boost::math::normal(), x));
}

double StatFunc::dnorm(double x) {
    return boost::math::pdf(boost::math::normal(), x);
}

// qnorm(p, upper=true): returns Phi^{-1}(p), the lower-tail normal quantile.
// The "upper" parameter controls how p is interpreted when upper=false:
//   upper=false: finds x such that P(Z > x) = p  (upper-tail quantile)
double StatFunc::qnorm(double p, bool upper) {
    if (upper)
        return boost::math::quantile(boost::math::normal(), p);
    return boost::math::quantile(
        boost::math::complement(boost::math::normal(), p));
}


    int i = 0, n = p_value.size();

    double c = 0.0, min_val = 1.0;
    vector<double> fdr(n);
    vector<pair<double, int>> pval_buf(n);
    vector<pair<int, int>> indx_buf(n);

    for( i = 0; i < n; i++ ) pval_buf[i] = make_pair(p_value[i], i);
    stable_sort(pval_buf.begin(), pval_buf.end(), [](const pair<double,int> a, const pair<double,int> b) {return a.first > b.first; });

    for( i = 0; i < n; i++ ) {
        c = (double) n / (double) (n-i) * pval_buf[i].first;
        if(c < min_val) min_val = c;
        fdr[pval_buf[i].second] = CommFunc::Min(1.0, min_val);
    }

    return (fdr);
}

// Default: upper-tail P(Z > x)
double StatFunc::pnorm(double x) {
    return boost::math::cdf(
        boost::math::complement(boost::math::normal(), x));
}

double StatFunc::dnorm(double x) {
    return boost::math::pdf(boost::math::normal(), x);
}

// qnorm(p, upper=true): returns Phi^{-1}(p), the lower-tail normal quantile.
// When upper=false: finds x such that P(Z > x) = p (upper-tail quantile).
double StatFunc::qnorm(double p, bool upper) {
    if (upper)
        return boost::math::quantile(boost::math::normal(), p);
    return boost::math::quantile(
        boost::math::complement(boost::math::normal(), p));
}

double StatFunc::pchisq(double x, double df, bool log_p) {
    if (x < 0) return log_p ? std::numeric_limits<double>::quiet_NaN() : -9;
    return chi_prob(df, x, log_p);
}

double StatFunc::qchisq(double q, double df) {
    if (q < 0) return -9;
    if (q >= 1) return 0;
    return boost::math::quantile(
        boost::math::complement(boost::math::chi_squared(df), q));
}

//#################
// functions to calculate pchisqsum 

double StatFunc::pchisqsum(double x, VectorXd lambda) {
    double pval = psadd(x, lambda);
    if (pval > 1.0) pval = psatt(x, lambda);
    return pval;
}

double StatFunc::psadd(double x, VectorXd lambda) {
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

double StatFunc::psatt(double x, VectorXd lambda) {
    double sum=lambda.sum();
    if(CommFunc::FloatEqual(sum,0.0)) return 2.0;
  
    double sq_sum=lambda.dot(lambda);
    double sum_sq=sum*sum;
    double a=sq_sum/sum;
    double b=sum_sq/sq_sum;
    
    return pchisq(x/a, b);
}

double StatFunc::K(double zeta, VectorXd &lambda) {
    return -0.5*(1.0 - 2.0 * zeta * lambda.array()).log().sum();
}

double StatFunc::Kp(double zeta, VectorXd &lambda) {
    return (lambda.array() / (1.0 - 2.0 * zeta * lambda.array())).sum();
}

double StatFunc::Kpp(double zeta, VectorXd &lambda) {
    return 2.0 * (lambda.array().square() / (1.0 - 2.0 * zeta * lambda.array()).array().square()).sum();
}

double StatFunc::Kp_min_x(double zeta, VectorXd &lambda, double x) {
    return Kp(zeta, lambda) - x;
}

double StatFunc::Brents_Kp_min_x(VectorXd &lambda, double x, double lowerLimit, double upperLimit, double errorTol) {
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


vector<size_t> StatFunc::sort_re_index(const vector<double> &x) {
  vector<size_t> x_index(x.size());
  std::iota(x_index.begin(), x_index.end(), 0);

  sort(x_index.begin(), x_index.end(), [&x](size_t index1, size_t index2) {
          return x[index1] > x[index2];});

  return x_index;
}
