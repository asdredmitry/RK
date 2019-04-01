#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <time.h>
const double EPS = 1e-308;
double c[13] = {0., 1./18., 1./12., 1./8., 5./16., 3./8., 59./400., 93./200., 5490023248./9719169821., 13./20., 1201146811./1299019798., 1., 1.};
double a[13][12] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                    1./18., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                    1./48., 1./16., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                    1./32., 0., 3./32., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                    5./16., 0., -75./64, 75./64., 0., 0., 0., 0., 0., 0., 0., 0.,
                    3./80., 0., 0., 3./16., 3./20., 0., 0., 0., 0., 0., 0., 0.,
                    29443841./614563906., 0., 0., 77736538./692538347., -28693883./1125000000., 23124283./1800000000., 0., 0., 0., 0., 0., 0.,
                    16016141./946692911., 0., 0., 61564180./158732637., 22789713./633445777., 545815736./2771057229., -180193667./1043307555., 0., 0., 0., 0., 0.,
                    39632708./573591083., 0., 0., -433636366./683701615., -421739975./2616292301., 100302831./723423059., 790204164./839813087., 800635310./3783071287., 0., 0., 0., 0.,
                    246121993./1340847787., 0., 0., -37695042795./15268766246., -309121744./1061227803., -12992083./490766935., 6005943493./2108947869., 393006217./1396673457., 123872331./1001029789., 0., 0., 0.,
                    -1028468189./846180014., 0., 0., 8478235783./508512852., 1311729495./1432422823., -10304129995./1701304382., -48777925059./3047939560., 15336726248./1032824649., -45442868181./3398467696., 3065993473./597172653., 0., 0.,
                    185892177./718116043., 0., 0., -3185094517./667107341., -477755414./1098053517., -703635378./230739211., 5731566787./1027545527, 5232866602./850066563., -4093664535./808688257., 3962137247./1805957418., 65686358./487910083., 0,
                    403863854./491063109., 0., 0., -5068492393./434740067., -411421997./543043805., 652783627./914296604., 11173962825./925320556., -13158990841./6184727034., 3936647629./1978049680., -160528059./685178525., 248638103./1413531060., 0.};
double b[13] = {14005451./335480064., 0., 0., 0., 0., -59238493./1068277825., 181606767./758867731., 561292985./797845732., -1041891430./1371343529., 760417239./1151165299., 118820643./751138087., -528747749./2220607170., 1./4.};
double b_[13] = {13451932./455176623., 0., 0., 0., 0., -808719846./976000145., 1757004468./5645159321., 656045339./265891186., -3867574721./1518517206., 465885868./322736535., 53011238./667516719., 2./45., 0.};

double c1[7] = {0., 0.5, 2.0/3.0,1.0/3.0, 5.0/6.0, 1.0/6.0, 1.0};
double b1[7] = {13.0/200.0, 0.0, 11.0/40.0, 11.0/40.0, 4.0/25.0 ,4.0/25.0 ,13.0/200.0};
double a1[7][6] = {0., 0., 0., 0., 0., 0.,
                   1./2., 0., 0., 0., 0., 0.,
                   2./9., 4./9., 0., 0., 0., 0.,
                   7./36., 2./9., -1./12, 0., 0., 0.,
                   -35./144., -55./36., 35./48., 15./8., 0., 0., 
                   -1./360, -11./36., -1./8., 1./2., 1./10., 0.,
                   -41./260., 22./13., 43./156., -118./39., 32./195., 80./39.};
typedef struct 
{
    double * data;
    int alloc;
    int used;
}vector;

void init(vector * v, int n)
{
    v->data = (double *)malloc(sizeof(double)*n);
    v->used = 0;
    v->alloc = n;
}
void reall(vector * v)
{
    v->alloc *= 2;
    v->data = (double *)realloc(v->data, sizeof(double)*(v->alloc));
}
void push_back(vector * v, double  val)
{
    if(v->alloc == v->used)
        reall(v);
    v->data[v->used] = val;
    v->used += 1;
}
void free_vec(vector * v)
{
    if(v->data != NULL)
        free(v->data);
    v->used = 0;
    v->alloc = 0;
    v->data = NULL;
}

double f1(double t, double y1, double y2)
{
    return y2;
}
double f2(double t, double y1, double y2)
{
    return -4*y1;
    //return -(1 + 0.2*y1*y1)*y1 + cos(t);
}
double check_sol(double t)
{
    return sin(2*t);
}


double max(double x1, double x2)
{
    return (x1 > x2) ? x1 : x2;
}
double min(double x1, double x2)
{
    return (x1 > x2) ? x2 : x1;
}
double fmax(double x1, double x2)
{
    return max(fabs(x1), fabs(x2));
}
double fmin(double x1, double x2)
{
    return min(fabs(x1), fabs(x2));
}


double get_h(double h, double err, double tol)
{
    double fac, facmax, facmin;
    fac = 0.8;
    facmax = 1.5;
    facmin = 0.0000001;
    if(err < EPS)
        err = 2*EPS;
    h *= min(facmax, max(facmin, fac*pow(tol/err, 1./7.)));
    if(h < EPS)
        h = 2*EPS;
    return h;
}
double norm(double x1, double y1, double x2, double y2)
{
    return fmax(x1 - x2, y1 - y2);
}


void dorman_prince(double h, double y1, double y2, double t, double * resy1, double * resy2, double * resy1_, double * resy2_)
{
    double tmpt, tmpy1, tmpy2;
    double k1[13];
    double k2[13];
    int i, j;
    k1[0] = f1(t, y1, y2);
    k2[0] = f2(t, y1, y2);
    for(i = 1; i < 13; i++)
    {
        tmpy1 = y1;
        tmpy2 = y2;
        tmpt = t + c[i]*h;
        for(j = 0; j < i; j++)
        {
            tmpy1 += a[i][j]*h*k1[j];
            tmpy2 += a[i][j]*h*k2[j];
        }
        k1[i] = f1(tmpt, tmpy1, tmpy2);
        k2[i] = f2(tmpt, tmpy1, tmpy2);
    }
    if(resy1 != NULL && resy2 != NULL)
    {
        tmpy1 = 0;
        tmpy2 = 0;
        for(i = 0; i < 13; i++)
        {
            tmpy1 += b[i]*k1[i];
            tmpy2 += b[i]*k2[i];
        }   
        tmpy1 *= h;
        tmpy2 *= h;
        tmpy1 += y1;
        tmpy2 += y2;
        *resy1 = tmpy1;
        *resy2 = tmpy2;
    }

    if(resy1_ != NULL && resy2_ != NULL)
    {
        tmpy1 = 0;
        tmpy2 = 0;
        for(i = 0; i < 13; i++)
        {
            tmpy1 += b_[i]*k1[i];
            tmpy2 += b_[i]*k2[i];
        }
        tmpy1 *= h;
        tmpy2 *= h;
        tmpy1 += y1;
        tmpy2 += y2;
        *resy1_ = tmpy1;
        *resy2_ = tmpy2;
    }
}
void runge_kutta(double h, double y1, double y2,double t, double * resy1, double * resy2)
{
    double tmpt, tmpy1, tmpy2;
    double k1[7];
    double k2[7];
    int i, j;
    k1[0] = f1(t, y1, y2);
    k2[0] = f2(t, y1, y2);
    for(i = 1; i < 7; i++)
    {
        tmpy1 = y1;
        tmpy2 = y2;
        tmpt = t + c1[i]*h;
        for(j = 0; j < i; j++)
        {
            tmpy1 += a1[i][j]*h*k1[j];
            tmpy2 += a1[i][j]*h*k2[j];
        }
        k1[i] = f1(tmpt, tmpy1, tmpy2);
        k2[i] = f2(tmpt, tmpy1, tmpy2);
    }
    tmpy1 = 0.;
    tmpy2 = 0.;
    for(i = 0;i < 7; i++)
    {
        tmpy1 += b1[i]*k1[i];
        tmpy2 += b1[i]*k2[i];
    }
    tmpy1 *= h;
    tmpy2 *= h;
    tmpy1 += y1;
    tmpy2 += y2;
    *resy1 = tmpy1;
    *resy2 = tmpy2;
}
void solve_dp(double l, double r, double y1, double y2, double tol, vector * t, vector * s1, vector * s2)
{
    double h, err;
    double resy1, resy2, resy1_, resy2_;
    init(t, 100);
    init(s1, 100);
    init(s2, 100);
    h = 0.1;
    push_back(t, l);
    push_back(s1, y1);
    push_back(s2, y2);
    while(l < r)
    {
        if(l + h > r)
            h = r - l;
        dorman_prince(h, y1, y2, l, &resy1, &resy2, &resy1_, &resy2_);
        //printf("%17g %17g \n", resy1 - resy1_, resy2 - resy2_);
        err = norm(resy1, resy2, resy1_, resy2_);
        if(err < tol)
        {
            //printf("%17g \n", h);
            l += h;
            push_back(t, l);
            push_back(s1, resy1_);
            push_back(s2, resy2_);
            y1 = resy1;
            y2 = resy2;
            h = get_h(h, err, tol);
        }
        else 
        {
            h = get_h(h, err, tol);
        }
    }
}
void solve_rk(double l, double r, double y1, double y2, double tol, vector * t, vector * s1, vector * s2, int accuracy)
{
    double h, err;
    double resy1, resy2;
    double resy1_, resy2_;
    double resy1_2h, resy2_2h;
    init(t, 100);
    init(s1, 100);
    init(s2, 100);
    h = 0.1;
    push_back(t, l);
    push_back(s1, y1);
    push_back(s2, y2);
    while(l < r)
    {
        if(l + h > r)
            h = r - l;
        if(accuracy == 6)
        {
            runge_kutta(h, y1, y2, l, &resy1, &resy2);
            runge_kutta(h, resy1, resy2, l + h, &resy1_, &resy2_);
            runge_kutta(2*h, y1, y2, l, &resy1_2h, &resy2_2h);
        }
        else if(accuracy == 8)
        {
            dorman_prince(h, y1, y2, l, &resy1, &resy2, NULL, NULL);
            dorman_prince(h, resy1, resy2, l + h, &resy1_, &resy2_, NULL, NULL);
            dorman_prince(2*h, y1, y2, l, &resy1_2h, &resy2_2h, NULL, NULL);
        }
        else if(accuracy == 7)
        {
            dorman_prince(h, y1, y2, l, NULL, NULL, &resy1, &resy2);
            dorman_prince(h, resy1, resy2, l + h, NULL, NULL, &resy1_, &resy2_);
            dorman_prince(2*h, y1, y2, l, NULL, NULL, &resy1_2h, &resy2_2h);
        }
        err = norm(resy1_2h, resy2_2h, resy1_, resy2_)/(pow(2, accuracy + 1) - 1);
        if(err < tol)
        {
            l += h;
            push_back(t, l);
            push_back(s1, resy1);
            push_back(s2, resy2);
            l += h;
            push_back(t, l);
            push_back(s1, resy1_);
            push_back(s2, resy2_);
            h = get_h(h, err, tol);
            y1 = resy1_;
            y2 = resy2_;
        }
        else 
        {
            h = get_h(h, err, tol);
        }
    }
}
void swap(double * x, double * y)
{
    double tmp = x[0];
    x[0] = y[0];
    y[0] = tmp;
}
double norm_full(vector * t, vector * s1)
{
    int i;
    double max = 0.;
    for(i = 0; i < t->used; i++)
    {
        max = fmax(s1->data[i] - check_sol(t->data[i]), max);
    }
    return max;
}
void write_data(vector * t, vector * s1, vector * s2, const char * name)
{
    int i;
    FILE * output = fopen(name, "w");
    for(i = 0; i < t->used; i++)
    {
        fprintf(output, "%17g %17g %17g \n", t->data[i], s1->data[i], s2->data[i]);
    }
    fclose(output);
}
void find_zero(double t1, double t2, double y11, double y12, double y21, double y22, double eps,double * rest)
{
    while(fabs(y22) > eps)
    {
        t2 = t2 - y22*(t2 - t1)/(y22 - y12);
        dorman_prince(t2 - t1, y11, y12, t1, &y21, &y22, NULL, NULL);
    }
    *rest = t2;
}
double find_period(double l, double r_max, double y1, double y2)
{
    double h, err, tol;
    double resy1, resy2, resy1_, resy2_;
    double zt, zy1, zy2;
    vector zero_t;
    vector zero_y1;
    vector zero_y2;
    init(&zero_t, 100);
    init(&zero_y1, 100);
    init(&zero_y2, 100);
    h = 0.1;
    tol = 1e-9;
    while(l < r_max)
    {
        //if(counter % 1000000 == 0)
        //    printf("%lf \n", l);
        dorman_prince(h, y1, y2, l, &resy1, &resy2, &resy1_, &resy2_);
        //printf("%17g %17g \n", resy1 - resy1_, resy2 - resy2_);
        err = norm(resy1, resy2, resy1_, resy2_);
        if(err < tol)
        {
            l += h;
            if(y2*resy2 < 0 && y2 > resy2)
            {
                find_zero(l - h, l, y1, y2, resy1, resy2, 1e-8, &zt);
                dorman_prince(zt - (l - h), y1, y2, (l - h), &zy1, &zy2, NULL, NULL);
                push_back(&zero_t, zt);
                push_back(&zero_y1, zy1);
                push_back(&zero_y2, zy2);
                if(zero_t.used > 1 && norm(zero_y1.data[0], zero_y2.data[0], zero_y1.data[zero_y1.used - 1], zero_y2.data[zero_y2.used - 1]) < 1e-3)
                {
                    goto out;
                    //printf(" \n %lf %lf \n", zero_t.data[0], zero_t.data[zero_t.used - 1]);
                    //return zero_t.data[0] - zero_t.data[zero_t.used - 1];
                }
            }
            h = get_h(h, err, tol);
            h = min(h, 0.1);
            y1 = resy1;
            y2 = resy2;
        }
        else 
        {
            h = get_h(h, err, tol);
            h = min(h, 0.1);
        }
    }
    write_data(&zero_t, &zero_y1, &zero_y2, "zeroes.dat");
    return -1;
    out:
    write_data(&zero_t, &zero_y1, &zero_y2, "zeroes.dat");
    return (-zero_t.data[0] + zero_t.data[zero_t.used - 1]);
}
void print_statistics(vector * t, vector * s1, vector * s2, double tol)
{
    double h;
    int i;
    double max_h = 0;
    double min_h = __DBL_MAX__;
    double e_h = 0;
    for(i = 1; i < t->used; i++)
    {  
        h = t->data[i] - t->data[i - 1];  
        e_h += h;
        min_h = min(min_h, h);
        max_h = max(max_h, h);        
    }
    e_h /= (t->used);
    printf("Method Statistics: \n");
    printf("[%lf, %lf] tol - %17g \n", t->data[0], t->data[t->used - 1], tol);
    printf("Min h - %17g ; Max h - %17g ; Mid - %17g \n", min_h, max_h, e_h);
    printf("Amount - %d \n", t->used);
}
void compare_methods(double l, double r, double tol, double y1, double y2)
{
    vector t_rk6, s1_rk6, s2_rk6;
    vector t_rk8, s1_rk8, s2_rk8;
    vector t_rk7, s1_rk7, s2_rk7;
    vector t_dp, s1_dp, s2_dp;
    r *= M_PI;
    solve_rk(l, r, y1, y2, tol, &t_rk8, &s1_rk8, &s2_rk8, 8);
    solve_dp(l, r, y1, y2, tol, &t_dp, &s1_dp, &s2_dp);
    solve_rk(l, r, y1, y2, tol, &t_rk6, &s1_rk6, &s2_rk6, 6);
    solve_rk(l, r, y1, y2, tol, &t_rk7, &s1_rk7, &s2_rk7, 7);
    write_data(&t_dp, &s1_dp, &s2_dp, "dormance_prince.dat");
    write_data(&t_rk6, &s1_rk6, &s2_rk6, "runge_kutta6.dat");
    write_data(&t_rk7, &s1_rk7, &s2_rk7, "runge_kutta7.dat");
    write_data(&t_rk8, &s1_rk8, &s2_rk8, "runge_kutta8.dat");
    printf("---------------------------------------------------------------- \n");
    printf("Runge Kutta 6 order: \n");
    print_statistics(&t_rk6, &s1_rk6, &s2_rk6, tol);
    printf("Global mistake : %17g \n", norm_full(&t_rk6, &s1_rk6));
    printf("---------------------------------------------------------------- \n");
    printf("Runge Kutta 7 order: \n");
    print_statistics(&t_rk7, &s1_rk7, &s2_rk7, tol);
    printf("Global mistake : %17g \n", norm_full(&t_rk7, &s1_rk7));
    printf("----------------------------------------------------------------\n");
    printf("Runge Kutta 8 order: \n");
    print_statistics(&t_rk8, &s1_rk8, &s2_rk8, tol);
    printf("Global mistake : %17g \n", norm_full(&t_rk8, &s1_rk8));
    printf("---------------------------------------------------------------- \n");
    printf("Dorman Prince 7 order: \n");
    print_statistics(&t_dp, &s1_dp, &s2_dp, tol);
    printf("Global mistake : %17g \n", norm_full(&t_dp, &s1_dp));
    printf("----------------------------------------------------------------\n");
    free_vec(&t_dp);
    free_vec(&s1_dp);
    free_vec(&s2_dp);
    free_vec(&t_rk6);
    free_vec(&s1_rk6);
    free_vec(&s2_rk6);
    free_vec(&t_rk7);
    free_vec(&s1_rk7);
    free_vec(&s2_rk7);
    free_vec(&t_rk8);
    free_vec(&s1_rk8);
    free_vec(&s2_rk8);
}
int main(void)
{
    double l, r, y1, y2;
    double tol;
    tol = 1e-9;
    printf("Enter l, r, y1, y2");
    scanf("%lf %lf %lf %lf", &l, &r, &y1, &y2);
    //compare_methods(l, r, tol, y1, y2);
    //printf("%lf ", find_period(l, r, y1, y2));
    return 0;
}