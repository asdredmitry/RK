#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
const double EPS = 1e-50;
const double H_MIN =  __DBL_MAX__;

const unsigned int FULL_NORM = 1;
const unsigned int STEPS_COUNT = 1 << 1;
const unsigned int LAST_NORM = 1 << 2;
const unsigned int NO_WRITE = 1 << 3;
const unsigned int GLOB_ERR = 1 << 4;

const double fac = 0.8;
const double facmin = 0.0;
const double facmax = 1.5;

const double c[13] = {0., 1./18., 1./12., 1./8., 5./16., 3./8., 59./400., 93./200., 5490023248./9719169821., 13./20., 1201146811./1299019798., 1., 1.};
const double a[13][12] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
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
const double b[13] = {14005451./335480064., 0., 0., 0., 0., -59238493./1068277825., 181606767./758867731., 561292985./797845732., -1041891430./1371343529., 760417239./1151165299., 118820643./751138087., -528747749./2220607170., 1./4.};
const double b_[13] = {13451932./455176623., 0., 0., 0., 0., -808719846./976000145., 1757004468./5645159321., 656045339./265891186., -3867574721./1518517206., 465885868./322736535., 53011238./667516719., 2./45., 0.};

typedef struct 
{
    double* data;
    int alloc;
    int used;
}vector;

void init(vector* v, int n);
void reall(vector* v);
void push_back(vector* v, double val);
void free_vec(vector* v);

double max(double x1, double x2);
double min(double x1, double x2);
double fmax(double x1, double x2);
double fmin(double x1, double x2);

double f1(double t, double y1, double y2);
double f2(double t, double y1, double y2);

double check_sol(double t);

double get_h(double h, double err, double tol);
double norm(double x1, double y1);

void dorman_prince(double h, double y1, double y2, double t, double* res);
void solve_dp(double l, double r, double y1, double y2, double tol, unsigned int attr, const char* name);

int main(void)
{
    solve_dp(0, 10, 0, 1, 1e-10, STEPS_COUNT | FULL_NORM | LAST_NORM | NO_WRITE, "data.dat");
    return EXIT_SUCCESS;
}

void init(vector* v, int n)
{
    if(v == NULL)
    {
        printf("invalid pointer");
        exit(EXIT_FAILURE);
    }
    v->data = NULL;
    v->data = (double*)malloc(n*sizeof(double));
    if(v->data == NULL)
    {
        printf("Cannot alloc mem");
        exit(EXIT_FAILURE);
    }
    v->used = 0;
    v->alloc = n;
}
void reall(vector* v)
{
    if(v == NULL)
    {
        printf("invalid pointer");
        exit(EXIT_FAILURE);
    }
    v->alloc *= 2;
    v->data = (double*)realloc(v->data, sizeof(double)*(v->alloc*2));
    if(v->data == NULL)
    {
        printf("Cannot realloc mem");
        exit(EXIT_FAILURE);
    }
}

void push_back(vector* v, double val)
{
    if(v == NULL)
    {
        printf("invalid pointer");
        exit(EXIT_FAILURE);
    }
    if(v->alloc == v->used)
        reall(v);
    v->data[v->used] = val;
    v->used += 1;
}

void free_vec(vector* v)
{
    if(v == NULL)
    {
        printf("invalid pointer");
        exit(EXIT_FAILURE);
    }
    if(v->data != NULL)
        free(v->data);
    v->data = NULL;
    v->used = 0;
    v->alloc = 0;    
}

double max(double x1, double x2)
{
    return (x1 > x2) ? x1 : x2;
}

double min(double x1, double x2)
{
    return (x1 < x2) ? x1 : x2;
}

double fmin(double x1, double x2)
{
    return min(fabs(x1), fabs(x2));
}

double fmax(double x1, double x2)
{
    return max(fabs(x1), fabs(x2));
}

double f1(double t, double y1, double y2)
{
    return y2;
    return y1;
    return t;
}

double f2(double t, double y1, double y2)
{
    return -y1;
    return y2;
    return t;
}

double check_sol(double t)
{
    return sin(t);
}

double norm(double x1, double x2)
{
    return fmax(x1, x2);
}

double get_h(double h, double err, double tol)
{
    if(err < EPS)
        err = 2*EPS;
    h *= min(facmax, max(facmin, fac*pow(tol/err, 1./7.)));
    if(h < EPS)
        h = 2*EPS;
    return h;
} 

void dorman_prince(double h, double y1, double y2, double t, double* res)
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
    res[0] = tmpy1;
    res[1] = tmpy2;

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
    res[2] = tmpy1;
    res[3] = tmpy2;
}

void solve_dp(double l, double r, double y1, double y2, double tol, unsigned int attr, const char* name)
{
    double h, err;
    FILE * output;
    double res[4];
    int steps = 0;
    double last = 0.;
    double max_norm = 0;
    h = 0.1;
    if(!(attr & NO_WRITE))
    {
        output = fopen(name, "w");
        if(output == NULL)
        {
            printf("Cannot open file to write");
            exit(EXIT_FAILURE);
        }
        fprintf(output, "%e %e %e \n", l, y1, y2);
    }
    while(l < r)
    {
        if(l + h > r)
            h = r - l;
        dorman_prince(h, y1, y2, l, res);
        err = norm(res[0] - res[2], res[1] - res[3]);
        if(err < tol)
        {
            l += h;
            if(!(attr & NO_WRITE))
            {
                fprintf(output, "%e %e %e \n", l, res[2], res[3]);
            }
            h = get_h(h, err, tol);
            y1 = res[2];
            y2 = res[3];
            if(attr & STEPS_COUNT)
                steps ++;
            if(attr & FULL_NORM)
                max_norm = max(max_norm, fabs(check_sol(l) - y1));
            last = y1;
        }
        else 
            h = get_h(h, err, tol);
    }
    if(!(attr & NO_WRITE))
        fclose(output);
    if(attr & FULL_NORM)
        printf("Norm - %e \n", max_norm);
    if(attr & STEPS_COUNT)
        printf("Steps - %d \n", steps);
    if(attr & LAST_NORM)
        printf("Last acc - %e \n", fabs(last - check_sol(l)));
}


