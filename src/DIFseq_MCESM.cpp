#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h> // pow, sqrt, lgamma
#include <cmath> //
#include "omprng.h"
#include <chrono> // time
#include <unistd.h> // getopt 
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <sys/stat.h>  //mkdir
#include <algorithm>    // sort
#include <vector>  //
#include <gsl/gsl_sf_psi.h> // trigamma and digamma
#include <gsl/gsl_sort.h>

using namespace std;

// Adam tuning parameters
// const double ALPHA = 0.001;
// const double BETA1 = 0.9;
// const double BETA2 = 0.999;
// const double EPSILON = 1e-8;

int rand_cate(double* _prop, omprng _rng){
    int res = 0;
    double u = _rng.runif();
    while(u > _prop[res]){
        u = u - _prop[res];
        res++;
    }
    return res;
}

int rand_Ber(double _prob, omprng _rng){
    int res = 1;
    double u = _rng.runif();
    if(u > _prob){
        res = 0;
    }
    return res;
}


int rand_NB(double _r, double _mu, omprng _rng){
    
    double lambda, nu, p, u;
    double STEP = 500;// steps for large lambda in Poisson distribution
    int res = 0;

    // sample lambda from gamma distribution
    lambda = _rng.rgamma(_r, _mu/_r);
        
    // sample x from Poisson distribution
    nu = lambda;
    p = 1.0;

    do{
        res ++;
        u = _rng.runif();
        p = p * u;
        while(p < 1 && nu > 0){
            if(nu > STEP){
                p = p * exp(STEP);
                nu = nu - STEP;
            }else{
                p = p * exp(nu);
                nu = 0;
            }
        }
    }while(p > 1);

    res--;
    return res;
}

void rand_Dir(double *_xi, int _K, omprng _rng, double *_pi){

    double *rn_gam = new double[_K];
    double sum_gam = 0.0;
    for(int k = 0; k < _K; k++){
        rn_gam[k] = _rng.rgamma(_xi[k], 1.0);
        sum_gam += rn_gam[k];
    }
    for(int k = 0; k < _K; k++){
        _pi[k] = rn_gam[k] / sum_gam;
    }

    delete [] rn_gam;
}

void _update_logmu(int _N, int _K,
    int* _t_infor, int* _b_infor,
    int* _W, double _alpha, double* _beta, double* _eta, double* _nu, double* _delta,//parameter
    double* _logmu) {

    int k, t, b;
    for (int i = 0; i < _N; i++) {
        k = _W[i];
        t = _t_infor[i];
        b = _b_infor[i];
        _logmu[i] = _alpha + _beta[k] + _eta[t * _K + k] + _nu[b] + _delta[i];
    }
}

void _update_zx(int _N, int* _b_infor,
    double** _gamma, double* _phi, double* _logmu,
    int* _Y, omprng _rng,
    int* _X, int* _Z) {

    int b;
    double log_rat, u, temp_max, mu_temp, prop_mean, log_acc_rat, logmu_ratio;
    int temp_x;


    for (int i = 0; i < _N; i++) {

        b = _b_infor[i];

        if (_Y[i] == 0) {
            //update z_{big}
            if (_X[i] == 0) {

                // printf("Sampling for Z.\n");

                log_rat = _gamma[b][0];
                _Z[i] = rand_Ber(1 / (1 + exp(-log_rat)), _rng);

            }
            else {

                _Z[i] = 1;

            }// end of if X == 0

            //update x_{big}
            if (_Z[i] == 1) {// Dropout event

                // printf("Sampling for X.\n");
                // MH sampling
                mu_temp = exp(_logmu[i]);
                prop_mean = _phi[b] * mu_temp * exp(_gamma[b][1]) / (mu_temp * (1 - exp(_gamma[b][1])) + _phi[b]);

                if (log(prop_mean) > 10.0) {
                    temp_x = rand_NB(_phi[b], exp(10.0), _rng);
                    if (temp_x > exp(15.0) || temp_x < 0) { // avoid the overflow of temp_x
                        temp_x = exp(15.0);
                    }
                    
                    logmu_ratio = _logmu[i] + _gamma[b][1] - log(mu_temp + _phi[b]) + log(exp(10.0) + _phi[b]) - 10.0;
                    log_acc_rat = (temp_x - _X[i]) * logmu_ratio;

                }
                else {
                    temp_x = rand_NB(_phi[b], prop_mean, _rng);
                    if (temp_x > exp(15.0) || temp_x < 0) { // avoid the overflow of temp_x
                        temp_x = exp(15.0);
                    }
                    log_acc_rat = 0.0;
                }

                log_acc_rat = log_acc_rat + log(1.0 + exp(_gamma[b][0] + _gamma[b][1] * _X[i])) - log(1.0 + exp(_gamma[b][0] + _gamma[b][1] * temp_x));



                // printf("Sampling for X.\n");
                u = _rng.runif();
                //make the exponential be compuational

                
                
                if (log(u) < log_acc_rat) {
                    _X[i] = temp_x;
                }
            }
            else { // No dropout
                _X[i] = 0;
            }// end of if Z == 0
        }// end of if Y == 0
    }// end of i for
}

void _E_w(int _b, int _t, int _p, int _G, int _K, 
    double* _pi, double* _alpha, double** _beta, double** _eta, double** _nu, double _delta, double** _phi,
    int* _Y,
    double* _u) {

    double* log_ratio = new double[_K];
    double max_logr, sum_explogr;

    for (int k = 0; k < _K; k++) {
        log_ratio[k] = log(_pi[k]);

        for (int g = 0; g < _G; g++) {
            log_ratio[k] = log_ratio[k] + (_beta[g][k] + _eta[g][_t * _K + k]) * _Y[g];
            log_ratio[k] = log_ratio[k] - log(exp(_alpha[g] + _beta[g][k] + _eta[g][_t * _K + k] + _nu[g][_b] + _delta) + _phi[g][_b]) * (_Y[g] + _phi[g][_b]);
        }

        if (k == 0) {
            max_logr = log_ratio[k];
        }
        else if(max_logr < log_ratio[k]){
            max_logr = log_ratio[k];
        }
    }

    sum_explogr = 0.0;
    for (int k = 0; k < _K; k++) {

        log_ratio[k] = log_ratio[k] - max_logr;
        sum_explogr = sum_explogr + exp(log_ratio[k]);

    }
  
    for (int k = 0; k < _K; k++) {
        _u[k] = exp(log_ratio[k]) / sum_explogr;
    }

    delete[] log_ratio;
}

void _E_L(int _K, // int _g,
    double* _beta, double _p, double _tau0, double _tau1,
    double* _PrL) {

    double ratio_L1, ratio_L0;
    for (int k = 1; k < _K; k++) {

        ratio_L1 = _p / sqrt(_tau1) * exp(-pow(_beta[k], 2.0) / 2.0 / _tau1);
        ratio_L0 = (1.0 - _p) / sqrt(_tau0) * exp(-pow(_beta[k], 2.0) / 2.0 / _tau0);
        _PrL[k] = ratio_L1 / (ratio_L1 + ratio_L0);

        /*
        if (_g == 0) {
            
            cout << "_tau1 = " << _tau1 << endl;
            cout << "_tau0 = " << _tau0 << endl;
            cout << "_beta[" << k << "] = " << _beta[k] << endl;
            cout << "ratio_L1 = " << ratio_L1 << endl;
            cout << "ratio_L0 = " << ratio_L0 << endl;
            cout << "_PrL["<< k << "] = " << _PrL[k] << endl;
        }
        */
    }
}


void _E_J(int _K, int _T,
    double* _eta, double _p, double _tau0, double _tau1,
    double* _PrJ) {

    double ratio_J1, ratio_J0;
    for (int t = 1; t < _T; t++) {
        for (int k = 0; k < _K; k++) {

            ratio_J1 = _p / sqrt(_tau1) * exp(-pow(_eta[t * _K + k], 2.0) / 2.0 / _tau1);
            ratio_J0 = (1 - _p) / sqrt(_tau0) * exp(-pow(_eta[t * _K + k], 2.0) / 2.0 / _tau0);
            _PrJ[t * _K + k] = ratio_J1 / (ratio_J1 + ratio_J0);

        }
    }

}

double _percentile_q(int _N, int _T, int _K, int _R, int* _b_infor, int* _t_infor,
    double _alpha, double* _beta, double* _eta, double* _nu, double* _delta, double* _phi,
    int** _w_mc, int** _x_mc,
    double* _PrJ, double _tau0, double _tau1, // prior for eta
    double _max_step, double _tol_nt, int _max_iter) {

    double err_nt;
    int iter_nt;
    double p_temp, fir_der, sec_der;
    int i, b, t, k, r;
    // update eta        
    double** sum_ux_eta = new double* [_T];
    for (t = 0; t < _T; t++) {
        sum_ux_eta[t] = new double[_K];
        for (k = 0; k < _K; k++) {
            sum_ux_eta[t][k] = 0.0;
        }
    }

    double eta_temp, eta_iter;

    // auxiliary variables
    double** mu_nocur = new double* [_N];
    for (i = 0; i < _N; i++) {

        mu_nocur[i] = new double[_K];

        // b = _b_infor[i];

        for (k = 0; k < _K; k++) {
            mu_nocur[i][k] = 0.0;
        }
    }

    // update mu_nocur
    for (i = 0; i < _N; i++) {
        b = _b_infor[i];
        for (k = 0; k < _K; k++) {
            mu_nocur[i][k] = exp(_alpha + _beta[k] + _nu[b] + _delta[i]); // not include eta!
        }
    }

    // calculate the auxiliary variable sum_{i} u_btik * y_btig    
    for (i = 0; i < _N; i++) {
        t = _t_infor[i];
        for (r = 0; r < _R; r++) {
            k = _w_mc[i][r];
            sum_ux_eta[t][k] += 1.0 * _x_mc[i][r] / _R;
        }
    }

    for (t = 1; t < _T; t++) {
        for (k = 0; k < _K; k++) {

            eta_iter = _eta[t * _K + k];

            // start the Newton's method
            err_nt = 1.0;
            iter_nt = 0;

            while (err_nt > _tol_nt && iter_nt < _max_iter) {

                // prior term
                fir_der = -eta_iter / _tau1 * _PrJ[t * _K + k] - eta_iter / _tau0 * (1.0 - _PrJ[t * _K + k]);
                sec_der = -1.0 / _tau1 * _PrJ[t * _K + k] - 1.0 / _tau0 * (1.0 - _PrJ[t * _K + k]);

                // log-likelihood term
                fir_der += sum_ux_eta[t][k];

                for (i = 0; i < _N; i++) {
                    if (_t_infor[i] == t) {

                        b = _b_infor[i];
                        for (r = 0; r < _R; r++) {
                            if (_w_mc[i][r] == k) {
                                p_temp = mu_nocur[i][k] / (mu_nocur[i][k] + _phi[b] / exp(eta_iter));
                                fir_der += -(_x_mc[i][r] + _phi[b]) * p_temp / _R;
                                sec_der += -_phi[b] * (_x_mc[i][r] + _phi[b]) * pow(p_temp, 2.0) / mu_nocur[i][k] / exp(eta_iter) / _R;
                            }
                        }
                    }
                }

                //if((_g == 1949 | _g == 2013) & (k == 2 | k == 3) ){
                //    cout << "The first derivative of eta["<<_g<< "]["<< t * _K + k << "] is " << fir_der << endl;
                //    cout << "The second derivative of eta["<<_g<< "]["<< t * _K + k << "] is " << sec_der << endl;
                //}

                double shift_eta = -fir_der / sec_der;

                if (shift_eta > _max_step) {
                    eta_temp = eta_iter + _max_step;
                }
                else if (shift_eta < -_max_step) {
                    eta_temp = eta_iter - _max_step;
                }
                else {
                    eta_temp = eta_iter + shift_eta;
                }
                // eta_temp = eta_iter - fir_der / sec_der;
                err_nt = abs(eta_temp - eta_iter);
                eta_iter = eta_temp;

                if (std::isnan(eta_iter)) {
                    cout << "t = " << t << endl;
                    cout << "k = " << k << endl;
                    cout << "p_temp = " << p_temp << endl;
                    cout << "fir_der = " << fir_der << endl;
                    cout << "sec_der = " << sec_der << endl;
                    exit(2);
                }

                iter_nt++;

                //if((_g == 1949 | _g == 2013) & (k == 2 | k == 3) ){
                //    cout << "At the "<< iter_nt  << "-th iteration, eta["<<_g<< "]["<< t * _K + k << "] = " << eta_iter << endl;
                //}
            }

            _eta[t * _K + k] = eta_iter;
        }
    }
    
    // Calculate the quasi Bayes factor
    double QBFg = 0.0;

    for (t = 1; t < _T; t++) {
        for (k = 0; k < _K; k++) {
            QBFg += -pow(_eta[t * _K + k], 2.0) / _tau1 / 2.0 * _PrJ[t * _K + k];
            QBFg += -pow(_eta[t * _K + k], 2.0) / _tau0 / 2.0 * (1.0 - _PrJ[t * _K + k]);
        }
    }

    for (i = 0; i < _N; i++) {
        t = _t_infor[i];
        b = _b_infor[i];
        for (r = 0; r < _R; r++) {

            k = _w_mc[i][r];
            QBFg += _eta[t * _K + k] * _x_mc[i][r] / _R;
            QBFg += -log(mu_nocur[i][k] * exp(_eta[t * _K + k]) + _phi[b]) * (_x_mc[i][r] + _phi[b])/_R;
            QBFg += log(mu_nocur[i][k] + _phi[b]) * (_x_mc[i][r] + _phi[b])/_R;
        }
    }

    // free memory
    for (i = 0; i < _N; i++) {

        delete[] mu_nocur[i];
    }

    delete[] mu_nocur;

    for (t = 0; t < _T; t++) {
        delete[] sum_ux_eta[t];
    }
    delete[] sum_ux_eta;

    return QBFg;
}

double _Stoc_Newton_abenp_GD(int _N, int _B, int _T, int _K, int _R, int* _b_infor, int* _t_infor,// int _g,
    double _alpha, double* _beta, double* _eta, double* _nu, double* _delta, double* _phi,
    int** _w_mc, int** _x_mc, int _control_ind,
    double* _nb, double* _nt, double* _nb_part, double* _nt_part,
    double _m_a, double _sigmasq_a,  // prior for alpha
    double _tau0, double _tau1, // prior for beta
    double* _PrL, double* _PrJ, // prior for eta
    double* _m_c, double _sigmasq_c, // prior for nu
    double* _phi_prior, // prior for phi
    int _update_stage, double _lr, double _max_step) {

    double p_temp, fir_der, sec_der;
    int i, b, t, k, r;
    // auxiliary variables
    double** mu = new double* [_N];
    for (i = 0; i < _N; i++) {
        b = _b_infor[i];
        t = _t_infor[i];
        mu[i] = new double[_K];
        for (k = 0; k < _K; k++) {
            mu[i][k] = exp(_alpha + _beta[k] + _eta[t * _K + k] + _nu[b] + _delta[i]);
        }
    }

    double N_tol = 0;
    for (int b = 0; b < _B; b++) {
        N_tol += _nb[b];
    }

    // update alpha
    double alpha_iter;

    // calculate the auxiliary variable sum_{i} y_btig
    alpha_iter = _alpha;

    // prior term
    fir_der = -(alpha_iter - _m_a) / _sigmasq_a * _N / N_tol;
    sec_der = -1.0 / _sigmasq_a * _N / N_tol;

    //if (_update_stage == 1 && _g == 0) {
    //   cout << "fir_der = " << fir_der << endl;
    //   cout << "sec_der = " << sec_der << endl;
    //}


    // log-likelihood term
    for (i = 0; i < _N; i++) {
        b = _b_infor[i];

        for (r = 0; r < _R; r++) {
            k = _w_mc[i][r];
            fir_der += 1.0 * _x_mc[i][r] / _R;
            p_temp = mu[i][k] / (mu[i][k] + _phi[b]);
            fir_der += -(_x_mc[i][r] + _phi[b]) * p_temp / _R;
            sec_der += -_phi[b] * (_x_mc[i][r] + _phi[b]) * pow(p_temp, 2.0) / mu[i][k] / _R;
        }

        /*
        if (std::isnan(fir_der) || std::isnan(sec_der)) {
            cout << "Updating the alpha of " << i << "-th cell" << endl;
            t = _t_infor[i];
            cout << "g = " << _g << endl;
            cout << "t = " << t << endl;
            cout << "alpha = " << _alpha << endl;
            cout << "beta = " << _beta[k] << endl;
            cout << "eta = " << _eta[t * _K + k] << endl;
            cout << "nu = " << _nu[b] << endl;
            cout << "delta = " << _delta[i] << endl;
            cout << "mu[" << i << "][" << k << "] = " << mu[i][k] << endl;
            for (r = 0; r < _R; r++) {
                cout << "w_ir = " << _w_mc[i][r] << endl;
                cout << "X_ir = " << _x_mc[i][r] << endl;
            }
            cout << "p_temp = " << p_temp << endl;
            cout << "fir_der = " << fir_der << endl;
            cout << "sec_der = " << sec_der << endl;
            exit(1);
        }
        */
    }


    double shift_alpha = -fir_der / sec_der;
    if (shift_alpha > _max_step) {
        alpha_iter = _alpha + _lr * _max_step;
    }
    else if (shift_alpha < -_max_step) {
        alpha_iter = _alpha - _lr * _max_step;
    }
    else {
        alpha_iter = _alpha + _lr * shift_alpha;
    }

    /*
    if (_update_stage == 1 && _g == 0) {
        cout << "Updating alpha..." << endl;
        cout << "fir_der = " << fir_der << endl;
        cout << "sec_der = " << sec_der << endl;
        cout << "alpha_iter = " << alpha_iter << endl;
    }
    */

    // update beta
    // calculate the auxiliary variable sum_{i} u_btik * y_btig
    for (k = 1; k < _K; k++) {
        //if (_g == 0) {
        //    cout << "Updating beta[" << k << "]..." << endl;
        //}
        // prior term
        fir_der = (-_beta[k] / _tau1 * _PrL[k] - _beta[k] / _tau0 * (1.0 - _PrL[k])) * _N / N_tol;
        sec_der = (-1.0 / _tau1 * _PrL[k] - 1.0 / _tau0 * (1.0 - _PrL[k])) * _N / N_tol;

        //if (_g == 0) {
        //    cout << "_PrL[" << k << "] = " << _PrL[k]  << endl;
        //    cout << "fir_der = " << fir_der << endl;
        //    cout << "sec_der = " << sec_der << endl;
        //}

        // log-likelihood term
        for (i = 0; i < _N; i++) {
            b = _b_infor[i];

            for (r = 0; r < _R; r++) {
                if (_w_mc[i][r] == k) {
                    fir_der += 1.0 * _x_mc[i][r] / _R;
                    p_temp = mu[i][k] / (mu[i][k] + _phi[b]);
                    fir_der += -(_x_mc[i][r] + _phi[b]) * p_temp / _R;
                    sec_der += -_phi[b] * (_x_mc[i][r] + _phi[b]) * pow(p_temp, 2.0) / mu[i][k] / _R;
                }
            }

            /*
            if (std::isnan(fir_der) || std::isnan(sec_der)) {
                cout << "Updating the beta by " << i << "-th cell" << endl;
                t = _t_infor[i];
                cout << "g = " << _g << endl;
                cout << "t = " << t << endl;
                cout << "alpha = " << _alpha << endl;
                cout << "beta = " << _beta[k] << endl;
                cout << "eta = " << _eta[t * _K + k] << endl;
                cout << "nu = " << _nu[b] << endl;
                cout << "delta = " << _delta[i] << endl;
                cout << "mu[" << i << "][" << k << "] = " << mu[i][k] << endl;
                for (r = 0; r < _R; r++) {
                    cout << "w_ir = " << _w_mc[i][r] << endl;
                    cout << "X_ir = " << _x_mc[i][r] << endl;
                }
                cout << "p_temp = " << p_temp << endl;
                cout << "fir_der = " << fir_der << endl;
                cout << "sec_der = " << sec_der << endl;
                exit(1);
            }
            */
        }

        //if (_g == 0) {
        //    cout << "fir_der = " << fir_der << endl;
        //    cout << "sec_der = " << sec_der << endl;
        //}

        double shift_beta = -fir_der / sec_der;
        if (shift_beta > _max_step) {
            _beta[k] = _beta[k] + _lr * _max_step;
        }
        else if (shift_beta < -_max_step) {
            _beta[k] = _beta[k] - _lr * _max_step;
        }
        else {
            _beta[k] = _beta[k] + _lr * shift_beta;
        }

    }

    // update eta
    if (_update_stage > 0) {

        if (_control_ind != 1) {

            for (t = 1; t < _T; t++) {
                for (k = 0; k < _K; k++) {

                    //if (_g == 0) {
                    //    cout << "Before updating eta" << endl;
                    //    cout << "eta = " << _eta[t * _K + k] << endl;
                    //}

                    // prior term
                    fir_der = (-_eta[t * _K + k] / _tau1 * _PrJ[t * _K + k] - _eta[t * _K + k] / _tau0 * (1.0 - _PrJ[t * _K + k])) * _nt_part[t] / _nt[t];
                    sec_der = (-1.0 / _tau1 * _PrJ[t * _K + k] - 1.0 / _tau0 * (1.0 - _PrJ[t * _K + k])) * _nt_part[t] / _nt[t];


                    //if (_g == 0) {
                    //    cout << "nt_part = " << _nt_part[t] << endl;
                    //    cout << "nt = " << _nt[t] << endl;
                    //    cout << "fir_der = " << fir_der << endl;
                    //    cout << "sec_der = " << sec_der << endl;
                    //}

                    // log-likelihood term
                    for (i = 0; i < _N; i++) {
                        if (_t_infor[i] == t) {
                            b = _b_infor[i];
                            for (r = 0; r < _R; r++) {
                                if (_w_mc[i][r] == k) {
                                    fir_der += 1.0 * _x_mc[i][r] / _R;
                                    p_temp = mu[i][k] / (mu[i][k] + _phi[b]);
                                    fir_der += -(_x_mc[i][r] + _phi[b]) * p_temp / _R;
                                    sec_der += -_phi[b] * (_x_mc[i][r] + _phi[b]) * pow(p_temp, 2.0) / mu[i][k] / _R;
                                }
                            }
                        }

                    }

                    double shift_eta = -fir_der / sec_der;
                    if (shift_eta > _max_step) {
                        _eta[t * _K + k] = _eta[t * _K + k] + _lr * _max_step;
                    }
                    else if (shift_eta < -_max_step) {
                        _eta[t * _K + k] = _eta[t * _K + k] - _lr * _max_step;
                    }
                    else {
                        _eta[t * _K + k] = _eta[t * _K + k] + _lr * shift_eta;
                    }
                    /*
                    if (_g == 0) {
                        cout << "After updating eta" << endl;
                        cout << "eta = " << _eta[t * _K + k] << endl;
                        cout << "fir_der = " << fir_der << endl;
                        cout << "sec_der = " << sec_der << endl;
                    }
                    */
                } // end of k
            } // end of t
        } // end of if
    }

    // update nu
    // calculate the auxiliary variable sum_{i}  y_btig

    for (b = 1; b < _B; b++) {

        // prior term
        fir_der = -(_nu[b] - _m_c[b]) / _sigmasq_c * _nb_part[b] / _nb[b];
        sec_der = -1.0 / _sigmasq_c * _nb_part[b] / _nb[b];

        // log-likelihood term
        for (i = 0; i < _N; i++) {
            if (_b_infor[i] == b) {
                for (r = 0; r < _R; r++) {
                    fir_der += 1.0 * _x_mc[i][r] / _R;
                    k = _w_mc[i][r];
                    p_temp = mu[i][k] / (mu[i][k] + _phi[b]);
                    fir_der += -(_x_mc[i][r] + _phi[b]) * p_temp / _R;
                    sec_der += -_phi[b] * (_x_mc[i][r] + _phi[b]) * pow(p_temp, 2.0) / mu[i][k] / _R;
                }
            }  // end of if

        } // end of i

        double shift_nu = -fir_der / sec_der;
        if (shift_nu > _max_step) {
            _nu[b] = _nu[b] + _lr * _max_step;
        }
        else if (shift_nu < -_max_step) {
            _nu[b] = _nu[b] - _lr * _max_step;
        }
        else {
            _nu[b] = _nu[b] + _lr * shift_nu;
        }
    }

    // update phi
    for (b = 0; b < _B; b++) {
        // prior term
        fir_der = ((_phi_prior[0] - 1) / _phi[b] - _phi_prior[1]) * _nb_part[b] / _nb[b];
        sec_der = -(_phi_prior[0] - 1) / pow(_phi[b], 2.0) * _nb_part[b] / _nb[b];

        // log-likelihood term
        for (i = 0; i < _N; i++) {

            if (_b_infor[i] == b) {

                t = _t_infor[i];

                fir_der += -gsl_sf_psi(_phi[b]) + log(_phi[b]) + 1.0;
                sec_der += -gsl_sf_psi_1(_phi[b]) + 1.0 / _phi[b];

                for (r = 0; r < _R; r++) {
                    fir_der += gsl_sf_psi(_phi[b] + _x_mc[i][r]) / _R;
                    sec_der += gsl_sf_psi_1(_phi[b] + _x_mc[i][r]) / _R;

                    k = _w_mc[i][r];
                    fir_der += -(log(mu[i][k] + _phi[b]) + (_x_mc[i][r] + _phi[b]) / (mu[i][k] + _phi[b])) / _R;
                    sec_der += -(1.0 / (mu[i][k] + _phi[b]) + (mu[i][k] - _x_mc[i][r]) / pow(mu[i][k] + _phi[b], 2.0)) / _R;

                }
            }  // end of if
        } // end of i

        double shift_phi = fir_der / sec_der;
        if (sec_der > 0) { // if the second derivative is greater than 0, then update phi for gradient ascent
            if (fir_der > 0) {
                _phi[b] = _phi[b] + _lr * _phi[b] * (exp(1.0) - 1.0);
            }
            else {
                _phi[b] = _phi[b] - _lr * _phi[b] * (1.0 - exp(-1.0));
            }
        }
        else {
            if ((_phi[b] - shift_phi) / _phi[b] < exp(-1.0)) {
                _phi[b] = _phi[b] - _lr * _phi[b] * (1.0 - exp(-1.0));
            }
            else {
                _phi[b] = _phi[b] - _lr * shift_phi;
            }
        }

        if (_phi[b] > 100.0) {
            _phi[b] = 100.0;
        }
    }
    

    for (i = 0; i < _N; i++) {
        delete[] mu[i];
    }
    delete[] mu;
    return alpha_iter;

}

double _Stoc_Newton_abenp_CGD(int _N, int _B, int _T, int _K, int _R, int* _b_infor, int* _t_infor,// int _g,
    double _alpha, double* _beta, double* _eta, double* _nu, double* _delta, double* _phi,
    int** _w_mc, int** _x_mc, int _control_ind,
    double* _nb, double* _nt, double* _nb_part, double* _nt_part,
    double _m_a, double _sigmasq_a,  // prior for alpha
    double _tau0, double _tau1, // prior for beta
    double* _PrL, double* _PrJ, // prior for eta
    double* _m_c, double _sigmasq_c, // prior for nu
    double* _phi_prior, // prior for phi
    int _update_stage, double _lr, double _max_step) {

    double p_temp, fir_der, sec_der;
    int i, b, t, k, r;
    // auxiliary variables
    double** mu = new double* [_N];
    for (i = 0; i < _N; i++) {
        b = _b_infor[i];
        t = _t_infor[i];
        mu[i] = new double[_K];
        for (k = 0; k < _K; k++) {
            mu[i][k] = exp(_alpha + _beta[k] + _eta[t * _K + k] + _nu[b] + _delta[i]);
        }
    }

    double N_tol = 0;
    for (int b = 0; b < _B; b++) {
        N_tol += _nb[b];
    }

    // update alpha
    // cout << "Update alpha" << endl;
    double alpha_iter;

    // calculate the auxiliary variable sum_{i} y_btig
    alpha_iter = _alpha;

    // prior term
    fir_der = -(alpha_iter - _m_a) / _sigmasq_a * _N / N_tol;
    sec_der = -1.0 / _sigmasq_a * _N / N_tol;

    //if (_update_stage == 1 && _g == 0) {
    //   cout << "fir_der = " << fir_der << endl;
    //   cout << "sec_der = " << sec_der << endl;
    //}


    // log-likelihood term
    for (i = 0; i < _N; i++) {
        b = _b_infor[i];

        for (r = 0; r < _R; r++) {
            k = _w_mc[i][r];
            fir_der += 1.0 * _x_mc[i][r] / _R;
            p_temp = mu[i][k] / (mu[i][k] + _phi[b]);
            fir_der += -(_x_mc[i][r] + _phi[b]) * p_temp / _R;
            sec_der += -_phi[b] * (_x_mc[i][r] + _phi[b]) * pow(p_temp, 2.0) / mu[i][k] / _R;
        }

    }


    double shift_alpha = -fir_der / sec_der;
    if (shift_alpha > _max_step) {
        alpha_iter = _alpha + _lr * _max_step;
    }
    else if (shift_alpha < -_max_step) {
        alpha_iter = _alpha - _lr * _max_step;
    }
    else {
        alpha_iter = _alpha + _lr * shift_alpha;
    }

    for (i = 0; i < _N; i++) {
        b = _b_infor[i];
        t = _t_infor[i];
        for (k = 0; k < _K; k++) {
            mu[i][k] = exp(alpha_iter + _beta[k] + _eta[t * _K + k] + _nu[b] + _delta[i]);
        }
    }

    // update beta
    // calculate the auxiliary variable sum_{i} u_btik * y_btig
    // cout << "Update beta" << endl;
    double* fir_der_beta = new double[_K];
    double* sec_der_beta = new double[_K];

    // prior term
    for (k = 0; k < _K; k++) {
        // prior term
        fir_der_beta[k] = (-_beta[k] / _tau1 * _PrL[k] - _beta[k] / _tau0 * (1.0 - _PrL[k])) * _N / N_tol;
        sec_der_beta[k] = (-1.0 / _tau1 * _PrL[k] - 1.0 / _tau0 * (1.0 - _PrL[k])) * _N / N_tol;
    }

    // log-likelihood term
    for (i = 0; i < _N; i++) {

        b = _b_infor[i];
        for (r = 0; r < _R; r++) {

            k = _w_mc[i][r];

            fir_der_beta[k] += 1.0 * _x_mc[i][r] / _R;
            p_temp = mu[i][k] / (mu[i][k] + _phi[b]);
            fir_der_beta[k] += -(_x_mc[i][r] + _phi[b]) * p_temp / _R;
            sec_der_beta[k] += -_phi[b] * (_x_mc[i][r] + _phi[b]) * pow(p_temp, 2.0) / mu[i][k] / _R;

        }
    }

    //if (_g == 0) {
    //    cout << "fir_der = " << fir_der << endl;
    //    cout << "sec_der = " << sec_der << endl;
    //}

    for (k = 0; k < _K; k++) {
        double shift_beta = -fir_der_beta[k] / sec_der_beta[k];
        if (shift_beta > _max_step) {
            _beta[k] = _beta[k] + _lr * _max_step;
        }
        else if (shift_beta < -_max_step) {
            _beta[k] = _beta[k] - _lr * _max_step;
        }
        else {
            _beta[k] = _beta[k] + _lr * shift_beta;
        }
    }

    // parameter adjustment
    alpha_iter = alpha_iter + _beta[0];
    for (k = 1; k < _K; k++) {
        _beta[k] = _beta[k] - _beta[0];
    }
    _beta[0] = 0.0;

    // update eta
    // cout << "Update eta" << endl;
    double** fir_der_eta = new double* [_T];
    double** sec_der_eta = new double* [_T];

    for (t = 0; t < _T;t++) {
        fir_der_eta[t] = new double[_K];
        sec_der_eta[t] = new double[_K];
    }

    if (_update_stage > 0) {

        if (_control_ind != 1) {

            for (i = 0; i < _N; i++) {
                b = _b_infor[i];
                t = _t_infor[i];
                for (k = 0; k < _K; k++) {
                    mu[i][k] = exp(alpha_iter + _beta[k] + _eta[t * _K + k] + _nu[b] + _delta[i]);
                }
            }

            // prior term
            for (t = 0; t < _T; t++) {
                for (k = 0; k < _K; k++) {

                    fir_der_eta[t][k] = (-_eta[t * _K + k] / _tau1 * _PrJ[t * _K + k] - _eta[t * _K + k] / _tau0 * (1.0 - _PrJ[t * _K + k])) * _nt_part[t] / _nt[t];
                    sec_der_eta[t][k] = (-1.0 / _tau1 * _PrJ[t * _K + k] - 1.0 / _tau0 * (1.0 - _PrJ[t * _K + k])) * _nt_part[t] / _nt[t];

                } // end of k
            } // end of t


            // log-likelihood term
            for (i = 0; i < _N; i++) {

                t = _t_infor[i];
                b = _b_infor[i];

                for (r = 0; r < _R; r++) {

                    k = _w_mc[i][r];
                    fir_der_eta[t][k] += 1.0 * _x_mc[i][r] / _R;
                    p_temp = mu[i][k] / (mu[i][k] + _phi[b]);
                    fir_der_eta[t][k] += -(_x_mc[i][r] + _phi[b]) * p_temp / _R;
                    sec_der_eta[t][k] += -_phi[b] * (_x_mc[i][r] + _phi[b]) * pow(p_temp, 2.0) / mu[i][k] / _R;
                }

            }


            for (t = 0; t < _T; t++) {

                for (k = 0; k < _K; k++) {

                    double shift_eta = -fir_der_eta[t][k] / sec_der_eta[t][k];
                    if (shift_eta > _max_step) {
                        _eta[t * _K + k] = _eta[t * _K + k] + _lr * _max_step;
                    }
                    else if (shift_eta < -_max_step) {
                        _eta[t * _K + k] = _eta[t * _K + k] - _lr * _max_step;
                    }
                    else {
                        _eta[t * _K + k] = _eta[t * _K + k] + _lr * shift_eta;
                    }
                }
            }

            // parameter adjustment
            alpha_iter = alpha_iter + _eta[0];
            for (k = 1; k < _K; k++) {
                _beta[k] = _beta[k] + _eta[k] - _eta[0];
            }

            for (t = 1; t < _T; t++) {
                for (k = 0; k < _K; k++) {
                    _eta[t * _K + k] = _eta[t * _K + k] - _eta[k];
                }
            }

            for (k = 0; k < _K; k++) {
                _eta[k] = 0.0;
            }

        } // end of if
    }

    for (i = 0; i < _N; i++) {
        b = _b_infor[i];
        t = _t_infor[i];
        for (k = 0; k < _K; k++) {
            mu[i][k] = exp(alpha_iter + _beta[k] + _eta[t * _K + k] + _nu[b] + _delta[i]);
        }
    }


    // update nu
    // calculate the auxiliary variable sum_{i}  y_btig
    // cout << "Update nu" << endl;
    double* fir_der_nu = new double[_B];
    double* sec_der_nu = new double[_B];

    // prior term
    for (b = 0; b < _B; b++) {
        fir_der_nu[b] = -(_nu[b] - _m_c[b]) / _sigmasq_c * _nb_part[b] / _nb[b];
        sec_der_nu[b] = -1.0 / _sigmasq_c * _nb_part[b] / _nb[b];
    }

    // log-likelihood term
    for (i = 0; i < _N; i++) {
        b = _b_infor[i];
        for (r = 0; r < _R; r++) {
            fir_der_nu[b] += 1.0 * _x_mc[i][r] / _R;
            k = _w_mc[i][r];
            p_temp = mu[i][k] / (mu[i][k] + _phi[b]);
            fir_der_nu[b] += -(_x_mc[i][r] + _phi[b]) * p_temp / _R;
            sec_der_nu[b] += -_phi[b] * (_x_mc[i][r] + _phi[b]) * pow(p_temp, 2.0) / mu[i][k] / _R;
        }  // end of if
    } // end of i

    for (b = 0; b < _B; b++) {
        double shift_nu = -fir_der_nu[b] / sec_der_nu[b];
        if (shift_nu > _max_step) {
            _nu[b] = _nu[b] + _lr * _max_step;
        }
        else if (shift_nu < -_max_step) {
            _nu[b] = _nu[b] - _lr * _max_step;
        }
        else {
            _nu[b] = _nu[b] + _lr * shift_nu;
        }
    }

    // parameter adjustment
    alpha_iter = alpha_iter + _nu[0];
    for (b = 1; b < _B; b++) {
        _nu[b] = _nu[b] - _nu[0];
    }
    _nu[0] = 0.0;

    // update phi
    // cout << "Update phi" << endl;
    double* fir_der_phi = new double[_B];
    double* sec_der_phi = new double[_B];


    for (i = 0; i < _N; i++) {
        b = _b_infor[i];
        t = _t_infor[i];
        for (k = 0; k < _K; k++) {
            mu[i][k] = exp(alpha_iter + _beta[k] + _eta[t * _K + k] + _nu[b] + _delta[i]);
        }
    }

    // prior term
    for (b = 0; b < _B; b++) {
        fir_der_phi[b] = ((_phi_prior[0] - 1) / _phi[b] - _phi_prior[1]) * _nb_part[b] / _nb[b];
        sec_der_phi[b] = -(_phi_prior[0] - 1) / pow(_phi[b], 2.0) * _nb_part[b] / _nb[b];
    }

    // log-likelihood term
    for (i = 0; i < _N; i++) {

        b = _b_infor[i];
        t = _t_infor[i];

        fir_der_phi[b] += -gsl_sf_psi(_phi[b]) + log(_phi[b]) + 1.0;
        sec_der_phi[b] += -gsl_sf_psi_1(_phi[b]) + 1.0 / _phi[b];

        for (r = 0; r < _R; r++) {

            fir_der_phi[b] += gsl_sf_psi(_phi[b] + _x_mc[i][r]) / _R;
            sec_der_phi[b] += gsl_sf_psi_1(_phi[b] + _x_mc[i][r]) / _R;

            k = _w_mc[i][r];
            fir_der_phi[b] += -(log(mu[i][k] + _phi[b]) + (_x_mc[i][r] + _phi[b]) / (mu[i][k] + _phi[b])) / _R;
            sec_der_phi[b] += -(1.0 / (mu[i][k] + _phi[b]) + (mu[i][k] - _x_mc[i][r]) / pow(mu[i][k] + _phi[b], 2.0)) / _R;

        }

    } // end of i


    for (b = 0; b < _B; b++) {
        double shift_phi = fir_der_phi[b] / sec_der_phi[b];
        if (sec_der_phi[b] > 0) { // if the second derivative is greater than 0, then update phi for gradient ascent
            if (fir_der_phi[b] > 0) {
                _phi[b] = _phi[b] + _lr * _phi[b] * (exp(1.0) - 1.0);
            }
            else {
                _phi[b] = _phi[b] - _lr * _phi[b] * (1.0 - exp(-1.0));
            }
        }
        else {
            if ((_phi[b] - shift_phi) / _phi[b] < exp(-1.0)) {
                _phi[b] = _phi[b] - _lr * _phi[b] * (1.0 - exp(-1.0));
            }
            else {
                _phi[b] = _phi[b] - _lr * shift_phi;
            }
        }

        if (_phi[b] > 100.0) {
            _phi[b] = 100.0;
        }
    }



    // cout << "Free memory for mu, beta and nu" << endl;
    for (i = 0; i < _N; i++) {
        delete[] mu[i];
    }
    delete[] mu;

    delete[] fir_der_phi;
    delete[] sec_der_phi;

    for (t = 0; t < _T; t++) {
        delete[] fir_der_eta[t];
        delete[] sec_der_eta[t];
    }

    delete[] fir_der_eta;
    delete[] sec_der_eta;

    delete[] fir_der_nu;
    delete[] sec_der_nu;

    delete[] fir_der_beta;
    delete[] sec_der_beta;

    return alpha_iter;

}

double _Stoc_Newton_d(int _b, int _t, int _G, int _K, int _R, // int _i, int _p,
    double* _alpha, double** _beta, double** _eta, double** _nu, double** _phi,
    int* _w_mc, int** _x_mc, double _m_d, double _sigmasq_d,
    double _lr, int _max_step,
    double _delta) {


    double fir_der, sec_der;
    double delta_iter;
    int g, k, r;


    double mu_temp, p_temp;

    // prior term
    fir_der = -(_delta - _m_d) / _sigmasq_d;
    sec_der = -1.0 / _sigmasq_d;

    // log-likelihood term


    for (g = 0; g < _G; g++) {

        for (r = 0; r < _R; r++) {

            fir_der += 1.0 * _x_mc[g][r] / _R;
            k = _w_mc[r];
            mu_temp = exp(_alpha[g] + _beta[g][k] + _eta[g][_t * _K + k] + _nu[g][_b] + _delta);
            p_temp = mu_temp / (mu_temp + _phi[g][_b]);
            fir_der += -(_phi[g][_b] + _x_mc[g][r]) * p_temp / _R;
            sec_der += -_phi[g][_b] * (_phi[g][_b] + _x_mc[g][r]) * pow(p_temp, 2.0) / mu_temp / _R;

        }
        /*
        if (std::isnan(fir_der) || std::isnan(sec_der)) {
            cout << "Updating the eta" << endl;
            cout << "g = " << g << endl;
            cout << "t = " << _t << endl;
            cout << "alpha = " << _alpha << endl;
            cout << "beta = " << _beta[k] << endl;
            cout << "eta = " << _eta[_t * _K + k] << endl;
            cout << "nu = " << _nu[_b] << endl;
            cout << "delta = " << _delta << endl;
            cout << "mu[" << k << "] = " << mu_temp << endl;
            cout << "p_temp = " << p_temp << endl;
            cout << "fir_der = " << fir_der << endl;
            cout << "sec_der = " << sec_der << endl;
            exit(1);
        }
        */
    } // end of g

    double shift_delta = -fir_der / sec_der;

    if (shift_delta > _max_step) {
        delta_iter = _delta + _lr * _max_step;
    }
    else if (shift_delta < -_max_step) {
        delta_iter = _delta - _lr * _max_step;
    }
    else {
        delta_iter = _delta + _lr * shift_delta;
    }

    return delta_iter;

}

double _Cal_log_complete_like(int _b, int _t, int _G, int _K,
    double* _pi, double* _alpha, double** _beta, double** _eta, double** _nu, double _delta, double** _phi, double* _u,
    int* _Y) {

    double log_comlike, log_mu;
    // double* log_ratio = new double[_K];

    log_comlike = 0.0;
    for (int k = 0; k < _K; k++) {

        // u_btik * log(pi_btk)
        log_comlike += _u[k] * log(_pi[k]);

        
        for (int g = 0; g < _G; g++) {

            // u_btik * (log(mu_btig) * x_btig - log(mu_btig + phi_bg) * (x_btig + phi_bg))
            log_mu = _alpha[g] + _beta[g][k] + _eta[g][_t * _K + k] + _nu[g][_b] + _delta;
            log_comlike += _u[k] * ( log_mu * _Y[g] - log(exp(log_mu) + _phi[g][_b]) * (_phi[g][_b] + _Y[g]));
        }
    }

    // log(Gamma(phi_bg + x_btig)) - log(Gamma(phi_bg)) + phi_bg * log(phi_bg)
    for (int g = 0; g < _G; g++) {
        log_comlike += lgamma(_Y[g] + _phi[g][_b]) - lgamma(_phi[g][_b]) + _phi[g][_b] * log(_phi[g][_b]);
    }

    return log_comlike;

}

double _Cal_loglike(int _b, int _t, int _G, int _K, int _R, // int _i,
    double* _pi, double** _gamma, double* _alpha, double** _beta, double** _eta, double** _nu, double _delta, double** _phi,
    int* _w_MC, int** _x_MC, int* _y) {

    double loglike, log_mu, gamma_temp;

    loglike = 0;

    for (int r = 0; r < _R; r++) {
        int k = _w_MC[r];
        loglike += log(_pi[k]);
        for (int g = 0; g < _G; g++) {

            // log(Sum_{z=0,1}Pr(z|x_{sig},y_{sig},gamma_b))
            gamma_temp = _gamma[_b][0] + _gamma[_b][1] * _x_MC[g][r];
            if (_y[g] > 0) {
                if (gamma_temp > 0) {
                    loglike += -gamma_temp - log(1.0 + exp(-gamma_temp));
                }
                else {
                    loglike += -log(1.0 + exp(gamma_temp));
                }
            }
            else if (_x_MC[g][r] > 0) {
                if (gamma_temp > 0) {
                    loglike += - log(1.0 + exp(-gamma_temp));
                }
                else {
                    loglike += gamma_temp  -log(1.0 + exp(gamma_temp));
                }
            }
            
            // log(x_{btig}|alpha, beta, eta, nu, delta, phi)
            if (_x_MC[g][r] > 0) {
                loglike += lgamma(_phi[g][_b] + _x_MC[g][r]) - lgamma(_phi[g][_b]);
            }
            loglike += -lgamma(_x_MC[g][r] + 1.0);

            //if (_i == 0 && g < 10) {
            //    cout << "After calculating phi + x + 1 choose phi, loglike = " << loglike << endl;
            //}

            log_mu = _alpha[g] + _beta[g][k] + _eta[g][_t * _K + k] + _nu[g][_b] + _delta;
            loglike += log_mu * _x_MC[g][r] + _phi[g][_b] * log(_phi[g][_b]);
            if (exp(log_mu) > _phi[g][_b]) {
                loglike += -(log_mu + log(1.0 + _phi[g][_b] * exp(-log_mu))) * (_x_MC[g][r] + _phi[g][_b]);
            }
            else {
                loglike += -log(exp(log_mu) + _phi[g][_b]) * (_x_MC[g][r] + _phi[g][_b]);
            }

        }  
    }

    loglike = loglike / _R;

    return loglike;

}


double _Cal_obs_loglike(int _b, int _t, int _G, int _K,
    double* _pi, double** _gamma, double* _alpha, double** _beta, double** _eta, double** _nu, double _delta, double** _phi,
    int* _Y, int _times) {

    double loglike, log_mu, p_temp, logp, log1mp, max_logr, sum_explogr, lr0_temp, sum_lr0, linear_temp, log_temp;
    double* log_ratio = new double[_K];
    
    //int x_max;

    for (int k = 0; k < _K; k++) {
        log_ratio[k] = log(_pi[k]);

        for (int g = 0; g < _G; g++) {
            log_mu = _alpha[g] + _beta[g][k] + _eta[g][_t * _K + k] + _nu[g][_b] + _delta;
            p_temp = exp(log_mu) / (exp(log_mu) + _phi[g][_b]);
            if (p_temp < pow(exp(-1.0), 100)) {
                logp = -100;
                log1mp = log(1 - p_temp);
            }
            else if (1 - p_temp < pow(exp(-1.0), 100)) {
                logp = log(p_temp);
                log1mp = -100;
            }
            else {
                logp = log(p_temp);
                log1mp = log(1 - p_temp);
            }

            if (_Y[g] > 0) {

                linear_temp = _gamma[_b][0] + _gamma[_b][1] * _Y[g];
                if (linear_temp > -100) {
                    log_ratio[k] += -log(1 + exp(linear_temp));
                }
                log_ratio[k] += lgamma(_Y[g] + _phi[g][_b]) - lgamma(_Y[g] + 1) - lgamma(_phi[g][_b]);
                log_ratio[k] += logp * _Y[g] + _phi[g][_b] * log1mp;
            }
            else {

                int t = 0;
                // x_max = (int)3 * exp(log_mu + _gamma[_b][1]);
                lr0_temp = _phi[g][_b] * log1mp; //x=0
                sum_lr0 = lr0_temp;

                double logq, logHt, log0t;

                while (t < _times) {
                    logq = log(_phi[g][_b]) - log(exp(log_mu) * (1 - exp((t + 1) * _gamma[_b][1])) + _phi[g][_b]);

                    logHt = (t + 1) * _gamma[_b][0] + _phi[g][_b] * logq;
                    log0t = lr0_temp + (t + 1) * _gamma[_b][0];

                    log_temp = log(1.0 + pow(-1.0, t) * exp(logHt - sum_lr0) + pow(-1.0, t + 1) * exp(log0t - sum_lr0));
                    sum_lr0 = sum_lr0 + log_temp;
                    t++;
                }

                sum_lr0 = sum_lr0 - log_temp / 2.0; // take the average between the lower bound and upper bound
                log_ratio[k] += sum_lr0;
                
            }
        }

        if (k == 0) {
            max_logr = log_ratio[k];
        }
        else if (max_logr < log_ratio[k]) {
            max_logr = log_ratio[k];
        }
    }

    sum_explogr = 0.0;
    for (int k = 0; k < _K; k++) {
        sum_explogr += exp(log_ratio[k] - max_logr);
    }

    loglike = max_logr + log(sum_explogr);

    delete[] log_ratio;

    return loglike;

}

double fdrDEindicator(double** _PPI, double _kappa, int _G, int _K) {

    double xi, fdr;
    double sum_xi = 0.0;
    int count_intri = 0;
    if (_kappa > 0) {
        for (int g = 0; g < _G; g++) {
            for (int k = 0; k < _K; k++) {
                xi = 1 - _PPI[g][k];
                if (xi <= _kappa) {
                    sum_xi += xi;
                    count_intri++;
                }
            }
        }
        fdr = sum_xi / count_intri;
    }
    else {
        fdr = 0.0;
    }

    return(fdr);
}

bool decreasing(double i, double j) { return (i > j); }

// calculate the DE posterior probability threshold
double postprob_DE_thr_fun(double** _L_est, double** _J_est, double _fdr_threshold, int _G, int _T, int _K) {

    double kappa;
    double postprob_DE_thr = 0.5;
    double fdr = 0.0;
    vector<double> vec_PPI;
    for (int g = 0; g < _G; g++) {
        copy(_L_est[g] + 1, _L_est[g] + _K, back_inserter(vec_PPI));
        copy(_J_est[g] + _K, _J_est[g] + _T * _K, back_inserter(vec_PPI));
    }
    // sort all PPIs decreasingly
    // cout << "The length of vec_PPI is " << vec_PPI.size() << endl;
    sort(vec_PPI.begin(), vec_PPI.end(), decreasing);

    // unique values in the PPI vector
    // cout << "The length of vec_PPI is " << vec_PPI.size() << endl;
    vector<double>::iterator it;
    it = unique(vec_PPI.begin(), vec_PPI.end());                                                           //                ^
    vec_PPI.resize(distance(vec_PPI.begin(), it));
    //cout << "The length of vec_PPI is " << vec_PPI.size() << endl;

    double** PPI_est = new double*[_G];
    for (int g = 0; g < _G; g++) {
        PPI_est[g] = new double[_K * _T - 1];
        for (int k = 0; k < _K - 1;k++) {
            PPI_est[g][k] = _L_est[g][k + 1];
        }
        for (int k = _K; k < _K * _T; k++) {
            PPI_est[g][k - 1] = _J_est[g][k];
        }
    }


    int i = 0;
    kappa = 1 - vec_PPI[0];
    fdr = fdrDEindicator(PPI_est, kappa, _G, _K * _T - 1);

    while (fdr <= _fdr_threshold && kappa <= postprob_DE_thr && i < vec_PPI.size()) {
        i++;
        kappa = 1 - vec_PPI[i];
        fdr = fdrDEindicator(PPI_est, kappa, _G, _K * _T - 1);
        // cout << "kappa = " << kappa << ", fdr = " << fdr << "." << endl;
    }

    if (i == 0) {
        kappa = 0.0;
    }
    else {
        kappa = 1 - vec_PPI[i - 1];
    }

    for (int g = 0; g < _G;g++) {
        delete[] PPI_est[g];
    }
    delete[] PPI_est;

    return kappa;
}

int IG_index(double* _PPI, double _kappa, int _K) {

    int IG_indicator = 0;
    for (int k = 1; k < _K; k++) {
        if (1 - _PPI[k] <= _kappa) {
            IG_indicator = 1;
        }
        // cout << "PPI_" << k << " = " << _PPI[k] << ", ";
    }
    // cout << "IG_indicator = " << IG_indicator << "." << endl;

    return IG_indicator;
}

int main(int argc, char** argv) {

    //////////////////////////////////////////////////////////
    //  0. Input the project and version from command line  //
    //////////////////////////////////////////////////////////
    string proj, wk_dir, readcount_dir;
    int ver, K, rep, R;
    // int iter_max; // the maximum iteration number and the tolerant
    int seed; // seed for random number generator 
    int nc; // number of codes for parallel
    int Stoc_N; // the number of cells in each partition of the stochastic M step
    int max_iter_steps, min_iter_steps, Iter_early_stop, check_per_iter, load_W, output_X;
    size_t Num_control_genes;
    int opt;

    // for beta and eta spike and slab prior
    double tau0, max_step_Newton, b_tau1, xi, hyper_sig;

    // learning rate
    double rate_coef = 1;
    int decay_coef = 1;
    
    // Size of the validation set
    double valid_size = 2.0;

    // default setting
    wk_dir.assign("./");
    readcount_dir.assign("./");
    // tol_para_order = 3;
    // tol_hyper_order = 4;
    // tol_Newton_order = 6;
    max_iter_steps = 1000;
    min_iter_steps = 100;
    check_per_iter = 5;
    Stoc_N = 1000;
    Iter_early_stop = 10;
    max_step_Newton = 1.0;
    Num_control_genes = 0;
    R = 10;
    load_W = 1;
    output_X = 1;

    // hyperparameters
    xi = 2.0;
    hyper_sig = 5.0;
    tau0 = 0.1;
    b_tau1 = 500;

    //while ((opt = getopt(argc, argv, "p:v:N:G:B:K:r:i:s:")) != -1) {
    while ((opt = getopt(argc, argv, ":d::r::p:v:K:L::D::M::N::O::P::Q::S::R::E::F::C::V::W::X::x::h::s:c:")) != -1) {
        switch (opt) {
        case 'd':
            if (optarg) {
                wk_dir.assign(optarg);
            }
            break;
        case 'r':
            if (optarg) {
                readcount_dir.assign(optarg);
            }
            break;
        case 'p':
            proj.assign(optarg);
            break;
        case 'v':
            ver = atoi(optarg);
            break;
        case 'K':
            K = atoi(optarg);
            break;
        case 'L':
            if (optarg) {
                rate_coef = atof(optarg);
            }
            break;
        case 'D':
            if (optarg) {
                decay_coef = atoi(optarg);
            }
            break;
        case 'M':
            if (optarg) {
                Stoc_N = atoi(optarg);
            }
            break;
        case 'N':
            if (optarg) {
                Iter_early_stop = atoi(optarg);
            }
            break;
        case 'O':
            if (optarg) {
                max_iter_steps = atoi(optarg);
            }
            break;
        case 'P':
            if (optarg) {
                min_iter_steps = atoi(optarg);
            }
            break;
        case 'Q':
            if (optarg) {
                check_per_iter = atoi(optarg);
            }
            break;
        case 'S':
            if (optarg) {
                max_step_Newton = atof(optarg);
            }
            break;
        case 'R':
            if (optarg) {
                R = atoi(optarg);
            }
            break;
        case 'E':
            if (optarg) {
                tau0 = atof(optarg);
            }
            break;
        case 'F':
            if (optarg) {
                b_tau1 = atof(optarg);
            }
            break;
        case 'C':
            if (optarg) {
                Num_control_genes = atoi(optarg);
            }
            break;
	    case 'V':
            if (optarg) {
                valid_size = atof(optarg);
            }
            break;
        case 'W':
            if (optarg) {
                load_W = atoi(optarg);
            }
            break;
        case 'X':
            if (optarg) {
                output_X = atoi(optarg);
            }
            break;
        case 'x':
            if (optarg) {
                xi = atof(optarg);
            }
            break;
        case 'h':
            if (optarg) {
                hyper_sig = atof(optarg);
            }
            break;
        case 's':
            seed = atoi(optarg);
            break;
        case 'c':
            nc = atoi(optarg);
            break;
        case '?':
            cout << "Unknown option " << optarg << endl;
            break;
        case ':':
            cout << "Missing option for " << optarg << endl;
            break;
        case 1:
            cout << "Non-option arg: " << optarg << endl;
            break;
        default: /* '?' */
            cout << "Error Usage: " << argv[0] << " [-d WorkingDir] [-p ProjectName] [-v Version] [-K CellTypeNumber] [-r ReplicationIndex] [-i MaxIterNum] [-P TolFirstStep] [-Q TolSecondStep] [-N TolNewton] [-M MaxIterNoNewton] [-S MaxStepLength] [-B SlabVarBeta] [-E SlabVarEta] [-s Seed] [-c Cores] name" << endl;
            exit(EXIT_FAILURE);
        }
    }


    // Check the input argument
    cout << "The working directory is " << wk_dir << endl;
    cout << "The directory to load read count data is " << readcount_dir << endl;
    cout << "The project name is " << proj << endl;
    cout << "The version is " << ver << endl;
    cout << "The total number of cell types is " << K << endl;
    cout << "Each partition contains around " << Stoc_N << " cells for stochastic updating;" << endl;
    cout << "The maximum step length of the Newton's method is " << max_step_Newton << endl;
    cout << "In the MC-E steps, " << R << " samples are generated." << endl;
    cout << "The specified spike variance of the prior is " << tau0 << endl;
    cout << "The seed is " << seed << endl;
    cout << "The number of cores for parallel computing is " << nc << endl;

    if (min_iter_steps * 2 > max_iter_steps) {
        min_iter_steps = max_iter_steps / 2;
        cout << "When the minimum iteration number is larger than half of the maximum iteration number, ";
        cout << "then the minimum iteration number is changed as half of the maximum iteration number." << endl;
    }

    // set count data matrix and the number of cells in each batch
    string fcount, fdim, fcontrol, finitial, ftbinfor;
    fcount = readcount_dir + "count_data_";
    fcount = fcount + proj; 
    fcount = fcount + "_v";
    fcount = fcount + to_string(ver);
    fcount = fcount + ".txt";

    fdim = readcount_dir + "dim_";
    fdim = fdim + proj;
    fdim = fdim + "_v";
    fdim = fdim + to_string(ver);
    fdim = fdim + ".txt";

	fcontrol = readcount_dir + "control_genes_";
	fcontrol = fcontrol + proj;
    fcontrol = fcontrol + "_v";
    fcontrol = fcontrol + to_string(ver);
    fcontrol = fcontrol + ".txt";

    finitial = readcount_dir + "w_initial_K";
    finitial = finitial + to_string(K);
    finitial = finitial + "_";
    finitial = finitial + proj;
    finitial = finitial + "_v";
    finitial = finitial + to_string(ver);
    finitial = finitial + ".txt";

	ftbinfor = readcount_dir + "tbinfor_";
	ftbinfor = ftbinfor + proj;
    ftbinfor = ftbinfor + "_v";
    ftbinfor = ftbinfor + to_string(ver);
    ftbinfor = ftbinfor + ".txt";

	auto start_overall = chrono::system_clock::now();

	// Set the number of cores for parallel
    omp_set_num_threads(nc);

    ////////////////////////
    // 1. Load count data //
    ////////////////////////
    // allocate memory for the observation
    int N, S, G, B, T, P;
	

    ifstream countFile, dimFile, controlFile, labelFile, tbinforFile;
    dimFile.open(fdim.c_str());
    if(!dimFile){
        cout << "Unable to open the file of dimensions from " << fdim << endl;
        exit(1); // terminate with error
    }
    dimFile >> N;
    cout << "There are totally " << N << " cells of ";
    dimFile >> S;
    cout << S << " samples and ";
    dimFile >> G;
    cout << G << " genes from ";
    dimFile >> B;
    cout << B << " batches under ";
    dimFile >> T;
    cout << T << " treatments." << endl;
    
    dimFile >> P;
    cout << "At the same time, there are totally " << P << " batch-treatment pairs." << endl;

    int **BT_pair = new int* [P];

    for (int p = 0; p < P; p++) {
        BT_pair[p] = new int[3];
        dimFile >> BT_pair[p][0];
        dimFile >> BT_pair[p][1];
        dimFile >> BT_pair[p][2];
        cout << BT_pair[p][2] << " cells are measured in batch " << BT_pair[p][0] << " under treatment " << BT_pair[p][1] << "." << endl;

        BT_pair[p][0] = BT_pair[p][0] - 1;
        BT_pair[p][1] = BT_pair[p][1] - 1;
    }

    dimFile.close();

    
    int **Y = new int*[G];
    for(int g = 0; g < G; g++){
        Y[g] = new int[N];// for parallel G
    }

    cout << "Loading the count data matrix from " << fcount << endl;
    
    countFile.open(fcount.c_str());
    if (!countFile) {
        cout << "Unable to open count data file" << endl;
        exit(1); // terminate with error
    }
    for(int g = 0; g < G; g++){
        for(int i = 0; i < N; i++){
            countFile >> Y[g][i];
        }
    }
    countFile.close();

	
	int *control_genes = new int[G];
    int num_control = 0;
    int exist_control = 1;

	cout << "Loading control genes from " << fcontrol << endl;
	controlFile.open(fcontrol.c_str());
    if (!controlFile) {
        cout << "Unable to open control gene file" << endl;
        exit(1); // terminate with error
    }
    for(int g = 0; g < G; g++){
        controlFile >> control_genes[g];
        num_control += control_genes[g];
    }
    controlFile.close();

    if (num_control == 0) {
        exist_control = 0;
    }

    int* b_infor = new int[N];
    int* t_infor = new int[N];
    int* p_infor = new int[N];
    int* s_infor = new int[N];
	
	cout << "Loading the information of batches and treatments for all cells from " << fcontrol << endl;
	tbinforFile.open(ftbinfor.c_str());
    if (!tbinforFile) {
        cout << "Unable to open tbinfor file" << endl;
        exit(1); // terminate with error
    }
    for(int i = 0; i < N; i++){

        tbinforFile >> b_infor[i]; // load batch index of i-th cell
        tbinforFile >> t_infor[i]; // load treatment index of i-th cell
        tbinforFile >> p_infor[i]; // load pair index of i-th cell
        tbinforFile >> s_infor[i]; // load sample index of i-th cell

		// adjust index for C++
		t_infor[i] = t_infor[i] - 1;
		b_infor[i] = b_infor[i] - 1;
        p_infor[i] = p_infor[i] - 1;
        s_infor[i] = s_infor[i] - 1;
    }
    tbinforFile.close();

    // count the total number of cells in each batch or under each treatment
    double* nb = new double[B];
    double* nt = new double[T];
    for (int b = 0; b < B; b++) {
        nb[b] = 0.0;
    }
    for (int t = 0; t < T; t++) {
        nt[t] = 0.0;
    }
    for (int p = 0; p < P; p++) {
        nb[BT_pair[p][0]] += BT_pair[p][2];
        nt[BT_pair[p][1]] += BT_pair[p][2];
    }

    int* ref_cell = new int[N];
    int* ind_ref_cell = new int[B];
    int cell_index = 0;
    for (int b = 0; b < B; b++) {
        ref_cell[cell_index] = 1;

        for (int i = 1; i < nb[b]; i++) {
            ref_cell[cell_index + i] = 0;
        }
        ind_ref_cell[b] = cell_index;
        cell_index += nb[b];
        // cout << "cell_index = " << cell_index << endl;
    }


    ///////////////////////////////////////////////////////
    // 2. Set hyper-parameters and initialize parameters //
    ///////////////////////////////////////////////////////
    cout << "Set initial values." << endl;
    omprng EM_Rng;
	// cout << "Initializing rng." << endl;
    EM_Rng.fixedSeed(seed);

    // intrinsic gene proportion and variance of non-intrinsic genes
    double p_beta, p_eta, tau1; // p and tau1 as well as tau1

    p_beta = EM_Rng.runif(0, 0.5);
    p_eta = EM_Rng.runif(0, 0.5);
    tau1 = b_tau1;
    
    /////////////////////////////////////
    // Initialize cell type proportion //
    double temp_sum = (K + 1) * K / 2;// sum up 1 to K
    double **prop = new double*[S];
    for(int s = 0; s < S; s++){
        prop[s] = new double[K];
        //if (b == 0) cout << "Take a look at pi[0]: ";
        for(int k = 0; k < K; k ++){
            prop[s][k] = (K - k)/temp_sum;
        }
    }

    /////////////////////////////////////
    // Initialize cell-specific effect //
    double* delta = new double[N];
    double* mu_d = new double[N];
    double* pointer_del = delta;

    double sum_y0;
    cell_index = 0;
    for (int p = 0; p < P; p++) {

        for (int i = 0; i < BT_pair[p][2]; i++) {

            if (i == 0) {
                *pointer_del = 0.0;// the first cell in each batch as reference
                sum_y0 = 0.0;// sum_g log(1+Y_b0g)
                for (int g = 0; g < G; g++) {
                    sum_y0 += Y[g][cell_index];
                }
                pointer_del++;

                delta[cell_index] = 0.0;
                mu_d[cell_index] = 0.0;
            }
            else {
                double sum_y;
                sum_y = 0.0;// sum_g log(1+Y_big)
                for (int g = 0; g < G; g++) {
                    sum_y += Y[g][cell_index];
                }
                delta[cell_index] = log(sum_y + 1) - log(sum_y0 + 1);
                mu_d[cell_index] = log(sum_y +1) - log(sum_y0 + 1) ;
                // cout << "mu_d[" << cell_index + 1 << "] = " << mu_d[cell_index] << endl;

            }
            cell_index++;
        }
    }

    ///////////////////////////////////////////////////////
    // Initialize baseline effects and cell type effects //
    int* W = new int[N];
    int* pointer_w = W;

    if (load_W != 0) {

        cout << "Loading initial cell type labels from " << finitial << endl;
        labelFile.open(finitial.c_str());
        if (!labelFile) {
            cout << "Unable to open the initial label file." << endl;
            exit(1); // terminate with error
        }

        for (int i = 0; i < N; i++) {
            labelFile >> W[i];
            W[i] = W[i] - 1;
        }
        labelFile.close();

    }
    else {

        // Randomly guess cell type indicator    
        for (int i = 0; i < N; i++) {
            int s = s_infor[i];
            *pointer_w = rand_cate(prop[s], EM_Rng);
            pointer_w++;
        }
    }
    
    cout << "Finish the initialization of W." << endl;
    // baseline effects
    double *alpha = new double[G];
	double *mu_a = new double[G];

    // cell-type specific effects
    double **beta = new double*[G];
    for(int g = 0; g < G; g++){
        beta[g] = new double[K];
    }
    double sum_betasq = 0.0;

    for (int g = 0; g < G; g++) {

        double* sum_logy = new double[K];
        int* count_type = new int[K];
        for (int k = 0; k < K; k++) {
            sum_logy[k] = 0.0;
            count_type[k] = 0;
        }

        cell_index = 0;
        // only in the first batch
        for (int i = 0; i < BT_pair[0][2]; i++) {
            int k = W[cell_index];
            sum_logy[k] += log(1 + Y[g][cell_index] / exp(delta[cell_index]));
            count_type[k] ++;
            cell_index++;
        }

        if (count_type[0] == 0) {
            alpha[g] = 0.0;
        }
        else {
            alpha[g] = sum_logy[0] / count_type[0];
        }

        mu_a[g] = alpha[g];
        // cout << "mu_a[" << g + 1 << "] = " << mu_a[g] << endl;
        beta[g][0] = 0.0;
        for (int k = 1; k < K; k++) {
            if (count_type[k] == 0) {
                beta[g][k] = 0.0;
            }
            else {
                beta[g][k] = sum_logy[k] / count_type[k] - alpha[g];
            }

            sum_betasq += pow(beta[g][k], 2.0);
        }

        

        delete[] sum_logy;
        delete[] count_type;
    }

    // roughly estimate the scale of spike variance 
    /*
    tau0_beta = sum_betasq / G / (K - 1);

    cout << "The sum of beta_gk is " << sum_betasq;
    cout << ", and the initial value of tau0sq_eta is " << tau0_beta << "." << endl;

    if (tau0_beta < 0.05) {
        tau0_beta = 0.05;
    }
    */

    //////////////////////////////
    // Initialize batch effects //
    double **nu = new double*[G];
	double **mu_c = new double*[G];
    for(int g = 0; g < G; g++){
        nu[g] = new double[B];
        mu_c[g] = new double[B];
    }
    // cout << "The empirical mean of nu: " << endl;
	int ind_nu = 0;

	// First case: regard the first batch under each treatment as the reference batch with no batch effect
	double *sum_logy = new double[B];
	int *count_batch = new int[B];

	// count the number of cells in each batch
	for(int b = 0; b < B; b++){
		count_batch[b] = 0;
	}

	for(int i = 0; i < N; i ++){
		int b = b_infor[i];
		count_batch[b]++;
	}

	// empirically estimated nu as prior
    for(int g = 0; g < G; g++){
        
		nu[g][0] = 0.0; 
        mu_c[g][0] = 0.0;

		for(int b = 0; b < B; b++){
			sum_logy[b] = 0.0;
		}

		for(int i = 0; i < N; i ++){
			int b = b_infor[i];
			sum_logy[b] += log(1+Y[g][i]/exp(delta[i]));
		}

		for(int b = 1; b < B; b++){
			nu[g][b] = sum_logy[b] / count_batch[b] - sum_logy[0] / count_batch[0];
			mu_c[g][b] = nu[g][b];

            // cout << "mu_c[" << g + 1 << "]["<< b + 1 <<"] = " << mu_c[g][b] << endl;
		}
    }

    ////////////////////////////////
    // Initialize over-dispersion //
    double** phi = new double* [G];
    for (int g = 0; g < G; g++) {
        phi[g] = new double[B];
        for (int b = 0; b < B; b++) {
            phi[g][b] = 5.0;
        }
    }

    //////////////////////////////////
    // Initialize treatment effects //
    double** eta = new double* [G];
    int ind_eta = 0;
    for (int g = 0; g < G; g++) {
        eta[g] = new double[T * K];
        for (int t = 0; t < T; t++) {
            for (int k = 0; k < K; k++) {
                ind_eta = t * K + k;
                eta[g][ind_eta] = 0.0;
            }
        }
    }
    
    ///////////////////////////////////////////////////////////
    // Initialize intercept and odds ratio of dropout events //
	cout << "Initializing gamma." << endl;
    double **gamma = new double*[B];
    for(int b = 0; b < B; b++){
        gamma[b] = new double[2];
        gamma[b][0] = 0.0;
        gamma[b][1] = -0.1;
    } 

    // underlying true read count
    int** X = new int* [G];
    for (int g = 0; g < G; g++) {
        X[g] = new int[N];
        for (int i = 0; i < N; i++) {
            X[g][i] = Y[g][i];
        }
    }

    // dropout indicator
    int** Z = new int* [G];
    for (int g = 0; g < G; g++) {
        Z[g] = new int[N];
        for (int i = 0; i < N; i++) {
            Z[g][i] = 0;
        }
    }

    // allocate memory to record the Monte Carlo samples
    int*** X_MC = new int**[G];
    int*** Z_MC = new int**[G];
    int** w_MC = new int* [N];
    for(int g = 0; g < G; g++){
        X_MC[g] = new int*[N];
        Z_MC[g] = new int*[N];

        for (int i = 0; i < N; i++) {
            X_MC[g][i] = new int[R];
            Z_MC[g][i] = new int[R];
            for (int r = 0; r < R; r++) {
                X_MC[g][i][r] = X[g][i];
                Z_MC[g][i][r] = Z[g][i];
            }
        }
    }

    for (int i = 0; i < N; i++) {
        w_MC[i] = new int[R];
        for (int r = 0; r < R;r++) {
            w_MC[i][r] = W[i];
        }
    }

    // Pr(L_{gk} = 1)
    double **PrL = new double*[G];
    for(int g = 0; g < G; g++){
        PrL[g] = new double[K];
    }

	// Pr(J_{tgk} = 1)
    double **PrJ = new double*[G];
	
    for(int g = 0; g < G; g++){
        PrJ[g] = new double[K * T];
        for (int t = 0; t < T; t++) {
            for(int k = 0; k < K; k++) {
                PrJ[g][t * K + k] = 0.0;
            }
        }
    }
    int num_noncontrol = 0;
    for (int g = 0; g < G; g++) {
        if (control_genes[g] != 1) { // non-control genes    
            num_noncontrol++;
        }
    }

    // hyper-parameters
    // for prop Dir prior
    // double xi = 2.0;

    // for alpha Normal prior
    double sigmasq_a = hyper_sig;

    // for nu Normal prior
    double sigmasq_c = hyper_sig;

    // for delta Normal prior
    double sigmasq_d = hyper_sig;

    // for gamma_b0 Normal prior
    double sigmasq_z = hyper_sig;
    
    // for gamma_b1 Gamma prior
    double gamma_prior[2];
    gamma_prior[0] = 0.001;
    gamma_prior[1] = 0.01;

    // for phi Gamma prior
    double phi_prior[2];
    phi_prior[0] = 1.0;
    phi_prior[1] = 0.1;
    

    // for pbeta  Beta prior
    double p_beta_prior[2];
    p_beta_prior[0] = 1.0;
    p_beta_prior[1] = 4.0;

    double p_eta_prior[2];
    p_eta_prior[0] = 1.0;
    p_eta_prior[1] = 9.0;

    // for tau0beta Gamma prior
    double tau1_prior[2];
    tau1_prior[0] = 2;
    tau1_prior[1] = b_tau1;

    //////////////////////
    // 3. MCESM algorithm //
    //////////////////////
    // auxiliary variable
    cout << "Running MCESM algorithm..." << endl;   
    int iter_MCESM = 0;
    int Iter_not_increase = 0;
    double max_loglike;

    // file variable to output the parameter values at each iteration
    ofstream post_File;
    

    // MCESM_iter_K5
    string output_dir = wk_dir + "MCESM_iter";
    output_dir = output_dir + "_K";
    output_dir = output_dir + to_string(K);
    output_dir = output_dir + "/";
    int check = mkdir(output_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (!check)
        cout << "Directory " << output_dir << " is created." << endl;
    else {
        cout << "Directory " << output_dir << " already exists." << endl;
    }

    string infer_dir = wk_dir + "Inference";
    infer_dir = infer_dir + "_K";
    infer_dir = infer_dir + to_string(K);
    infer_dir = infer_dir + "/";
    check = mkdir(infer_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (!check)
        cout << "Directory " << infer_dir << " is created." << endl;
    else {
        cout << "Directory " << infer_dir << " already exists." << endl;
    }

    ///////////////////////////////////////
    // record initial values into files. //
    ///////////////////////////////////////
    string out_file;

    // auxiliary variable to update prop
    double** count_w = new double*[S];
    for (int s = 0; s < S; s++) {
        count_w[s] = new double[K];
    }
    
    double** logmu = new double*[G];

   // record initial values into files.
    out_file = output_dir + "alpha_iter.txt";
    // cout << "Writing alpha_iter into " << out_file << endl;
    post_File.open(out_file.c_str(), ios::out | ios::app);
    for (int g = 0; g < G; g++) {
        post_File << alpha[g];
        post_File << " ";
    }
    post_File << endl;
    post_File.close();

    out_file = output_dir + "beta_iter.txt";
    // cout << "Writing beta_iter into " << out_file << endl;
    post_File.open(out_file.c_str(), ios::out | ios::app);
    for (int g = 0; g < G; g++) {
        for (int k = 0; k < K; k++) {
            post_File << beta[g][k];
            post_File << " ";
        }
    }
    post_File << endl;
    post_File.close();

    out_file = output_dir + "eta_iter.txt";
    // cout << "Writing eta_iter into " << out_file << endl;
    post_File.open(out_file.c_str(), ios::out | ios::app);
    for (int g = 0; g < G; g++) {
        for (int t = 0; t < T; t++) {
            for (int k = 0; k < K; k++) {
                post_File << eta[g][t * K + k];
                post_File << " ";
            }
        }
    }
    post_File << endl;
    post_File.close();

    out_file = output_dir + "nu_iter.txt";
    // cout << "Writing nu_iter into " << out_file << endl;
    post_File.open(out_file.c_str(), ios::out | ios::app);
    for (int g = 0; g < G; g++) {
        for (int b = 0; b < B; b++) {
            post_File << nu[g][b];
            post_File << " ";
        }
    }
    post_File << endl;
    post_File.close();

    out_file = output_dir + "delta_iter.txt";
    // cout << "Writing delta_iter into " << out_file << endl;
    post_File.open(out_file.c_str(), ios::out | ios::app);

    for (int i = 0; i < N; i++) {
        post_File << delta[i];
        post_File << " ";
    }
    post_File << endl;
    post_File.close();


    out_file = output_dir + "gamma_iter.txt";
    // cout << "Writing gamma_iter into " << out_file << endl;
    post_File.open(out_file.c_str(), ios::out | ios::app);

    for (int b = 0; b < B; b++) {
        post_File << gamma[b][0];
        post_File << " ";
        post_File << gamma[b][1];
        post_File << " ";
    }
    post_File << endl;
    post_File.close();

    out_file = output_dir + "phi_iter.txt";
    // cout << "Writing phi_iter into " << out_file << endl;
    post_File.open(out_file.c_str(), ios::out | ios::app);
    for (int g = 0; g < G; g++) {
        for (int b = 0; b < B; b++) {
            post_File << phi[g][b];
            post_File << " ";
        }
    }
    post_File << endl;
    post_File.close();


    out_file = output_dir + "pi_iter.txt";
    // cout << "Writing pi_iter into " << out_file << endl;
    post_File.open(out_file.c_str(), ios::out | ios::app);
    for (int s = 0; s < S; s++) {
        for (int k = 0; k < K; k++) {
            post_File << prop[s][k];
            post_File << " ";
        }
    }
    post_File << endl;
    post_File.close();

    out_file = output_dir + "ptau1_iter.txt";
    // cout << "Writing p_iter into " << out_file << endl;
    post_File.open(out_file.c_str(), ios::out | ios::app);
    post_File << tau1;
    post_File << " ";
    post_File << p_beta;
    post_File << " ";
    post_File << p_eta;
    post_File << endl;
    post_File.close();

    /////////////////////////////////
    // 4. Stochatstic EM algorithm //
    /////////////////////////////////
    // Record the running time
    auto start_MCESM = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds_MCESM;

    // Do not update eta, pbeta, tau0beta, peta, tau0eta until other paremers almost converge
    int update_stage = 0;
    int iter_control = 0;
    int* iter_stage = new int[3];

    // double tol_beta = pow(0.1, tol_beta_order);
    // When stage = 0, do not update eta, pbeta, tau0beta, peta and tau0eta
    // When stage = 1, do not update peta and tau0eta

    // auxiliary variable to update w
    double* proposal_pi = new double[K];
    for (int k = 0; k < K; k++) {
        proposal_pi[k] = 1.0 / K;
    }
    int w_proposal, w_current;
    double log_proposal, log_current;

    // Record the running time of each steps
    double* time_consumption = new double[8];

    // Stochastic updating control
    int Stoc_part = (N - 1) / Stoc_N + 1;
    double* assign_prob = new double[Stoc_part];
    int m;
    for (m = 0; m < Stoc_part; m++) {
        assign_prob[m] = 1.0 / Stoc_part;
    }
    int* group_ind = new int[N];
    int* N_group = new int[Stoc_part];

    // Generate a validation set for stopping rule
    int* valid_index = new int[N];
    double* valid_prob = new double[2];
    valid_prob[0] = valid_size / Stoc_part;
    if (valid_prob[0] > 1.0) {
        valid_prob[0] = 1.0;
    }
    valid_prob[1] = 1.0 - valid_prob[0];
    // int N_valid = 0; 
    for (int i = 0; i < N; i++) {
        valid_index[i] = rand_cate(valid_prob, EM_Rng);
    }
    
    // calculate the observed loglikelihood of initial values
    double obs_loglike;
    double obs_loglike_prev;

    obs_loglike = 0.0;

#pragma omp parallel
    {
        // auxiliary variable
        int* Y_i = new int[G];
        double loglike_thread = 0.0;

        int b, t, s, g;
#pragma omp for
        for (int i = 0; i < N; i++) {
            if (valid_index[i] == 0) {
                b = b_infor[i];
                t = t_infor[i];
                // p = p_infor[i];
                s = s_infor[i];

                for (g = 0; g < G; g++) {
                    Y_i[g] = Y[g][i];
                }

                loglike_thread += _Cal_obs_loglike(b, t, G, K,
                    prop[s], gamma, alpha, beta, eta, nu, delta[i], phi,
                    Y_i, 5);
            }
        }
#pragma omp critical
        {
            obs_loglike += loglike_thread;
            delete[] Y_i;
        }
    }// end of omp parallel


    max_loglike = obs_loglike;
    // cout << "The observed log-likelihood function of initial values is " << obs_loglike << endl;

    while ((update_stage < 2 || Iter_not_increase < Iter_early_stop) && iter_control < max_iter_steps) {

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // 4.0 Randomly partition N cells into several minibatches such that each block contains about Stoc_N cells for stochastic updating //
        if (iter_MCESM % Stoc_part == 0) {
            for (m = 0; m < Stoc_part; m++) {
                N_group[m] = 0;
            }

            for (int i = 0; i < N; i++) {
                group_ind[i] = rand_cate(assign_prob, EM_Rng);
                N_group[group_ind[i]]++;
            }

            // save out the mini-batch labels to file
            /*
            out_file = output_dir + "minibatch_label_iter.txt";
            // cout << "Writing W into " << out_file << endl;
            post_File.open(out_file.c_str(), ios::out | ios::app);
            for (int i = 0; i < N; i++) {
                post_File << group_ind[i];
                post_File << " ";
            }
            post_File << endl;
            post_File.close();
            */

            m = 0;
        }
        else {

            // use the next block for stochastic updating
            m++;

        }

        /////////////////////////
        // 4.1 Regular E steps //
        // cout << "Take the expectation of intrinsic gene indicators." << endl;

        // Pr(L_gk = 1)
        // cout << "Update Pr(L_gk = 1)." << endl;
        auto start_LJ = chrono::system_clock::now();

#pragma omp parallel for
        for (int g = 0; g < G; g++) {

            _E_L(K, // g,
                beta[g], p_beta, tau0, tau1,
                PrL[g]);
        }

        // cout << "p_beta = " << p_beta << endl;
        // cout << "tau0_beta = " << tau0_beta << endl;
        // cout << "tau1_beta = " << tau1_beta << endl;

        // Pr(J_tgk = 1)
        // cout << "Update Pr(J_tgk = 1)." << endl;
#pragma omp parallel for
        for (int g = 0; g < G; g++) {
            if (control_genes[g] != 1) {
                _E_J(K, T,
                    eta[g], p_eta, tau0, tau1,
                    PrJ[g]);
            }
        }

        //cout << "PrJ = " << endl;
        //for(int g = 0; g < G; g++){
        //    for(int q = 0; q < K * T; q ++){
        //            cout <<  PrJ[g][q] << " ";           
        //    }
        //    cout << endl;
        //}

        auto end_LJ = chrono::system_clock::now();
        chrono::duration<double> elapsed_seconds_LJ = end_LJ - start_LJ;
        time_consumption[0] = elapsed_seconds_LJ.count();

        /////////////////
        // 4.2 M steps //
        // cout << "Update pi." << endl;
        auto start_Estep = chrono::system_clock::now();
        cell_index = 0;
        double sum_w;

        for (int s = 0; s < S; s++) {
            for (int k = 0; k < K; k++) {
                count_w[s][k] = xi - 1.0;
            }
        }

        for (int i = 0; i < N; i++) {
            int s = s_infor[i];
            for (int r = 0; r < R; r++) {
                int k = w_MC[i][r];
                count_w[s][k] = count_w[s][k] + 1.0 / R;
            }
        }

        for (int s = 0; s < S; s++) {
            sum_w = 0.0;
            for (int k = 0; k < K; k++) {
                sum_w += count_w[s][k];
            }

            for (int k = 0; k < K; k++) {
                prop[s][k] = count_w[s][k] / sum_w;
            }
        }

        auto end_Estep = chrono::system_clock::now();
        chrono::duration<double> elapsed_seconds_Estep = end_Estep - start_Estep;
        time_consumption[1] = elapsed_seconds_Estep.count();

        ////////////////////////////
        // 4.2 Monte Carlo E step //
        int N_sub = N_group[m];

        // cout << "Building the variables of subset data." << endl;
        // Generate auxilliary variables, including Ind_sub, Y_sub, u_sub, bt_infor_sub
        // 
        // count the number of cells in each batch or in each treatment for the given set of cells
        double* nb_part = new double[B];
        double* nt_part = new double[T];

        for (int b = 0; b < B; b++) {
            nb_part[b] = 0.0;
        }
        for (int t = 0; t < T; t++) {
            nt_part[t] = 0.0;
        }

        // Ind_sub
        int* Ind_sub = new int[N_sub];
        int i_sub = 0;
        for (int i = 0; i < N; i++) {
            if (group_ind[i] == m) {
                Ind_sub[i_sub] = i;
                int b = b_infor[i];
                int t = t_infor[i];
                nb_part[b]++;
                nt_part[t]++;
                i_sub++;
                //if(i_sub < 5){
                //    cout << "The " << i_sub << "-th selected cells in the " << m << "-th iterations is cell " << i << endl;
                //}
            }
        }

        // logmu = alpha + beta + eta + nu + delta
        for (int g = 0; g < G; g++) {
            logmu[g] = new double[N_sub];
        }

        // X_sub
        int** X_sub = new int* [G];
        int** Z_sub = new int* [G];
        int** Y_sub = new int* [G];
        int* w_sub = new int [N_sub];
        for (int g = 0; g < G; g++) {
            X_sub[g] = new int[N_sub];
            Z_sub[g] = new int[N_sub];
            Y_sub[g] = new int[N_sub];
            for (int j = 0; j < N_sub; j++) {
                X_sub[g][j] = X[g][Ind_sub[j]];
                Z_sub[g][j] = Z[g][Ind_sub[j]];
                Y_sub[g][j] = Y[g][Ind_sub[j]];
            }
        }

        for (int j = 0; j < N_sub; j++) {
            w_sub[j] = W[Ind_sub[j]];
        }

        // X_MC_sub
        int*** X_MC_sub = new int** [G];
        int*** Z_MC_sub = new int** [G];
        // int** w_MC_sub = new int* [N_sub];
        for (int g = 0; g < G; g++) {
            X_MC_sub[g] = new int* [N_sub];
            Z_MC_sub[g] = new int* [N_sub];
            for (int j = 0; j < N_sub; j++) {
                X_MC_sub[g][j] = X_MC[g][Ind_sub[j]];
                Z_MC_sub[g][j] = Z_MC[g][Ind_sub[j]];
                //for (int r = 0; r < R; r++) {
                //    X_sub[g][j][r] = X_MC[g][Ind_sub[j]][r];
                //    Z_sub[g][j][r] = Z_MC[g][Ind_sub[j]][r];
                //}
            }
        }
        int** w_MC_sub = new int* [N_sub];
        for (int j = 0; j < N_sub; j++) {
            w_MC_sub[j] = w_MC[Ind_sub[j]];
        }

        // u_sub, bt_infor, delta_sub, ref_cell_sub
        int* b_infor_sub = new int[N_sub];
        int* t_infor_sub = new int[N_sub];
        int* p_infor_sub = new int[N_sub];
        int* s_infor_sub = new int[N_sub];
        double* delta_sub = new double[N_sub];
        int* ref_cell_sub = new int[N_sub]; // the indicators of whether a cell is a reference cell

        for (int j = 0; j < N_sub; j++) {
            b_infor_sub[j] = b_infor[Ind_sub[j]];
            t_infor_sub[j] = t_infor[Ind_sub[j]];
            p_infor_sub[j] = p_infor[Ind_sub[j]];
            s_infor_sub[j] = s_infor[Ind_sub[j]];
            delta_sub[j] = delta[Ind_sub[j]];
            ref_cell_sub[j] = ref_cell[Ind_sub[j]];
        }

        // Sample w, X, Z from the posterior distribution
        // cout << "Generate Monte Calro samples of w, X and Z." << endl;
        auto start_MCE = chrono::system_clock::now();

        for (int r = 0; r < R; r++) {

            double logr_w;
            for (int j = 0; j < N_sub; j++) {

                int b = b_infor_sub[j];
                int t = t_infor_sub[j];
                int p = p_infor_sub[j];
                int s = s_infor_sub[j];

                w_proposal = rand_cate(proposal_pi, EM_Rng);
                w_current = w_sub[j];

                if (w_proposal != w_current) {


                    log_proposal = log(prop[s][w_proposal]);
                    log_current = log(prop[s][w_current]);

                    //calculate the posterior ratio in log scale
#pragma omp parallel
                    {
                        //double logr_thread = 0.0;
                        double log_proposal_thread, log_current_thread;
                        log_proposal_thread = 0.0;
                        log_current_thread = 0.0;
                        double temp_logmu;
                        int* X_thread;
                        double* beta_thread;
                        double* eta_thread;
                        double* nu_thread;
                        double* phi_thread;
#pragma omp for
                        for (int g = 0; g < G; g++) {
                            X_thread = X_sub[g];
                            beta_thread = beta[g];
                            eta_thread = eta[g];
                            nu_thread = nu[g];
                            phi_thread = phi[g];

                            temp_logmu = alpha[g] + beta_thread[w_proposal] + eta_thread[t * K + w_proposal] + nu_thread[b] + delta_sub[j];
                            log_proposal_thread += (beta_thread[w_proposal] + eta_thread[t * K + w_proposal]) * X_thread[j];
                            log_proposal_thread += -(phi_thread[b] + X_thread[j]) * log(phi_thread[b] + exp(temp_logmu));

                            //ind_beta = j + w_current * _G;
                            temp_logmu = alpha[g] + beta_thread[w_current] + eta_thread[t * K + w_current] + nu_thread[b] + delta_sub[j];
                            log_current_thread += (beta_thread[w_current] + eta_thread[t * K + w_current]) * X_thread[j];
                            log_current_thread += -(phi_thread[b] + X_thread[j]) * log(phi_thread[b] + exp(temp_logmu));

                        }
#pragma omp critical
                        {
                            log_proposal += log_proposal_thread;
                            log_current += log_current_thread;
                        }
                    }
                    logr_w = log_proposal - log_current;
                    if (logr_w > log(EM_Rng.runif())) {  
                        w_sub[j] = w_proposal;
                    }
                }// end of if(w_proposal != w_current)
            }// end of i

#pragma omp parallel for
            for (int g = 0; g < G; g++) {
                _update_logmu(N_sub, K, t_infor_sub, b_infor_sub,
                    w_sub, alpha[g], beta[g], eta[g], nu[g], delta_sub,//parameter
                    logmu[g]);
            }

#pragma omp parallel for
            for (int g = 0; g < G; g++) {
                _update_zx(N_sub, b_infor_sub,
                    gamma, phi[g], logmu[g],
                    Y_sub[g], EM_Rng,
                    X_sub[g], Z_sub[g]);
            }

            // Record these values
            for (int g = 0; g < G; g++) {
                for (int j = 0; j < N_sub; j++) {

                    X[g][Ind_sub[j]] = X_sub[g][j];
                    Z[g][Ind_sub[j]] = Z_sub[g][j];

                    X_MC_sub[g][j][r] = X_sub[g][j]; // also update X_MC
                    Z_MC_sub[g][j][r] = Z_sub[g][j]; // also update Z_MC
                    /*
                    if (X_MC_sub[g][j][r] < 0) {
                        int b = b_infor_sub[j];
                        int t = t_infor_sub[j];

                        cout << "At the " << iter_MCESM << "-th iteration," << endl;
                        cout << "X_MC_sub[" << g << "][" << j << "][" << r << "] = " << X_MC_sub[g][j][r] << ";" << endl;
                        cout << "X_sub[" << g << "][" << j << "] = " << X_sub[g][j] << ";" << endl;
                        cout << "b = " << b_infor_sub[j] << endl;
                        cout << "gamma[0] = " << gamma[0] << ", gamma[1] = " << gamma[1] << endl;
                        cout << "alpha[" << g << "] = " << alpha[g] << endl;
                        cout << "beta[" << g << "][" << w_sub[j] << "] = " << beta[g][w_sub[j]] << endl;
                        cout << "eta[" << g << "][" << t * K + w_sub[j] << "] = " << eta[g][t * K + w_sub[j]] << endl;
                        cout << "nu[" << g << "][" << b << "] = " << nu[g][b] << endl;
                        cout << "delta_sub[" << j << "] = " << delta_sub[j] << endl;
                        cout << "phi[" << g << "][" << b << "] = " << phi[g][b] << endl;
                        cout << "logmu[" << g << "][" << j << "] = " << logmu[g][j] << endl;
                        cout << "Y_sub[" << g << "][" << j << "] = " << Y_sub[g][j] << endl;
                        cout << "Z_sub[" << g << "][" << j << "] = " << Z_sub[g][j] << endl;
                        exit(3);
                    }
                    */
                    //if (r == 9 && g < 5 && j < 5) {
                    //    cout << "X_sub[" << g << "][" << j << "] = " << X_sub[g][j] << ", and ";
                    //    cout << "Z_sub[" << g << "][" << j << "] = " << Z_sub[g][j] << endl;
                    //}
                }
            }
            for (int j = 0; j < N_sub; j++) {
                
                W[Ind_sub[j]] = w_sub[j];
                w_MC[Ind_sub[j]][r] = w_sub[j];

                //if (r == 9 && j < 5) {
                //    cout << "w_sub[" << j << "] = " << w_sub[j] << endl;
                //}
            }
        }

        auto end_MCE = chrono::system_clock::now();
        chrono::duration<double> elapsed_seconds_MCE = end_MCE - start_MCE;
        time_consumption[2] = elapsed_seconds_MCE.count();

        //////////////////////////////////////////////
        // 4.4 stochastic M steps by Adam algorithm //
        // cout << "Start stochastic updating." << endl;
        auto start_update_gamma = chrono::system_clock::now();
        double learn_rate;
        if (update_stage == 0) {
            learn_rate = rate_coef / log(iter_MCESM / Stoc_part / decay_coef + 3.0);
        }
        else if(update_stage == 1) {
            learn_rate = rate_coef / log((iter_MCESM - iter_stage[0]) / Stoc_part / decay_coef + 3.0);
        }
       // else if(update_stage == 2) {
       //     learn_rate = rate_coef / log((iter_MCESM - iter_stage[1] - iter_stage[0]) / Stoc_part / decay_coef + 3.0);
       //}
        else {
            learn_rate = rate_coef / log((iter_MCESM - iter_stage[1] - iter_stage[0]) / Stoc_part / decay_coef + 3.0);
        }
        
        // Record the parameter values of the previous iteration
        // for (int m = 0; m < Stoc_part; m++) {

        // Stochastically update gamma
        // Adam algorithm
        // cout << "Stochastically update gamma." << endl;
        double exp_temp, linear_term, first_der, second_der;
        for (int b = 0; b < B; b++) {
            // double gamma0_prev = gamma[b][0];
            // double gamma1_prev = gamma[b][1];

            // update gamma_b0
            // prior term
            first_der = -gamma[b][0] / sigmasq_z * nb_part[b] / nb[b];
            second_der = -1.0 / sigmasq_z * nb_part[b] / nb[b];

            // log-likelihood term
            for (int j = 0; j < N_sub; j++) {
                if (b_infor_sub[j] == b) {
                    for (int g = 0; g < G; g++) {
                        for (int r = 0; r < R; r++) {

                            first_der += 1.0 * Z_MC_sub[g][j][r] / R;
                            linear_term = gamma[b][0] + gamma[b][1] * X_MC_sub[g][j][r];
                            if (linear_term > -100) {
                                exp_temp = exp(linear_term);
                                first_der += -exp_temp / (1.0 + exp_temp) / R;
                                second_der += -exp_temp / pow(1.0 + exp_temp, 2.0) / R;
                            }

                        }// end of r
                    }// end of g
                }// end of if
            } // end of j

            double shift_gamma = -first_der / second_der;
            if (shift_gamma > max_step_Newton) {
                gamma[b][0] = gamma[b][0] + learn_rate * max_step_Newton;
            }
            else if (shift_gamma < -max_step_Newton) {
                gamma[b][0] = gamma[b][0] - learn_rate * max_step_Newton;
            }
            else {
                gamma[b][0] = gamma[b][0] + learn_rate * shift_gamma;
            }


            // update gamma_b1
            first_der = ((gamma_prior[0] - 1.0) / gamma[b][1] + gamma_prior[1]) * nb_part[b] / nb[b];
            second_der = -(gamma_prior[0] - 1.0) / pow(gamma[b][1], 2.0) * nb_part[b] / nb[b];

            // log-likelihood term
            for (int j = 0; j < N_sub; j++) {
                if (b_infor_sub[j] == b) {
                    for (int g = 0; g < G; g++) {
                        for (int r = 0; r < R; r++) {

                            first_der += 1.0 * X_MC_sub[g][j][r] * Z_MC_sub[g][j][r] / R;
                            linear_term = gamma[b][0] + gamma[b][1] * X_MC_sub[g][j][r];
                            if (linear_term > -100) {
                                exp_temp = exp(linear_term);
                                first_der += -exp_temp * X_MC_sub[g][j][r] / (1.0 + exp_temp) / R;
                                second_der += -exp_temp * pow(X_MC_sub[g][j][r], 2.0) / pow(1.0 + exp_temp, 2.0) / R;
                            }
                        }
                    }
                }
            }
            
            shift_gamma = -first_der / second_der;
            if ((gamma[b][1] + shift_gamma) / gamma[b][1] < exp(-1.0)) {
                gamma[b][1] = gamma[b][1] - learn_rate * gamma[b][1] * (1.0 - exp(-1.0));
            }
            else if ((gamma[b][1] + shift_gamma) / gamma[b][1] > exp(1.0)) {
                gamma[b][1] = gamma[b][1] + learn_rate * gamma[b][1] * (exp(1.0) - 1.0);
            }
            else {
                gamma[b][1] = gamma[b][1] + learn_rate * shift_gamma;
            }
        }
        auto end_update_gamma = chrono::system_clock::now();
        chrono::duration<double> elapsed_seconds_update_gamma = end_update_gamma - start_update_gamma;
        time_consumption[3] = elapsed_seconds_update_gamma.count();


        // Stochastically update alpha, beta, eta, nu, phi
        // cout << "Stochastically update alpha, beta, eta, nu and phi." << endl;
        auto start_update_abenp = chrono::system_clock::now();

        if (update_stage == 0) {
#pragma omp parallel for
            for (int g = 0; g < G; g++) {

                alpha[g] = _Stoc_Newton_abenp_GD(N_sub, B, T, K, R, b_infor_sub, t_infor_sub,
                    // parameters and latent variables
                    alpha[g], beta[g], eta[g], nu[g], delta_sub, phi[g],
                    // observed data and latent variables
                    w_MC_sub, X_MC_sub[g], control_genes[g],
                    // prior
                    nb, nt, nb_part, nt_part,
                    mu_a[g], sigmasq_a,
                    tau0, tau1,
                    PrL[g], PrJ[g],
                    mu_c[g], sigmasq_c,
                    phi_prior,
                    // the setting of Newton's method
                    update_stage, learn_rate, max_step_Newton);
            }
        }
        else {

#pragma omp parallel for
            for (int g = 0; g < G; g++) {

                alpha[g] = _Stoc_Newton_abenp_CGD(N_sub, B, T, K, R, b_infor_sub, t_infor_sub,
                    // parameters and latent variables
                    alpha[g], beta[g], eta[g], nu[g], delta_sub, phi[g],
                    // observed data and latent variables
                    w_MC_sub, X_MC_sub[g], control_genes[g],
                    // prior
                    nb, nt, nb_part, nt_part,
                    mu_a[g], sigmasq_a,
                    tau0, tau1,
                    PrL[g], PrJ[g],
                    mu_c[g], sigmasq_c,
                    phi_prior,
                    // the setting of Newton's method
                    update_stage, learn_rate, max_step_Newton);
            }

        }
 
        auto end_update_abenp = chrono::system_clock::now();
        chrono::duration<double> elapsed_seconds_update_abenp = end_update_abenp - start_update_abenp;
        time_consumption[4] = elapsed_seconds_update_abenp.count();

        // Stochastically update delta
        // cout << "Stochastically update delta." << endl;
        auto start_update_delta = chrono::system_clock::now();
#pragma omp parallel
        {
            // auxiliary variable to update gamma
            int** X_i = new int* [G];
            int b, t, g, r;
            for (g = 0; g < G; g++) {
                X_i[g] = new int[R];
            }
#pragma omp for
            for (int j = 0; j < N_sub; j++) {

                b = b_infor_sub[j];
                t = t_infor_sub[j];

                for (g = 0; g < G; g++) {
                    for (r = 0; r < R; r++) {
                        X_i[g][r] = X_MC_sub[g][j][r];
                    }
                }

                delta_sub[j] = _Stoc_Newton_d(b, t, G, K, R, // i, p,
                    alpha, beta, eta, nu, phi,
                    w_MC_sub[j], X_i, mu_d[Ind_sub[j]], sigmasq_d,
                    learn_rate, max_step_Newton,
                    delta_sub[j]);

            }
#pragma omp critical
            {
                for (g = 0; g < G; g++) {
                    delete[] X_i[g];
                }
                delete[] X_i;
            }
        }// end of omp parallel

        for (int j = 0; j < N_sub; j++) {
            delta[Ind_sub[j]] = delta_sub[j];
        }

        // Adjust for alpha and nu such that alpha_g + nu_bg + delta_sig = alpha_g' + nu_bg' + delta_sig'
        for (int g = 0; g < G; g++) {

            alpha[g] = alpha[g] + delta[ind_ref_cell[0]];
            for (int b = 1;b < B;b++) {
                nu[g][b] = nu[g][b] + delta[ind_ref_cell[b]] - delta[ind_ref_cell[0]];
            }

        }

        // Adjust for delta
        for (int i = 0; i < N; i++) {
            int b = b_infor[i];
            int ind_ref = ind_ref_cell[b];
            if (ref_cell[i] == 0) {
                delta[i] = delta[i] - delta[ind_ref];
            }
        }

        for (int b = 0; b < B; b++) {
            delta[ind_ref_cell[b]] = 0.0;
        }

        // free memory
        for (int g = 0; g < G; g++) {
            delete[] X_sub[g];
            delete[] Z_sub[g];
            delete[] Y_sub[g];
            delete[] X_MC_sub[g];
            delete[] Z_MC_sub[g];
        }
        delete[] X_sub;
        delete[] Z_sub;
        delete[] Y_sub;
        delete[] w_sub;
        delete[] X_MC_sub;
        delete[] Z_MC_sub;
        delete[] w_MC_sub;
        
        for (int g = 0; g < G; g++) {
            delete[] logmu[g];
        }
        delete[] b_infor_sub;
        delete[] t_infor_sub;
        delete[] delta_sub;
        delete[] ref_cell_sub;
        delete[] Ind_sub;
        delete[] nb_part;
        delete[] nt_part;
        auto end_update_delta = chrono::system_clock::now();
        chrono::duration<double> elapsed_seconds_update_delta = end_update_delta - start_update_delta;
        time_consumption[5] = elapsed_seconds_update_delta.count();


        //////////////////////////////////////
        // 4.2.3 update pbeta and peta //
        auto start_ptau = chrono::system_clock::now();
        // if (update_stage == 1) {
        //    
        //}
        // pbeta_iter = p_beta;
        if (update_stage > 0) {
            // cout << "Update p_beta." << endl;
            double cbeta_temp = p_beta_prior[0] - 1;
            for (int g = 0; g < G; g++) {
                for (int k = 1; k < K; k++) {
                    cbeta_temp += PrL[g][k];
                }
            }
            p_beta = cbeta_temp / (G * (K - 1) + p_beta_prior[0] + p_beta_prior[1] - 2);
        }

        // peta_iter = p_eta;
        if (update_stage > 1) {
            // cout << "Update p_eta." << endl;
            double ceta_temp = p_eta_prior[0] - 1;
            for (int g = 0; g < G; g++) {
                if (control_genes[g] != 1) {
                    for (int t = 1; t < T; t++) {
                        for (int k = 0; k < K; k++) {
                            ceta_temp += PrJ[g][t * K + k];
                        }
                    }
                }
            }
            p_eta = ceta_temp / (num_noncontrol * (T - 1) * K + p_eta_prior[0] + p_eta_prior[1] - 2);
        }

        /////////////////////////////////////////////////////
        // 4.2.4 update sigmasq_{beta,0} and sigma_{eta,0} //
        //if (update_stage == 1) {
        //    cout << "Update beta0sq." << endl;
        //}
        if (update_stage > 0) {
            double a_temp = tau1_prior[0] - 1;
            double b_temp = tau1_prior[1];

            for (int g = 0; g < G; g++) {
                for (int k = 1; k < K; k++) {
                    a_temp += PrL[g][k] / 2.0;
                    b_temp += PrL[g][k] * pow(beta[g][k], 2.0) / 2.0;
                }
            }

            for (int g = 0; g < G; g++) {
                if (control_genes[g] != 1) {
                    for (int t = 1; t < T; t++) {
                        for (int k = 0; k < K; k++) {
                            a_temp += PrJ[g][t * K + k] / 2.0;
                            b_temp += PrJ[g][t * K + k] * pow(eta[g][t * K + k], 2.0) / 2.0;
                        }
                    }
                }
            }

            tau1 = b_temp / a_temp;
        }

        auto end_ptau = chrono::system_clock::now();
        chrono::duration<double> elapsed_seconds_ptau = end_ptau - start_ptau;
        time_consumption[6] = elapsed_seconds_ptau.count();


        ////////////////////////////////
        // 4.3 Early stop criterion //
        auto start_loglike = chrono::system_clock::now();
        if ((iter_MCESM + 1) % check_per_iter == 0) {


            obs_loglike_prev = obs_loglike;
            obs_loglike = 0.0;

#pragma omp parallel
            {
                // auxiliary variable
                int* Y_i = new int[G];
                double loglike_thread = 0.0;

                int b, t, s, g;
#pragma omp for
                for (int i = 0; i < N; i++) {
                    if (valid_index[i] == 0) {
                        b = b_infor[i];
                        t = t_infor[i];
                        // p = p_infor[i];
                        s = s_infor[i];

                        for (g = 0; g < G; g++) {
                            Y_i[g] = Y[g][i];
                        }

                        loglike_thread += _Cal_obs_loglike(b, t, G, K,
                            prop[s], gamma, alpha, beta, eta, nu, delta[i], phi,
                            Y_i, 5);
                    }
                }
#pragma omp critical
                {
                    obs_loglike += loglike_thread;
                    delete[] Y_i;
                }
            }// end of omp parallel


            if (iter_MCESM == 0 || obs_loglike > max_loglike) {
                max_loglike = obs_loglike;
                Iter_not_increase = 0;
            }
            else {
                Iter_not_increase++;
            }

            cout << endl;
            cout << "At the " << iter_MCESM + 1 << "-th iteration with the updating stage " << update_stage + 1 << ", ";
            cout << "the difference of loglikelihood is " << obs_loglike - obs_loglike_prev << ", ";
            cout << "and the observed log-likelihood does not decrease in the past " << Iter_not_increase << " iterations. ";
            cout << "Currently, the maximum observed loglikelihood is " << max_loglike << "." << endl;
            cout << endl;

            // output w after calculating the approximated log-likelihood
            out_file = output_dir + "w_iter.txt";
            // cout << "Writing W into " << out_file << endl;
            post_File.open(out_file.c_str(), ios::out | ios::app);
            for (int i = 0; i < N; i++) {
                post_File << W[i];
                post_File << " ";
            }
            post_File << endl;
            post_File.close();

            // record parameter values at each iteration into files.
            // cout << "Writing the parameter values at "<< iter_MCESM + 1 << "-th iteration into the directory " << output_dir << endl;
            out_file = output_dir + "alpha_iter.txt";
            // cout << "Writing alpha_iter into " << out_file << endl;
            post_File.open(out_file.c_str(), ios::out | ios::app);
            for (int g = 0; g < G; g++) {
                post_File << alpha[g];
                post_File << " ";
            }
            post_File << endl;
            post_File.close();

            out_file = output_dir + "beta_iter.txt";
            // cout << "Writing beta_iter into " << out_file << endl;
            post_File.open(out_file.c_str(), ios::out | ios::app);
            for (int g = 0; g < G; g++) {
                for (int k = 0; k < K; k++) {
                    post_File << beta[g][k];
                    post_File << " ";
                }
            }
            post_File << endl;
            post_File.close();

            out_file = output_dir + "eta_iter.txt";
            // cout << "Writing eta_iter into " << out_file << endl;
            post_File.open(out_file.c_str(), ios::out | ios::app);
            for (int g = 0; g < G; g++) {
                for (int t = 0; t < T; t++) {
                    for (int k = 0; k < K; k++) {
                        post_File << eta[g][t * K + k];
                        post_File << " ";
                    }
                }
            }
            post_File << endl;
            post_File.close();

            out_file = output_dir + "nu_iter.txt";
            // cout << "Writing nu_iter into " << out_file << endl;
            post_File.open(out_file.c_str(), ios::out | ios::app);
            for (int g = 0; g < G; g++) {
                for (int b = 0; b < B; b++) {
                    post_File << nu[g][b];
                    post_File << " ";
                }
            }
            post_File << endl;
            post_File.close();

            out_file = output_dir + "delta_iter.txt";
            // cout << "Writing delta_iter into " << out_file << endl;
            post_File.open(out_file.c_str(), ios::out | ios::app);

            for (int i = 0; i < N; i++) {
                post_File << delta[i];
                post_File << " ";
            }
            post_File << endl;
            post_File.close();


            out_file = output_dir + "gamma_iter.txt";
            // cout << "Writing gamma into " << out_file << endl;
            post_File.open(out_file.c_str(), ios::out | ios::app);

            for (int b = 0; b < B; b++) {
                post_File << gamma[b][0];
                post_File << " ";
                post_File << gamma[b][1];
                post_File << " ";
            }
            post_File << endl;
            post_File.close();

            out_file = output_dir + "phi_iter.txt";
            // cout << "Writing phi_iter into " << out_file << endl;
            post_File.open(out_file.c_str(), ios::out | ios::app);
            for (int g = 0; g < G; g++) {
                for (int b = 0; b < B; b++) {
                    post_File << phi[g][b];
                    post_File << " ";
                }
            }
            post_File << endl;
            post_File.close();


            out_file = output_dir + "pi_iter.txt";
            // cout << "Writing pi_iter into " << out_file << endl;
            post_File.open(out_file.c_str(), ios::out | ios::app);
            for (int s = 0; s < S; s++) {
                for (int k = 0; k < K; k++) {
                    post_File << prop[s][k];
                    post_File << " ";
                }
            }
            post_File << endl;
            post_File.close();

            out_file = output_dir + "ptau1_iter.txt";
            // cout << "Writing p_iter into " << out_file << endl;
            post_File.open(out_file.c_str(), ios::out | ios::app);
            post_File << tau1;
            post_File << " ";
            post_File << p_beta;
            post_File << " ";
            post_File << p_eta;
            post_File << endl;
            post_File.close();

        }

        auto end_loglike = chrono::system_clock::now();
        chrono::duration<double> elapsed_seconds_loglike = end_loglike - start_loglike;
        time_consumption[7] = elapsed_seconds_loglike.count();

        // save the running time of the current iteration into file 
        out_file = infer_dir + "time_consumption.txt";
        post_File.open(out_file.c_str(), ios::out | ios::app);
        for (int step = 0; step < 8; step++) {
            post_File << time_consumption[step];
            post_File << " ";
        }
        post_File << endl;
        post_File.close();

        iter_control++;

        if (update_stage == 0 && (Iter_not_increase == Iter_early_stop || iter_control >= max_iter_steps) && iter_control > min_iter_steps) {

            update_stage = 1;
            iter_stage[0] = iter_MCESM + 1;

            double* qvalue = new double[G];
#pragma omp parallel for
            for (int g = 0; g < G; g++) {
                qvalue[g] = _percentile_q(N, T, K, R, b_infor, t_infor, // g,
                    // parameters and latent variables
                    alpha[g], beta[g], eta[g], nu[g], delta, phi[g],
                    // observed data
                    w_MC, X_MC[g],
                    // prior
                    PrJ[g], tau0, tau1,
                    // the setting of Newton's method
                    max_step_Newton, 0.001, 10);
            }

            double sum_etasq = 0.0;
            for (int g = 0;g < G; g++) {
                for (int t = 1; t < T; t++) {
                    for (int k = 0; k < K; k++) {
                        sum_etasq += pow(eta[g][t * K + k], 2.0);
                    }
                }
            }

            // automatically select control genes
            if (exist_control == 0) {

                // Extract the genes with the first 100 smallest QBF as control genes 
                size_t* smallest_index = new size_t[Num_control_genes];
                int stage_sort;
                stage_sort = gsl_sort_smallest_index(smallest_index, Num_control_genes, qvalue, 1, G);

                for (int i = 0; i < Num_control_genes; i++) {
                    int g = smallest_index[i];
                    control_genes[g] = 1;
                    for (int t = 0; t < T; t++) {
                        for (int k = 0; k < K; k++) {
                            eta[g][t * K + k] = 0.0;
                        }
                    }
                    // cout << "Gene " << smallest_index[i] << " is selected as a control gene!" << endl;
                }

                // record the initial value of eta
                out_file = output_dir + "eta_initial.txt";
                // cout << "Writing eta_iter into " << out_file << endl;
                post_File.open(out_file.c_str(), ios::out | ios::app);
                for (int g = 0; g < G; g++) {
                    for (int t = 0; t < T; t++) {
                        for (int k = 0; k < K; k++) {
                            post_File << eta[g][t * K + k];
                            post_File << " ";
                        }
                    }
                }
                post_File << endl;
                post_File.close();

                // record QBF
                out_file = output_dir + "qvalue.txt";
                // cout << "Writing QBF into " << out_file << endl;
                post_File.open(out_file.c_str(), ios::out | ios::app);
                for (int g = 0; g < G; g++) {
                    post_File << qvalue[g];
                    post_File << " ";
                }
                post_File << endl;
                post_File.close();


                // record control genes
                out_file = infer_dir + "control_genes.txt";
                // cout << "Writing QBF into " << out_file << endl;
                post_File.open(out_file.c_str(), ios::out | ios::app);
                for (int g = 0; g < G; g++) {
                    post_File << control_genes[g];
                    post_File << endl;
                }
                post_File.close();
            }

            // Switch the major cell type with the reference cell type
            int* cum_u = new int[K];
            int max_cum = 0.0;
            int max_celltype = 0;

            for (int k = 0; k < K;k++) {
                cum_u[k] = 0;
            }


            for (int i = 0; i < N; i++) {
                int k = W[i];
                cum_u[k] += 1;
            }

            for (int k = 0; k < K; k++) {

                if (cum_u[k] > max_cum) {
                    max_cum = cum_u[k];
                    max_celltype = k;
                }

                // cout << "cum_u[" << k << "] = " << cum_u[k] << endl;
                // cout << "max_cum = " << max_cum << endl;
                // cout << "max_celltype = " << max_celltype << endl;
            }
            cout << "Set the " << max_celltype << "-th cell type as the reference cell type." << endl;

            // modify W[i]
            for (int i = 0; i < N; i++) {
                int temp_W;
                if (W[i] == max_celltype) {
                    temp_W = 0;
                }
                else if (W[i] == 0) {
                    temp_W = max_celltype;
                }
                else {
                    temp_W = W[i];
                }
                W[i] = temp_W;
            }

            if (max_celltype != 0) {
                double alpha_temp;
                double eta_temp;
                for (int g = 0; g < G; g++) {
                    alpha_temp = alpha[g] + beta[g][max_celltype];

                    for (int k = 1; k < K; k++) {
                        if (k == max_celltype) {
                            // beta'_gk* = - beta_gk*  
                            beta[g][k] = alpha[g] - alpha_temp;
                        }
                        else {
                            // beta'_gk = alpha_g + beta_gk - alpha'_g
                            beta[g][k] = alpha[g] + beta[g][k] - alpha_temp;
                        }
                    }
                    alpha[g] = alpha_temp;
                }
            }

            Iter_not_increase = 0;
            iter_control = 0;
            delete[] cum_u;
            delete[] qvalue;

        }

        if (update_stage == 1 && (Iter_not_increase == Iter_early_stop || iter_control >= max_iter_steps) && iter_control > min_iter_steps) {

            update_stage = 2;
            iter_stage[1] = iter_MCESM + 1 - iter_stage[0];

            // Switch the major cell type with the reference cell type
            int* cum_u = new int[K];
            int max_cum = 0.0;
            int max_celltype = 0;

            for (int k = 0; k < K;k++) {
                cum_u[k] = 0;
            }


            for (int i = 0; i < N; i++) {
                int k = W[i];
                cum_u[k] += 1;
            }

            for (int k = 0; k < K; k++) {

                if (cum_u[k] > max_cum) {
                    max_cum = cum_u[k];
                    max_celltype = k;
                }

                // cout << "cum_u[" << k << "] = " << cum_u[k] << endl;
                // cout << "max_cum = " << max_cum << endl;
                // cout << "max_celltype = " << max_celltype << endl;
            }
            cout << "Set the " << max_celltype << "-th cell type as the reference cell type." << endl;
            // modify W[i]
            for (int i = 0; i < N; i++) {
                int temp_W;
                if (W[i] == max_celltype) {
                    temp_W = 0;
                }
                else if (W[i] == 0) {
                    temp_W = max_celltype;
                }
                else {
                    temp_W = W[i];
                }
                W[i] = temp_W;
            }

            if (max_celltype != 0) {
                double alpha_temp;
                double eta_temp;
                for (int g = 0; g < G; g++) {
                    alpha_temp = alpha[g] + beta[g][max_celltype];

                    // switch eta_tg1 and eta_tgk*
                    for (int t = 1; t < T; t++) {
                        for (int k = 0; k < K; k++) {
                            eta_temp = eta[g][t * K + max_celltype];
                            eta[g][t * K + max_celltype] = eta[g][t * K];
                            eta[g][t * K] = eta_temp;
                        }
                    }

                    for (int k = 1; k < K; k++) {
                        if (k == max_celltype) {
                            // beta'_gk* = - beta_gk*  
                            beta[g][k] = alpha[g] - alpha_temp;
                        }
                        else {
                            // beta'_gk = alpha_g + beta_gk - alpha'_g
                            beta[g][k] = alpha[g] + beta[g][k] - alpha_temp;
                        }
                    }
                    alpha[g] = alpha_temp;
                }
            }

            Iter_not_increase = 0;
            iter_control = 0;
            delete[] cum_u;

        }
        iter_MCESM++;

    }

    delete[] group_ind;
    delete[] assign_prob;
    delete[] N_group;
    
    iter_stage[2] = iter_MCESM - iter_stage[1] - iter_stage[0];
    auto end_MCESM = chrono::system_clock::now();
    elapsed_seconds_MCESM = end_MCESM - start_MCESM;
    cout << "In total, it has taken " << elapsed_seconds_MCESM.count() << "s to finish " << iter_MCESM << " iterations." << endl;

    auto start_inf = chrono::system_clock::now();

    ofstream est_File;
    cout << "Storing the parameter values in the last iteration..." << infer_dir << endl;
    // write the last iteration as the estimated values in "Inference_K5" folder
    // write out the iteration numbers and running times
    out_file = infer_dir + "iter_infor.txt";
    // cout << "Writing alpha_est into " << out_file << endl;
    est_File.open(out_file.c_str(), ios::out | ios::app);
    
    est_File << iter_MCESM;
    est_File << endl;
    est_File << elapsed_seconds_MCESM.count();
    est_File << endl;
    for (int st = 0; st < 3; st++) {
        est_File << iter_stage[st];
        est_File << endl;
    }
    est_File.close();
    
    // write out latent variables
    out_file = infer_dir + "alpha_est.txt";
    // cout << "Writing alpha_est into " << out_file << endl;
    est_File.open(out_file.c_str(), ios::out | ios::app);
    for (int g = 0; g < G; g++) {
        est_File << alpha[g];
        est_File << endl;
    }
    est_File.close();

    out_file = infer_dir + "beta_est.txt";
    // cout << "Writing beta_est into " << out_file << endl;
    est_File.open(out_file.c_str(), ios::out | ios::app);
    for (int g = 0; g < G; g++) {
        for (int k = 0; k < K; k++) {
            est_File << beta[g][k];
            est_File << " ";
        }
        est_File << endl;
    }
    est_File.close();

    out_file = infer_dir + "eta_est.txt";
    // cout << "Writing eta_est into " << out_file << endl;
    est_File.open(out_file.c_str(), ios::out | ios::app);
    for (int g = 0; g < G; g++) {
        for (int t = 0; t < T; t++) {
            for (int k = 0; k < K; k++) {
                est_File << eta[g][t * K + k];
                est_File << " ";
            }
        }
        est_File << endl;
    }
    est_File.close();

    out_file = infer_dir + "nu_est.txt";
    // cout << "Writing nu_est into " << out_file << endl;
    est_File.open(out_file.c_str(), ios::out | ios::app);
    for (int g = 0; g < G; g++) {
        for (int b = 0; b < B; b++) {
            est_File << nu[g][b];
            est_File << " ";
        }
        est_File << endl;
    }
    est_File.close();

    out_file = infer_dir + "delta_est.txt";
    // cout << "Writing delta_est into " << out_file << endl;
    est_File.open(out_file.c_str(), ios::out | ios::app);
    for (int i = 0; i < N; i++) {
        est_File << delta[i];
        est_File << " ";
    }
    est_File << endl;
    est_File.close();

    out_file = infer_dir + "gamma_est.txt";
    // cout << "Writing gamma_est into " << out_file << endl;
    est_File.open(out_file.c_str(), ios::out | ios::app);
    for (int b = 0; b < B; b++) {
        est_File << gamma[b][0];
        est_File << " ";
        est_File << gamma[b][1];
        est_File << endl;
    }
    est_File.close();

    out_file = infer_dir + "phi_est.txt";
    // cout << "Writing phi_est into " << out_file << endl;
    est_File.open(out_file.c_str(), ios::out | ios::app);
    for (int g = 0; g < G; g++) {
        for (int b = 0; b < B; b++) {
            est_File << phi[g][b];
            est_File << " ";
        }
        est_File << endl;
    }
    est_File.close();

    out_file = infer_dir + "pi_est.txt";
    // cout << "Writing pi_est into " << out_file << endl;
    est_File.open(out_file.c_str(), ios::out | ios::app);
    for (int s = 0; s < S; s++) {
        for (int k = 0; k < K; k++) {
            est_File << prop[s][k];
            est_File << " ";
        }
        est_File << endl;
    }
    est_File.close();

    out_file = infer_dir + "ptau1_est.txt";
    // cout << "Writing p_est into " << out_file << endl;
    est_File.open(out_file.c_str(), ios::out | ios::app);
    est_File << tau1;
    est_File << endl;
    est_File << p_beta;
    est_File << endl;
    est_File << p_eta;
    est_File << endl;
    est_File.close();

    out_file = infer_dir + "w_est.txt";
    // cout << "Writing w_post into " << out_file << endl;
    est_File.open(out_file.c_str(), ios::out | ios::app);
    for (int i = 0; i < N; i++) {
        est_File << W[i] + 1;
        est_File << endl;
    }
    est_File.close();

    out_file = infer_dir + "w_MC.txt";
    // cout << "Writing w_post into " << out_file << endl;
    est_File.open(out_file.c_str(), ios::out | ios::app);
    for (int i = 0; i < N; i++) {
        for (int r = 0; r < R; r++) {
            est_File << w_MC[i][r] + 1;  
            est_File << " ";
        }
        est_File << endl;
    }
    est_File.close();

    out_file = infer_dir + "PrL_est.txt";
    // cout << "Writing w_post into " << out_file << endl;
    est_File.open(out_file.c_str(), ios::out | ios::app);
    for (int g = 0; g < G; g++) {
        for (int k = 0; k < K; k++) {
            est_File << log(1.0 - PrL[g][k]);
            est_File << " ";
        }
        est_File << endl;
    }
    est_File.close();

    out_file = infer_dir + "PrJ_est.txt";
    // cout << "Writing w_post into " << out_file << endl;
    est_File.open(out_file.c_str(), ios::out | ios::app);
    for (int g = 0; g < G; g++) {
        for (int t = 0; t < T; t++) {
            for (int k = 0; k < K; k++) {
                est_File << log(1.0 - PrJ[g][t * K + k]);
                est_File << " ";
            }
        }
        est_File << endl;
    }
    est_File.close();

    auto end_inf = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds_inf = end_inf - start_inf;
    cout << "It takes " << elapsed_seconds_inf.count() << "s to store the estimated parameter values." << endl;

    // Identify intrinsic genes
    cout << "Identifying control genes by controlling FDR...";
    auto start_intri = chrono::system_clock::now();

    double fdr_threshold = 0.05, kappa;
    int count_gene = 0;
    int* D = new int[G];// intrinsic gene indicator


    int** J_est = new int* [G]; // cell-type-spcific treatment-specific DE gene indicators
    for (int g = 0; g < G; g++) {
        J_est[g] = new int[T * K];
        for (int t = 0;t < T; t++) {
            for (int k = 0; k < K; k++) {
                J_est[g][t * K + k] = 0;
            }
        }
    }

    kappa = postprob_DE_thr_fun(PrL, PrJ, fdr_threshold, G, T, K);
    cout << "The threshold to identify intrinsic genes is " << kappa << " to control FDR at level " << fdr_threshold << "." << endl;

    for (int g = 0; g < G; g++) {
        D[g] = IG_index(PrL[g], kappa, K);
        count_gene += D[g];
    }
    cout << "There are " << count_gene << " identified intrisic genes." << endl;
    for (int g = 0; g < G; g++) {
        if (control_genes[g] != 1) {
            for (int t = 1; t < T;t++) {
                for (int k = 0; k < K;k++) {
                    if (1 - PrJ[g][t * K + k] <= kappa) {
                        J_est[g][t * K + k] = 1;
                    }
                }
            }
        }
    }

    out_file = infer_dir + "D_est.txt";
    // cout << "Writing w_post into " << out_file << endl;
    est_File.open(out_file.c_str(), ios::out | ios::app);
    for (int g = 0; g < G; g++) {
        est_File << D[g];
        est_File << endl;
    }
    est_File.close();

    out_file = infer_dir + "E_est.txt";
    // cout << "Writing w_post into " << out_file << endl;
    est_File.open(out_file.c_str(), ios::out | ios::app);
    for (int g = 0; g < G; g++) {
        for (int t = 0; t < T;t++) {
            for (int k = 0; k < K; k++) {
                est_File << J_est[g][t * K + k];
                est_File << " ";
            }
        }

        est_File << endl;
    }
    est_File.close();

    auto end_intri = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds_intri = end_intri - start_intri;
    cout << "It takes " << elapsed_seconds_intri.count() << "s to identify intrinsic genes." << endl;


    

    ////////////////////////////////////////////////////////
    // Calculate BIC values for model comparison: Style I //
    ////////////////////////////////////////////////////////
    cout << "Calculating BIC values for all cells...";
    auto start_BIC = chrono::system_clock::now();

    obs_loglike = 0.0;

#pragma omp parallel
    {
        // auxiliary variable
        int* Y_i = new int[G];
        double loglike_thread = 0.0;

        int b, t, s, g;
#pragma omp for
        for (int i = 0; i < N; i++) {

            b = b_infor[i];
            t = t_infor[i];
            // p = p_infor[i];
            s = s_infor[i];

            for (g = 0; g < G; g++) {
                Y_i[g] = Y[g][i];
            }

            loglike_thread += _Cal_obs_loglike(b, t, G, K,
                prop[s], gamma, alpha, beta, eta, nu, delta[i], phi,
                Y_i, 5);

        }
#pragma omp critical
        {
            obs_loglike += loglike_thread;
            delete[] Y_i;
        }
    }// end of omp parallel

    out_file = infer_dir + "loglike.txt";
    // cout << "Writing p_est into " << out_file << endl;
    est_File.open(out_file.c_str(), ios::out | ios::app);
    est_File << obs_loglike;
    est_File << endl;
    est_File.close();

    double BIC;
    BIC = -2 * obs_loglike + log(G * N) * (S * (K - 1) + G * K + num_noncontrol * (T - 1) * K + G * (2 * B - 1) + N - P + B * 2);
    out_file = infer_dir + "BIC.txt";
    // cout << "Writing p_est into " << out_file << endl;
    est_File.open(out_file.c_str(), ios::out | ios::app);
    est_File << BIC;
    est_File << endl;
    est_File.close();

    auto end_BIC = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds_BIC = end_BIC - start_BIC;
    cout << "It takes " << elapsed_seconds_BIC.count() << "s to calculate BIC." << endl;

    if (output_X == 1) {
        auto start_X = chrono::system_clock::now();
        out_file = infer_dir + "imputed_count.txt";
        // cout << "Writing w_post into " << out_file << endl;
        est_File.open(out_file.c_str(), ios::out | ios::app);
        for (int g = 0;g < G; g++) {
            for (int i = 0; i < N; i++) {
                est_File << X[g][i];
                est_File << " ";
            }
            est_File << endl;
        }
        est_File.close();
        auto end_X = chrono::system_clock::now();
        elapsed_seconds_inf = end_X - start_X;
        cout << "It takes " << elapsed_seconds_inf.count() << "s to store X." << endl;
    }
    
    //free the memory
    //cout << "Free the memory." << endl;
    for (int g = 0; g < G; g++) {
        delete[] J_est[g];
    }
    delete[] J_est;
    // delete[] count_gene_eta;

    // intrinsic gene indicators
    for (int g = 0; g < G; g++) {
        delete[] PrJ[g];
    }
    delete[] PrJ;

    // for (int g = 0; g < num_noncontrol; g++) {
    //    delete[] PrJ_temp[g];
    // }
    // delete[] PrJ_temp;

    // intrinsic gene indicators
    for (int g = 0; g < G; g++) {
        delete[] PrL[g];
    }
    delete[] PrL;
    delete[] D;

    for (int s = 0; s < S;s++) {
        delete[] count_w[s];
    }
    delete[] count_w;

    //cout << "Free the memory of iteration variables." << endl;
    delete[] time_consumption;
    delete[] iter_stage;

    // cell-specific effect
    delete[] delta;
    delete[] mu_d;
    delete[] ref_cell;
    delete[] ind_ref_cell;
    
    //cout << "Free the memory of parameters." << endl;
    // treatment effects
    for (int g = 0; g < G; g++) {
        delete[] eta[g];
        // delete[] eta_prev[g];
    }
    delete[] eta;
    // delete[] eta_prev;


    // overdispersion
    for (int g = 0; g < G; g++) {
        delete[] phi[g];
        // delete[] phi_prev[g];
    }
    delete[] phi;
    // delete[] phi_prev;
    

    // batch effects
    for (int g = 0; g < G; g++) {
        delete[] nu[g];
        delete[] mu_c[g];
        // delete[] nu_prev[g];
    }
    delete[] nu;
    // delete[] nu_prev;
    delete[] mu_c;
    delete[] sum_logy;
    delete[] count_batch;

    // cell-type specific effects
    for (int g = 0; g < G; g++) {
        delete[] beta[g];
        // delete[] beta_prev[g];
    }
    delete[] beta;
    // delete[] beta_prev;

   
    //cout << "Free the memory of alpha." << endl;
    // baseline effects
    delete[] alpha;
    delete[] mu_a;
    // delete[] alpha_prev;
    
    for(int b = 0; b < B; b++){
        delete[] gamma[b];
    }
    delete[] gamma;


    //cout << "Free the memory of w." << endl;
    for(int g= 0; g < G; g++){
        delete[] X[g];
        delete[] Z[g];
        for(int i = 0 ; i < N; i++){
             delete[] X_MC[g][i]; 
             delete[] Z_MC[g][i];
        }
        delete[] X_MC[g];
        delete[] Z_MC[g];
    }
    delete[] X;
    delete[] Z;
    delete[] X_MC;
    delete[] Z_MC;


    for(int i = 0; i < N; i++){
        delete[] w_MC[i];
    }
    delete[] w_MC;
    delete[] W;
    delete[] logmu;
    


    // cell type proportion
    //cout << "Free the memory of prop." << endl;
    for (int s = 0; s < S; s++) {
        delete[] prop[s];
    }
    delete[] prop;

    delete[] nb;
    delete[] nt;

    delete[] p_infor;
    delete[] t_infor;
    delete[] b_infor;
    delete[] control_genes;

    for (int g = 0; g < G; g++) {
        delete[] Y[g];
    }
    delete[] Y;

    for (int p = 0; p < P; p++) {
        delete[] BT_pair[p];
    }
    delete[] BT_pair;


    return 0;
}