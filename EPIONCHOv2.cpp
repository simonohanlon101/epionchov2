#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <iostream>
//#include <globalvars.h>

/* *************************************************************
* EPIONCHOv2
************************************************************* */

/* *************************************************************
 * written by Martin Walker 2017
************************************************************* */

/* *************************************************************
show debugging information
0 = no information
1 = show info on which steps are reached
2 = output counter during time steps
*************************************************************** */
int DEBUG=1;

/* *************************************************************
declare population dynamics parameters
*************************************************************** */

/* the annual biting rate */
double ABR0, ABR;
double m0, m, BR;

/* within host population parameter */
double deltaH0, deltaHinfty, cH;
double alpha0, beta0, lambda0, lambda0adj;

/* distrubution parameters */
double kW0,kW1, kM0, kM1;

/* skin snipping parameters */
double wtss;
int nss;

/* within vector population paarameters */
double deltaV0, cV, aV, alphaV, aH, muL0, sigma1, sigma2, h, g, muV;

/* exposure function parameters */
double q, E0, Q, alphaF, alphaM;

/* ivermectin or moxidectin treatment parameters */
int mox = 0;
double beta1Max, gammabeta, epsilonbeta, nu, omega, mu1Max, epsilonmu, zeta;

/* default treatment program parameters */
double startTreat, startVC;
int ntrt1, ntrt2, ntrt3, ntrt4, ntrt5;
double ftrt1, ftrt2, ftrt3, ftrt4, ftrt5, cov1, cov2, cov3, cov4, cov5, noncmp;
double contra = 5;

/* vector control parameters */
double effVC, durVC;

/* globally fixed age & sex structure parameters */
double muH = 0.04;
int na = 80;
int amax = 80;
double da = 1;
double dda = 1/da;
int ns = 1;
double rhoF = 0.45;
double gammaM, gammaF, EM, EF;
int fcp = 2;
int cnst = 100;
int ncp = fcp+2;

/* total number of events */
int tnev;
 /* total number of equations */
int neq;

/* number of states modelled */

/* no L4 states */
int nL4;
/* mean prepatent period */
double prepat;
/* rate of progression */
 double gammaL;

/* no adult worm states */
int nW;
/* life expecancy */
double LE;
/* mortality rate */
 double mu0;

/* no mf states  */
int nM;
/* life expectancy mf */
double LEMf;
/* mf mortality rate */
double muM0;

int nst;
int nlst = 3;

/* indicator varibles for states*/
 std::vector<int> L4id(1), Nid(1), Fid(1), Mid(1);

/* integers for jumping vector elements */
int jumpj, jumpk, jumpl, jump, jumptoL1, jumptoL2, jumptoL3, nex, jpast;

/* dynamically allocated constants & arrays */
 std::vector<double> age(na), pa(na), pa5(na), pa20(na), ucontracum(na+1), pacumF(na+1), pacumM(na+1),psex(ns);

 /* resize at run time */
std::vector<int> ntr(1), ntr1(1), ntr2(1), ntr3(1), ntr4(1), ntr5(1), nev(1), nevcum(1);
std::vector<double> ftr1(1), ftr2(1), ftr3(1), ftr4(1), ftr5(1);


 /* compliance & coverage structure - resized at runtime  */
std::vector<double> pcp1(1), pcp2(1), pcp3(1), pcp4(1), pcp5(1), therpcp1(1), therpcp2(1), therpcp3(1), therpcp4(1), therpcp5(1);


 std::vector<double> therpcpcum1(1), therpcpcum2(1), therpcpcum3(1), therpcpcum4(1),therpcpcum5(1);
double pcmp1, pcmp2, pcmp3, pcmp4, pcmp5, psemi1, psemi2, psemi3, psemi4, psemi5, ucontra, actcov1, actcov2, actcov3, actcov4, actcov5;
int sumto;

 /* event vector  */
 std::vector<double> event(1);

/* exposure function vector */
std::vector<double> OmegaM(na), OmegaF(na), Omega(ns*na);

/* integration variables */
double t, stepbi, stepdi, stepai, switcht1, switcht2, switcht3, switcht4;
/* this is for looking one year after last treatment for pOTTIS */
double survperiod = 1;
 /* this is for checking the breakpoint */
double elimperiod = 18;

/* step functions resized depending on parameters
 and then dynamically updated through time */
std::vector<double> stepivm(1), stepmox(1), expo(1);

/* temporally dynamic parameters whihc are resized
based on parameter values and then updated through time */
std::vector<double> beta1(1), mu1(1), muM1(1), psi(1);

/* population size at each time step - resized */
std::vector<double> y, ydot;

/* storage vectors for summing over worm age groups */
std::vector<double> L4iicum(1), Niicum(1), Fiicum(1), Miicum(1);

/* storage vector for summing over net mf contribution from
adult worms in different nominal age states */
std::vector<double> lambdaiicum(1);

/* storage vectors for summming over exposure groups */
std::vector<double> L4icum(1), Nicum(1), Ni(1), Ficum(1), Fi(1), Micum(1), netlambdacum(1);

 /* storage vector for summing over net mf contribution from
 differentially exposed adult worms */
 std::vector<double> lambdaicum(1);

/* storage vectors for age dependent model output */
std::vector<double> Nj(1);
std::vector<double> Fj(1);
std::vector<double> Wj(1);
std::vector<double> WTj(1);
std::vector<double> Mj(1);
std::vector<double> Mpj(1);
std::vector<double> pFj(1);
std::vector<double> Phij(1);
std::vector<double> lambdaj(1);
std::vector<double> PiV1j(1);
std::vector<double> PiV2j(1);
std::vector<double> L1j(1);
std::vector<double> L2j(1);
std::vector<double> L3j(1);
std::vector<double> kM(1);
std::vector<double> kW(1);
std::vector<double> A(1);
std::vector<double> B(1);
std::vector<double> Mjsq(1);
std::vector<double> varMj(1);

/* storage vector for states averaged over age */
std::vector<double> Nk(1);
std::vector<double> Fk(1);
std::vector<double> Wk(1);
std::vector<double> Mk(1);
std::vector<double> M5k(1);
std::vector<double> M20k(1);
std::vector<double> Mpk(1);
std::vector<double> Mp5k(1);
std::vector<double> Mp20k(1);
std::vector<double> Phik(1);
std::vector<double> lambdak(1);
std::vector<double> L1k(1);
std::vector<double> L2k(1);
std::vector<double> L3k(1);

/* storage vector for states averaged over age & sex */
std::vector<double> Nl(1);
std::vector<double> Fl(1);
std::vector<double> Wl(1);
std::vector<double> Ml(1);
std::vector<double> M5l(1);
std::vector<double> M20l(1);
std::vector<double> Mpl(1);
std::vector<double> Mp5l(1);
std::vector<double> Mp20l(1);
std::vector<double> Phil(1);
std::vector<double> lambdal(1);
std::vector<double> L1l(1);
std::vector<double> L2l(1);
std::vector<double> L3l(1);

/* population averages*/
double N, F, W, M, M5, M20, Mp, Mp5, Mp20, Phi, L1, L2, L3, lambda, varM, Msq;

/* *************************************************************
prototype functions
*************************************************************** */

    void fillAgeStruc(double da, double muH, double amax, double contra);

    void fillTreatStruc(int ntrt1, int ntrt2, int ntrt3, int ntrt4, int ntrt5, double ftrt1, double ftrt2, double ftrt3, double ftrt4, double ftrt5, int fcp);

    void calcNormFacs();

    void fillExpoStruc();

    void calcEvents(std::vector<int> ntr);

    void fillCompStruc();

    void fillEventStruc();

    void calcCov();

    double checkSum(std::vector<double> x);

    void updateSteps();

    void updateParms();

    void updateABR(double t, double ABR0);

    void updateL3();

    void updateHostStages();

    void setStateIndicators();

    void fillInits();

/* *************************************************************
 main dirivative function
 *************************************************************** */

void Diff(std::vector<double> y)
    {
            int l, k, j, i, ii;
            int jumplcum[ncp+1];
            jumplcum[0] = 0;


            for (l=0; l<ncp; l++)
            {
              int jumpitl = ns*(nst*nev[l]*na+na*nlst);
              jumplcum[l+1] = jumplcum[l] + jumpitl;

              /* tmp vectors for summing over exposure groups to
              * calculate adult worm and mf totals by age group */
                Nicum.resize(nev[l]+1);
                Ni.resize(nev[l]+1);
                Ficum.resize(nev[l]+1);
                Fi.resize(nev[l]+1);
                Micum.resize(nev[l]+1);
                netlambdacum.resize(nev[l] + 1);

                lambdaicum.resize(nev[l]+1);
                Nicum[0] = 0;
                Ni[0] = 0;
                Fi[0] = 0;
                Ficum[0] = 0;
                Micum[0] = 0;
                lambdaicum[0] = 0;
                netlambdacum[0] = 0;

              /* number of exposure groups within each age group
              * within current compliance group*/
                nex = nst*nev[l];

                /* number of places to skip in current compliance group */
                jumpl = jumplcum[l];

                for (k=0; k<ns; k++)
                  {
                    jumpk = k*(nex*na+na*nlst);
                    for (j=0;j<na;j++)
                    {
                      jumpj = j*nex;
                      jump = jumpj+jumpk+jumpl;
                      jumptoL1 = nex*na+jumpk+jumpl;
                      jumptoL2 = nex*na+na+jumpk+jumpl;
                      jumptoL3 = nex*na+2*na+jumpk+jumpl;

                      /* by worm exposure to differnt numbers of treatments */
                        for (i=0;i<nev[l];i++)
                        {

                          for (ii = 0; ii<nL4; ii++)
                          {
                            if (j < 1) {
                              if (ii<1) {
                                ydot[i+L4id[ii]*nev[l]+jump] =  0.5*expo[i+nevcum[l]]*m*BR*Omega[j]*L3*( (deltaH0 + deltaHinfty*cH*m*BR*L3)*pow(1 + cH*m*BR*L3, -1) )
                                - (gammaL+dda)*y[i+L4id[ii]*nev[l]+jump];
                              } else {
                                ydot[i+L4id[ii]*nev[l]+jump] = gammaL*y[i+L4id[ii-1]*nev[l]+jump] - (gammaL+dda)*y[i+L4id[ii]*nev[l]+jump];
                              }
                            } else if (age[j]<contra) {
                              if (ii<1) {
                                ydot[i+L4id[ii]*nev[l]+jump] =  0.5*expo[i+nevcum[l]]*m*BR*Omega[j]*L3*( (deltaH0 + deltaHinfty*cH*m*BR*L3)*pow(1 + cH*m*BR*L3, -1) )
                                - (gammaL +  dda)*y[i+L4id[ii]*nev[l]+jump] + dda*y[i+L4id[ii]*nev[l]-nex+jump];
                              } else {
                                ydot[i+L4id[ii]*nev[l]+jump] = gammaL*y[i+L4id[ii-1]*nev[l]+jump] - (gammaL + dda)*y[i+L4id[ii]*nev[l]+jump] + dda*y[i+L4id[ii]*nev[l]-nex+jump];
                              }
                            } else {
                              if (ii<1) {
                                ydot[i+L4id[ii]*nev[l]+jump] =  0.5*expo[i+nevcum[l]]*m*BR*Omega[j]*L3*( (deltaH0 + deltaHinfty*cH*m*BR*L3)*pow(1 + cH*m*BR*L3, -1) )
                                - (gammaL +  dda)*y[i+L4id[ii]*nev[l]+jump] + dda*y[i+L4id[ii]*nev[l]-nex+jump];
                              } else {
                                ydot[i+L4id[ii]*nev[l]+jump] = gammaL*y[i+L4id[ii-1]*nev[l]+jump] - (gammaL + dda)*y[i+L4id[ii]*nev[l]+jump] + dda*y[i+L4id[ii]*nev[l]-nex+jump];
                              }
                            }
                          }

                          for (ii=0; ii<nW; ii++) {
                            /* first age group -  no incoming from aging and no treatment  */

                              if(j<1)
                              {
                                if (ii<1) {
                                  /* non-fertile worms in each exposure group  */
                                    ydot[i+Nid[ii]*nev[l]+jump] = 1*gammaL*y[i+L4id[nL4-1]*nev[l]+jump] +
                                      beta0*y[i+Fid[ii]*nev[l]+jump] - (alpha0 + mu0 + dda)*y[i+Nid[ii]*nev[l]+jump];
                                    /* fertile worms in each exposure group */
                                      ydot[i+Fid[ii]*nev[l]+jump] = 0*gammaL*y[i+L4id[nL4-1]*nev[l]+jump] +
                                        alpha0*y[i+Nid[ii]*nev[l]+jump] - (beta0 + mu0 + dda )*y[i+Fid[ii]*nev[l]+jump];
                                } else {
                                  /* non-fertile worms in each exposure group */
                                    ydot[i+Nid[ii]*nev[l]+jump] = mu0*y[i+Nid[ii-1]*nev[l]+jump] +
                                      beta0*y[i+Fid[ii]*nev[l]+jump] - (alpha0 + mu0 + dda )*y[i+Nid[ii]*nev[l]+jump];
                                    /* fertile worms in each exposure group - half of incoming larvae are female */
                                      ydot[i+Fid[ii]*nev[l]+jump] = mu0*y[i+Fid[ii-1]*nev[l]+jump] +
                                        alpha0*y[i+Nid[ii]*nev[l]+jump] - (beta0 + mu0 + dda )*y[i+Fid[ii]*nev[l]+jump];
                                }

                              }
                            else if (age[j]<contra)
                            {
                              if (ii<1) {
                                /* non-fertile worms in each exposure group  */
                                  ydot[i+Nid[ii]*nev[l]+jump] = 1*gammaL*y[i+L4id[nL4-1]*nev[l]+jump] +
                                    (beta0)*y[i+Fid[ii]*nev[l]+jump] - (alpha0 + mu0 + dda)*y[i+Nid[ii]*nev[l]+jump] + dda*y[i+Nid[ii]*nev[l]-nex+jump];
                                  /* fertile worms in each exposure group */
                                    ydot[i+Fid[ii]*nev[l]+jump] = 0*gammaL*y[i+L4id[nL4-1]*nev[l]+jump] +
                                      alpha0*y[i+Nid[ii]*nev[l]+jump] - (beta0 + mu0 + dda)*y[i+Fid[ii]*nev[l]+jump] + dda*y[i+Fid[ii]*nev[l]-nex+jump];
                              } else {
                                /* non-fertile worms in each exposure group */
                                  ydot[i+Nid[ii]*nev[l]+jump] = mu0*y[i+Nid[ii-1]*nev[l]+jump]+
                                    (beta0)*y[i+Fid[ii]*nev[l]+jump] - (alpha0 + mu0 + dda)*y[i+Nid[ii]*nev[l]+jump] + dda*y[i+Nid[ii]*nev[l]-nex+jump];
                                  /* fertile worms in each exposure group */
                                    ydot[i+Fid[ii]*nev[l]+jump] = mu0*y[i+Fid[ii-1]*nev[l]+jump] +
                                      alpha0*y[i+Nid[ii]*nev[l]+jump] - (beta0 + mu0 + dda)*y[i+Fid[ii]*nev[l]+jump] + dda*y[i+Fid[ii]*nev[l]-nex+jump];
                              }
                            }
                            else
                            {
                              if (ii<1) {
                                /* non-fertile worms in each exposure group */
                                  ydot[i+Nid[ii]*nev[l]+jump] = 1*gammaL*y[i+L4id[nL4-1]*nev[l]+jump] +
                                    (beta0 + beta1[i+nevcum[l]])*y[i+Fid[ii]*nev[l]+jump] - (alpha0 + mu0 + mu1[i+nevcum[l]]+ dda)*y[i+Nid[ii]*nev[l]+jump] + dda*y[i+Nid[ii]*nev[l]-nex+jump];
                                  /* fertile worms in each exposure group */
                                    ydot[i+Fid[ii]*nev[l]+jump] =  0*gammaL*y[i+L4id[nL4-1]*nev[l]+jump] +
                                      alpha0*y[i+Nid[ii]*nev[l]+jump] - (beta0 + beta1[i+nevcum[l]] + mu0 + mu1[i+nevcum[l]] + dda )*y[i+Fid[ii]*nev[l]+jump] + dda*y[i+Fid[ii]*nev[l]-nex+jump];
                              } else {
                                /* non-fertile worms in each exposure group */
                                  ydot[i+Nid[ii]*nev[l]+jump] = mu0*y[i+Nid[ii-1]*nev[l]+jump] +
                                    (beta0 + beta1[i+nevcum[l]])*y[i+Fid[ii]*nev[l]+jump] - (alpha0 + mu0 + mu1[i+nevcum[l]]+ dda )*y[i+Nid[ii]*nev[l]+jump] + dda*y[i+Nid[ii]*nev[l]-nex+jump];
                                  /* fertile worms in each exposure group */
                                    ydot[i+Fid[ii]*nev[l]+jump] =  mu0*y[i+Fid[ii-1]*nev[l]+jump] +
                                      alpha0*y[i+Nid[ii]*nev[l]+jump] - (beta0 + beta1[i+nevcum[l]] + mu0 + mu1[i+nevcum[l]] + dda )*y[i+Fid[ii]*nev[l]+jump] + dda*y[i+Fid[ii]*nev[l]-nex+jump];
                              }
                            }
                          } // end ii (worm age) loop
                        } //end i loop

                      /* sum worms age groups and
                      and then over exposure groups
                      within each age class */

                        for (i=0; i<nev[l]; i++)
                        {
                          Niicum[0] = 0;
                          Fiicum[0] = 0;
                          for (ii=0; ii<nW; ii++) {
                            Niicum[ii+1] = Niicum[ii] + y[i+Nid[ii]*nev[l]+jump];
                            Fiicum[ii+1] = Fiicum[ii] + y[i+Fid[ii]*nev[l]+jump];
                          }
                          Ni[i] = Niicum[nW];
                          Fi[i] = Fiicum[nW];
                          Nicum[i+1] = Nicum[i] + Niicum[nW];
                          Ficum[i+1] = Ficum[i] + Fiicum[nW];

                        }


                      Nj[j+k*na+l*ns*na] = Nicum[nev[l]];
                      Fj[j+k*na+l*ns*na] = Ficum[nev[l]];
                      Wj[j+k*na+l*ns*na] = Nj[j+k*na+l*ns*na] + Fj[j+k*na+l*ns*na];
                      WTj[j+k*na+l*ns*na] = 2*Wj[j+k*na+l*ns*na];
                      /* mating probability in each age group assuming a linear kW */
                        kW[j+k*na+l*ns*na] = kW0 + Wj[j+k*na+l*ns*na]*kW1;

                      /* proportion of total worm population that are fertile females */
                        pFj[j+k*na+l*ns*na] = Fj[j+k*na+l*ns*na]/(WTj[j+k*na+l*ns*na]);

                      /* mating probability */
                      Phij[j+k*na+l*ns*na] = 1 - pow(1 +  Wj[j+k*na+l*ns*na]*pow(kW[j+k*na+l*ns*na], -1), -(kW[j+k*na+l*ns*na]+1));

                      /* loop again over the Mf exposure groups using all-worm mating probability
                      and net incoming from all nominal adul worm age groups*/
                        for (i=0; i < nev[l]; i++)
                        {
                          /* loop over adult worm nominal age groups to calculate net mf production */
                            for (ii=0; ii<nW; ii++)
                            {
                              lambdaiicum[0] = 0;
                              if (j<contra) {
                                lambdaiicum[ii+1] = lambdaiicum[ii] +  Phij[j+k*na+l*ns*na]*1*lambda0adj*(y[i+Fid[ii]*nev[l]+jump] );
                              } else {
                                lambdaiicum[ii+1] = lambdaiicum[ii] +  Phij[j+k*na+l*ns*na]*psi[i+nevcum[l]]*lambda0adj*(y[i+Fid[ii]*nev[l]+jump]);
                              }
                            }

                          lambdaicum[i] = lambdaiicum[nW];


                          for (int ii=0; ii<nM; ii++) {
                            /* microfilarial populations from different exposure groups */
                              /* first age group -  no incoming from aging and no treatment */
                              if ( j < 1 ) {
                                if (ii<1) {
                                  ydot[i+Mid[ii]*nev[l]+jump] =  lambdaicum[i] - (muM0 + dda)*y[i+Mid[ii]*nev[l]+jump];
                                } else {
                                  ydot[i+Mid[ii]*nev[l]+jump] =  muM0*y[i+Mid[ii-1]*nev[l]+jump] - (muM0 + dda)*y[i+Mid[ii]*nev[l]+jump];
                                }
                              } else if (age[j] < contra){
                                if (ii<1) {
                                  ydot[i+Mid[ii]*nev[l]+jump] =  lambdaicum[i]  - (muM0 +  dda)*y[i+Mid[ii]*nev[l]+jump] + dda*y[i+Mid[ii]*nev[l]-nex+jump];
                                } else {
                                  ydot[i+Mid[ii]*nev[l]+jump] =  muM0*y[i+Mid[ii-1]*nev[l]+jump]  - (muM0 +  dda)*y[i+Mid[ii]*nev[l]+jump] + dda*y[i+Mid[ii]*nev[l]-nex+jump];
                                }
                              } else {
                                if (ii<1) {
                                  ydot[i+Mid[ii]*nev[l]+jump] =  lambdaicum[i] - (muM0 + muM1[i+nevcum[l]] + dda)*y[i+Mid[ii]*nev[l]+jump] + dda*y[i+Mid[ii]*nev[l]-nex+jump];
                                } else {
                                  ydot[i+Mid[ii]*nev[l]+jump] =  muM0*y[i+Mid[ii-1]*nev[l]+jump] - (muM0 + muM1[i+nevcum[l]] + dda)*y[i+Mid[ii]*nev[l]+jump] + dda*y[i+Mid[ii]*nev[l]-nex+jump];
                                }
                              }
                          } // end ii loop
                        } // end i loop



                      /* sum all the mf contributions from differntially exposed adults */
                        for (i=0; i<nev[l]; i++)
                        {
                          Miicum[ii] = 0;
                          /* first sum over nominal age categories */
                            for (ii=0; ii<nM; ii++) {
                              Miicum[ii+1] = Miicum[ii] + y[i+Mid[ii]*nev[l]+jump];
                            }
                          Micum[i+1] = Micum[i] + Miicum[nM];
                        }

                      Mj[j+k*na+l*ns*na] = Micum[nev[l]];
                      /* overdispersion mf in the skin function */
                        kM[j+k*na+l*ns*na] = kM0*exp(Mj[j+k*na+l*ns*na]*kM1);
                      /* prevalence model  */
                        Mpj[j+k*na+l*ns*na] = (1 - pow(1+wtss* Mj[j+k*na+l*ns*na]/kM[j+k*na+l*ns*na], -kM[j+k*na+l*ns*na]*nss))*(1 - pow(1+Fj[j+k*na+l*ns*na]/kW[j+k*na+l*ns*na], -kW[j+k*na+l*ns*na]));

                      /* calculate the net proportion fecundity parameter acerageing over exposure groups */
                        for (i=0; i<nev[l]; i++)
                        {
                          netlambdacum[i+1] = netlambdacum[i] +  psi[i+nevcum[l]]*lambda0adj*(Ni[i] + Fi[i])/(Fj[j+k*na+l*ns*na]+Nj[j+k*na+l*ns*na]);
                        }

                      lambdaj[j+k*na+l*ns*na] = netlambdacum[nev[l]];

                        /* compute necessary components of variance function */
                        A[j+k*na+l*ns*na] = (1/(nss*wtss) + lambda0/muM0*(1+1/(nss*kM[j+k*na+l*ns*na])));
                        B[j+k*na+l*ns*na] = ( (1+ 1/(nss*kM[j+k*na+l*ns*na]))*(1 + 1/kW[j+k*na+l*ns*na]) - 1  );
                        Mjsq[j+k*na+l*ns*na] = Mj[j+k*na+l*ns*na]*Mj[j+k*na+l*ns*na];

                        varMj[j+k*na+l*ns*na]=A[j+k*na+l*ns*na]*Mj[j+k*na+l*ns*na] + B[j+k*na+l*ns*na]*Mjsq[j+k*na+l*ns*na];


                      /* calculate within blackfly larval population dynamics, biting on each age group
                      - set to equilibrium */
                        ydot[j+jumptoL1] = 0*y[j+jumptoL1];  /* L1 */
                        ydot[j+jumptoL2] = 0*y[j+jumptoL2]; /* L2 */
                        ydot[j+jumptoL3] = 0*y[j+jumptoL3]; /* L3 */

                        /* L1 density dependence incorporates overdispersion in adult worm burden as per Churcher et al (2006, 2008) but uses Filipe et al notatation (piV) and includes excess mortality on fly as a separate density dependence (i.e piV1 & piV2) */

                        double alpha;
                        alpha = Mj[j+k*na+l*ns*na]*aV/(Mj[j+k*na+l*ns*na]*aV+kW[j+k*na+l*ns*na]);


                        /* density-dependent establishment in blackfly */
                        PiV1j[j+k*na+l*ns*na] = pow((1-alpha)/(1-alpha*exp(-cV)), kW[j+k*na+l*ns*na])*(kW[j+k*na+l*ns*na]/(kW[j+k*na+l*ns*na]+(1-exp(-cV))*Mj[j+k*na+l*ns*na]*aV));

                        /* density-dependent excess mortality of blackly */
                        PiV2j[j+k*na+l*ns*na] = alphaV*Mj[j+k*na+l*ns*na];

                      if (k==0)
                      {

                        L1j[j+k*na+l*ns*na] =  y[j+jumptoL1] - y[j+jumptoL1] +
                          BR*deltaV0*PiV1j[j+k*na+l*ns*na]*Mj[j+k*na+l*ns*na]*pow(muV + muH + alphaF +PiV2j[j+k*na+l*ns*na]+sigma1, -1);
                      } else if (k==1)
                      {

                        L1j[j+k*na+l*ns*na] =  y[j+jumptoL1] - y[j+jumptoL1] +
                          BR*deltaV0*PiV1j[j+k*na+l*ns*na]*Mj[j+k*na+l*ns*na]*pow(muV + muH + alphaM + PiV2j[j+k*na+l*ns*na]+sigma1, -1);

                      }

                      L2j[j+k*na+l*ns*na] = y[j+jumptoL2] - y[j+jumptoL2] +
                        L1j[j+k*na+l*ns*na]*sigma1*pow( muV + sigma2, -1);

                      L3j[j+k*na+l*ns*na] = y[j+jumptoL3] - y[j+jumptoL3] +
                        L2j[j+k*na+l*ns*na]*sigma2*pow(aH*pow(g, -1) + muV + muL0,-1);

                    } //end j loop

                  }//end k loop

            } // end l loop

            return;

              }//end function

            void RungeKutta(double);

            /* *************************************************************
              main program
            *************************************************************** */

              // [[Rcpp::export]]
            Rcpp::List runEPIONCHO(std::vector<double> theta, int itervtn)
            {

            /* ********************************
              read in parameter values
            ******************************** */
              ABR0 = theta[0];
            effVC = theta[1];
            durVC = theta[2];
            ntrt1 = (int)theta[3];
            ntrt2 = (int)theta[4];
            ntrt3 = (int)theta[5];
            ntrt4 = (int)theta[6];
            ntrt5 = (int)theta[7];
            ftrt1 = theta[8];
            ftrt2 = theta[9];
            ftrt3 = theta[10];
            ftrt4 = theta[11];
            ftrt5 = theta[12];
            cov1 = theta[13];
            cov2 = theta[14];
            cov3 = theta[15];
            cov4 = theta[16];
            cov5 = theta[17];
            noncmp = theta[18];
            kW0 = theta[19];
            kW1 = theta[20];
            kM0 = theta[21];
            kM1 = theta[22];
            prepat = theta[23];
            deltaH0 = theta[24];
            deltaHinfty = theta[25];
            cH = theta[26];
            alpha0 = theta[27];
            beta0 = theta[28];
            lambda0 = theta[29];
            deltaV0 = theta[30];
            cV = theta[31];
            aV = theta[32];
            alphaV = theta[33];
            aH = theta[34];
            muL0 = theta[35];
            sigma1  = theta[36];
            sigma2 =  theta[37];
            h = theta[38];
            g = theta[39];
            muV = theta[40];

            E0 = theta[41];
            alphaF = theta[42];
            beta1Max = theta[43];
            gammabeta = theta[44];
            epsilonbeta = theta[45];
            nu = theta[46];
            omega = theta[47];
            mu1Max = theta[48];
            epsilonmu = theta[49];
            zeta = theta[50];

            LE = theta[51];
            LEMf = theta[52];

            nL4 = (int)theta[53];
            nW = (int)theta[54];
            nM = (int)theta[55];

            wtss = theta[56];
            nss = (int)theta[57];

                if (itervtn==0) {
                    ntrt1 = 0;
                    ntrt2 = 0;
                    ntrt3 = 0;
                    ntrt4 = 0;
                    ntrt5 = 0;
                    fcp = 0;
                    ncp = 1;
                } else { fcp = 2;
                    ncp = fcp+2;}

                /* resize vectors now that ncp and fcp are redefined */

                ntr.resize(fcp+2);
                ntr1.resize(fcp+2);
                ntr2.resize(fcp+2);
                ntr3.resize(fcp+2);
                ntr4.resize(fcp+2);
                ntr5.resize(fcp+2);
                nev.resize(fcp+2);
                nevcum.resize(fcp+3);
                ftr1.resize(fcp+2);
                ftr2.resize(fcp+2);
                ftr3.resize(fcp+2);
                ftr4.resize(fcp+2);
                ftr5.resize(fcp+2);

                pcp1.resize(ncp);
                pcp2.resize(ncp);
                pcp3.resize(ncp);
                pcp4.resize(ncp);
                pcp5.resize(ncp);

                therpcp1.resize(ncp);
                therpcp2.resize(ncp);
                therpcp3.resize(ncp);
                therpcp4.resize(ncp);
                therpcp5.resize(ncp);

            /* reszize storage vectors for model averages  now that ncp and ns are defined */
                Nj.resize(ncp*ns*na);
                Fj.resize(ncp*ns*na);
                Wj.resize(ncp*ns*na);
                WTj.resize(ncp*ns*na);
                Mj.resize(ncp*ns*na);
                Mpj.resize(ncp*ns*na);
                pFj.resize(ncp*ns*na);
                Phij.resize(ncp*ns*na);
                lambdaj.resize(ncp*ns*na);
                PiV1j.resize(ncp*ns*na);
                PiV2j.resize(ncp*ns*na);
                L1j.resize(ncp*ns*na);
                L2j.resize(ncp*ns*na);
                L3j.resize(ncp*ns*na);
                kM.resize(ncp*ns*na);
                kW.resize(ncp*ns*na);
                A.resize(ncp*ns*na);
                B.resize(ncp*ns*na);
                Mjsq.resize(ncp*ns*na);
                varMj.resize(ncp*ns*na);

                /* storage vector for states averaged over age */
                Nk.resize(ncp*ns);
                Fk.resize(ncp*ns);
                Wk.resize(ncp*ns);
                Mk.resize(ncp*ns);
                M5k.resize(ncp*ns);
                M20k.resize(ncp*ns);

                Mpk.resize(ncp*ns);
                Mp5k.resize(ncp*ns);
                Mp20k.resize(ncp*ns);
                Phik.resize(ncp*ns);
                lambdak.resize(ncp*ns);
                L1k.resize(ncp*ns);
                L2k.resize(ncp*ns);
                L3k.resize(ncp*ns);

                /* storage vector for states averaged over age & sex */
                Nl.resize(ncp);
                Fl.resize(ncp);
                Wl.resize(ncp);
                Ml.resize(ncp);
                M5l.resize(ncp);
                M20l.resize(ncp);

                Mpl.resize(ncp);
                Mp5l.resize(ncp);
                Mp20l.resize(ncp);
                Phil.resize(ncp);
                lambdal.resize(ncp);
                L1l.resize(ncp);
                L2l.resize(ncp);
                L3l.resize(ncp);

            /* derived parameter values */
            mu0 = nW*(1/LE);
            gammaL = nL4*(1/prepat);
            muM0 = nM*(1/LEMf);

            nst = nL4 + 2*nW + nM;
            BR = h/g;
            m0 = ABR0/BR;

            lambda0adj = lambda0*(alpha0+beta0+1/LE)/alpha0;

            L4iicum.resize(nL4+1);
            Niicum.resize(nW+1);
            Fiicum.resize(nW+1);
            lambdaiicum.resize(nW+1);
            Miicum.resize(nM+1);

            L4id.resize(nL4);
            Nid.resize(nW);
            Fid.resize(nW);
            Mid.resize(nM);

            startTreat=150;
            startVC = startTreat;


            /* ********************************
              initialize global array structures &
              and other globally defined variables
            ******************************** */

              /* sex structure - females (0) and males (1) */
              if (ns<2) {
                /* removes sex structuring to speed up sims */
                  Q = 1;
                  psex[0] = 1;
                  psex[1] = 0;
                  alphaM = alphaF;
                  q = 0;
              } else {
                psex[0] = rhoF;
                psex[1] = 1 - rhoF;
              }


            /*relative exposure of males versus females */
              EM = 1/(rhoF/Q + 1 - rhoF);
            EF = EM/Q;

            /* create age structure */
              fillAgeStruc(da, muH, amax, contra);

            /* create treatment structure */
              fillTreatStruc(ntrt1, ntrt2, ntrt3, ntrt4, ntrt5, ftrt1, ftrt2, ftrt3, ftrt4, ftrt5, fcp);

            /* calculate number events */
              calcEvents(ntr);
            tnev = nevcum[ncp];

            /* resize event vector with length (available at run time)
            tnev (and which requires an additional "end" constant
                  term for each compliance group)
            */
              event.resize(tnev+ncp);
            fillEventStruc();

            /* caluclate exposure structure
            * - first calculating normalization
            * factors */
              calcNormFacs();
            fillExpoStruc();

            /* create compliance structure */
                if (itervtn==0) {
                    pcmp1 = 1;
                } else {
            pcmp1 = std::fmax(cov1 - (1/(1-fcp))*(cov1+noncmp-1), 0);
            pcmp2 = std::fmax(cov2 -  (1/(1-fcp))*(cov2+noncmp-1), 0);
            pcmp3 = std::fmax(cov3 - (1/(1-fcp))*(cov3+noncmp-1), 0);
            pcmp4 = std::fmax(cov4 - (1/(1-fcp))*(cov4+noncmp-1), 0);
            pcmp5 = std::fmax(cov5 - (1/(1-fcp))*(cov5+noncmp-1), 0);

            psemi1 =  cov1 - std::fmax(cov1 - (1/(1-fcp))*(cov1+noncmp-1), 0);
            psemi2 =  cov2 - std::fmax(cov2 -  (1/(1-fcp))*(cov2+noncmp-1), 0);
            psemi3 =  cov3 - std::fmax(cov3 - (1/(1-fcp))*(cov3+noncmp-1), 0);
            psemi4 =  cov4 - std::fmax(cov4 - (1/(1-fcp))*(cov4+noncmp-1), 0);
            psemi5 =  cov5 - std::fmax(cov5 - (1/(1-fcp))*(cov5+noncmp-1), 0);
                }

            /* fill the treatment coverage of each
             compliance group */
              fillCompStruc();

            /* caluclate the acheived population-level
            * coverage (accounting for contraindicated age group) */
              if (fcp==0)
              {
                sumto = 2;
              }
            else
            {
              sumto = 3;
            }
            therpcpcum1.resize(sumto);
            therpcpcum2.resize(sumto);
            therpcpcum3.resize(sumto);
            therpcpcum4.resize(sumto);
            therpcpcum5.resize(sumto);

            calcCov();

            /* resize step functions depending on parameter values */
              stepivm.resize(tnev+ncp);
            stepmox.resize(tnev+ncp);
            expo.resize(tnev);

            /* resize temporally dynamic rate parameters */
              beta1.resize(tnev);
            mu1.resize(tnev);
            muM1.resize(tnev);
            psi.resize(tnev);

            /* ********************************
              debugging information
            ******************************** */
              if(DEBUG>0) {
                Rcpp::Rcout<< "running debugging checks... " << std::endl;
                Rcpp::Rcout<< "age structure summation = " << checkSum(pa) << std::endl;
                Rcpp::Rcout<< "total events = " << tnev << std::endl;
                for (int l=0; l<ncp; l++) {
                  if (l==0) {
                    Rcpp::Rcout<< "number of treatments in fully-adherent group " << l+1  << " = " << ntr[l] << std::endl;
                  } else if (l>0 && l <- (ncp-1)) {
                    Rcpp::Rcout<< "number of treatments in semi-adherent group " << l+1  << " = " << ntr[l] << std::endl;
                  } else
                  {
                    Rcpp::Rcout<< "number of treatments in non-adherent group " << l+1  << " = " << ntr[l] << std::endl;
                  }
                }
                double cumsumM[na+1], cumsumF[na+1];
                cumsumM[0] = 0;
                cumsumF[0] = 0;
                for (int i=0; i<na; i++) {
                  cumsumM[i+1] = cumsumM[i] + pa[i]*(OmegaM[i]/EM);
                  cumsumF[i+1] = cumsumF[i] + pa[i]*(OmegaF[i]/EF);
                }
                Rcpp::Rcout << "male exposure structure summation = " << cumsumM[na] << std::endl;
                Rcpp::Rcout << "female exposure structure summation = " << cumsumF[na] << std::endl;
                  if (itervtn > 0) {

                Rcpp::Rcout << "coverage summation before switch = " << checkSum(pcp1) << std::endl;
                Rcpp::Rcout << "coverage summation after first switch = " << checkSum(pcp2) << std::endl;
                Rcpp::Rcout << "coverage summation after second switch = " << checkSum(pcp3) << std::endl;
                Rcpp::Rcout << "coverage summation after third switch = " << checkSum(pcp4) << std::endl;
                Rcpp::Rcout << "coverage summation after fourth switch = " << checkSum(pcp5) << std::endl;
                Rcpp::Rcout << "population coverage before switch = " << actcov1 << std::endl;
                Rcpp::Rcout << "population coverage after first switch = " << actcov2 << std::endl;
                Rcpp::Rcout << "population coverage after second switch = " << actcov3 << std::endl;
                Rcpp::Rcout << "population coverage after third switch = " << actcov4 << std::endl;
                Rcpp::Rcout << "population coverage after fourth switch = " << actcov5 << std::endl;


                for (int i=0; i<(tnev+ncp); i++)
                {
                  if (i==0) {
                    Rcpp::Rcout << "event times = " << event[i];
                  } else if (i==(tnev+ncp-1)) {
                    Rcpp::Rcout << ", " << event[i] << std::endl;
                  } else{
                    Rcpp::Rcout << ", " << event[i];
                  }

                }
                  }

              }

            /* ***************************************************************
              initialize time variables, step size & number of equations in model
            ***************************************************************** */

              t=0;
            /* step size before intervention (every 1/2 year) */
              stepbi=.25;
            /* step size during intervention (every day: 1/365.25) */
            /* This controls the size of step for the algorithmn */
              stepdi=0.002737850787;
            /* step size after intervention (approx. every 28 days: 1/365.25*28) */
            /* This controls how fine-grained the output is */
              stepai=stepdi*28;
            double Everydi = stepai; //(approx. every month)
            double Everyai = stepai;
            /* number of model equations */
              std::vector<double> neqcum(ncp+1);
            neqcum[0] = 0;
            for (int l = 0; l<ncp; l++)
            {
              neqcum[l+1] = neqcum[l] + ns*na*(nst*nev[l]+nlst);
            }
            neq = neqcum[ncp];

            switcht1 = startTreat + (ntr1[0])*ftr1[0];
            switcht2 = startTreat + (ntr1[0])*ftr1[0] + (ntr2[0]*ftr2[0]);
            switcht3 = startTreat + (ntr1[0])*ftr1[0] + (ntr2[0]*ftr2[0]) + (ntr3[0]*ftr3[0]);
            switcht4 = startTreat + (ntr1[0])*ftr1[0] + (ntr2[0]*ftr2[0]) + (ntr3[0]*ftr3[0]) +
              (ntr4[0]*ftr4[0]);

            double maxt0, maxt1, maxt, treatLength;

                if (itervtn>0)
                {
                  treatLength = ntr1[0]*ftr1[0] + ntr2[0]*ftr2[0] + ntr3[0]*ftr3[0] +
                    ntr4[0]*ftr4[0] + ntr5[0]*ftr5[0];
            maxt0 = startTreat + ntr1[0]*ftr1[0] + ntr2[0]*ftr2[0] + ntr3[0]*ftr3[0] +
              ntr4[0]*ftr4[0] + ntr5[0]*ftr5[0] + survperiod;
            maxt1 = maxt0 + elimperiod;
            maxt = maxt0 + maxt1;
                } else {
                    maxt = startTreat;
                    maxt0 = maxt;
                    maxt1 = maxt;
                }

            if(DEBUG>0) {
              Rcpp::Rcout << "number of states in model = " << neq << std::endl;
            }

            /* resize y to be of approrpirate length */
              y.resize(neq), ydot.resize(neq);

            if(DEBUG>0) {
              Rcpp::Rcout << "initialising model equations..." << std::endl;
            }

            setStateIndicators();
            fillInits();
            Diff(y);

            /* create output array for string specified
            population averages */
              double pretreat = 0.5;
            double startOutput = (startTreat - pretreat);
                int testnrow, nrow;
                if (itervtn>0){
                    testnrow =  floor( (maxt0 - Everydi - startOutput) / Everydi )  + floor( (maxt1 - Everyai - maxt0) / Everyai );
                } else {
                    testnrow =  floor( (maxt0 - Everydi - startOutput) / Everydi ) + 1;
                }
            if (testnrow>0) {
              nrow = testnrow;
              if(DEBUG>0) {
                Rcpp::Rcout << "Length of output vector = " << nrow << std::endl;
              }
            }
            else { nrow = 1; }
            std::vector<double> tout(nrow);
            std::vector<double> L3out(nrow);
            std::vector<double> Mout(nrow);
            std::vector<double> M5out(nrow);
            std::vector<double> M20out(nrow);

            std::vector<double> Mpout(nrow);
            std::vector<double> Mp5out(nrow);
            std::vector<double> Mp20out(nrow);
            std::vector<double> Nout(nrow);
            std::vector<double> Fout(nrow);
            std::vector<double> Wout(nrow);
            std::vector<double> ABRout(nrow);
            std::vector<double> ATPout(nrow);
            std::vector<double> lambdaout(nrow);
            std::vector<double> Mfull(nrow);
            std::vector<double> Msemi1(nrow);
            std::vector<double> Msemi2(nrow);
            std::vector<double> Mnoncmp(nrow);
            std::vector<double> Mp5full(nrow);
            std::vector<double> Mp5semi1(nrow);
            std::vector<double> Mp5semi2(nrow);
            std::vector<double> Mp5noncmp(nrow);

            if(DEBUG>0) {
              Rcpp::Rcout << "running model..." << std::endl;
            }


            do {
              /*update step functions */
                updateSteps();
              /* update dynamic treatment parameters */
                updateParms();
              /*update ABR (for vector control) */
                updateABR(t, ABR0);
              /* update mean number of L3 */
                updateL3();
              /* update mean number of host parasite stages */
                updateHostStages();
              /* integrate to determine states at next time step */
                RungeKutta(stepbi);
              t+=stepbi;
            } while (t<(startTreat-pretreat));

            if(DEBUG>0) {
              Rcpp::Rcout << "Completed running to equilibrium..." << std::endl;
            }

            int counter = 0;
            /*  EDIT: @SOHANLON //CHANGING TIME DURING INTERVENTION NOT ALLOWED// */
            /*t = startTreat-pretreat;*/
            if(DEBUG>0) {
              Rcpp::Rcout << "Now recording outputs..." << std::endl;
              Rcpp::Rcout << "Time is currently:\t" << t << std::endl;
              Rcpp::Rcout << "Interventions will start at time:\t" << startTreat << std::endl;
            }

            do {
              /*update step functions */
                updateSteps();
              /* update dynamic treatment parameters */
                updateParms();
              /*update ABR (for vector control) */
                updateABR(t, ABR0);
              /* update mean number of L3 */
                updateL3();
              /* update mean number of host parasite stages */
                updateHostStages();
              /* integrate to determine states at next time step */
                RungeKutta(stepdi);
                t+=stepdi;

                if(DEBUG>1) {
                  if(counter % 1 == 0){
                    Rcpp::Rcout << "\rValue of time is:\t" << t;
                    Rcpp::Rcout << "\tValue of counter is:\t" << counter;
                  }
                }

              /* fill storeage matrix with current values */
                  if( floor((t-(startTreat-pretreat))/Everydi) > floor(((t-(startTreat-pretreat))-stepdi)/Everydi) )

                  {
                    tout[counter]=t-startTreat;
                    L3out[counter]=L3;
                    Mout[counter]=M;
                    M5out[counter]=M5;
                    M20out[counter]=M20;

                    Mpout[counter]=Mp;
                    Mp5out[counter]=Mp5;
                    Mp20out[counter]=Mp20;
                    Nout[counter] = N;
                    Fout[counter]=F;
                    Wout[counter] = W;
                    ABRout[counter] = ABR;
                    ATPout[counter] = ABR*L3;
                    lambdaout[counter] = lambda;
                    Mfull[counter] = Ml[0];
                    Msemi1[counter] = Ml[1];
                    Msemi2[counter] = Ml[2];
                    Mnoncmp[counter] = Ml[3];
                    Mp5full[counter] = Mp5l[0];
                    Mp5semi1[counter] = Mp5l[1];
                    Mp5semi2[counter] = Mp5l[2];
                    Mp5noncmp[counter] = Mp5l[3];
                    counter+=1;
                  }
            } while (t < maxt0);

            if(DEBUG>0) {
              Rcpp::Rcout << "\nCompleted intervention scenario + 1 year surviellance looking for pOTTIS." << std::endl;
              Rcpp::Rcout << "Interventions ran over:\t" << treatLength;
              Rcpp::Rcout << " years.\nTime is currently:\t" << t << std::endl;
            }

            if (itervtn>0) {
              if(DEBUG>0) {
                Rcpp::Rcout << "\nRunning post-treatment timeline to see if elimination acheived.\nTime:" << std::endl;
              }
              do {
                Rcpp::Rcout << "\t" << t;
                /*update step functions */
                updateSteps();
                /* update dynamic treatment parameters */
                updateParms();
                /*update ABR (for vector control) */
                updateABR(t, ABR0);
                /* update mean number of L3 */
                updateL3();
                /* update mean number of host parasite stages */
                updateHostStages();
                /* integrate to determine states at next time step */
                RungeKutta(stepai);
                /* print temp dyamic par */
                t+=stepai;

                if( floor((t-(startTreat-pretreat))/Everyai) > floor(((t-(startTreat-pretreat))-stepai)/Everyai) )

                {
                  if( counter < nrow )
                  {
                    tout[counter]=t-startTreat;
                    L3out[counter]=L3;
                    Mout[counter]=M;
                    M5out[counter]=M5;
                    M20out[counter]=M20;

                    Mpout[counter]=Mp;
                    Mp5out[counter]=Mp5;
                    Mp20out[counter]=Mp20;
                    Nout[counter] = N;
                    Fout[counter]=F;
                    Wout[counter] = W;
                    ABRout[counter] = ABR;
                    ATPout[counter] = ABR*L3;
                    lambdaout[counter] = lambda;
                    Mfull[counter] = Ml[0];
                    Msemi1[counter] = Ml[1];
                    Msemi2[counter] = Ml[2];
                    Mnoncmp[counter] = Ml[3];
                    Mp5full[counter] = Mp5l[0];
                    Mp5semi1[counter] = Mp5l[1];
                    Mp5semi2[counter] = Mp5l[2];
                    Mp5noncmp[counter] = Mp5l[3];
                    counter+=1;
                  }
                }

                if(DEBUG>1) {
                  if(counter % 1 == 0){
                    Rcpp::Rcout << "\rValue of time is:\t" << t;
                    Rcpp::Rcout << "\tValue of counter is:\t" << counter;
                  }
                }
              } while ( t<maxt1 );
              }


              Rcpp::Rcout << "\nCompleted integration." << std::endl;
              Rcpp::Rcout << "Final length of output vectors was " << nrow;
              Rcpp::Rcout << " whilst the program had " << tout.size();
              Rcpp::Rcout << " time steps." << std::endl;
              Rcpp::Rcout << "Time is currently:\t" << t << std::endl;

            return  Rcpp::List::create(Rcpp::Named("time") = tout,
                                       Rcpp::Named("L3") = L3out,
                                       Rcpp::Named("M") = Mout,
                                       Rcpp::Named("M5") = M5out,
                                       Rcpp::Named("M20") = M20out,
                                       Rcpp::Named("Mp") = Mpout,
                                       Rcpp::Named("Mp5") = Mp5out,
                                       Rcpp::Named("Mp20") = Mp20out,
                                       Rcpp::Named("N") = Nout,
                                       Rcpp::Named("F") = Fout,
                                       Rcpp::Named("W") = Wout,
                                       Rcpp::Named("ABR") = ABRout);

            }


            /* *************************************************************
              utility functions...
            *************************************************************** */

              /* *************************************************************
              function to create age structure
            *************************************************************** */
              void fillAgeStruc(double da, double muH, double amax, double contra)
            {
            int j;
            ucontracum[0] = 0;
            for (j = 0; j<na; j++)
            {
              age[j] = j*da+da*pow(2,-1);
              pa[j] = da*muH*exp(-muH*age[j])*pow(1-exp(-muH*amax), -1);
              if (age[j]<5 && age[j]<20)
              {
                pa5[j] = 0;
                pa20[j] = 0;
              }
              else if (age[j]>=5 && age[j]<20)
              {
                pa20[j] = 0;
                pa5[j] = da*muH*exp(-muH*age[j])*pow( (1-exp(-muH*amax)) - (1-exp(-muH*5)), -1);
              } else if (age[j]>=5 && age[j]>=20)
              {
                pa5[j] = da*muH*exp(-muH*age[j])*pow( (1-exp(-muH*amax)) - (1-exp(-muH*5)), -1);
                pa20[j] = da*muH*exp(-muH*age[j])*pow( (1-exp(-muH*amax)) - (1-exp(-muH*20)), -1);
              }
              if (age[j]<contra)
              {
                ucontracum[j+1] = ucontracum[j]+pa[j];
              }
              else
              {
                ucontracum[j+1] = ucontracum[j];
              }
            }
            ucontra = ucontracum[na];
              }

            /* *************************************************************
              function to fill number & frequency of treatments
            in each compliance group
            *************************************************************** */
              void fillTreatStruc(int ntrt1, int ntrt2, int ntrt3, int ntrt4, int ntrt5, double ftrt1, double ftrt2, double ftrt3, double ftrt4, double ftrt5, int fcp)
            {
            int l;
            for (l = 0; l<ncp; l++)
            {
              /* set the number of treatments for the complicance group 0
              that is treated every round */
                if (l == 0)
                {
                  ntr1[l] =  ntrt1;
                  ntr2[l] =  ntrt2;
                  ntr3[l] = ntrt3;
                  ntr4[l] = ntrt4;
                  ntr5[l] = ntrt5;
                  ntr[l] = ntr1[l] + ntr2[l] + ntr3[l] + ntr4[l] + ntr5[l];
                  ftr1[l] = ( ftrt1 );
                  ftr2[l] = ( ftrt2 );
                  ftr3[l] = ( ftrt3 );
                  ftr4[l] = ( ftrt4 );
                  ftr5[l] = ( ftrt5 );
                }
              /* set the number of treatments for the complicance group 1
              that is treated once every fcp rounds and begins treatment
              at the same time as the fully compliant group 0*/
                else if (l > 0 && l < (ncp-1) && (l % 2) > 0 ) /* odd compliance groups */
                {
                  ntr1[l] = (int) ceil( double(ntrt1)/fcp );
                  /* if preveuois number is even */
                    if ( ntrt1 % 2 == 0) {
                      ntr2[l] =  (int) ceil( double(ntrt2)/fcp );
                    } else {
                      ntr2[l] = (int) floor( double(ntrt2)/fcp );
                    }
                  if ( (ntrt1 + ntrt2) % 2 == 0) {
                    ntr3[l] = (int) ceil( double(ntrt3)/fcp);
                  } else {
                    ntr3[l] = (int) floor( double(ntrt3)/fcp);
                  }
                  if ( (ntrt1 + ntrt2 + ntrt3) % 2 == 0) {
                    ntr4[l] = (int) ceil( double(ntrt4)/fcp);
                  } else {
                    ntr4[l] = (int) floor( double(ntrt4)/fcp);
                  }
                  if ( (ntrt1 + ntrt2 + ntrt3 + ntrt4) % 2 == 0) {
                    ntr5[l] = (int) ceil( double(ntrt5)/fcp);
                  } else {
                    ntr5[l] = (int) floor( double(ntrt5)/fcp);
                  }
                  ntr[l] = ntr1[l] + ntr2[l] + ntr3[l] + ntr4[l] + ntr5[l];
                  ftr1[l] = ( fcp*ftrt1 );
                  ftr2[l] = ( fcp*ftrt2 );
                  ftr3[l] = ( fcp*ftrt3 );
                  ftr4[l] = ( fcp*ftrt4 );
                  ftr5[l] = ( fcp*ftrt5 );
                }
              /* set the number of treatments for the other complicance groups >1
              that are treated once every fcp rounds and begins treatment
              shifted forward by l-1 events */
                else if (l > 1 && l < (ncp-1) &&  (l % 2) == 0) /* even compliance groups */
                {
                  ntr1[l] = (int) floor( double(ntrt1)/fcp );
                  /* if preveuois number is even */
                    if ( (ntrt1) % 2 == 0) {
                      ntr2[l] =  (int) floor( double(ntrt2)/fcp );
                    } else {
                      ntr2[l] = (int) ceil( double(ntrt2)/fcp );
                    }
                  if ( (ntrt1 + ntrt2) % 2 == 0) {
                    ntr3[l] = (int) floor( double(ntrt3)/fcp);
                  } else {
                    ntr3[l] = (int) ceil( double(ntrt3)/fcp);
                  }
                  if ( (ntrt1 + ntrt2 + ntrt3) % 2 == 0) {
                    ntr4[l] = (int) floor( double(ntrt4)/fcp);
                  } else {
                    ntr4[l] = (int) ceil( double(ntrt4)/fcp);
                  }
                  if ( (ntrt1 + ntrt2 + ntrt3 + ntrt4) % 2 == 0) {
                    ntr5[l] = (int) floor( double(ntrt5) /fcp);
                  } else {
                    ntr5[l] = (int) ceil( double(ntrt5)/fcp);
                  }
                  ntr[l] = ntr1[l] + ntr2[l] + ntr3[l] + ntr4[l] + ntr5[l];
                  ftr1[l] = ( fcp*ftrt1 );
                  ftr2[l] = ( fcp*ftrt2 );
                  ftr3[l] = ( fcp*ftrt3 );
                  ftr4[l] = ( fcp*ftrt4 );
                  ftr5[l] = ( fcp*ftrt5 );
                }
              /* set the number and frequency of treatments for the non-compliant group) */
  else if (l == (ncp-1))
  {
    ntr1[l] = 0;
    ntr2[l] = 0;
    ntr3[l] = 0;
    ntr4[l] = 0;
    ntr5[l] = 0;
    ntr[l] = 0;
    ftr1[l] = fcp*ftrt1; /* arbitrary constants */
      ftr2[l] = fcp*ftrt2;
      ftr3[l] = fcp*ftrt3;
      ftr4[l] = fcp*ftrt4;
      ftr5[l] = fcp*ftrt5;
  }
            }
              }

            void calcNormFacs()
            {

            pacumM[0] = 0;
            pacumF[0] = 0;
            int j;
            for (j = 0; j<na; j++)
            {
              if (age[j] < q)
              {
                pacumF[j+1] = pacumF[j] + pa[j]*E0;
                pacumM[j+1] = pacumM[j] + pa[j]*E0;
              }
              else
              {
                pacumF[j+1] = pacumF[j] + pa[j]*exp(-alphaF*(age[j]-q));
                pacumM[j+1] = pacumM[j] + pa[j]*exp(-alphaM*(age[j]-q));
              }
            }
            gammaF = pow(pacumF[na], -1);
            gammaM = pow(pacumM[na], -1);
            }

            void fillExpoStruc()
            {
            int j;
            for (j = 0; j<na; j++)
            {
              if (age[j] < q)
              {
                OmegaF[j] = EF*gammaF*E0;
                OmegaM[j] = EM*gammaM*E0;
              }
              else
              {
                OmegaF[j] = EF*gammaF*exp(-alphaF*(age[j]-q));
                OmegaM[j] = EM*gammaM*exp(-alphaM*(age[j]-q));
              }
            }

            /* join exposure functions together into a single vector
            * females first folowed by males*/
              for (j = 0; j< (ns*na); j++)
              {
                if (j<na)
                {
                  Omega[j] = OmegaF[j];
                }
                else
                {
                  Omega[j] = OmegaM[j-na];
                }
              }
            }

            /* *************************************************************
              calculate number of events in eatch comoliance group
            and cumulative number of events - depends on ntr vector
            *************************************************************** */

              void calcEvents(std::vector<int> ntr)
            {
            /* events vector with a dummy event
            at time 0 and at a time far in the future */
              int l;
            nev[0] = 0;
            nevcum[0] = 0;
            for (l = 0; l<ncp; l++)
            {
              nev[l] = ntr[l] + 1;
              nevcum[l+1] = nevcum[l] + nev[l];
            }
              }

            /* *************************************************************
              function to create event structure
            *************************************************************** */
              void fillEventStruc()
            {
            for (int l = 0; l<ncp; l++)
            {
              int jumpi = nevcum[l]+l;

              /* set the time of events (treatments + dummy events)
              for the complicance group 0 that is treated every round */
                if (l==0)
                {
                  int ii;
                  for (int i = 0; i<(nev[l]+1); i++)
                  {
                    if (i==0) /* dummy event zero at time 0 */
                    {
                      event[i+jumpi] = 0;
                    }
                    else if (i==nev[l]) /* dummy events some time in the distant future  chenged for EPIONCHOv8 */
                    {
                      event[i+jumpi] = startTreat + ntr1[l]*ftr1[l] + ntr2[l]*ftr2[l] + ntr3[l]*ftr3[l] + ntr4[l]*ftr4[l] + ntr5[l]*ftr5[l] + cnst;
                    }
                    else if (i < (ntr1[l]+1)) /* events before switch */
                    {
                      ii = i - 1;
                      event[i + jumpi] = startTreat + ii*ftr1[l];
                    }
                    else if (i > ntr1[l] && i < (ntr1[l]+ntr2[l]+1)) /* events after the first switch */
                    {
                      ii = i - ntr1[l] - 1;
                      event[i + jumpi] = startTreat + (ntr1[l])*ftr1[l] + ii*ftr2[l];
                    }
                    else if (i > (ntr1[l]+ntr2[l]) && i< (ntr1[l]+ntr2[l]+ntr3[l]+1))
                    {
                      ii = i - ntr1[l] - ntr2[l] - 1;
                      event[i + jumpi] = startTreat + (ntr1[l])*ftr1[l] + ntr2[l]*ftr2[l] + ii*ftr3[l];
                    }
                    else if (i > (ntr1[l]+ntr2[l]+ntr3[l]) && i< (ntr1[l]+ntr2[l]+ntr3[l]+ntr4[l]+1))
                    {
                      ii = i - ntr1[l] - ntr2[l] - ntr3[l] - 1;
                      event[i + jumpi] = startTreat + (ntr1[l])*ftr1[l] + ntr2[l]*ftr2[l] + ntr3[l]*ftr3[l] + ii*ftr4[l];
                    }
                    else if (i > (ntr1[l]+ntr2[l]+ntr3[l]+ntr4[l]))
                    {
                      ii = i - ntr1[l] - ntr2[l] - ntr3[l] - ntr4[l] - 1;
                      event[i + jumpi] = startTreat + (ntr1[l])*ftr1[l] + ntr2[l]*ftr2[l] + ntr3[l]*ftr3[l] + ntr4[l]*ftr4[l] + ii*ftr5[l];
                                        }
                  }

                }

              /* set the times of events (treatments + dummy events) for
              the complicance groups that are treated once every fcp rounds -
                this is done by picking appropriate elements from the fully
              compliant event vector */

                else if (l>0)
                {
                  for (int i = 0; i<(nev[l]+1); i++)
                  {
                    if ( i == 0)
                    {
                      event[i+jumpi] = 0;
                    }
                    else if ((i > 0) & (i < nev[l]))
                    {
                      event[i+jumpi] = event[i + (l-1) + (i-1)*(fcp-1)];
                    }
                    else if (i == nev[l])
                    {
                      event[i+jumpi] = startTreat +  ntr1[l]*ftr1[l] + ntr2[l]*ftr2[l] + ntr3[l]*ftr3[l] +
                        ntr4[l]*ftr4[l] + ntr5[l]*ftr5[l] + cnst;
                    }

                  }
                }
            }
              }


            /* *************************************************************
              function to create compliance structure
            *************************************************************** */
              void fillCompStruc()
            {
            int l;
            for (l=0; l<ncp; l++)
            {
              if (fcp==0 && l==0)
              {
                pcp1[l] = pcmp1;
                therpcp1[l] = pcp1[l]*(1-ucontra);
                pcp2[l] = pcmp2;
                therpcp2[l] = pcp2[l]*(1-ucontra);
                pcp3[l] = pcmp3;
                therpcp3[l] = pcp3[l]*(1-ucontra);
                pcp4[l] = pcmp4;
                therpcp4[l] = pcp4[l]*(1-ucontra);
                pcp5[l] = pcmp5;
                therpcp5[l] = pcp5[l]*(1-ucontra);

              }
              else if (fcp==0 && l>0)
              {
                pcp1[l] = noncmp;
                therpcp1[l] = 0;
                pcp2[l] = noncmp;
                therpcp2[l] = 0;
                pcp3[l] = noncmp;
                therpcp3[l] = 0;
                pcp4[l] = noncmp;
                therpcp4[l] = 0;
                pcp5[l] = noncmp;
                therpcp5[l] = 0;
              }
              else if (fcp>0 && l==0)
              {
                pcp1[l] = pcmp1;
                therpcp1[l] = pcp1[l]*(1-ucontra);
                pcp2[l] = pcmp2;
                therpcp2[l] = pcp2[l]*(1-ucontra);
                pcp3[l] = pcmp3;
                therpcp3[l] = pcp3[l]*(1-ucontra);
                pcp4[l] = pcmp4;
                therpcp4[l] = pcp4[l]*(1-ucontra);
                pcp5[l] = pcmp5;
                therpcp5[l] = pcp5[l]*(1-ucontra);
              }
              else if (fcp>0 && l>0 && l<(ncp-1))
              {
                pcp1[l] = psemi1;
                therpcp1[l] = pcp1[l]*(1-ucontra);
                pcp2[l] = psemi2;
                therpcp2[l] = pcp2[l]*(1-ucontra);
                pcp3[l] = psemi3;
                therpcp3[l] = pcp3[l]*(1-ucontra);
                pcp4[l] = psemi4;
                therpcp4[l] = pcp4[l]*(1-ucontra);
                pcp5[l] = psemi5;
                therpcp5[l] = pcp5[l]*(1-ucontra);
              }
              else
              {
                pcp1[l] = noncmp;
                therpcp1[l] = 0;
                pcp2[l] = noncmp;
                therpcp2[l] = 0;
                pcp3[l] = noncmp;
                therpcp3[l] = 0;
                pcp4[l] = noncmp;
                therpcp4[l] = 0;
                pcp5[l] = noncmp;
                therpcp5[l] = 0;
              }
            }

              }

            void calcCov()
            {
            therpcpcum1[0] = 0;
            therpcpcum2[0] = 0;
            therpcpcum3[0] = 0;
            therpcpcum4[0] = 0;
            therpcpcum5[0] = 0;
            int l;
            for (l=0; l<(sumto-1); l++)
            {
              therpcpcum1[l+1]=therpcpcum1[l]+therpcp1[l];
              therpcpcum2[l+1]=therpcpcum2[l]+therpcp2[l];
              therpcpcum3[l+1]=therpcpcum3[l]+therpcp3[l];
              therpcpcum4[l+1]=therpcpcum4[l]+therpcp4[l];
              therpcpcum5[l+1]=therpcpcum5[l]+therpcp5[l];
            }

            actcov1 = therpcpcum1[sumto-1];
            actcov2 = therpcpcum2[sumto-1];
            actcov3 = therpcpcum3[sumto-1];
            actcov4 = therpcpcum4[sumto-1];
            actcov5 = therpcpcum5[sumto-1];
            }

            /* *************************************************************
              function to check summation to 1
            *************************************************************** */
              double checkSum(std::vector<double> x) {
                int n = x.size();
                double cumsum[n+1];
                cumsum[0] = 0;
                for (int i=0; i<n; i++)
                {
                  cumsum[i+1] = cumsum[i]+x[i];
                }
                return cumsum[n];
              }

            /* *************************************************************
              update step functions at each time step
            *************************************************************** */
              void updateSteps()
            {
            int l, jumpi1, jumpi2, i;
            for (l = 0; l < ncp; l++)
            {
              jumpi1 = nevcum[l]+l;
              jumpi2 = nevcum[l];
              for (i = 0; i < nev[l]; i++)
              {
                if (t > event[i+jumpi1] && t < event[i+1+jumpi1])
                {
                  expo[i+jumpi2] = 1;
                }
                else
                {
                  expo[i+jumpi2] = 0;
                }
              }
            }

            for (l = 0; l < ncp; l++)
            {
              int jumpi = nevcum[l]+l;
              for (i = 0; i < (nev[l]+1); i++)
              {
                if (i==0)
                {
                  stepivm[i+jumpi] = 0;
                }
                else if (i==nev[l])
                {
                  stepivm[i+jumpi] = 0;
                }
                else if (t > event[i+jumpi])
                {
                  stepivm[i+jumpi] = 1;
                }
                else
                {
                  stepivm[i+jumpi] = 0;
                }
              }
            }

            for (l = 0; l < ncp; l++)
            {
              int jumpi = nevcum[l]+l;
              for (i = 0; i < (nev[l]+1); i++)
              {
                if (i==0)
                {
                  stepmox[i+jumpi] = 0;
                }
                else if (i==nev[l])
                {
                  stepmox[i+jumpi] = 0;
                }
                else if (t > event[i+jumpi] && t < event[i+jumpi+1])
                {
                  stepmox[i+jumpi] = 1;
                }
                else
                {
                  stepmox[i+jumpi] = 0;
                }
              }
            }
              }


            /* *************************************************************
              update ABR  at each time step (for vector control, rectangular distribution)
            *************************************************************** */

              /* to make changes to ABR if it is called as an argument, ut must be called by
            reference not by value */
              void updateABR(double t, double ABR0)
            {
            /* save value at address ABR
            before the update */
              double ABRold = ABR0;
            double ABRnew;
            if (t<startVC)
            {
              ABRnew = ABRold;

            } else if (t>(startVC+durVC))
            {
              ABRnew = ABRold;
            } else
            {
              ABRnew = (1-effVC)*ABRold;
            }
            /* update the value of ABR */
              ABR = ABRnew;
            m = ABR/BR;
              }



            /* *************************************************************
              update average L3  at each time step
            *************************************************************** */
              void updateL3()
            {
            int l, k, j;
            std::vector<double> L3jcum(na+1);
            std::vector<double> L3kcum(ns+1);
            std::vector<double> L3lcum(ncp+1);
            L3jcum[0] = 0;
            L3kcum[0] = 0;
            L3lcum[0] = 0;
            for (l=0; l<ncp; l++)
            {
              for (k=0; k<ns; k++)
              {
                for (j=0; j<na; j++)
                {
                  L3jcum[j+1] = L3jcum[j] + L3j[j+k*na+l*ns*na]*pa[j]*Omega[j+k*na]*Omega[j+k*na];
                }
                L3k[k+l*ns] = L3jcum[na];
                L3kcum[k+1] = L3kcum[k] + L3k[k+l*ns]*psex[k];
              }
              L3l[l] = L3kcum[ns];
              if (t <= switcht1)
              {
                L3lcum[l+1] = L3lcum[l] + L3l[l]*pcp1[l];
              }
              else if (t>switcht1 && t<=switcht2)
              {
                L3lcum[l+1] = L3lcum[l] + L3l[l]*pcp2[l];
              } else if (t>switcht2 && t<=switcht3)
              {
                L3lcum[l+1] = L3lcum[l] + L3l[l]*pcp3[l];
              } else if (t>switcht3 && t<=switcht4)
              {
                L3lcum[l+1] = L3lcum[l] + L3l[l]*pcp4[l];
              } else
              {
                L3lcum[l+1] = L3lcum[l] + L3l[l]*pcp5[l];
              }
            }
            L3 = L3lcum[ncp];
              }

            /* *************************************************************
              update human parasite stage averages  at each time step
            *************************************************************** */
              void updateHostStages()
            {
            /* these are getting redefined with every step !*/
            int l,k,j;
            std::vector<double> Mjcum(ns*na*ncp+1);
            std::vector<double> Mkcum(ns*ncp+1);
            std::vector<double> Mlcum(ncp+1);
            Mjcum[0] = 0;
            Mkcum[0] = 0;
            Mlcum[0] = 0;
            std::vector<double> M5jcum(ns*na*ncp+1);
            std::vector<double> M5kcum(ns*ncp+1);
            std::vector<double> M5lcum(ncp+1);
            M5jcum[0] = 0;
            M5kcum[0] = 0;
            M5lcum[0] = 0;
            std::vector<double> M20jcum(ns*na*ncp+1);
            std::vector<double> M20kcum(ns*ncp+1);
            std::vector<double> M20lcum(ncp+1);
            M20jcum[0] = 0;
            M20kcum[0] = 0;
            M20lcum[0] = 0;
            std::vector<double> Mpjcum(ns*na*ncp+1);
            std::vector<double> Mpkcum(ns*ncp+1);
            std::vector<double> Mplcum(ncp+1);
            Mpjcum[0] = 0;
            Mpkcum[0] = 0;
            Mplcum[0] = 0;
            std::vector<double> Mp5jcum(ns*na*ncp+1);
            std::vector<double> Mp5kcum(ns*ncp+1);
            std::vector<double> Mp5lcum(ncp+1);
            Mp5jcum[0] = 0;
            Mp5kcum[0] = 0;
            Mp5lcum[0] = 0;
            std::vector<double> Mp20jcum(ns*na*ncp+1);
            std::vector<double> Mp20kcum(ns*ncp+1);
            std::vector<double> Mp20lcum(ncp+1);
            Mp20jcum[0] = 0;
            Mp20kcum[0] = 0;
            Mp20lcum[0] = 0;
            std::vector<double> Fjcum(ns*na*ncp+1);
            std::vector<double> Fkcum(ns*ncp+1);
            std::vector<double> Flcum(ncp+1);
            Fjcum[0] = 0;
            Fkcum[0] = 0;
            Flcum[0] = 0;
            std::vector<double> Njcum(ns*na*ncp+1);
            std::vector<double> Nkcum(ns*ncp+1);
            std::vector<double> Nlcum(ncp+1);
            Njcum[0] = 0;
            Nkcum[0] = 0;
            Nlcum[0] = 0;
            std::vector<double> Wjcum(ns*na*ncp+1);
            std::vector<double> Wkcum(ns*ncp+1);
            std::vector<double> Wlcum(ncp+1);
            Wjcum[0] = 0;
            Wkcum[0] = 0;
            Wlcum[0] = 0;
            std::vector<double> lambdajcum(ns*na*ncp+1);
            std::vector<double> lambdakcum(ns*ncp+1);
            std::vector<double> lambdalcum(ncp+1);
            lambdajcum[0] = 0;
            lambdakcum[0] = 0;
            lambdalcum[0] = 0;


                std::vector<double> Mjsqcum(ns*na*ncp+1);
                Mjsqcum[0] = 0;
                std::vector<double> varMjcum(ns*na*ncp+1);
                varMjcum[0] = 0;

            for (l = 0; l<ncp; l++)
            {
              for (k=0; k<ns; k++)
              {
                for (j=0; j<na; j++)
                {
                  Mjcum[j+1] = Mjcum[j] + Mj[j+k*na+l*ns*na]*pa[j];
                  M5jcum[j+1] = M5jcum[j] + Mj[j+k*na+l*ns*na]*pa5[j];
                  M20jcum[j+1] = M20jcum[j] + Mj[j+k*na+l*ns*na]*pa20[j];
                  Mpjcum[j+1] = Mpjcum[j] + Mpj[j+k*na+l*ns*na]*pa[j];
                  Mp5jcum[j+1] = Mp5jcum[j] + Mpj[j+k*na+l*ns*na]*pa5[j];
                  Mp20jcum[j+1] = Mp20jcum[j] + Mpj[j+k*na+l*ns*na]*pa20[j];
                  Njcum[j+1] = Njcum[j] + Nj[j+k*na+l*ns*na]*pa[j];
                  Fjcum[j+1] = Fjcum[j] + Fj[j+k*na+l*ns*na]*pa[j];
                  Wjcum[j+1] = Wjcum[j] + Wj[j+k*na+l*ns*na]*pa[j];

                    /* component s to construct variance assume to come from 5+pop */
                    Mjsqcum[j+1] = Mjsqcum[j] + Mjsq[j+k*na+l*ns*na]*pa5[j];
                    varMjcum[j+1] = varMjcum[j] + varMj[j+k*na+l*ns*na]*pa5[j];

                  lambdajcum[j+1] = lambdajcum[j] + lambdaj[l+k*na+l*ns*na]*pa[j];

                }
                Mk[k+l*ns] = Mjcum[na];
                Mkcum[k+1] = Mkcum[k] + Mk[k+l*ns]*psex[k];
                M5k[k+l*ns] = M5jcum[na];
                M5kcum[k+1] = M5kcum[k] + M5k[k+l*ns]*psex[k];
                M20k[k+l*ns] = M20jcum[na];
                M20kcum[k+1] = M20kcum[k] + M20k[k+l*ns]*psex[k];
                Mpk[k+l*ns] = Mpjcum[na];
                Mpkcum[k+1] = Mpkcum[k] + Mpk[k+l*ns]*psex[k];
                Mp5k[k+l*ns] = Mp5jcum[na];
                Mp5kcum[k+1] = Mp5kcum[k] + Mp5k[k+l*ns]*psex[k];
                Mp20k[k+l*ns] = Mp20jcum[na];
                Mp20kcum[k+1] = Mp20kcum[k] + Mp20k[k+l*ns]*psex[k];
                Nk[k+l*ns] = Njcum[na];
                Nkcum[k+1] = Nkcum[k] + Nk[k+l*ns]*psex[k];
                Fk[k+l*ns] = Fjcum[na];
                Fkcum[k+1] = Fkcum[k] + Fk[k+l*ns]*psex[k];
                Wk[k+l*ns] = Wjcum[na];
                Wkcum[k+1] = Wkcum[k] + Wk[k+l*ns]*psex[k];
                lambdak[k+l*ns] = lambdajcum[na];
                lambdakcum[k+1] = lambdakcum[k] + lambdak[k+l*ns]*psex[k];

              }
              Ml[l] = Mkcum[ns];
              M5l[l] = M5kcum[ns];
              M20l[l] = M20kcum[ns];
              Mpl[l] = Mpkcum[ns];
              Mp5l[l] = Mp5kcum[ns];
              Mp20l[l] = Mp20kcum[ns];
              Nl[l] = Nkcum[ns];
              Fl[l] = Fkcum[ns];
              Wl[l] = Wkcum[ns];
              lambdal[l] = lambdakcum[ns];

              if (t <= switcht1)
              {
                Mlcum[l+1] = Mlcum[l] + Ml[l]*pcp1[l];
                M5lcum[l+1] = M5lcum[l] + M5l[l]*pcp1[l];
                M20lcum[l+1] = M20lcum[l] + M20l[l]*pcp1[l];
                Mplcum[l+1] = Mplcum[l] + Mpl[l]*pcp1[l];
                Mp5lcum[l+1] = Mp5lcum[l] + Mp5l[l]*pcp1[l];
                Mp20lcum[l+1] = Mp20lcum[l] + Mp20l[l]*pcp1[l];
                Nlcum[l+1] = Nlcum[l] + Nl[l]*pcp1[l];
                Flcum[l+1] = Flcum[l] + Fl[l]*pcp1[l];
                Wlcum[l+1] = Wlcum[l] + Wl[l]*pcp1[l];
                lambdalcum[l+1] = lambdalcum[l]+lambdal[l]*pcp1[l];

              }
              else if (t>switcht1 && t<=switcht2)
              {
                Mlcum[l+1] = Mlcum[l] + Ml[l]*pcp2[l];
                M5lcum[l+1] = M5lcum[l] + M5l[l]*pcp2[l];
                M20lcum[l+1] = M20lcum[l] + M20l[l]*pcp2[l];

                Mplcum[l+1] = Mplcum[l] + Mpl[l]*pcp2[l];
                Mp5lcum[l+1] = Mp5lcum[l] + Mp5l[l]*pcp2[l];
                Mp20lcum[l+1] = Mp20lcum[l] + Mp20l[l]*pcp2[l];
                Nlcum[l+1] = Nlcum[l] + Nl[l]*pcp2[l];
                Flcum[l+1] = Flcum[l] + Fl[l]*pcp2[l];
                Wlcum[l+1] = Wlcum[l] + Wl[l]*pcp2[l];
                lambdalcum[l+1] = lambdalcum[l]+lambdal[l]*pcp2[l];

              }
              else if (t>switcht2 && t<=switcht3)
              {
                Mlcum[l+1] = Mlcum[l] + Ml[l]*pcp3[l];
                M5lcum[l+1] = M5lcum[l] + M5l[l]*pcp3[l];
                M20lcum[l+1] = M20lcum[l] + M20l[l]*pcp3[l];

                Mplcum[l+1] = Mplcum[l] + Mpl[l]*pcp3[l];
                Mp5lcum[l+1] = Mp5lcum[l] + Mp5l[l]*pcp3[l];
                Mp20lcum[l+1] = Mp20lcum[l] + Mp20l[l]*pcp3[l];
                Nlcum[l+1] = Nlcum[l] + Nl[l]*pcp3[l];
                Flcum[l+1] = Flcum[l] + Fl[l]*pcp3[l];
                Wlcum[l+1] = Wlcum[l] + Wl[l]*pcp3[l];
                lambdalcum[l+1] = lambdalcum[l]+lambdal[l]*pcp3[l];

              }    else if (t>switcht3 && t<=switcht4)
              {
                Mlcum[l+1] = Mlcum[l] + Ml[l]*pcp4[l];
                M5lcum[l+1] = M5lcum[l] + M5l[l]*pcp4[l];
                M20lcum[l+1] = M20lcum[l] + M20l[l]*pcp4[l];

                Mplcum[l+1] = Mplcum[l] + Mpl[l]*pcp4[l];
                Mp5lcum[l+1] = Mp5lcum[l] + Mp5l[l]*pcp4[l];
                Mp20lcum[l+1] = Mp20lcum[l] + Mp20l[l]*pcp4[l];
                Nlcum[l+1] = Nlcum[l] + Nl[l]*pcp4[l];
                Flcum[l+1] = Flcum[l] + Fl[l]*pcp4[l];
                Wlcum[l+1] = Wlcum[l] + Wl[l]*pcp4[l];
                lambdalcum[l+1] = lambdalcum[l]+lambdal[l]*pcp4[l];

              }
              else
              {
                Mlcum[l+1] = Mlcum[l] + Ml[l]*pcp5[l];
                M5lcum[l+1] = M5lcum[l] + M5l[l]*pcp5[l];
                M20lcum[l+1] = M20lcum[l] + M20l[l]*pcp5[l];

                Mplcum[l+1] = Mplcum[l] + Mpl[l]*pcp5[l];
                Mp5lcum[l+1] = Mp5lcum[l] + Mp5l[l]*pcp5[l];
                Mp20lcum[l+1] = Mp20lcum[l] + Mp20l[l]*pcp5[l];
                Nlcum[l+1] = Nlcum[l] + Nl[l]*pcp5[l];
                Flcum[l+1] = Flcum[l] + Fl[l]*pcp5[l];
                Wlcum[l+1] = Wlcum[l] + Wl[l]*pcp5[l];
                lambdalcum[l+1] = lambdalcum[l]+lambdal[l]*pcp5[l];

              }
            }
            M = Mlcum[ncp];
            M5 = M5lcum[ncp];
            M20 = M20lcum[ncp];

            Mp = Mplcum[ncp];
            Mp5 = Mp5lcum[ncp];
            Mp20 = Mp20lcum[ncp];
            N = Nlcum[ncp];
            F = Flcum[ncp];
            W = Wlcum[ncp];
            lambda = lambdalcum[ncp];

                Msq = Mjsqcum[na];
                varM = varMjcum[na] + Msq - M*M;

              }

            /* *************************************************************
              update parameter values at each time step
            *************************************************************** */
              void updateParms()
            {
            int l, i, ii;
            for (l = 0; l<ncp; l++)
            {
              int jump = nevcum[l];
              for (i=0; i<nev[l]; i++)
              {
                int len = nev[l]-i;
                double beta1i[len];
                double mu1i[len];
                double muM1i[len];
                double psii[len];
                double beta1cum[len+1];
                double mu1cum[len+1];
                double muM1cum[len+1];
                double psicum[len+1];
                beta1cum[0] = 0;
                mu1cum[0] = 0;
                muM1cum[0] = 0;
                psicum[0] = 0;
                int start = i+1;
                for (ii=start; ii<(nev[l]+1); ii++)
                {
                  if (mox<1)
                  {
                    if ( stepivm[ii+nevcum[l]+l] > 0 )
                    {
                      muM1i[ii-start] = pow( (t-event[ii+jump+l]) + nu, -omega);
                      beta1i[ii-start] = beta1Max*exp(-gammabeta*(t-event[ii+jump+l]) + epsilonbeta*(ii-start));
                      mu1i[ii-start] = mu1Max*exp(epsilonmu*(ii-start));
                      if (ii==start)
                      {
                        psii[ii-start] = 0;
                      }
                      else
                      {
                        psii[ii-start] = log(1-zeta);
                      }
                    }
                    else
                    {
                      muM1i[ii-start] = 0;
                      beta1i[ii-start] = 0;
                      mu1i[ii-start] = 0;
                      psii[ii-start] = 0;
                    }
                  }
                  else if (mox>0)
                  {
                    if ( stepmox[ii+nevcum[l]+l] > 0 )
                    {
                      muM1i[ii-start] = pow( (t-event[ii+jump+l]) + nu, -omega);
                      beta1i[ii-start] = beta1Max*exp(-gammabeta*(t-event[ii+jump+l]) + epsilonbeta*(ii-start));
                      mu1i[ii-start] = mu1Max*exp(epsilonmu*(ii-start));
                      if (ii==start)
                      {
                        psii[ii-start] = 0;
                      }
                      else
                      {
                        psii[ii-start] = log(1-zeta);
                      }
                    }
                    else
                    {
                      muM1i[ii-start] = 0;
                      beta1i[ii-start] = 0;
                      mu1i[ii-start] = 0;
                      psii[ii-start] = 0;
                    }
                  }
                  beta1cum[ii-start+1] = beta1cum[ii-start] + beta1i[ii-start];
                  mu1cum[ii-start+1] = mu1cum[ii-start] + mu1i[ii-start];
                  muM1cum[ii-start+1] = muM1cum[ii-start] + muM1i[ii-start];
                  psicum[ii-start+1] = psicum[ii-start] + psii[ii-start];
                }
                beta1[i+jump] = beta1cum[len];
                mu1[i+jump] = mu1cum[len];
                muM1[i+jump]= muM1cum[len];
                psi[i+jump] = exp(psicum[len]);
              }
            }
              }

            /* *************************************************************
              set state indicators for different worm age groups
            *************************************************************** */
              void setStateIndicators()
            {
            int ii;
            for (ii=0;ii<nL4;ii++)
            {
              L4id[ii] = ii;
            }

            for (ii=0;ii<nW;ii++)
            {
              Nid[ii] = nL4 + 0*nW + ii;
              Fid[ii] = nL4 + 1*nW + ii;
            }
            for (ii=0; ii<nM;ii++)
            {
              Mid[ii] = nL4+2*nW + ii;
            }
              }


            /* *************************************************************
              initial value generator
            *************************************************************** */
              void fillInits()
            {
            int l,k,j,i,ii;
            double jumplcum[ncp+1];
            jumplcum[0] = 0;
            for (l = 0; l<ncp; l++)
            {
              jumpl = ns*(nst*nev[l]*na+na*nlst);
              jumplcum[l+1] = jumplcum[l] + jumpl;
              for (k = 0; k < ns; k++)
              {
                jumpk = k*(nst*nev[l]*na+na*nlst);
                jumptoL1 = nst*nev[l]*na + jumplcum[l] + jumpk;
                jumptoL2 = nst*nev[l]*na + na + jumplcum[l] + jumpk;
                jumptoL3 = nst*nev[l]*na + 2*na + jumplcum[l] + jumpk;
                for (j = 0; j < na; j++)
                {
                  jumpj = j*nst*nev[l];
                  jump = jumplcum[l] + jumpk + jumpj;
                  for (i = 0; i < nev[l]; i++)
                  {
                    for (ii=0;ii<nL4;ii++)
                    {
                      y[i+L4id[ii]*nev[l]+jump] = 0;
                    }
                    for (ii=0;ii<nW;ii++)
                    {
                      if (i<1)
                      {
                        y[i+Nid[ii]*nev[l]+jump] = 1;
                        y[i+Fid[ii]*nev[l]+jump] = 1;
                      }
                      else
                      {
                        y[i+Nid[ii]*nev[l]+jump] = 0;
                        y[i+Fid[ii]*nev[l]+jump] = 0;
                      }

                    }// end 2nd ii loop
                    for (ii=0;ii<nM;ii++)
                    {
                      if (i<1) {
                        y[i+Mid[ii]*nev[l]+jump] = 1;
                      } else {
                        y[i+Mid[ii]*nev[l]+jump] = 0;
                      }
                    } // end 3rd ii loop

                  }// end i loop
                  y[j+jumptoL1] = 0.02;
                  y[j+jumptoL2] = 0.02;
                  y[j+jumptoL3] = 0.02;
                } // end j loop
              } // end k loop
            } // end l loop
              }




            /* *************************************************************
              Runga Kutta integrator
            *************************************************************** */
              void RungeKutta(double step)
            {
            int i;
            std::vector<double> ydot1(neq), ydot2(neq), ydot3(neq), ydot4(neq);
            std::vector<double> tmpy(neq), inity(neq);
            for (i=0; i<neq; i++) {
              inity[i] = y[i];
            }

            Diff(inity);
            for (i=0;i<neq;i++)
            {
              ydot1[i]=ydot[i];
              tmpy[i]=inity[i]+step*ydot1[i]/2;
            }

            Diff(tmpy);
            for (i=0;i<neq;i++)
            {
              ydot2[i]=ydot[i];
              tmpy[i]=inity[i]+step*ydot2[i]/2;
            }

            Diff(tmpy);
            for (i=0;i<neq;i++)
            {
              ydot3[i]=ydot[i];
              tmpy[i]=inity[i]+step*ydot3[i]/2;
            }

            Diff(tmpy);
            for (i=0;i<neq;i++)
            {
              ydot4[i]=ydot[i];
              tmpy[i]=inity[i]+(ydot1[i]/6 + ydot2[i]/3 + ydot3[i]/3 + ydot4[i]/6)*step;
            }
            for (i=0;i<neq;i++)
            {
              y[i] = tmpy[i];
            }
            return;
              }
