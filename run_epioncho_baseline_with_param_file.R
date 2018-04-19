#  Two arguments. 1) argument to set ABR, 2) present WD script was called from
args <- commandArgs(TRUE)

cat( args , sep = "\n" )

#  Function to install any missing packaeges upon load
pkg_load_or_install <- function( ... ){
  pkgs <- list(...)
  for( p in pkgs ){
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( p , character.only = TRUE ) ){
      #  If package was not able to be loaded then re-install
      install.packages( p , dependencies = TRUE )
      #  Load package after installing
      suppressPackageStartupMessages( library( p , character.only = TRUE ) )
    }
  }
}

#  Load or install pacakges
pkg_load_or_install("Rcpp","data.table")

# change to WD of script
setwd( args[1] )
#  does an output directory exist in the current WD?
if (! dir.exists("output") ) {
  dir.create("output")
}

#  Need to enable C++11 features in the compiler on the cluster
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

#  Compile model and load default parameter values
sourceCpp(file="EPIONCHOv2.cpp")

#  Construct list of default parameters
theta <- structure(list(ABR = 5000,
                        effVC = 0L,
                        durVC = 0L,

                        #  Treatment
                        ntrt1 = 30L,
                        ntrt2 = 0L,
                        ntrt3 = 0L,
                        ntrt4 = 0L,
                        ntrt5 = 0L,

                        ftrt1 = 1L,
                        ftrt2 = 1L,
                        ftrt3 = 1L,
                        ftrt4 = 1L,
                        ftrt5 = 1L,

                        cov1 = 0.8,
                        cov2 = 0.8,
                        cov3 = 0.8,
                        cov4 = 0.8,
                        cov5 = 0.8,

                        noncmp = 0.05,

                        kW0 = 0.27,
                        kW1 = 0.000803011838339735,
                        kM0 = 0.427173079725123,
                        kM1 = 0.0209161328822218,
                        prepat = 1.55997585928652,
                        deltaH0 = 0.119432572107205,
                        deltaHinfty = 0.00457074548385125,
                        cH = 0.0119049694560962,
                        alpha0 = 0.571604787155046,
                        beta0 = 0.308478333109165,
                        lambda0 = 0.75425667188659,
                        deltaV0 = 0.0169204898681346,
                        cV = 0.00453100323141645,
                        aV = 0.564059195143618,
                        alphaV = 0.435659080103379,
                        aH = 0.815583904282093,
                        muL0 = 62.4920639576994,
                        sigma1 = 73.7061712354941,
                        sigma2 = 135.895802924147,
                        h = 0.672584581065538,
                        g = 0.00992092335828835,
                        muV = 23.0282152589951,
                        E0 = 0.1,
                        alphaF.M = 0.00582443050496265,
                        beta1Max = 32.4,
                        gammabeta = 19.6,
                        epsilonbeta = 0L,
                        nu = 0.0096,
                        omega = 1.25,
                        mu1Max = 0L,
                        epsilonmu = 0L,
                        zeta = 0.35,
                        LE = 9.18559223539608,
                        LEMf = 1.47455868351041,
                        nL4 = 1L,
                        nW = 18L,
                        nM = 1L,
                        wtss = 2,
                        nss = 2),
                   .Names = c("ABR","effVC", "durVC", "ntrt1", "ntrt2", "ntrt3", "ntrt4", "ntrt5",
                              "ftrt1", "ftrt2", "ftrt3", "ftrt4", "ftrt5", "cov1", "cov2",
                              "cov3", "cov4", "cov5", "noncmp", "kW0", "kW1", "kM0", "kM1",
                              "prepat", "deltaH0", "deltaHinfty", "cH", "alpha0", "beta0",
                              "lambda0", "deltaV0", "cV", "aV", "alphaV", "aH", "muL0", "sigma1",
                              "sigma2", "h", "g", "muV", "E0", "alphaF.M", "beta1Max", "gammabeta",
                              "epsilonbeta", "nu", "omega", "mu1Max", "epsilonmu", "zeta",
                              "LE", "LEMf", "nL4", "nW", "nM", "wtss", "nss")
                   )
#  // END PARAMETER SET //


#  Parameter values to alter in this simulation
theta$ABR <- args[2]
theta$kW0 <- args[3]


#  Output data.table
dt <- data.table( "time" = numeric(),
                  "L3" = numeric(),
                  "M" = numeric(),
                  "M5" = numeric(),
                  "M20" = numeric(),
                  "Mp" = numeric(),
                  "Mp5" = numeric(),
                  "Mp20" = numeric(),
                  "N" = numeric(),
                  "F" = numeric(),
                  "W" = numeric(),
                  "ABR" = numeric(),
                  "k" = numeric() )

#  Run model
cat( paste0("Running ABR=" , theta$ABR, " and k=" , theta$kW0 , "\t" ,  system("echo $PARALLEL_SEQ") ) )
out <- runEPIONCHO(theta = as.double(theta), itervtn = 0)
dt <- rbind( dt , as.list( out ) )


#  Output file
write.table( dt ,
  file = paste0( "output/epioncho_renata_params_baseline___ABR_" , args[2] , "___k_" , args[3] , "___job_number_" , args[4] , ".txt" ) ,
  quote = FALSE ,
  sep = "\t" ,
  row.names = FALSE
  )
