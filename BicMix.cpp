#include <stdlib.h>
#include <ctime>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/LU>
#include "myHeader.cpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>
//#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;
using namespace std;

using Eigen::MatrixXd;

//const double PI  =3.141592653589793238462;
/*
Chuan Gao C++, modified by Bryan Quach
*/

/* 
usage
./BicMix --nf 50 --y ./sim_data/Bicluster/Y_noise1_2.txt --out result --sep space

Args:
    --nf: The max number of latent factors to estimate.
    --y: A gene (rows) x sample (columns) expression matrix with delimiter specified by --sep.
    --out: Output directory for results files. Assumes the directory already exists.
    --sep: Either `space` or `tab` to denote space or tab delimited values in --y.
    --a: Sparse factor analysis model parameter for local shrinkage.
    --b: Sparse factor analysis model parameter for local shrinkage.
    --c: Sparse factor analysis model parameter for factor-specific shrinkage.
    --d: Sparse factor analysis model parameter for factor-specific shrinkage.
    --e: Sparse factor analysis model parameter for global shrinkage.
    --f: Sparse factor analysis model parameter for global shrinkage.
    --alpha: Sparse factor analysis model parameter.
    --beta: Sparse factor analysis model parameter.
    --seed (deprecated)
    --itr: The max number of variational EM iterations to run.
    --interval: The number of warm start iterations before estimates begin saving to file.
    --output-interval: The iteration interval between output files of model estimates. 
        For example, if --interval=200 and --output-interval=50, model estimates would be output 
	to file every 50th iteration after 200.
*/

int main(int argc,char *argv[]){

    //IOFormat HeavyFmt(FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
    //IOFormat HeavyFmt(FullPrecision);
    cout.precision(10);
    
    // declare variables 
    int nf=50,nf2=0,s_n=0,d_y=0,i_seed=0,n_itr=10001,nswitch1=0,nswitch2=0,out_interval=100;
    double aall=200,ball=200;
    double a=aall,b=ball,c=aall,d=ball,g=aall,h=ball,alpha=1,beta=1;
    //double a=0.5,b=0.5,c=0.5,d=0.5,g=0.5,h=0.5,alpha=1,beta=1;
    
    string file_y,dir_out,sep;
    stringstream ss;

    sep="space";
    
    int interval=200;
    
    // read in argument
    string s_df="--nf",s_y="--y",s_out="--out",s_sep="--sep",s_a="--a",s_b="--b",s_c="--c",s_d="--d",s_e="--e",s_f="--f",s_alpha="--alpha",s_beta="--beta",s_seed="--seed",s_itr="--itr",s_interval="--interval",s_output_interval="--output-interval";
    for(int i=0;i<argc;i++){
        if(s_df.compare(argv[i])==0){nf=atoi(argv[i+1]);}
        if(s_y.compare(argv[i])==0){file_y=argv[i+1];}
        if(s_out.compare(argv[i])==0){dir_out=argv[i+1];}
        if(s_sep.compare(argv[i])==0){sep=argv[i+1];}
        if(s_a.compare(argv[i])==0){a=atof(argv[i+1]);}
        if(s_b.compare(argv[i])==0){b=atof(argv[i+1]);}
        if(s_c.compare(argv[i])==0){c=atof(argv[i+1]);}
        if(s_d.compare(argv[i])==0){d=atof(argv[i+1]);}
        if(s_e.compare(argv[i])==0){g=atof(argv[i+1]);}
        if(s_f.compare(argv[i])==0){h=atof(argv[i+1]);}
        if(s_alpha.compare(argv[i])==0){alpha=atof(argv[i+1]);}
        if(s_beta.compare(argv[i])==0){beta=atof(argv[i+1]);}
        if(s_itr.compare(argv[i])==0){n_itr=atoi(argv[i+1])+1;}
        if(s_interval.compare(argv[i])==0){interval=atoi(argv[i+1]);}
        if(s_output_interval.compare(argv[i])==0){out_interval=atoi(argv[i+1]);}
    }
     
    // calculate the sample size and the gene numbers 
    string line;
    string field;
    ifstream f;
    
    f.open(file_y.c_str());
    if (! f.is_open())
    {
        printf("Gene expression file open failed\n");exit(0);
    }
    getline(f,line);
    s_n++;
    istringstream iss(line);
    if(sep.compare("space")==0){
        while(getline(iss,field,' ')){d_y++;}
    }else if(sep.compare("tab")==0){
        while(getline(iss,field,'\t')){d_y++;}
    }else{
        cout << "Please specify a valid separator." << endl << endl;
    }
    while(getline(f,line)){s_n++;}


    // write command into file for later references, also write the dimension of the gene expression matrix 
    ss.str("");
    ss.clear();
    ss << dir_out << "/command.txt";
    ofstream f_com (ss.str().c_str());
    if (f_com.is_open()){
        for(int i=0;i<argc;i++){
            f_com << argv[i] << " ";
        }
        f_com << endl;
    }
    f_com << endl << "Y_TMP has dimension of " << s_n << " by " << d_y << endl << endl;
    
    // read in the Y matrix 
    MatrixXd Y_TMP=MatrixXd::Constant(s_n,d_y,0);
   
  
  
    f.clear();
    f.seekg (0, ios_base::beg);
    int i=0,j=0;
    while(getline(f,line)){
        istringstream iss(line);
        j=0;
        if(sep.compare("space")==0){
            while(getline(iss,field,' ')){
                Y_TMP(i,j)=atof(field.c_str());
                j++;
            }
            i++;
        }else{
            while(getline(iss,field,'\t')){
                Y_TMP(i,j)=atof(field.c_str());
                j++;
            }
            i++;
        }
    }
    f.close();
 
    
    // chunk Y to a smaller sub matrix to test if necessary
    //s_n=1000;
    MatrixXd Y=MatrixXd::Constant(s_n,d_y,0);
    //Y=(Y_TMP.block(0,0,s_n,d_y)).transpose();
    Y=Y_TMP.block(0,0,s_n,d_y);
    f_com << "Submatrix of Y has dimension of " << s_n << " by " << d_y << endl;
    f_com.close();
    
    //int snbak=s_n;
    //s_n=d_y;
    //d_y=snbak;
    VectorXd mu = VectorXd::Constant(s_n,0);
    for(int i=0;i<d_y;i++){
        mu += Y.col(i);
    }
    mu *= double(1)/d_y;

    for(int i=0;i<d_y;i++){
        Y.col(i) -= mu;
    }
    
    // mu
    ss.str("");
    ss.clear();
    ss << dir_out << "/mu";
    ofstream f_mu (ss.str().c_str());
    f_mu << mu << endl;
    f_mu.close();
       
    /*
    MatrixXd MatrixXd::Constant(d_y,nf,0); 
    MatrixXd VLM=MatrixXd::Constant(s_n,nf,0);
    */
    // Declare variables independent of factor number to prepare for the EM algorithm
    
    VectorXd psi_v = VectorXd::Constant(s_n,1);
    MatrixXd PSI=MatrixXd::Constant(s_n,s_n,0);
    MatrixXd PSI_INV=MatrixXd::Constant(s_n,s_n,0);
    MatrixXd LX=MatrixXd::Constant(s_n,d_y,0);
    
    // Declare variables that are dependent of the factor number
    nf2 = nf;
    int nt=nf;
    
    MatrixXd EX=MatrixXd::Constant(nt,d_y,0);
    MatrixXd TEX=MatrixXd::Constant(d_y,nt,0);

    MatrixXd EXX=MatrixXd::Constant(nt,nt,0);
    MatrixXd ELL=MatrixXd::Constant(nt,nt,0);

    MatrixXd LAM=MatrixXd::Constant(s_n,nt,0);
    MatrixXd THETA=MatrixXd::Constant(s_n,nf,1);
    MatrixXd DELTA=MatrixXd::Constant(s_n,nf,1);
    VectorXd PHI = VectorXd::Constant(nf,1);
    VectorXd TAU = VectorXd::Constant(nf,1);
    double nu = 1;
    double ETA = 1;
    double GAMMA = 1;
    
    VectorXd count_lam = VectorXd::Constant(nf,0);
    MatrixXd count_lam_M = MatrixXd::Constant(n_itr,nf,0);
    VectorXd index = VectorXd::Constant(nf,0);
    
    double nmix=2;
    double zi = double(1)/nmix;
    MatrixXd Z = MatrixXd::Constant(nmix,nf,zi);
    MatrixXd logZ = MatrixXd::Constant(nmix,nf,log(zi));
    MatrixXd LOGV = MatrixXd::Constant(nmix,1,log(zi));
    
    
    long seed;
    
    gsl_rng *r;  // random number generator
    r=gsl_rng_alloc(gsl_rng_mt19937);

    seed = time (NULL) * getpid();    
    //gsl_rng_set (r,2250814569120);
    gsl_rng_set (r, seed);
    //gsl_rng_set (r, 22713859596316);
    
    // set seed
    //srand (time(NULL));

    ss.str("");
    ss.clear();
    ss << dir_out << "/seed";
    ofstream f_seed (ss.str().c_str());
    f_seed << seed << endl;
    f_seed.close();
       
    
    PSI.diagonal() = psi_v;
    inv_psi(PSI,PSI_INV,s_n);

    // fill in the lambda matrix  
    for (int i=0; i<s_n; i++) {
        for(int j=0;j<nt;j++){
            LAM(i,j)=gsl_ran_gaussian(r,1);
        }
    }

    VectorXd lam_count_v = VectorXd::Constant(n_itr,0);
    VectorXd log_det = VectorXd::Constant(n_itr,0);
    VectorXd bic = VectorXd::Constant(n_itr,0);

    //declare and initialize parameters related to X
    double XI=1,VARSIG=1,OMEGA=1;
    VectorXd KAPPA = VectorXd::Constant(nf,1);
    VectorXd LAMX = VectorXd::Constant(nf,1);
    MatrixXd RHO=MatrixXd::Constant(nf,d_y,1);
    MatrixXd SIGMA=MatrixXd::Constant(nf,d_y,1);
    
    VectorXd count_x = VectorXd::Constant(nf,0);
    MatrixXd count_x_M = MatrixXd::Constant(n_itr,nf,0);

    MatrixXd O = MatrixXd::Constant(nmix,nf,zi);
    MatrixXd logO = MatrixXd::Constant(nmix,nf,log(zi));
    MatrixXd LOGVO = MatrixXd::Constant(nmix,1,log(zi));
    
    MatrixXd LPL = MatrixXd::Constant(nf,nf,0);
    VectorXd vLXL = VectorXd::Constant(s_n,0);
    //VectorXd VLVX = VectorXd::Constant(s_n,0);
    MatrixXd partR = MatrixXd::Constant(nf,d_y,0);
    MatrixXd partL = MatrixXd::Constant(s_n,nf,0);
    // fill in the EX matrix
    for (int i=0; i<nt; i++) {
        for(int j=0;j<d_y;j++){
            EX(i,j)=gsl_ran_gaussian(r,1);
        }
    }

    MatrixXd vx_m = MatrixXd::Constant(2,nf,1);
    
   
    EXX=EX*EX.transpose();
    ELL=LAM*LAM.transpose();
   
    VectorXd x_count_v = VectorXd::Constant(n_itr,0);
    LPL=LAM.transpose()*PSI_INV*LAM;  

    for(int itr=0;itr<n_itr;itr++){  

        if(itr<0){
            aall=1000,ball=1000;
            a=aall,b=ball,c=aall,d=ball,g=aall,h=ball,alpha=1,beta=1;
        }
        //}else{
        //    aall=0.5,ball=0.5;
        //    a=aall,b=ball,c=aall,d=ball,g=aall,h=ball,alpha=1,beta=1;
        //}


       int lam_count=0;
        for(int i=0;i<s_n;i++){
            for(int j=0;j<nf;j++){
                if(LAM(i,j)!=0){
                    lam_count++;
                }
            }
        }
        
        int x_count=0;
        for(int i=0;i<nf;i++){
            for(int j=0;j<d_y;j++){
                if(EX(i,j)!=0){
                    x_count++;
                }
            }
        }

        
        lam_count_v(itr)=lam_count;
        x_count_v(itr)=x_count;
        

        VARSIG=double(g+h)/(OMEGA+XI);
        OMEGA=double((d*nf+g))/(VARSIG+KAPPA.sum());
        for(int i=0;i<nf;i++){
            KAPPA(i)=double(c+d)/(LAMX(i)+OMEGA);
        }
                
        for(int i=0;i<nf;i++){
            double x_sum=0;
            double sum_c=d_y*b*O(0,i)+c-1-0.5*d_y*O(1,i);
            double at = 2*(KAPPA(i)+O(0,i)*(RHO.row(i).sum()));
        
            //x_sum=EXX(i,i);
            
            x_sum=EX.row(i).dot(EX.row(i));
            //}
            double bt = O(1,i)*x_sum;
            LAMX(i)=double(sum_c+sqrt(sum_c*sum_c+at*bt))/at;
            if(LAMX(i)<1e-10){
                LAMX(i)=1e-10;
            }
        }
        //cout << "LAMX " << endl << LAMX.head(5) << endl;
        //cout << "good1" << endl;
        
        //////////// Mazimization step
        GAMMA=double(g+h)/(ETA+nu);
        ETA=double((d*nf+g))/(GAMMA+TAU.sum());
        
        for(int i=0;i<nf;i++){
            TAU(i)=double(c+d)/(PHI(i)+ETA);
        }
        for(int i=0;i<nf;i++){
            double lam_sum=0;
            double sum_c=s_n*b*Z(0,i)+c-1-0.5*s_n*Z(1,i);
            double at = 2*(TAU(i)+Z(0,i)*(DELTA.col(i).sum()));
            
            //if(itr<nswitch1){
            //lam_sum=ELL(i,i);
                //}else{
            
            lam_sum=LAM.col(i).dot(LAM.col(i));
                //}
            
            double bt = Z(1,i)*lam_sum;
            PHI(i)=double(sum_c+sqrt(sum_c*sum_c+at*bt))/at;
            if(PHI(i)<1e-10){
              PHI(i)=1e-10;
            }
        }
        //cout << "good2" << endl;
        //cout << "PHI " << endl << PHI.head(5) << endl;
        
        // count the number of x that are set to all zero
        count_x.setZero();
        for(int i=0;i<nf;i++){
            for(int j=0;j<d_y;j++){
                if(EX(i,j)!=0){
                    count_x(i) +=  1;
                }
            }
        }
        count_x_M.row(itr)=count_x;
        // count the number of loadings that are set to all zero
        count_lam.setZero();
        for(int i=0;i<nf;i++){
            for(int j=0;j<s_n;j++){
                if(LAM(j,i)!=0){
                    count_lam(i) +=  1;
                }
            }
        }
        count_lam_M.row(itr)=count_lam;

            
        // Count the number of loadings that are active, either all zero or PHI_k zero will kill
        int count_nonzero = 0;
        for(int i=0;i<nf;i++){
            if(count_lam(i)!=0&&PHI(i)!=0&&count_x(i)!=0&&LAMX(i)!=0){
                //if(count_lam(i)!=0&&PHI(i)!=0){
                index(count_nonzero)=i;
                count_nonzero++;
            }
        }

        // remove inactive factors, loadings etc. and assign to new matrix
        if(count_nonzero != nf){
 
            nf=count_nonzero;
            nt=nf;
            
            MatrixXd vx_m2=MatrixXd::Constant(2,nf,0);
            /*
            MatrixXd VLM2=MatrixXd::Constant(s_n,nf,0);
            MatrixXd VXM2=MatrixXd::Constant(d_y,nf,0);
            */
            MatrixXd EX2=MatrixXd::Constant(nt,d_y,0);
            MatrixXd TEX2=MatrixXd::Constant(d_y,nt,0);
            //MatrixXd VX2=MatrixXd::Constant(nt,nt,0);
            MatrixXd EXX2=MatrixXd::Constant(nt,nt,0);
            MatrixXd ELL2=MatrixXd::Constant(nt,nt,0);
            
            MatrixXd LAM2=MatrixXd::Constant(s_n,nt,0);
            
            MatrixXd THETA2=MatrixXd::Constant(s_n,nf,0);
            MatrixXd DELTA2=MatrixXd::Constant(s_n,nf,0);
            VectorXd PHI2 = VectorXd::Constant(nf,0);
            VectorXd TAU2 = VectorXd::Constant(nf,0);
            MatrixXd Z2 = MatrixXd::Constant(nmix,nf,zi);
            MatrixXd logZ2 = MatrixXd::Constant(nmix,nf,log(zi));
            VectorXd count_lam2 = VectorXd::Constant(nf,0);
            MatrixXd count_lam_M2 = MatrixXd::Constant(n_itr,nf,0);
            VectorXd index2 = VectorXd::Constant(nf,0);
            //ID2.diagonal()=id_v2;
            
            // for EX related
            VectorXd KAPPA2 = VectorXd::Constant(nf,0);
            VectorXd LAMX2 = VectorXd::Constant(nf,0);
            MatrixXd RHO2=MatrixXd::Constant(nf,d_y,0);
            MatrixXd SIGMA2=MatrixXd::Constant(nf,d_y,0);
            
            VectorXd count_x2 = VectorXd::Constant(nf,0);
            MatrixXd count_x_M2 = MatrixXd::Constant(n_itr,nf,0);
            MatrixXd O2 = MatrixXd::Constant(nmix,nf,zi);
            MatrixXd logO2 = MatrixXd::Constant(nmix,nf,log(zi));
            //MatrixXd LOGVO2 = MatrixXd::Constant(nmix,1,log(0.5));
            
            MatrixXd LPL2 = MatrixXd::Constant(nf,nf,0);
            MatrixXd partR2 = MatrixXd::Constant(nf,d_y,0);
            MatrixXd partL2 = MatrixXd::Constant(s_n,nf,0);
 
            for(int i=0;i<nf;i++){
                
                vx_m2.col(i)=vx_m.col(index(i));
                /*
                VLM2.col(i)=VLM.col(index(i));
                VXM2.col(i)=VXM.col(index(i));
                */
                EX2.row(i)=EX.row(index(i));
                TEX2=EX2.transpose();
                for(int j=0;j<nf;j++){
                    EXX2(i,j)=EXX(index(i),index(j));
                    ELL2(i,j)=ELL(index(i),index(j));
                    LPL2(i,j)=LPL(index(i),index(j));
                    
                }
                //LP2.row(i)=LP.row(index(i));
                LAM2.col(i)=LAM.col(index(i));
                
                THETA2.col(i)=THETA.col(index(i));
                DELTA2.col(i)=DELTA.col(index(i));
                PHI2(i)=PHI(index(i));
                TAU2(i)=TAU(index(i));
                Z2.col(i)=Z.col(index(i));
                logZ2.col(i)=logZ.col(index(i));
                
                count_lam2(i)=count_lam(index(i));
                count_lam_M2.col(i)=count_lam_M.col(index(i));
                index2(i)=index(i);
                
                // EX related
                KAPPA2(i)=KAPPA(index(i));
                LAMX2(i)=LAMX(index(i));
                RHO2.row(i)=RHO.row(index(i));
                SIGMA2.row(i)=SIGMA.row(index(i));
                count_x2(i)=count_x(index(i));
                count_x_M2.col(i)=count_x_M.col(index(i));
                O2.col(i)=O.col(index(i));
                logO2.col(i)=logO.col(index(i));
                
               partR2.row(i)=partR.row(index(i));
                partL2.col(i)=partL.col(index(i));
                
            }
         
            // Assign the new parameters back
            
            vx_m=vx_m2;
            /*
            VLM=VLM2;
            VXM=VXM2;
            */
            EX=EX2;
            TEX=TEX2;
            //VX=VX2;
            EXX=EXX2;
            //LP=LP2;
            //ID=ID2;
            LAM=LAM2;
            ELL=ELL2;
            //LAM_BAK=LAM_BAK2;
            THETA=THETA2;
            DELTA=DELTA2;
            PHI=PHI2;
            TAU=TAU2;
            Z=Z2;
            logZ=logZ2;
            //LAM_TOP=LAM_TOP2;
            //LAM_BOT=LAM_BOT2;
            count_lam=count_lam2;
            count_lam_M=count_lam_M2;
            index=index2;
            // EX related
            KAPPA=KAPPA2;
            LAMX=LAMX2;
            RHO=RHO2;
            SIGMA=SIGMA2;
            count_x=count_x2;
            count_x_M=count_x_M2;
            O=O2;
            logO=logO2;
            
            LPL=LPL2;
            partR=partR2;
            partL=partL2;
            //partV=partV2;

          }
        
      
        ////////////////////////////////////////  EX related
   
        for(int i=0;i<d_y;i++){
            for(int j=0;j<nf;j++){
                RHO(j,i)=double((a+b))/(SIGMA(j,i)+LAMX(j));
            }
        }
        //cout << "good3" << endl;
        //cout << "RHO " << endl << RHO.block(0,0,5,5) << endl;
        // this has been 
        /*
        if(itr<nswitch2){
            for(int i=0;i<d_y;i++){
                for(int j=0;j<nf;j++){
                    double a23=(2*a-3);
                    SIGMA(j,i)=double(a23+sqrt(a23*a23+8*(EX(j,i)*EX(j,i)+VXM(i,j))*RHO(j,i)))/4/RHO(j,i);
                }
            }
            
        }else{
        */



        for(int i=0;i<d_y;i++){
            for(int j=0;j<nf;j++){
                double a23=(2*a-3);
                SIGMA(j,i)=double(a23+sqrt(a23*a23+8*(EX(j,i)*EX(j,i))*RHO(j,i)))/4/RHO(j,i);
                //SIGMA(j,i)=double(a23+sqrt(a23*a23+8*(EX(j,i)*EX(j,i)+VXM(i,j))*RHO(j,i)))/4/RHO(j,i);
                // if(SIGMA(j,i)<1e-50){
                //     SIGMA(j,i)=0;
                // }
            }
        }
            //}
            
        //cout << "SIGMA " << endl << SIGMA.block(0,0,5,5) << endl;
        
        // continue updating parameters
        for(int i=0;i<s_n;i++){
            for(int j=0;j<nf;j++){
                DELTA(i,j)=double((a+b))/(THETA(i,j)+PHI(j));
            }
        }
        /*
        if(itr<nswitch2){
            for(int i=0;i<s_n;i++){
                for(int j=0;j<nf;j++){
                    double a23=(2*a-3);
                    THETA(i,j)=double(a23+sqrt(a23*a23+8*(LAM(i,j)*LAM(i,j)+VLM(i,j))*DELTA(i,j)))/4/DELTA(i,j);
                }
            }
            
        }else{
        */

        for(int i=0;i<s_n;i++){
            for(int j=0;j<nf;j++){
                double a23=(2*a-3);
                THETA(i,j)=double(a23+sqrt(a23*a23+8*(LAM(i,j)*LAM(i,j))*DELTA(i,j)))/4/DELTA(i,j);
                //THETA(i,j)=double(a23+sqrt(a23*a23+8*(LAM(i,j)*LAM(i,j)+VLM(i,j))*DELTA(i,j)))/4/DELTA(i,j);
                // if(THETA(i,j)<1e-50){
                //     THETA(i,j)=0;
                // }
            }
        }
            //}
        
        //cout << "good4" << endl;
    
        //cout << "THETA " << endl << THETA.block(0,0,5,5) << endl;
        //////////////////////////////////////////////
       
        partR=LAM.transpose()*PSI_INV*Y;
        //cout << "partR " << endl << partR.block(0,0,5,5) << endl;
        
        //cout << "PSI_INV " << PSI_INV.block(0,0,10,10) << endl;
       
        VectorXd indexALL = VectorXd::Constant(nf,0);
        EXX.setZero();
        //EX.setZero();
        
        MatrixXd partV=MatrixXd::Constant(nf,nf,0);
        
        //cout << "indexLAM " << endl << indexLAM << endl;
        for(int j=0;j<d_y;j++){
            int count_indexALL=0;
            
            for(int i=0;i<nf;i++){
                if(SIGMA(i,j)!=0&&LAMX(i)!=0){
                    indexALL(count_indexALL)=i;
                    count_indexALL++;
                }
            }
  
            //cout << "indexSIG " << endl << indexSIG << endl;
            //cout << "indexALL " << endl << indexALL << endl;
            
            if(count_indexALL==0){
                EX.col(j).setZero();
                continue;
            }
            
            partV.setZero();
            
            for(int i1=0;i1<count_indexALL;i1++){
                for(int i2=0;i2<count_indexALL;i2++){
                    partV(indexALL(i1),indexALL(i2)) = LPL(indexALL(i1),indexALL(i2));
                }
                partV(indexALL(i1),indexALL(i1)) += O(0,indexALL(i1))/SIGMA(indexALL(i1),j)+O(1,indexALL(i1))/LAMX(indexALL(i1));
            }
            
            // if(itr>1303){
            //     if(j==233){
            //         cout << "partV X for gene " << j << endl;
            //         cout << partV << endl;
            //     }
            // }


            MatrixXd partVI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
            for(int i1=0;i1<count_indexALL;i1++){
                for(int i2=0;i2<count_indexALL;i2++){
                    partVI(i1,i2)=partV(indexALL(i1),indexALL(i2));
                }
            }
            MatrixXd partRI=MatrixXd::Constant(count_indexALL,1,0);
            MatrixXd EXI=MatrixXd::Constant(count_indexALL,1,0);
            MatrixXd EXXI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
            cpy_col_matrix(partRI,partR,indexALL,count_indexALL,0,j);


            MatrixXd IDNF = MatrixXd::Identity(count_indexALL, count_indexALL);
            MatrixXd vx=MatrixXd::Constant(count_indexALL,count_indexALL,0);

            vx = partVI.lu().solve(IDNF);
            EXI=vx*partRI;
            /*
            for(int i=0;i<count_indexALL;i++){
                VXM(j,indexALL(i))=vx(i,i);
            }
            */
            
            for(int i=0;i<count_indexALL;i++){
                if(SIGMA(indexALL(i),j)==0){
                    EXI(i)=0;
                }
            }
            
            cpy_col_matrix_bak(EX,EXI,indexALL,count_indexALL,j,0);
            EXXI=EXI*EXI.transpose();
            for(int i1=0;i1<count_indexALL;i1++){
                for(int i2=0;i2<count_indexALL;i2++){
                    EXX(indexALL(i1),indexALL(i2)) += EXXI(i1,i2)+vx(i1,i2);
                    //EXX(indexALL(i1),indexALL(i2)) += EXXI(i1,i2);
                }
            }

        }
    
        for(int i=0;i<nf;i++){
            for(int j=0;j<d_y;j++){
                if(SIGMA(i,j)==0){
                    EX(i,j)=0;
                }
            }
        }
        //cout << "good5" << endl;

        ///////////////////////////////////////////////////// to calculate Lambda
        
      
        //VectorXd indexALL = VectorXd::Constant(nf,0);
        indexALL.setZero();
        partL=PSI_INV*Y*EX.transpose();
        //cout << "partL " << endl << partL.block(0,0,5,5) << endl;
        //LPL=LAM.transpose()*PSI_INV*LAM;
        LPL.setZero();
        vLXL.setZero();
        ELL.setZero();
        //LAM.setZero();
        
        //MatrixXd partV=MatrixXd::Constant(nf,nf,0);
        
        for(int j=0;j<s_n;j++){
            int count_indexALL=0;
            for(int i=0;i<nf;i++){              
                if(THETA(j,i)!=0 && PHI(i)!=0){
                    indexALL(count_indexALL)=i;
                    count_indexALL++;
                }
            }
            
            if(count_indexALL==0){
                LAM.row(j).setZero();
                continue;
            }
            
            partV.setZero();
            for(int i1=0;i1<count_indexALL;i1++){
                for(int i2=0;i2<count_indexALL;i2++){
                    partV(indexALL(i1),indexALL(i2)) = PSI_INV(j,j)*EXX(indexALL(i1),indexALL(i2));
                }
                partV(indexALL(i1),indexALL(i1)) += Z(0,indexALL(i1))/THETA(j,indexALL(i1))+Z(1,indexALL(i1))/PHI(indexALL(i1));
            }  
       
            MatrixXd partVI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
            for(int i1=0;i1<count_indexALL;i1++){
                for(int i2=0;i2<count_indexALL;i2++){
                    partVI(i1,i2) = partV(indexALL(i1),indexALL(i2));
                }
            }
            
    
            
            MatrixXd partLI=MatrixXd::Constant(1,count_indexALL,0);
            MatrixXd LAMI=MatrixXd::Constant(1,count_indexALL,0);
            MatrixXd LLI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
            cpy_row_matrix(partLI,partL,indexALL,count_indexALL,j);
            
            MatrixXd IDNF = MatrixXd::Identity(count_indexALL, count_indexALL);

            //LDLT<MatrixXd> ldltOfA(partVI);
            //MatrixXd vl=ldltOfA.solve(IDNF);
            
            MatrixXd vl = partVI.lu().solve(IDNF);

            // if(itr>1303){
            //     if(j==233){
            //         cout << "partV lam for gene " << j << endl;
            //         cout << partV << endl;
            //     }
            // }

            //MatrixXd vl = partV.inverse();
            LAMI=partLI*vl;


            
            /*
            for(int i=0;i<count_indexALL;i++){
                VLM(j,indexALL(i))=vl(i,i);
            }
            */
            
            for(int i=0;i<count_indexALL;i++){
                if(THETA(j,indexALL(i))==0){
                    LAMI(i)=0;
                }
            }
            

            cpy_row_matrix_bak(LAM,LAMI,indexALL,count_indexALL,j);
            LLI=LAMI.transpose()*LAMI;


            for(int i1=0;i1<count_indexALL;i1++){
                for(int i2=0;i2<count_indexALL;i2++){
                    LPL(indexALL(i1),indexALL(i2)) += PSI_INV(j,j)*(LLI(i1,i2)+vl(i1,i2));
                    ELL(indexALL(i1),indexALL(i2)) += LLI(i1,i2)+vl(i1,i2);
                    //LPL(indexALL(i1),indexALL(i2)) += PSI_INV(j,j)*(LLI(i1,i2));
                    //ELL(indexALL(i1),indexALL(i2)) += LLI(i1,i2);
                    vLXL(j) += vl(i1,i2)*EXX(indexALL(i1),indexALL(i2));
                    
                }
            }
        }

        //cout << "good6" << endl;
   
     //cout << "LAM " << endl << LAM.block(0,0,10,10) << endl;
     
         // If theta=0, then LAM=0
        
        for(int i=0;i<s_n;i++){
            for(int j=0;j<nf;j++){
                if(THETA(i,j)==0){
                    LAM(i,j)=0;
                }
            }
        }
        
        // logO

        
        for(int i=0;i<nf;i++){
            logO(0,i)=LOGVO(0,0);
            logO(1,i)=LOGVO(1,0);
            for(int j=0;j<d_y;j++){
                logO(0,i)=logO(0,i)+log_norm(EX(i,j),0,SIGMA(i,j))+log_gamma(SIGMA(i,j),a,RHO(i,j))+log_gamma(RHO(i,j),b,LAMX(i));
                logO(1,i)=logO(1,i)+log_norm(EX(i,j),0,LAMX(i));
            }
        }
        
        // O
        for(int i=0;i<nf;i++){
            O(0,i)=double(1)/(1+exp(logO(1,i)-logO(0,i)));
            O(1,i)=1-O(0,i);
        }
        
        double ps1=alpha;
        double ps2=beta;
        
        // Probability
        for(int i=0;i<nf;i++){
            ps1=ps1+O(0,i);
            ps2=ps2+O(1,i);
        }
        double dgama = gsl_sf_psi(ps1+ps2);
        LOGVO(0,0)=gsl_sf_psi(ps1)-dgama;
        LOGVO(1,0)=gsl_sf_psi(ps2)-dgama;
         
        //cout << "O " << endl << O.block(0,0,2,10) << endl;
        
        // logZ 
        for(int i=0;i<nf;i++){
            logZ(0,i)=LOGV(0,0);
            logZ(1,i)=LOGV(1,0);
            for(int j=0;j<s_n;j++){
                logZ(0,i)=logZ(0,i)+log_norm(LAM(j,i),0,THETA(j,i))+log_gamma(THETA(j,i),a,DELTA(j,i))+log_gamma(DELTA(j,i),b,PHI(i));
                logZ(1,i)=logZ(1,i)+log_norm(LAM(j,i),0,PHI(i));
            }
        }
        
        // Z 
        for(int i=0;i<nf;i++){
            Z(0,i)=double(1)/(1+exp(logZ(1,i)-logZ(0,i)));
            Z(1,i)=1-Z(0,i);
        }
        
        ps1=alpha;
        ps2=beta;
        
        // Probability 
        for(int i=0;i<nf;i++){
            ps1=ps1+Z(0,i);
            ps2=ps2+Z(1,i);
        }
        dgama = gsl_sf_psi(ps1+ps2);
        LOGV(0,0)=gsl_sf_psi(ps1)-dgama;
        LOGV(1,0)=gsl_sf_psi(ps2)-dgama;
        //LOGV(0,0)=log(ps1)-log(ps1+ps2);
        //LOGV(1,0)=log(ps2)-log(ps1+ps2);

        //cout << "Z " << endl << Z.block(0,0,2,10) << endl;
        //cout << "vLXL " << endl << vLXL.block(0,0,10,10) << endl;
        //cout << "EXX " << endl << EXX.block(0,0,10,10) << endl;
             
        // PSI 
        LX=LAM*EX;
        for(int i=0;i<s_n;i++){
            //PSI(i,i)=Y.row(i).dot(Y.row(i))-2*(LX.row(i)).dot(Y.row(i))+(LAM.row(i)*EXX).dot(LAM.row(i))+vLXL(i);
            PSI(i,i)=(0.5*(Y.row(i).dot(Y.row(i))-2*(LX.row(i)).dot(Y.row(i))+(LAM.row(i)*EXX).dot(LAM.row(i)))+1)/(double(d_y)/2+1);
        }
        //PSI=PSI/d_y;
        inv_psi(PSI,PSI_INV,s_n);
        /*
        cout << "PSI " << endl << PSI.diagonal().head(5) << endl;
        cout << "LPL " << endl << LPL.block(0,0,5,5) << endl;
         */
        //cout << "EXX " << endl << EXX.block(0,0,10,10) << endl;
        
        // LAM active 
 
        for(int i=0;i<s_n;i++){
            log_det(itr) += log(PSI(i,i));
        }
        log_det(itr)=-0.5*d_y*log_det(itr);
        bic(itr)=2*log_det(itr)+lam_count*log(2*PI)+lam_count*log(s_n);
        

    

        MatrixXd range_LAM = MatrixXd::Constant(2,nf,0);
        MatrixXd range_EX = MatrixXd::Constant(2,nf,0);    
         
        range_colwise(range_LAM,LAM,s_n,nf);
        range_rowwise(range_EX,EX,nf,d_y);


        if(itr%10==0){
            cout << "Iteration............................................................... " << itr << endl;
            cout << "Number of components................................................. " << nf << endl;
                    
            //cout << "range LAM " << LAM.minCoeff() << " " << LAM.maxCoeff() << endl;
            //cout << "range EX " << EX.minCoeff() << " " << EX.maxCoeff() << endl;
            cout << "Number of genes in each loading \n" << count_lam.transpose() << endl;
            cout << "Number of samples in each factor \n" << count_x.transpose() << endl;
            
            
            //cout << "range LAM ori\n" << std::setprecision(3) << range_LAM_ori << endl;
            //cout << "range LAM\n" << std::setprecision(3) << range_LAM << endl;
            //cout << "range EX ori\n" << std::setprecision(3) << range_EX_ori << endl;
            //cout << "range EX\n" << std::setprecision(3) << range_EX << endl;

            //cout << "before PHI update " << endl << PHI.transpose() << endl;
            //cout << "before LAMX update " << endl << LAMX.transpose() << endl;
            
            /*
            cout << "vx_m\n " << vx_m << endl;
            */
        }



        if(itr>(interval-1)){

            
            double count_lam_each=0,count_x_each=0;
            for(int i=0;i<nf;i++){
                if(count_lam_M(itr,i)-count_lam_M(itr-interval,i)>-1){
                    count_lam_each++;
                }
                if(count_x_M(itr,i)-count_x_M(itr-interval,i)>-1){
                    count_x_each++;
                }
            }

            if(itr%out_interval==0){
                    ss.str("");
                ss.clear();
                ss << dir_out << "/LAM_" << itr;
                ofstream f_lam (ss.str().c_str());
                if (f_lam.is_open()){
                    f_lam << LAM << endl;
                }
                f_lam.close();
      
               
                ss.str("");
                ss.clear();
                ss << dir_out << "/Z_" << itr;
                ofstream f_Z (ss.str().c_str());
                if (f_Z.is_open()){
                    f_Z << Z << endl;
                }
                f_Z.close();
                
                
                ss.str("");
                ss.clear();
                ss << dir_out << "/EX_" << itr;
                ofstream f_EX (ss.str().c_str());
                if (f_EX.is_open()){
                    f_EX << EX << endl;
                }
                f_EX.close();

                ss.str("");
                ss.clear();
                ss << dir_out << "/EXX_" << itr;
                ofstream f_EXX (ss.str().c_str());
                if (f_EXX.is_open()){
                    f_EXX << EXX << endl;
                }
                f_EXX.close();


                ss.str("");
                ss.clear();
                ss << dir_out << "/O_" << itr;
                ofstream f_O (ss.str().c_str());
                if (f_O.is_open()){
                    f_O << O << endl;
                }
                f_O.close();

                ss.str("");
                ss.clear();
                ss << dir_out << "/PSI_" << itr;
                ofstream f_PSI (ss.str().c_str());
                if (f_PSI.is_open()){
                    f_PSI << PSI.diagonal() << endl;
                }
                f_PSI.close();

            }
            if((count_lam_each==nf&&count_x_each==nf)||itr==n_itr){
                //if(abs(log_det(itr)-log_det(itr-50))<0.05){
                    //if(itr==n_itr){
            

                ss.str("");
                ss.clear();
                ss << dir_out << "/LAM";
                ofstream f_lam (ss.str().c_str());
                if (f_lam.is_open()){
                    f_lam << LAM << endl;
                }
                f_lam.close();
      
               
                ss.str("");
                ss.clear();
                ss << dir_out << "/Z";
                ofstream f_Z (ss.str().c_str());
                if (f_Z.is_open()){
                    f_Z << Z << endl;
                }
                f_Z.close();
                
                
                ss.str("");
                ss.clear();
                ss << dir_out << "/EX";
                ofstream f_EX (ss.str().c_str());
                if (f_EX.is_open()){
                    f_EX << EX << endl;
                }
                f_EX.close();

                ss.str("");
                ss.clear();
                ss << dir_out << "/EXX";
                ofstream f_EXX (ss.str().c_str());
                if (f_EXX.is_open()){
                    f_EXX << EXX << endl;
                }
                f_EXX.close();


                ss.str("");
                ss.clear();
                ss << dir_out << "/O";
                ofstream f_O (ss.str().c_str());
                if (f_O.is_open()){
                    f_O << O << endl;
                }
                f_O.close();

                ss.str("");
                ss.clear();
                ss << dir_out << "/PSI";
                ofstream f_PSI (ss.str().c_str());
                if (f_PSI.is_open()){
                    f_PSI << PSI.diagonal() << endl;
                }
                f_PSI.close();

                exit(0);
            }
        }
    }
}



