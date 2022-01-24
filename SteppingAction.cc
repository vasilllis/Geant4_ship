#include "RunAction.hh"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "G4SteppingManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VParticleChange.hh"
#include "G4VProcess.hh"
#include "G4Step.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4SystemOfUnits.hh"

//extern std::ofstream outbrick;
//extern  G4int NtrackKaon, NtrackPi0, Ntrackgamma1, Ntrackgamma2;
//extern  G4double Eneutron_summ;
//extern  G4double SumPartEnKilled;
//extern  G4int nsob_part, partident;

extern G4double Det_level[52],ll1[52],ll2[52];
extern G4double Ekin_rand;
//extern G4double EnDepSi1[1000], EnDepSi2[1000], EnDepCsI[1000], EnDepSi3[1000];


extern G4int  Nev;
extern G4int n_muMsc, n_muPP, n_muBr,n_muSS; 

extern G4int nscat;

extern FILE *fp2;
extern FILE *fp5;
extern FILE *fp55;
extern FILE *fp6;



SteppingAction::SteppingAction()
{
}

SteppingAction::~SteppingAction()
{}


//int PP[10];
//for (i=0;i<10;i++) PP[i]=0; 
void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
    
    


//    G4cout<<"STEPPING ACTION!"<<G4endl;
    G4Track* aTrack = aStep->GetTrack();
    G4int MyTrackID= aTrack->GetTrackID();
    
    int i, nip=5;
    int MyParentID;
    int ktaucharm=0;
    int kkaon=0, mpid=-3 ;
    int saveMyTrackID=-1;
    int pant=-888;
    int tr= MyTrackID;
    int kkk, kpee=0, iem=0;
    double dli;


  
    fp2=fopen( "mu-trec-info-opt0-exp-0-1cm","a");
    fp5=fopen("mu-scattering-info-opt0-exp-0-1cm","a");
    fp55=fopen("mu-coul-scatt-delta","a");    
        fp6=fopen("mu-end-point-opt0-exp-0-1cm","a");

//    fp2=fopen( "mu-trec-info-50GeV-SS-exp-0-1cm","a");
//    fp5=fopen("mu-scattering-info-50GeV-SS-exp-0-1cm","a");
//        fp6=fopen("mu-end-point-50GeV-WENT-SS-0-1cm","a");
  
    
    G4double x,y,z;
    G4double xx,yy,zz;
    G4double Ekin_thr=10.0*MeV;
    G4double MyVertKinEnerg=0.0,pmx_vert=0.0,pmy_vert=0.0,pmz_vert=0.0;
    int PPCH, mstnumb, PPCH_nu=0;
    G4double mx,my,mz,Partenergy;
    G4double pmx_pre,pmy_pre,pmz_pre,pmx_post,pmy_post,pmz_post,ptot_pre,ptot_post;
    G4double hall_x,hall_y,l;


    G4String myprocess="An";
    G4String myprocess1="An_pre";
    G4String myprocess2="An_post";        
    G4String part_proc1, part_proc2;
    
    G4StepPoint* MyPreStepPoint=aStep->GetPreStepPoint();
    G4StepPoint* MyPostStepPoint=aStep->GetPostStepPoint();
    G4VPhysicalVolume* MyPreStepVolume= MyPreStepPoint->GetPhysicalVolume();
    G4VPhysicalVolume* MyPostStepVolume= MyPostStepPoint->GetPhysicalVolume();


    G4double MyKineticEnergy = MyPreStepPoint->GetKineticEnergy();
    G4double MyKineticEnergyPre = MyPreStepPoint->GetKineticEnergy();
    G4double MyKineticEnergyPost = MyPostStepPoint->GetKineticEnergy();    
    const G4DynamicParticle* MyDynamicParticle = aTrack->GetDynamicParticle();


    if (MyDynamicParticle!=NULL)
    {
//	G4cout<<"MyDynamicParticle!=NULL"<<G4endl;
	const G4ParticleDefinition* MyParticleDefinition = MyDynamicParticle->GetDefinition();
	if (MyParticleDefinition!=NULL)
	{  
//	    G4cout<<"MyParticleDefinition!=NULL"<<G4endl;
	    const G4String MyParticleName = MyParticleDefinition->GetParticleName();
//	    G4cout<<MyParticleName<<G4endl;


	    G4double theta;
//	    const G4VProcess* PreStepProcess=MyPreStepPoint->GetProcessDefinedStep();
//	    const G4VProcess* PostStepProcess=MyPostStepPoint->GetProcessDefinedStep();
//	    const G4VProcess* MyCreatorProcess = aTrack->GetCreatorProcess();

            if(aStep->GetTrack()->GetCreatorProcess()){ 
	    myprocess = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();}

//	    G4cout<<myprocess<<G4endl;

             const G4VProcess* process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
             if(process) myprocess1 =process ->GetProcessName();
             
              process = aStep->GetPreStepPoint()->GetProcessDefinedStep();
             if(process) myprocess2 =process ->GetProcessName();


//	    G4cout<<myprocess1<<G4endl;	    
//	    G4cout<<myprocess2<<G4endl;	    	    

//	    const G4String MyProcessName = aTrack->GetCreatorProcess()->GetProcessName();	    
//            if(aTrack->GetCreatorProcess()->GetProcessName() =="Decay"){ppp=1 ;}
	    
	    x= MyPreStepPoint->GetPosition()[0];
	    y= MyPreStepPoint->GetPosition()[1];
	    z= MyPreStepPoint->GetPosition()[2];
	    xx= MyPostStepPoint->GetPosition()[0];
	    yy= MyPostStepPoint->GetPosition()[1];
	    zz= MyPostStepPoint->GetPosition()[2];

	    pmx_pre= MyPreStepPoint->GetMomentum()[0];
	    pmy_pre= MyPreStepPoint->GetMomentum()[1];
	    pmz_pre= MyPreStepPoint->GetMomentum()[2];
	    
	    ptot_pre=sqrt(pmx_pre*pmx_pre+pmy_pre*pmy_pre+pmz_pre*pmz_pre);

	    pmx_post= MyPostStepPoint->GetMomentum()[0];
	    pmy_post= MyPostStepPoint->GetMomentum()[1];
	    pmz_post= MyPostStepPoint->GetMomentum()[2];

	    ptot_post=sqrt(pmx_post*pmx_post+pmy_post*pmy_post+pmz_post*pmz_post);	    
	    
	    hall_x= 500.0*cm;
	    hall_y= 500.0*cm;
//	    l=6.0*mm;
	    l=60.0*m;
	    mx= MyPreStepPoint->GetMomentumDirection()[0];
	    my= MyPreStepPoint->GetMomentumDirection()[1];
	    mz= MyPreStepPoint->GetMomentumDirection()[2];
            theta= atan(sqrt(mx*mx+my*my)/mz);
	    if ((xx>= -hall_x)&&(xx<= hall_x)&&(yy>= -hall_y)&&(yy<= hall_y)&&(zz>=-1.0*cm)&&(zz<= l)&&(MyPreStepVolume!= NULL)&&(MyPostStepVolume!= NULL))
	    {
	        G4double R;
		G4double RR;
	        RR= sqrt(x*x+y*y+z*z);
		R= sqrt(xx*xx+yy*yy+zz*zz);
//a 21.02.10        if (


		    
                    PPCH=-1;
                    PPCH_nu=0;
                    ktaucharm=0;
                    kkaon=0;
                    
                    if (MyParticleName=="tau-") {PPCH=0;pant= 15;}
		    if (MyParticleName=="pi-") {PPCH=1;pant= -211;}
		    if (MyParticleName=="pi+") {PPCH=2;pant= 211;}
		    if (MyParticleName=="e+")  {PPCH=3;pant= -11;}
		    if (MyParticleName=="e-") {PPCH=4;pant= 11;}
		    if (MyParticleName=="mu-") {PPCH=5;pant= 13;}
		    if (MyParticleName=="mu+") {PPCH=6;pant= -13;}
		    if (MyParticleName=="proton"){PPCH=7;pant= 2212;}
		    if (MyParticleName=="D0") {PPCH=8; pant= 421;}
		    if (MyParticleName=="kaon-") {PPCH=9;pant= -321;}
		    if (MyParticleName=="anti_kaon0") {PPCH=10;pant=-311;}		    		    		    
		    if (MyParticleName=="pi0") {PPCH=11;pant= 111;}
    		    if (MyParticleName=="nu_e") PPCH=12;
    		    if (MyParticleName=="nu_mu") PPCH=13;    		    	            
		    if (MyParticleName=="rho+") PPCH=14;    		    	            
    		    if (MyParticleName=="rho0") PPCH=15;    		    	                		    
    		    if (MyParticleName=="rho-") PPCH=16;    		    	                		    		    
    		    if (MyParticleName=="nu_tau") PPCH=17;    		    	                		    
    		    if (MyParticleName=="anti_nu_tau") PPCH=18;    		    	                		    
    		    if (MyParticleName=="anti_nu_mu") PPCH=19;    		    	            
    		    if (MyParticleName=="anti_nu_e") PPCH=20;    		    
    		    if (MyParticleName=="gamma") {PPCH=21; pant= 22;}
    		    if (MyParticleName=="D+"){PPCH=22; pant= 411;}
    		    if (MyParticleName=="Ds+") {PPCH=23; pant= 431;} 
    		    if (MyParticleName=="neutron") PPCH=24;    
    		    if (MyParticleName=="kaon+") {PPCH=25;pant= 321;}    
    		    if (MyParticleName=="kaon0") {PPCH=26;pant= 311;}
    		    if (MyParticleName=="kaon0L") {PPCH=27;pant= 130;}    		        
    		    if (MyParticleName=="kaon0S") {PPCH=28;pant= 310;}
    		    if (MyParticleName=="eta") PPCH=29;
    		    if (MyParticleName=="anti_k_star-") PPCH=30;
		    
		    		   

             MyParentID= aTrack->GetParentID();
            mpid=MyParentID;  
             MyTrackID= aTrack->GetTrackID();

               Partenergy= MyPreStepPoint->GetKineticEnergy();

		            nip=1;
//		            fprintf(fp2,"IP: x=%7.6e y=%7.6e z=%7.6e  %d  %d\n", xx,yy,zz,PPCH,tr);
//outbrick << std::setw(20)<< " 1 "<< MyParticleName <<G4endl;	    

          mstnumb= aTrack-> GetCurrentStepNumber();


G4double edepStep = aStep->GetTotalEnergyDeposit(); 
//fprintf(fp2,"   mstnumb=%d   zz=   %7.6e cm,    edepStep=  %7.6e \n",  mstnumb,  zz/cm,  edepStep);           

//processName="muMsc"
//processName = "muPairProd"
//processName = "muBrems"

 if(myprocess== "mumsc" ) {n_muMsc++; G4cout<<"n_muMsc(mumsc)="<<n_muMsc<<G4endl;}
 if(myprocess== "msc" ) {n_muMsc++; G4cout<<"n_muMsc(msc)="<<n_muMsc<<G4endl;}
 if(myprocess== "muMsc" ) {n_muMsc++; G4cout<<"n_muMsc(muMsc)="<<n_muMsc<<G4endl;}
 
 if(myprocess== "muss" ) {n_muSS++; G4cout<<"n_muSS="<<n_muSS<<G4endl;} 
 
 if(myprocess== "WentzelVIUni") {n_muMsc++; G4cout<<"n_muMsc="<<n_muMsc<<G4endl;}  
 
 
// if(myprocess== "CoulombScat") {n_muPP++; G4cout<<"n_muCoulombScat="<<n_muPP<<G4endl;
// for nuclear recoil Fe!!!!:
//G4cout<<"*Coulomb  PART NAME "<<  MyParticleName  <<G4endl;
//                if (PPCH==5){ 
// fprintf(fp55,"track: pre  %d  %d  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e %d %d %d\n",
// Nev, pant,x/cm,y/cm,z/cm,MyKineticEnergyPre/GeV,pmx_pre/GeV,pmy_pre/GeV,pmz_pre/GeV,PPCH,tr,mpid);
// fprintf(fp55,"track: post %d  %d  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e %d %d %d\n",
// Nev, pant,xx/cm,yy/cm,zz/cm,MyKineticEnergyPost/GeV,pmx_post/GeV,pmy_post/GeV,pmz_post/GeV,PPCH,tr,mpid);
//                               }

 if(PPCH==5 && myprocess1== "CoulombScat") {G4cout<<"* myprocess1  Coulomb  PART NAME "<<  MyParticleName  <<G4endl;}
 if(PPCH==5 && myprocess2== "CoulombScat") {G4cout<<"* myprocess2  Coulomb  PART NAME "<<  MyParticleName  <<G4endl;}
   
int ilev=zz/cm;

 if(myprocess2== "CoulombScat") {
// if(myprocess2== "muBrems") { 
 
//n_muPP++; G4cout<<"n_muCoulombScat="<<n_muPP<<G4endl;
 
//G4cout<<"**1** Coulomb  PARTICLE NAME "<<  MyParticleName  <<G4endl;
//if (PPCH==5)
// fprintf(fp5,"%d  %d  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e %d %d %d\n",
// Nev, pant,x/cm,y/cm,z/cm,MyKineticEnergyPre/GeV,pmx_pre/GeV,pmy_pre/GeV,pmz_pre/GeV,PPCH,tr,mpid);

  nscat++;

if (PPCH==5){
 fprintf(fp5,"%d  %d  %d  %d 888 %7.6e %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e %d %d %d\n",
 Nev, pant, ilev, nscat ,Ekin_rand, xx/cm,yy/cm,zz/cm,MyKineticEnergyPost/GeV,pmx_post/GeV,pmy_post/GeV,pmz_post/GeV,PPCH,tr,mpid);
n_muPP++; G4cout<<"n_muCoulombScat="<<n_muPP<<G4endl;

// fprintf(fp55,"pre^ %d  %d  %d  %d 888 %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e %d %d %d\n",
// Nev, pant, ilev, nscat ,x/cm,y/cm,z/cm,MyKineticEnergyPre/GeV,pmx_pre/GeV,pmy_pre/GeV,pmz_pre/GeV,PPCH,tr,mpid);
// fprintf(fp55,"post^ %d  %d  %d  %d 888 %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e %d %d %d\n",
// Nev, pant, ilev, nscat ,xx/cm,yy/cm,zz/cm,MyKineticEnergyPost/GeV,pmx_post/GeV,pmy_post/GeV,pmz_post/GeV,PPCH,tr,mpid);

double p_pre=sqrt(pmx_pre*pmx_pre+pmy_pre*pmy_pre+pmz_pre*pmz_pre);
double p_post=sqrt(pmx_post*pmx_post+pmy_post*pmy_post+pmz_post*pmz_post);
double delta_pz=pmz_pre-pmz_post;

// fprintf(fp55,"    \n");  
// fprintf(fp55,"p_pre=%7.6e   p_post=%7.6e  pmz_pre=%7.6e  pmz_post=%7.6e \n",
//        p_pre/GeV, p_post/GeV, pmz_pre/GeV, pmz_post/GeV);


double teta_z_pre=acos(pmz_pre/p_pre); 
double teta_z_post=acos(pmz_post/p_post);
double delta_teta=teta_z_post-teta_z_pre;


// fprintf(fp55," ggggggggggg n_muCoulombScat=  \n");  
 fprintf(fp55,"%d  %d  %d  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e \n",
     Nev, pant, nscat,Ekin_rand, zz/cm, MyKineticEnergyPre/GeV,MyKineticEnergyPost/GeV,delta_pz/GeV,teta_z_pre, teta_z_post, delta_teta);
            }}

   
if (PPCH==5 && MyKineticEnergyPost==0.0 &&  MyKineticEnergyPre==0.0)
 fprintf(fp6,"%d  %d  %d  %d 777  %7.6e %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e %d %d %d\n",
 Nev,pant,ilev , nscat, Ekin_rand, xx/cm,yy/cm,zz/cm,MyKineticEnergyPost/GeV,pmx_post/GeV,pmy_post/GeV,pmz_post/GeV,PPCH,tr,mpid);


if ( PPCH==5   &&  zz == ll1[49] && nscat == 0 )
 fprintf(fp5,"%d  %d  %d  0  777 %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e %7.6e  %7.6e %d %d  %d\n",
 Nev,pant,ilev,Ekin_rand,  xx/cm,yy/cm,zz/cm,MyKineticEnergyPost/GeV,pmx_post/GeV,pmy_post/GeV,pmz_post/GeV,PPCH,tr,mpid);

if ( PPCH==5   &&  zz == ll1[49] && nscat > 0 )
 fprintf(fp5,"%d  %d  %d  %d  333 %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e %7.6e  %7.6e %d %d  %d\n",
 Nev,pant,ilev,nscat,Ekin_rand, xx/cm,yy/cm,zz/cm,MyKineticEnergyPost/GeV,pmx_post/GeV,pmy_post/GeV,pmz_post/GeV,PPCH,tr,mpid);


   
 if(myprocess== "muBrems" ) {n_muBr++; G4cout<<"n_muBr="<<n_muBr<<G4endl; 
// G4cout<<"*Brems  PRTICLE NAME "<<  MyParticleName  <<G4endl;
                             }
  

if (PPCH==5){
        for (i=0; i<50;i++) {
//	for (i=0; i<10;i++) {
	
// if  ((MyPreStepPoint->GetStepStatus() == fGeomBoundary || z>=0.0)
//                                && zz > lll[i] && zz<= lll[i+1])
 
// if  ( zz > ll1[i] && zz<= ll2[i])  Det_level[i]+=edepStep;

 if  ( zz > ll1[i] && zz<= ll1[i+1])  Det_level[i+1]+=edepStep;

if (zz == ll1[i+1])
if (abs(zz - ll1[i+1]) < 0.00001)
fprintf(fp2,"%d  00%d  %d  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e %7.6e  %7.6e %d %d  %d\n",
  i, Nev, pant,Ekin_rand,xx/cm,yy/cm,zz/cm,MyKineticEnergyPost/GeV,pmx_post/GeV,pmy_post/GeV,pmz_post/GeV,PPCH,tr,mpid);

                            }
           }
           
/*			    
if (PPCH==5 && (zz == ll1[1]))
fprintf(fp2,"%d  001 %d  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e %7.6e  %7.6e %d %d  %d\n",
  Nev, pant,xx/cm,yy/cm,zz/cm,MyKineticEnergyPost/GeV,pmx_post/GeV,pmy_post/GeV,pmz_post/GeV,PPCH,tr,mpid);

if (PPCH==5 && (zz == ll1[5]))
fprintf(fp2,"%d  002 %d  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e %7.6e  %7.6e %d %d  %d\n",
 Nev, pant,xx/cm,yy/cm,zz/cm,MyKineticEnergyPost/GeV,pmx_post/GeV,pmy_post/GeV,pmz_post/GeV,PPCH,tr,mpid);

 			    
if (PPCH==5 && (zz == ll1[10]))
fprintf(fp2,"%d  003 %d  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e %7.6e  %7.6e %d %d  %d\n",
 Nev, pant,xx/cm,yy/cm,zz/cm,MyKineticEnergyPost/GeV,pmx_post/GeV,pmy_post/GeV,pmz_post/GeV,PPCH,tr,mpid);

if (PPCH==5 && (zz == ll1[15]))
fprintf(fp2,"%d  004 %d  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e %7.6e  %7.6e %d %d  %d\n",
 Nev, pant,xx/cm,yy/cm,zz/cm,MyKineticEnergyPost/GeV,pmx_post/GeV,pmy_post/GeV,pmz_post/GeV,PPCH,tr,mpid);
			    

if (PPCH==5 && (zz == ll1[20]))
fprintf(fp2,"%d  005 %d  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e %7.6e  %7.6e %d %d  %d\n",
 Nev, pant,xx/cm,yy/cm,zz/cm,MyKineticEnergyPost/GeV,pmx_post/GeV,pmy_post/GeV,pmz_post/GeV,PPCH,tr,mpid);

if (PPCH==5 && (zz == ll1[25]))
fprintf(fp2,"%d  006 %d  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e %7.6e  %7.6e %d %d  %d\n",
 Nev, pant,xx/cm,yy/cm,zz/cm,MyKineticEnergyPost/GeV,pmx_post/GeV,pmy_post/GeV,pmz_post/GeV,PPCH,tr,mpid);


if (PPCH==5 && (zz == ll1[30]))
fprintf(fp2,"%d  007 %d  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e %7.6e  %7.6e %d %d  %d\n",
 Nev, pant,xx/cm,yy/cm,zz/cm,MyKineticEnergyPost/GeV,pmx_post/GeV,pmy_post/GeV,pmz_post/GeV,PPCH,tr,mpid);

if (PPCH==5 && (zz == ll1[35]))
fprintf(fp2,"%d  008 %d  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e %7.6e  %7.6e %d %d  %d\n",
 Nev, pant,xx/cm,yy/cm,zz/cm,MyKineticEnergyPost/GeV,pmx_post/GeV,pmy_post/GeV,pmz_post/GeV,PPCH,tr,mpid);


if (PPCH==5 && (zz == ll1[40]))
fprintf(fp2,"%d  009 %d  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e %7.6e  %7.6e %d %d  %d\n",
 Nev, pant,xx/cm,yy/cm,zz/cm,MyKineticEnergyPost/GeV,pmx_post/GeV,pmy_post/GeV,pmz_post/GeV,PPCH,tr,mpid);

if (PPCH==5 && (zz == ll1[45]))
fprintf(fp2,"%d  0010 %d  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e %7.6e  %7.6e %d %d  %d\n",
 Nev, pant,xx/cm,yy/cm,zz/cm,MyKineticEnergyPost/GeV,pmx_post/GeV,pmy_post/GeV,pmz_post/GeV,PPCH,tr,mpid);

if (PPCH==5 && MyKineticEnergyPost==0.0 && pmz_pre>0.0)fprintf(fp2,"%d  888 %d  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e %7.6e  %7.6e %d %d  %d\n",
 Nev, pant,x/cm,y/cm,z/cm,MyKineticEnergyPre/GeV,pmx_pre/GeV,pmy_pre/GeV,pmz_pre/GeV,PPCH,tr,mpid);
*/


//outbrick << std::setw(15)<< xx << std::setw(15) <<  yy << std::setw(15) << zz << std::setw(5)<< PPCH << std::setw(5)<< tr<<G4endl;
                        if(zz>0.0)       {
			     nip=0;                          
          mstnumb= aTrack-> GetCurrentStepNumber();
         
     if(mstnumb==1){
             MyVertKinEnerg= aTrack->GetVertexKineticEnergy();
           pmx_vert= aTrack->GetVertexMomentumDirection()[0];
           pmy_vert= aTrack->GetVertexMomentumDirection()[1];
           pmz_vert= aTrack->GetVertexMomentumDirection()[2];



//   fprintf(fp2,"    8  8  8.8  8.8  8.8  8.8  8.8  8.8  8.8   8  8  8 \n");
//                 fprintf(fp2,"   %d  %7.6e  %7.6e  %7.6e  %7.6e  %7.6e %7.6e  %7.6e %d %d  %d\n",
//                  pant,x,y,z,MyKineticEnergyPre/GeV,pmx_pre/GeV,pmy_pre/GeV,pmz_pre/GeV,PPCH,tr,mpid);
//    fprintf(fp2," %7.6e %7.6e %7.6e %d  %d\n", x,y,z,PPCH,tr);           
                   }

 	    		                  }       
                                    
//                fclose(fp2);            
                 
//G4cout<<"E="<<MyKineticEnergy<<" "<<theta<<" "<<R<<" z= " <<z<<"time "<<MyPreStepPoint->GetGlobalTime()<<"Track= "<<MyTrackID<<aTrack->GetVolume()->GetCopyNo()<<G4endl;

	    }
	    
 	}
 	
	
    }
            fclose(fp5);
        fclose(fp55);            
            fclose(fp2);
            fclose(fp6);            
                           
}
