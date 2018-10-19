/***********************************************************/
/* (3-D)  trianglur-surface model  Brownian dynamics

   repulsive potential 80exp(1/(x-0.85))/(x-0.67) 
   bond potential 80exp(1/(1.15-x))/(1.33-x ) 
   do not bond-flip if new bond share 3 faces or another bond
   Time-stamp: <2018-03-27 22:25:14 noguchi>  
	 last modified: Maike Jung, JGU Mainz: 25.04.2018
*/
/***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define  DT       0.0001         	 /* time interval of a unit step */
#define  NUMBER   2562           	 /* number of vertices  */
#define  KBT      0.00001            	 /* temprature k_BT */
#define  BETA    (-1.0/KBT)
#define  BMASSB   1.0            	 /* mass of the vertices */

#define  MAXBOND  20             	 /* maximum of bonds */ 
#define  NF       (2*(NUMBER-2 ))    /* number of faces */
#define  NE       (3*(NUMBER-2 ))    /* number of vertices */

/* parameters for U_bond and U_rep */
/* IMPORTANT: unit length of bond (a) is ONE */
#define  BL_MAX   1.33   /* maximum of bond length */
#define  BL_MAXIN  (BL_MAX-0.001)
#define  BL_MAXIN2  (BL_MAXIN*BL_MAXIN)
#define  BL_MIN   0.67     /* mimimum of bond length */
#define  BL_MINOU  (BL_MIN+0.001)
#define  BL_MINOU2  (BL_MINOU*BL_MINOU)
#define  RCUTO  1.15
#define  RCUTI  0.85
#define  RCUTIN2 0.72012   /* numerical cut-off length^2 */
#define  RCUTOU2 1.32572   /* outside */
#define  LN80    4.38202663467388 

/* parameters for area and volume potentials */
#define  KAPPA    20.
#define  ARMIN    0.00001
#define  KARE     2.  /* area potential */
#define  ARE0     (0.41*NF )
#define  KVOL     1.  /* volume potential */


#define  MAXBK   20      /* maximum number of Verlet booking */ 
#define  RBOOK2  1.69

#define  ITIME_D0  20

#define  NHV     200  
#define  VMAXH   2.0  
#define  PI     3.1415926535897932
#define  PIX2   6.2831853071795864
#define  ARD00   7.0898154036220639488

void read_init();
void write_cood();
void energy();
void record_draw();
void record_drawVMD();
void cal_force();
void calaf();
void calvof();
void calrn();
void bonflip();
void cal_vbook();
void calrn_new();
void calrn_old();
void init_para();
void cal_rg();
void check_V();
void check_bond();
void mean();
void distmov();
void sort_bd();
void wnoise();
double gasdev();
double ran2();
void set_gg();
double calc_cm();
void find_points();
void ext_const_force();
void ext_const_force_single();
double calc_max_dist();

int main(){
	double  bx[NUMBER],  by[NUMBER],  bz[NUMBER];
	double  bxv[NUMBER], byv[NUMBER], bzv[NUMBER];
	double  bxp[NUMBER], byp[NUMBER], bzp[NUMBER];
	double  bxq[NUMBER], byq[NUMBER], bzq[NUMBER];
	double fx[NUMBER], fy[NUMBER], fz[NUMBER];
	double wn_f[(NUMBER)*3+1],gg[3];
	int iface[NF][3],iedge[NE][2];
	int ibon_all[NUMBER][MAXBOND],nbon_all[NUMBER+1];
	int iedb[NE][2],iedf[NE][2],ifaed[NF][3];
	int ibook[NUMBER*MAXBK],ibooki[NUMBER+1];
	int index_pull[2];
	double pull[2][3];
	double en[9],rg[5],ar[5],vol0,vol00[2],vol_int,vol_fin,pf,fdist;
	double en_m[9],en_m2[9],cv[9],rgm[12],arm[10];
	double ari[NUMBER],rni[NUMBER][3],eni[NUMBER],gni[NUMBER][3];
	double extforce[NUMBER][3];
	double maxdist;
	double xcm, ycm, zcm;
	  
	int istep,maxstep,itime_m,itime_l,imstep,imsteph, runnumber;
	int i,iii;
	int icnt[3],ibdmami[3];
	long   iseed=-11;
	double gm,dif[2],dtm, gmpulse,gmminus,gmrate,dtymgmp;
	double bqmax[2],rbonl[4];
	FILE *finp,*dp,*arp,*enp,*rp,*mp,*forcedev,*dpVMD;

	/* read simulation parameters from file */
	finp=fopen("para","r"); 
	fscanf(finp,"%lf %lf",&vol_int,&vol_fin);
	fscanf(finp,"%d",&maxstep);
	fscanf(finp,"%d",&itime_m);
	fscanf(finp,"%d",&itime_l);
	fscanf(finp,"%d",&imstep);
	fscanf(finp,"%ld",&iseed);
	fscanf(finp,"%lf",&pf);
	fscanf(finp,"%d",&runnumber);
	fclose(finp);
	imsteph=imstep/2;

	printf("START SIMULATION \n f = %.2f; run %d \n",pf , runnumber);	
	
	for(i=0;i<2;i++)
		bqmax[i]=0.;
	ibdmami[0]=0;
	ibdmami[1]=MAXBOND;
	ibdmami[2]=0;
	for(i=0;i<3;i++)
		icnt[i]=0;

	/* set inital parameters for area and volume potential */
	ar[2]= sqrt(ARE0*0.25/PI);
	vol00[1]= (4.*PI*ar[2]*ar[2]*ar[2])/3.;
	vol00[0]= vol_int;
	vol0=vol00[0]*vol00[1];

	/* parameters for integration (leapfrog with Langevin thermostat) */
	//gm=1.;
	gm = 1.0;
	dtm=DT/BMASSB;
	dif[0]=sqrt(2.*gm*KBT/DT);
	gmpulse= 1. + gm*DT/(2.*BMASSB);
	gmminus= 1. - gm*DT/(2.*BMASSB);
	gmrate =gmminus/gmpulse;
	dtymgmp=dtm/gmpulse;

	/* output files for later */
	char buf_m[0x100], buf_draw[0x100], buf_ar[0x100], buf_en[0x100], buf_rad[0x100], buf_forcedev[0x100], buf_drawVMD[0x100];
	snprintf(buf_m, sizeof(buf_m), "Data/m_pf%.0f_%d.dat", pf,runnumber);
	snprintf(buf_draw, sizeof(buf_draw), "Data/draw_pf%.0f_%d.dat", pf,runnumber);
	snprintf(buf_ar, sizeof(buf_ar), "Data/ar_pf%.0f_%d.dat", pf,runnumber);
	snprintf(buf_en, sizeof(buf_en), "Data/en_pf%.0f_%d.dat", pf,runnumber);
	snprintf(buf_rad, sizeof(buf_rad), "Data/rad_pf%.0f_%d.dat", pf,runnumber);
	snprintf(buf_forcedev, sizeof(buf_forcedev), "Data/time_evolution_pf%.0f_%d.dat", pf,runnumber);
	snprintf(buf_drawVMD, sizeof(buf_drawVMD), "Data/drawVMD_pf%.0f_%d.xyz", pf,runnumber);

	mp=fopen(buf_m,"w");
	dp=fopen(buf_draw,"w"); 
	dpVMD=fopen(buf_drawVMD,"w"); 
	//fprintf(dp,"%lf %d\n", 0.0,0); // no meaning at all, required by the drawing program
	arp=fopen(buf_ar,"w"); 
	enp=fopen(buf_en,"w"); 
	rp =fopen(buf_rad,"w"); 
	forcedev = fopen(buf_forcedev,"w");

	/* read initial configuration */
	read_init(bx,by,bz,bxv,byv,bzv,bxq,byq,bzq,iface,iedge,ibon_all,nbon_all,iedb,iedf,ifaed,mp,index_pull);
	init_para(en_m,en_m2,rgm,arm);

	/* select two opposing beads, which will be used for pulling later */	
	find_points(index_pull, bx, by, bz, pull,pf);


	for(i=0;i<NUMBER;i++){
		bxp[i]  =bx[i];
		byp[i]  =by[i];
		bzp[i]  =bz[i];
	}

	rbonl[2]=-1.;
	rbonl[3]= 1.;

	cal_vbook(bx,by,bz,ibook,ibooki,ibdmami);
	cal_force(bx,by,bz,fx,fy,fz,ari,rni,eni,en,ar,vol0,iface,iedge,ibon_all,nbon_all,ibook,ibooki,icnt,index_pull, pull, 0);
	energy(bx,by,bz,bxv,byv,bzv,rni,gni,iface,iedge,en,ar,rbonl);
	cal_rg(bx,by,bz,rg,gg,ar);
  
	istep=0;

	record_draw(bx,by,bz,iface,0,dp);
	record_drawVMD(bx,by,bz,0,dpVMD);
	xcm, ycm, zcm = calc_cm(bx, by,bz);
	fprintf(forcedev,"%12lg %12lg %12lg %12lg %12lg %12lg\n", istep*DT, sqrt((bx[index_pull[0]] - bx[index_pull[1]])*(bx[index_pull[0]] - bx[index_pull[1]]) + (by[index_pull[0]] - by[index_pull[1]])*(by[index_pull[0]] - by[index_pull[1]]) + (bz[index_pull[0]] - bz[index_pull[1]])*(bz[index_pull[0]] - bz[index_pull[1]])),sqrt((bx[index_pull[0]] - bx[index_pull[1]])*(bx[index_pull[0]] - bx[index_pull[1]]) + (by[index_pull[0]] - by[index_pull[1]])*(by[index_pull[0]] - by[index_pull[1]]) + (bz[index_pull[0]] - bz[index_pull[1]])*(bz[index_pull[0]] - bz[index_pull[1]])),xcm, ycm, zcm);
	fprintf(enp,"%12lg %12lg %12lg %12lg %12lg %12lg %12lg %12lg \n", istep*DT, en[0],en[1]/NUMBER,en[2]/NUMBER,en[3],en[4],en[5],en[6]/NUMBER);
	fprintf(rp,"%12lg %12lg %12lg %12lg %12lg %12lg %12lg %12lg %12lg\n", istep*DT, rg[0],rg[2],rg[3],rg[4],rbonl[0],gg[0],gg[1],gg[2]);
	fprintf(arp,"%12lg %12lg %12lg %12lg %12lg %12lg\n", istep*DT, ar[0],ar[1],ar[4],ar[3]/ARD00,vol00[0]);



	/*----------------main loop----------------*/
	for(istep=1;istep<=maxstep;istep++){
		
		if(istep%10000==0){
			printf("step: %d; time: %f \n", istep, istep*DT);
		}

		if( (icnt[1])>0||(icnt[2]>0)){
			fprintf(mp,"error istep %d %d %d\n",istep,icnt[1],icnt[2]);
			fclose(mp);
			mp=fopen(buf_m,"a");
		}

		if(istep<imsteph){
			vol00[0]= vol_int+ ((vol_fin-vol_int)*istep)/imsteph;
		}else{
			vol00[0]= vol_fin;
		}
		vol0=vol00[0]*vol00[1];


		/* performe MD moves - with external pulling force */
		for(i=0;i<NUMBER;i++){
			bx[i] = bxp[i] + DT*bxv[i];
			by[i] = byp[i] + DT*byv[i];
			bz[i] = bzp[i] + DT*bzv[i];
		}

		/* calculate forces and velocities */
		cal_force(bx,by,bz,fx,fy,fz,ari,rni,eni,en,ar,vol0,iface,iedge,ibon_all,nbon_all,ibook,ibooki,icnt,index_pull, pull, istep);
		wnoise(&iseed,wn_f,dif);
		//ext_const_force(extforce, index_pull, pull, istep);
		ext_const_force_single(extforce, index_pull, pull, istep);
		for(i=0;i<NUMBER;i++){
			bxv[i]  = gmrate*bxv[i] + dtymgmp*(fx[i]+wn_f[3*i]+extforce[i][0]);
			byv[i]  = gmrate*byv[i] + dtymgmp*(fy[i]+wn_f[3*i+1]+extforce[i][1]);
			bzv[i]  = gmrate*bzv[i] + dtymgmp*(fz[i]+wn_f[3*i+2]+extforce[i][2]);
		}

		

		/* performe bond flips every ITIME_D0 steps */
		if( (istep/ITIME_D0)*ITIME_D0==istep){
			bonflip(bx,by,bz,ari,rni,eni,en,ar,vol0,iface,iedge,ibon_all,nbon_all,iedb,iedf,ifaed,icnt,ibdmami,&iseed);

			cal_vbook(bx,by,bz,ibook,ibooki,ibdmami);
			distmov(bx,by,bz,bxq,byq,bzq,bqmax);
			energy(bx,by,bz,bxv,byv,bzv,rni,gni,iface,iedge,en,ar,rbonl);
			cal_rg(bx,by,bz,rg,gg,ar);
		  
			if( istep>imstep ){
				mean(en,en_m,en_m2,rg,rgm,ar,arm,rbonl);
			}
			if( (istep/itime_m)*itime_m==istep){
				fprintf(enp,"%12lg %12lg %12lg %12lg %12lg %12lg %12lg %12lg \n", istep*DT, en[0],en[1]/NUMBER,en[2]/NUMBER,en[3],en[4],en[5],en[6]/NUMBER);
				fprintf(rp,"%12lg %12lg %12lg %12lg %12lg %12lg %12lg %12lg %12lg\n", istep*DT, rg[0],rg[2],rg[3],rg[4],rbonl[0],gg[0],gg[1],gg[2]);
				fprintf(arp,"%12lg %12lg %12lg %12lg %12lg %12lg\n", istep*DT, ar[0],ar[1],ar[4],ar[3]/ARD00,vol00[0]);
				if( (istep/itime_l)*itime_l==istep){
					record_draw(bx,by,bz,iface,istep,dp);
					record_drawVMD(bx,by,bz,istep,dpVMD);
					maxdist = calc_max_dist(bx, by, bz);
					xcm, ycm, zcm = calc_cm(bx, by,bz);
					fprintf(forcedev,"%12lg %12lg %12lg %12lg %12lg %12lg \n", istep*DT, sqrt((bx[index_pull[0]] - bx[index_pull[1]])*(bx[index_pull[0]] - bx[index_pull[1]]) + (by[index_pull[0]] - by[index_pull[1]])*(by[index_pull[0]] - by[index_pull[1]]) + (bz[index_pull[0]] - bz[index_pull[1]])*(bz[index_pull[0]] - bz[index_pull[1]])),maxdist,xcm, ycm, zcm);
					check_V(bx,by,bz,gg,ar,iface,iedge,mp);
					if( ( istep/(itime_l*10) )*(itime_l*10)==istep){
						// resort faces and edges for faster computation 
						sort_bd(bx, by, bz,bxv, byv, bzv,bxq,byq,bzq,iface,iedge,ibon_all,nbon_all,iedb,iedf,ifaed,mp,index_pull);
						cal_vbook(bx,by,bz,ibook,ibooki,ibdmami);
					}
				}
			}
		} 

		for(i=0;i<NUMBER;i++){
			bxp[i]  =bx[i];
			byp[i]  =by[i];
			bzp[i]  =bz[i];
		}

	}
	/*------------------------------------------*/ 
	iii = (maxstep-imstep)/ITIME_D0;

	write_cood(bx,by,bz,bxv,byv,bzv,iface,ibon_all,nbon_all,pf,runnumber);

	for(i=0;i<9;i++){
		en_m[i]=en_m[i]/(iii);
		en_m2[i]=en_m2[i]/(iii);
		en_m2[i]=  (en_m2[i] - en_m[i]*en_m[i])*(iii)/(iii-1.);
		cv[i]= en_m2[i]/(KBT*KBT);
		en_m2[i]= sqrt(en_m2[i]);
	}
	en_m[1]=en_m[1]/NUMBER;
	en_m2[1]=en_m2[1]/NUMBER;
	en_m[2]=en_m[2]/NUMBER;
	en_m2[2]=en_m2[2]/NUMBER;
	en_m[6]=en_m[6]/NUMBER;
	en_m2[6]=en_m2[6]/NUMBER;

	for(i=0;i<6;i++){
		rgm[2*i] = rgm[2*i]/iii;
		rgm[2*i+1]=rgm[2*i+1]/iii;
		rgm[2*i+1]=sqrt( (rgm[2*i+1] - rgm[2*i]*rgm[2*i])*(iii)/(iii-1.));
	}
	for(i=0;i<5;i++){
		arm[2*i] = arm[2*i]/iii;
		arm[2*i+1]=arm[2*i+1]/iii;
		arm[2*i+1]=sqrt( (arm[2*i+1] - arm[2*i]*arm[2*i])*(iii)/(iii-1.));
	}
	fprintf(mp,"en %lg  %lg  %lg  %lg  %lg  %lg  %lg\n",en_m[0],en_m[1],en_m[2],en_m[3],en_m[4],en_m[5],en_m[7]);
	fprintf(mp,"are %lg  %lg  %lg  %lg  %lg  %lg\n",arm[8],arm[9],arm[0],arm[1],arm[2],arm[3]);
	fprintf(mp,"ard %lg  %lg  %lg  %lg  %lg  %lg\n",arm[6]/ARD00,arm[7]/ARD00,arm[6],arm[7],arm[4],arm[5]);
	fprintf(mp,"rg %lg  %lg  %lg  %lg  %lg  %lg  %lg  %lg  %lg  %lg\n",rgm[0],rgm[1],rgm[2],rgm[3],rgm[4],rgm[5],rgm[6],rgm[7],rgm[8],rgm[9]);
	fprintf(mp,"en2 %lg  %lg  %lg  %lg  %lg  %lg  %lg\n",en_m2[0],en_m2[1],en_m2[2],en_m2[3],en_m2[4],en_m2[5],en_m2[7]);
	fprintf(mp,"cv %lg  %lg %lg  %lg, %lg %lg\n",cv[0],cv[1],cv[2],cv[3],cv[4],cv[5]);
	fprintf(mp,"kin en %lg  %lg q %lg  %lg\n",en_m[6],en_m2[6],en_m[8],en_m2[8]);
	fprintf(mp,"rbonl m %lg  %lg,max %lg min %lg\n",rgm[10],rgm[11],rbonl[2],rbonl[3]);
	fprintf(mp,"bond num, max %d min %d vbook %d\n",ibdmami[0],ibdmami[1],ibdmami[2]);
	fprintf(mp,"bond flip %d  %lg\n",icnt[0],1.0*icnt[0]/(1.*NE*ITIME_D0*maxstep));
	fprintf(mp,"number %d kappa %lg karea, %lg  ARE0 %lg\n",NUMBER,KAPPA,KARE,ARE0);
	fprintf(mp,"kvol, %lg  \n",KVOL);
	fprintf(mp,"vol %lg %lg  v_ar %lg %lg\n",vol_int,vol_fin,vol00[1],vol0);
	bqmax[0]=sqrt(bqmax[0]);
	fprintf(mp,"bqmax %lg\n",bqmax[0]);
	if(icnt[1]>0)
		fprintf(mp,"error too long bond %d\n",icnt[1]);
	if(icnt[2]>0)
		fprintf(mp,"error negative ari %d\n",icnt[2]);
	if(ibdmami[0]>MAXBOND)
		fprintf(mp,"error too many bonds %d %d\n",ibdmami[0],MAXBOND);
	if(ibdmami[2]>MAXBK*NUMBER)
		fprintf(mp,"error too many book %lg %d\n",ibdmami[2]*(1./NUMBER),MAXBK);
	if(sqrt(RBOOK2)<(2.*bqmax[0]+RCUTI) )
		fprintf(mp,"error too small rbook %lg %lg\n",2.*bqmax[0]+RCUTI,sqrt(RBOOK2) );

	fdist = sqrt((bx[index_pull[0]] - bx[index_pull[1]])*(bx[index_pull[0]] - bx[index_pull[1]]) + (by[index_pull[0]] - by[index_pull[1]])*(by[index_pull[0]] - by[index_pull[1]]) + (bz[index_pull[0]] - bz[index_pull[1]])*(bz[index_pull[0]] - bz[index_pull[1]]));
	printf("FINISHED: d = %lg ; pf = %lg \n", fdist,pf );
  
	fclose(mp);
	fclose(dp);
	fclose(dpVMD);
	fclose(arp);
	fclose(enp);
	fclose(rp);
	fclose(forcedev);

	return 0;
}


/*++++++++++++++++++++++++++++++++++++++++++++++*/
void cal_force(bx,by,bz,fx,fy,fz,ari,rni,eni,en,ar,vol0,iface,iedge,ibon_all,nbon_all,ibook,ibooki,icnt,index_pull, pull, istep)
     double  bx[],  by[],  bz[];
     double  fx[],  fy[],  fz[];
     double ari[],rni[NUMBER][3],eni[],en[],ar[],vol0;
     int iface[NF][3],iedge[NE][2];
     int ibon_all[NUMBER][MAXBOND],nbon_all[];
     int ibook[],ibooki[], icnt[]; 
	 int index_pull[2];
	 int istep;
	 double pull[2][3];
{
  int i,j,k,ib,jb,kb,j1,ial[2],jal[2],kal[2];
  double f0,r0,ftmp,fe,swct;
  double fcv0[NUMBER][MAXBOND+1][3];
  double fcv1[NUMBER][MAXBOND+1][3];
  double fcv2[NUMBER][MAXBOND+1][3];
  double fcvv[NUMBER][MAXBOND+1][3];
  double rsij[4],rsik[4],rsjk[4],rsij2[4],rsik2[4],rsjk2[4];
  double ar0,yar,gn[5];
  double fxa[NUMBER],fya[NUMBER],fza[NUMBER];
  double fxvo[NUMBER],fyvo[NUMBER],fzvo[NUMBER];
    
  for(i=0;i<NUMBER;i++){
    fx[i]=0.0;
    fy[i]=0.0;
    fz[i]=0.0;
    fxa[i]=0.0;
    fya[i]=0.0;
    fza[i]=0.0;
    fxvo[i]=0.0;
    fyvo[i]=0.0;
    fzvo[i]=0.0;
  }
  for(i=0;i<NUMBER;i++){
    eni[i]=0.;
    ari[i]=0.;
    for(j=0;j<3;j++){
      rni[i][j]=0.;
    }
    for(j=0;j<MAXBOND+1;j++){
      for(k=0;k<3;k++){
	fcv0[i][j][k]=0.;
	fcv1[i][j][k]=0.;
	fcv2[i][j][k]=0.;
	fcvv[i][j][k]=0.;
      }
    }
  }
  for(i=1;i<4;i++)
    en[i]=0.;
  ar[0]=0.;
  ar[1]=0.;
  ar0=0.;

	/* calculate U_bond */
  for(i=0;i<NE;i++){
    ib=iedge[i][0];
    jb=iedge[i][1];
    rsij[1]= bx[ib]-bx[jb];
    rsij[2]= by[ib]-by[jb];
    rsij[3]= bz[ib]-bz[jb];
    rsij[0]= rsij[1]*rsij[1]+rsij[2]*rsij[2]+rsij[3]*rsij[3];
    r0=sqrt(rsij[0]);
    if(rsij[0]<RCUTOU2 ){
      ftmp=0.;
      swct=0.;
    }else if(rsij[0]<BL_MAXIN2 ){
      ftmp= 1./(BL_MAX-r0);
      swct=1./(RCUTO-r0);
    }else{
      ftmp= 1./(BL_MAX-BL_MAXIN);
      swct=1./(RCUTO-BL_MAXIN);
      icnt[1]++;
    }
    fe=exp(swct+LN80)*ftmp;
    f0= (ftmp+swct*swct)*(-fe)/r0;
    en[1]=en[1]+fe;
    fx[ib]=fx[ib] +f0*rsij[1];
    fx[jb]=fx[jb]-f0*rsij[1];
    fy[ib]=fy[ib] +f0*rsij[2];
    fy[jb]=fy[jb]-f0*rsij[2];
    fz[ib]=fz[ib] +f0*rsij[3];
    fz[jb]=fz[jb]-f0*rsij[3];
  }

  for(i=0;i<NF;i++){
    ib=iface[i][0];
    jb=iface[i][1];
    kb=iface[i][2];

    rsij[1]=bx[ib]-bx[jb];
    rsij[2]=by[ib]-by[jb];
    rsij[3]=bz[ib]-bz[jb];
    rsij[0]= rsij[1]*rsij[1]+rsij[2]*rsij[2]+rsij[3]*rsij[3];

    for(k=0;k<nbon_all[ib];k++){
      if(jb==ibon_all[ib][k])
	ial[0]=k+1;
      if(kb==ibon_all[ib][k])
	ial[1]=k+1;
    }
    for(k=0;k<nbon_all[jb];k++){
      if(ib==ibon_all[jb][k])
	jal[0]=k+1;
      if(kb==ibon_all[jb][k])
	jal[1]=k+1;
    }
    for(k=0;k<nbon_all[kb];k++){
      if(ib==ibon_all[kb][k])
	kal[0]=k+1;
      if(jb==ibon_all[kb][k])
	kal[1]=k+1;
    }
    rsik[1]=bx[ib]-bx[kb];
    rsik[2]=by[ib]-by[kb];
    rsik[3]=bz[ib]-bz[kb];
    rsik[0]= rsik[1]*rsik[1]+rsik[2]*rsik[2]+rsik[3]*rsik[3];
    rsjk[1]=bx[jb]-bx[kb];
    rsjk[2]=by[jb]-by[kb];
    rsjk[3]=bz[jb]-bz[kb];
    rsjk[0]= rsjk[1]*rsjk[1]+rsjk[2]*rsjk[2]+rsjk[3]*rsjk[3];
    for(k=1;k<4;k++){
      rsij2[k]=rsij[k]*rsij[k];
      rsik2[k]=rsik[k]*rsik[k];
      rsjk2[k]=rsjk[k]*rsjk[k];
    }
    gn[1]= rsij[2]*rsik[3] -rsij[3]*rsik[2];
    gn[2]= rsij[3]*rsik[1] -rsij[1]*rsik[3];
    gn[3]= rsij[1]*rsik[2] -rsij[2]*rsik[1];
    gn[4]= gn[1]*gn[1]+gn[2]*gn[2]+gn[3]*gn[3];
    gn[0]= sqrt(gn[4]);
    
		/* area forece, volume force and curvature potential*/
    calaf(ib,jb,kb,rsij,rsik,rsjk,rsij2,rsik2,rsjk2,fxa,fya,fza,ar,gn);
    calvof(ib,jb,kb,bx,by,bz,fxvo,fyvo,fzvo,ar,gn);
    calrn(ib,jb,kb,ial,jal,kal,bx,by,bz,fcv0,fcv1,fcv2,fcvv,rsij,ari,rni,rsik,rsjk,gn);
  }
  
  for(i=0;i<NUMBER;i++){
    if(ari[i]<ARMIN){
      icnt[2]++;
    }
    yar=1./ari[i];

    eni[i]=(rni[i][0]*rni[i][0]+rni[i][1]*rni[i][1]+rni[i][2]*rni[i][2])*yar;
    for(j=0;j<nbon_all[i];j++){
      jb=ibon_all[i][j];
      j1=j+1;

      for(k=0;k<3;k++)
	fcv0[i][j1][k]= fcv0[i][j1][k]*rni[i][0]+fcv1[i][j1][k]*rni[i][1]+fcv2[i][j1][k]*rni[i][2];
      fx[jb] =fx[jb] + KAPPA*(eni[i]*fcvv[i][j1][0]+2.0*fcv0[i][j1][0])*yar;
      fy[jb] =fy[jb] + KAPPA*(eni[i]*fcvv[i][j1][1]+2.0*fcv0[i][j1][1])*yar;
      fz[jb] =fz[jb] + KAPPA*(eni[i]*fcvv[i][j1][2]+2.0*fcv0[i][j1][2])*yar;
    }
    for(k=0;k<3;k++)
      fcv0[i][0][k]= fcv0[i][0][k]*rni[i][0]+fcv1[i][0][k]*rni[i][1]+fcv2[i][0][k]*rni[i][2];
    fx[i] =fx[i] + KAPPA*(eni[i]*fcvv[i][0][0]+2.0*fcv0[i][0][0])*yar;
    fy[i] =fy[i] + KAPPA*(eni[i]*fcvv[i][0][1]+2.0*fcv0[i][0][1])*yar;
    fz[i] =fz[i] + KAPPA*(eni[i]*fcvv[i][0][2]+2.0*fcv0[i][0][2])*yar;

    en[3]=en[3] + eni[i];
    ar0=ar0+ari[i];
  }
  ar0=ar0*0.125;
  ar[0]=ar[0]*0.5;
  en[3]=en[3]*KAPPA;
	/*area force finish*/
  r0= 0.5*KARE*(ar[0]-ARE0);
  ar[1]=ar[1]/6.;
  ftmp= (KVOL/6.)*(ar[1]-vol0);
  for(i=0;i<NUMBER;i++){
    fx[i]=fx[i]+ r0*fxa[i] +ftmp*fxvo[i];
    fy[i]=fy[i]+ r0*fya[i] +ftmp*fyvo[i];
    fz[i]=fz[i]+ r0*fza[i] +ftmp*fzvo[i];
  }
  en[4]= 0.5*KARE*(ar[0]-ARE0)*(ar[0]-ARE0);
  en[5]= 0.5*KVOL*(ar[1]-vol0)*(ar[1]-vol0);
	//printf("test vol %lg %lg \n", ar[1], vol0);
  if( fabs(ar[0]-ar0)>0.0001)
    printf("ar %lg     %lg %lg area\n",ar[0]-ar0,ar[0],ar0);
   
  /*****      2-pair repulsive interaction *****/

  for(i=0;i<NUMBER;i++){
    for(j=ibooki[i];j<ibooki[i+1];j++){
      jb=ibook[j];

      rsij[1]=bx[jb]-bx[i];
      rsij[2]=by[jb]-by[i];
      rsij[3]=bz[jb]-bz[i];
      rsij[0]= rsij[1]*rsij[1]+rsij[2]*rsij[2]+rsij[3]*rsij[3];
      r0=sqrt(rsij[0]);
      if(rsij[0]>RCUTIN2 ){
	ftmp= 0.;
	swct= 0.;
      }else if(rsij[0]>BL_MINOU2 ){
	ftmp= 1./(r0-BL_MIN);
	swct=1./(r0-RCUTI);
      }else{
	ftmp= 1./(BL_MINOU-BL_MIN);
	swct=1./(BL_MINOU-RCUTI);
	icnt[1]++;
      }
      fe=exp(swct+LN80)*ftmp;
      f0= (ftmp+swct*swct)*fe/r0;
      en[2] =en[2] + fe;
      fx[i]=fx[i]-f0*rsij[1];
      fx[jb]=fx[jb]+f0*rsij[1];
      fy[i]=fy[i]-f0*rsij[2];
      fy[jb]=fy[jb]+f0*rsij[2];
      fz[i]=fz[i]-f0*rsij[3];
      fz[jb]=fz[jb]+f0*rsij[3];
    }
  }

  /* external force, pulling at two oposing beads  */

	/*fx[index_pull[0]] += pull[0][0];
	fy[index_pull[0]] += pull[0][1];
	fz[index_pull[0]] += pull[0][2];
	fx[index_pull[1]] += pull[1][0];
	fy[index_pull[1]] += pull[1][1];
	fz[index_pull[1]] += pull[1][2]; */  
	
}

void calaf(ib,jb,kb,rsij,rsik,rsjk,rsij2,rsik2,rsjk2,fxa,fya,fza,ar,gn)
     int ib,jb,kb;
     double rsij[],rsik[],rsjk[],rsij2[],rsik2[],rsjk2[];
     double  fxa[],  fya[],  fza[], ar[],gn[];
{
	/* area force calculation */
  double af0;
  // gn[0]: length of normal vector, af0: normalization factor
  af0= 1./gn[0];
  fxa[ib]=fxa[ib] +af0*( -rsij[1]*(rsjk2[2]+rsjk2[3])+rsjk[1]*(rsij[2]*rsjk[2]+rsij[3]*rsjk[3]));
  fya[ib]=fya[ib] +af0*( -rsij[2]*(rsjk2[1]+rsjk2[3])+rsjk[2]*(rsij[1]*rsjk[1]+rsij[3]*rsjk[3]));
  fza[ib]=fza[ib] +af0*( -rsij[3]*(rsjk2[1]+rsjk2[2])+rsjk[3]*(rsij[1]*rsjk[1]+rsij[2]*rsjk[2]));
  fxa[jb]=fxa[jb] +af0*(rsij[1]*(rsik2[2]+rsik2[3])-rsik[1]*(rsij[2]*rsik[2]+rsij[3]*rsik[3]));
  fya[jb]=fya[jb] +af0*(rsij[2]*(rsik2[1]+rsik2[3])-rsik[2]*(rsij[1]*rsik[1]+rsij[3]*rsik[3]));
  fza[jb]=fza[jb] +af0*(rsij[3]*(rsik2[1]+rsik2[2])-rsik[3]*(rsij[1]*rsik[1]+rsij[2]*rsik[2]));
  fxa[kb]=fxa[kb] +af0*(rsik[1]*(rsij2[2]+rsij2[3])-rsij[1]*(rsij[2]*rsik[2]+rsij[3]*rsik[3]));
  fya[kb]=fya[kb] +af0*(rsik[2]*(rsij2[1]+rsij2[3])-rsij[2]*(rsij[1]*rsik[1]+rsij[3]*rsik[3]));
  fza[kb]=fza[kb] +af0*(rsik[3]*(rsij2[1]+rsij2[2])-rsij[3]*(rsij[1]*rsik[1]+rsij[2]*rsik[2]));
  ar[0]= ar[0] + gn[0];
}

void calvof(ib,jb,kb,bx,by,bz,fxvo,fyvo,fzvo,ar,gn)
     int ib,jb,kb;
     double bx[],by[],bz[];
     double  fxvo[],  fyvo[],  fzvo[], ar[],gn[];
{
	/* helps calculate volume force */
  fxvo[ib]=fxvo[ib] - gn[1];
  fxvo[jb]=fxvo[jb] - gn[1];
  fxvo[kb]=fxvo[kb] - gn[1];
  fyvo[ib]=fyvo[ib] - gn[2];
  fyvo[jb]=fyvo[jb] - gn[2];
  fyvo[kb]=fyvo[kb] - gn[2];
  fzvo[ib]=fzvo[ib] - gn[3];
  fzvo[jb]=fzvo[jb] - gn[3];
  fzvo[kb]=fzvo[kb] - gn[3];

  ar[1]= ar[1] + (gn[1]*bx[ib]+gn[2]*by[ib]+gn[3]*bz[ib]);
}

void calrn(i,jb,kb,ial,jal,kal,bx,by,bz,fcv0,fcv1,fcv2,fcvv,rsij,ari,rni,rsik,rsjk,gn)
     int i,jb,kb,ial[],jal[],kal[];
     double bx[], by[], bz[];
     double fcv0[NUMBER][MAXBOND+1][3];
     double fcv1[NUMBER][MAXBOND+1][3];
     double fcv2[NUMBER][MAXBOND+1][3];
     double fcvv[NUMBER][MAXBOND+1][3];
     double rsij[],ari[],rni[NUMBER][3],rsik[],rsjk[],gn[];
{
	/* calculate curvature potential */
  int k;
  double sgij[4],sgik[4],sgjk[4],sar[2];
  double rpdi,rpdj,rpdk,sg0,stmp0,stmp1,stmp2,stmp3;
  
  rpdi= rsij[1]*rsik[1]+rsij[2]*rsik[2]+rsij[3]*rsik[3];
  rpdj= -(rsij[1]*rsjk[1]+rsij[2]*rsjk[2]+rsij[3]*rsjk[3]);
  rpdk= rsik[1]*rsjk[1]+rsik[2]*rsjk[2]+rsik[3]*rsjk[3];

  sar[0]=1./gn[0];
  sar[1]=1./gn[4];

  sg0= rpdi*sar[0];
  for(k=0;k<4;k++)
    sgjk[k]= sg0*rsjk[k];
  fcv0[jb][0][0]    = fcv0[jb][0][0]- sg0;
  fcv0[jb][jal[1]][0]=fcv0[jb][jal[1]][0]+ sg0;
  fcv0[kb][0][0]    = fcv0[kb][0][0]- sg0;
  fcv0[kb][kal[1]][0]=fcv0[kb][kal[1]][0]+ sg0;
  fcv1[jb][0][1]    = fcv1[jb][0][1]- sg0;
  fcv1[jb][jal[1]][1]=fcv1[jb][jal[1]][1]+ sg0;
  fcv1[kb][0][1]    = fcv1[kb][0][1]- sg0;
  fcv1[kb][kal[1]][1]=fcv1[kb][kal[1]][1]+ sg0;
  fcv2[jb][0][2]    = fcv2[jb][0][2]- sg0;
  fcv2[jb][jal[1]][2]=fcv2[jb][jal[1]][2]+ sg0;
  fcv2[kb][0][2]    = fcv2[kb][0][2]- sg0;
  fcv2[kb][kal[1]][2]=fcv2[kb][kal[1]][2]+ sg0;
  for(k=0;k<3;k++){
    stmp0=2.*sgjk[k+1];
    fcvv[jb][0][k]=fcvv[jb][0][k]+ stmp0;
    fcvv[jb][jal[1]][k]=fcvv[jb][jal[1]][k]- stmp0;
    fcvv[kb][0][k]=fcvv[kb][0][k]- stmp0;
    fcvv[kb][kal[1]][k]=fcvv[kb][kal[1]][k]+ stmp0;
  }
  sg0= rpdj*sar[0];
  for(k=0;k<4;k++)
    sgik[k]= sg0*rsik[k];
  fcv0[i][0][0]     = fcv0[i][0][0]- sg0;
  fcv0[i][ial[1]][0] =fcv0[i][ial[1]][0] + sg0;
  fcv0[kb][0][0]    = fcv0[kb][0][0]- sg0;
  fcv0[kb][kal[0]][0]=fcv0[kb][kal[0]][0]+ sg0;
  fcv1[i][0][1]     = fcv1[i][0][1]- sg0;
  fcv1[i][ial[1]][1] =fcv1[i][ial[1]][1] + sg0;
  fcv1[kb][0][1]    = fcv1[kb][0][1]- sg0;
  fcv1[kb][kal[0]][1]=fcv1[kb][kal[0]][1]+ sg0;
  fcv2[i][0][2]     = fcv2[i][0][2]- sg0;
  fcv2[i][ial[1]][2] =fcv2[i][ial[1]][2] + sg0;
  fcv2[kb][0][2]    = fcv2[kb][0][2]- sg0;
  fcv2[kb][kal[0]][2]=fcv2[kb][kal[0]][2]+ sg0;
  for(k=0;k<3;k++){
    stmp0=2.*sgik[k+1];
    fcvv[i][0][k]=fcvv[i][0][k]+ stmp0;
    fcvv[i][ial[1]][k]=fcvv[i][ial[1]][k]- stmp0;
    fcvv[kb][0][k]=fcvv[kb][0][k]- stmp0;
    fcvv[kb][kal[0]][k]=fcvv[kb][kal[0]][k]+ stmp0;
  }
  sg0= rpdk*sar[0];
  for(k=0;k<4;k++)
    sgij[k]= sg0*rsij[k];
  fcv0[i][0][0]     = fcv0[i][0][0] - sg0;
  fcv0[i][ial[0]][0] =fcv0[i][ial[0]][0] + sg0;
  fcv0[jb][0][0]    = fcv0[jb][0][0]- sg0;
  fcv0[jb][jal[0]][0]=fcv0[jb][jal[0]][0]+ sg0;
  fcv1[i][0][1]     = fcv1[i][0][1] - sg0;
  fcv1[i][ial[0]][1] =fcv1[i][ial[0]][1] + sg0;
  fcv1[jb][0][1]    = fcv1[jb][0][1]- sg0;
  fcv1[jb][jal[0]][1]=fcv1[jb][jal[0]][1]+ sg0;
  fcv2[i][0][2]     = fcv2[i][0][2] - sg0;
  fcv2[i][ial[0]][2] =fcv2[i][ial[0]][2] + sg0;
  fcv2[jb][0][2]    = fcv2[jb][0][2]- sg0;
  fcv2[jb][jal[0]][2]=fcv2[jb][jal[0]][2]+ sg0;
  for(k=0;k<3;k++){
    stmp0=2.*sgij[k+1];
    fcvv[i][0][k]=fcvv[i][0][k]+ stmp0;
    fcvv[i][ial[0]][k]=fcvv[i][ial[0]][k]- stmp0;
    fcvv[jb][0][k]=fcvv[jb][0][k]- stmp0;
    fcvv[jb][jal[0]][k]=fcvv[jb][jal[0]][k]+ stmp0;
  }

  ari[i]=ari[i]+ (sgij[0]+sgik[0]);
  rni[i][0]=rni[i][0] +(sgij[1]+sgik[1]);
  rni[i][1]=rni[i][1] +(sgij[2]+sgik[2]);
  rni[i][2]=rni[i][2] +(sgij[3]+sgik[3]);
  ari[jb]=ari[jb]+ (sgij[0]+ sgjk[0]);
  rni[jb][0]=rni[jb][0] +(-sgij[1]+sgjk[1]);
  rni[jb][1]=rni[jb][1] +(-sgij[2]+sgjk[2]);
  rni[jb][2]=rni[jb][2] +(-sgij[3]+sgjk[3]);
  ari[kb]=ari[kb]+ (sgik[0]+ sgjk[0]);
  rni[kb][0]=rni[kb][0] -(sgik[1]+sgjk[1]);
  rni[kb][1]=rni[kb][1] -(sgik[2]+sgjk[2]);
  rni[kb][2]=rni[kb][2] -(sgik[3]+sgjk[3]);

  for(k=0;k<3;k++){
    sg0=(rsjk[k+1]-rpdk*sar[1]*(rsjk[0]*rsik[k+1]-rpdk*rsjk[k+1]) )*sar[0];
    stmp0=sg0*rsij[0];
    stmp1=sg0*rsij[1];
    stmp2=sg0*rsij[2];
    stmp3=sg0*rsij[3];
    fcvv[i][0][k]   =  fcvv[i][0][k]     + stmp0;
    fcvv[i][ial[1]][k]=fcvv[i][ial[1]][k]- stmp0;
    fcvv[jb][jal[0]][k]=fcvv[jb][jal[0]][k]+ stmp0;
    fcvv[jb][jal[1]][k]=fcvv[jb][jal[1]][k]- stmp0;
    fcv0[i][0][k]   =  fcv0[i][0][k]     - stmp1;
    fcv0[i][ial[1]][k]=fcv0[i][ial[1]][k]+ stmp1;
    fcv0[jb][jal[0]][k]=fcv0[jb][jal[0]][k]+ stmp1;
    fcv0[jb][jal[1]][k]=fcv0[jb][jal[1]][k]- stmp1;
    fcv1[i][0][k]   =  fcv1[i][0][k]     - stmp2;
    fcv1[i][ial[1]][k]=fcv1[i][ial[1]][k]+ stmp2;
    fcv1[jb][jal[0]][k]=fcv1[jb][jal[0]][k]+ stmp2;
    fcv1[jb][jal[1]][k]=fcv1[jb][jal[1]][k]- stmp2;
    fcv2[i][0][k]   =  fcv2[i][0][k]     - stmp3;
    fcv2[i][ial[1]][k]=fcv2[i][ial[1]][k]+ stmp3;
    fcv2[jb][jal[0]][k]=fcv2[jb][jal[0]][k]+ stmp3;
    fcv2[jb][jal[1]][k]=fcv2[jb][jal[1]][k]- stmp3;
  }
  for(k=0;k<3;k++){
    sg0=(rsik[k+1]-rpdk*sar[1]*(rsik[0]*rsjk[k+1]-rpdk*rsik[k+1]) )*sar[0];
    stmp0=sg0*rsij[0];
    stmp1=sg0*rsij[1];
    stmp2=sg0*rsij[2];
    stmp3=sg0*rsij[3];
    fcvv[i][ial[0]][k]=fcvv[i][ial[0]][k]+ stmp0;
    fcvv[i][ial[1]][k]=fcvv[i][ial[1]][k]- stmp0;
    fcvv[jb][0][k]  =   fcvv[jb][0][k]     + stmp0;
    fcvv[jb][jal[1]][k]=fcvv[jb][jal[1]][k]- stmp0;
    fcv0[i][ial[0]][k]=fcv0[i][ial[0]][k]- stmp1;
    fcv0[i][ial[1]][k]=fcv0[i][ial[1]][k]+ stmp1;
    fcv0[jb][0][k]    = fcv0[jb][0][k]     + stmp1;
    fcv0[jb][jal[1]][k]=fcv0[jb][jal[1]][k]- stmp1;
    fcv1[i][ial[0]][k]=fcv1[i][ial[0]][k]- stmp2;
    fcv1[i][ial[1]][k]=fcv1[i][ial[1]][k]+ stmp2;
    fcv1[jb][0][k]    = fcv1[jb][0][k]     + stmp2;
    fcv1[jb][jal[1]][k]=fcv1[jb][jal[1]][k]- stmp2;
    fcv2[i][ial[0]][k]=fcv2[i][ial[0]][k]- stmp3;
    fcv2[i][ial[1]][k]=fcv2[i][ial[1]][k]+ stmp3;
    fcv2[jb][0][k]    = fcv2[jb][0][k]     + stmp3;
    fcv2[jb][jal[1]][k]=fcv2[jb][jal[1]][k]- stmp3;
  }

  for(k=0;k<3;k++){
    sg0=(-rsjk[k+1]-rpdj*sar[1]*(rsjk[0]*rsij[k+1]+rpdj*rsjk[k+1]) )*sar[0];
    stmp0=sg0*rsik[0];
    stmp1=sg0*rsik[1];
    stmp2=sg0*rsik[2];
    stmp3=sg0*rsik[3];
    fcvv[i][0][k]   =  fcvv[i][0][k]     + stmp0;
    fcvv[i][ial[0]][k]=fcvv[i][ial[0]][k]- stmp0;
    fcvv[kb][kal[0]][k]=fcvv[kb][kal[0]][k]+ stmp0;
    fcvv[kb][kal[1]][k]=fcvv[kb][kal[1]][k]- stmp0;
    fcv0[i][0][k]   =  fcv0[i][0][k]     - stmp1;
    fcv0[i][ial[0]][k]=fcv0[i][ial[0]][k]+ stmp1;
    fcv0[kb][kal[0]][k]=fcv0[kb][kal[0]][k]+ stmp1;
    fcv0[kb][kal[1]][k]=fcv0[kb][kal[1]][k]- stmp1;
    fcv1[i][0][k]   =  fcv1[i][0][k]     - stmp2;
    fcv1[i][ial[0]][k]=fcv1[i][ial[0]][k]+ stmp2;
    fcv1[kb][kal[0]][k]=fcv1[kb][kal[0]][k]+ stmp2;
    fcv1[kb][kal[1]][k]=fcv1[kb][kal[1]][k]- stmp2;
    fcv2[i][0][k]   =  fcv2[i][0][k]     - stmp3;
    fcv2[i][ial[0]][k]=fcv2[i][ial[0]][k]+ stmp3;
    fcv2[kb][kal[0]][k]=fcv2[kb][kal[0]][k]+ stmp3;
    fcv2[kb][kal[1]][k]=fcv2[kb][kal[1]][k]- stmp3;
  }
  for(k=0;k<3;k++){
    sg0=(rsij[k+1]-rpdj*sar[1]*(-rsij[0]*rsjk[k+1]-rpdj*rsij[k+1]) )*sar[0];
    stmp0=sg0*rsik[0];
    stmp1=sg0*rsik[1];
    stmp2=sg0*rsik[2];
    stmp3=sg0*rsik[3];
    fcvv[i][ial[1]][k]=fcvv[i][ial[1]][k]+ stmp0;
    fcvv[i][ial[0]][k]=fcvv[i][ial[0]][k]- stmp0;
    fcvv[kb][0][k]  =   fcvv[kb][0][k]     + stmp0;
    fcvv[kb][kal[1]][k]=fcvv[kb][kal[1]][k]- stmp0;
    fcv0[i][ial[1]][k]=fcv0[i][ial[1]][k]- stmp1;
    fcv0[i][ial[0]][k]=fcv0[i][ial[0]][k]+ stmp1;
    fcv0[kb][0][k]    = fcv0[kb][0][k]     + stmp1;
    fcv0[kb][kal[1]][k]=fcv0[kb][kal[1]][k]- stmp1;
    fcv1[i][ial[1]][k]=fcv1[i][ial[1]][k]- stmp2;
    fcv1[i][ial[0]][k]=fcv1[i][ial[0]][k]+ stmp2;
    fcv1[kb][0][k]    = fcv1[kb][0][k]     + stmp2;
    fcv1[kb][kal[1]][k]=fcv1[kb][kal[1]][k]- stmp2;
    fcv2[i][ial[1]][k]=fcv2[i][ial[1]][k]- stmp3;
    fcv2[i][ial[0]][k]=fcv2[i][ial[0]][k]+ stmp3;
    fcv2[kb][0][k]    = fcv2[kb][0][k]     + stmp3;
    fcv2[kb][kal[1]][k]=fcv2[kb][kal[1]][k]- stmp3;
  }

  for(k=0;k<3;k++){
    sg0=(-rsik[k+1]-rpdi*sar[1]*(-rsik[0]*rsij[k+1]+rpdi*rsik[k+1]) )*sar[0];
    stmp0=sg0*rsjk[0];
    stmp1=sg0*rsjk[1];
    stmp2=sg0*rsjk[2];
    stmp3=sg0*rsjk[3];
    fcvv[jb][0][k]   =  fcvv[jb][0][k]     + stmp0;
    fcvv[jb][jal[0]][k]=fcvv[jb][jal[0]][k]- stmp0;
    fcvv[kb][kal[1]][k]=fcvv[kb][kal[1]][k]+ stmp0;
    fcvv[kb][kal[0]][k]=fcvv[kb][kal[0]][k]- stmp0;
    fcv0[jb][0][k]   =  fcv0[jb][0][k]     - stmp1;
    fcv0[jb][jal[0]][k]=fcv0[jb][jal[0]][k]+ stmp1;
    fcv0[kb][kal[1]][k]=fcv0[kb][kal[1]][k]+ stmp1;
    fcv0[kb][kal[0]][k]=fcv0[kb][kal[0]][k]- stmp1;
    fcv1[jb][0][k]   =  fcv1[jb][0][k]     - stmp2;
    fcv1[jb][jal[0]][k]=fcv1[jb][jal[0]][k]+ stmp2;
    fcv1[kb][kal[1]][k]=fcv1[kb][kal[1]][k]+ stmp2;
    fcv1[kb][kal[0]][k]=fcv1[kb][kal[0]][k]- stmp2;
    fcv2[jb][0][k]   =  fcv2[jb][0][k]     - stmp3;
    fcv2[jb][jal[0]][k]=fcv2[jb][jal[0]][k]+ stmp3;
    fcv2[kb][kal[1]][k]=fcv2[kb][kal[1]][k]+ stmp3;
    fcv2[kb][kal[0]][k]=fcv2[kb][kal[0]][k]- stmp3;
  }
  for(k=0;k<3;k++){
    sg0=(-rsij[k+1]-rpdi*sar[1]*(-rsij[0]*rsik[k+1]+rpdi*rsij[k+1]) )*sar[0];
    stmp0=sg0*rsjk[0];
    stmp1=sg0*rsjk[1];
    stmp2=sg0*rsjk[2];
    stmp3=sg0*rsjk[3];
    fcvv[jb][jal[1]][k]=fcvv[jb][jal[1]][k]+ stmp0;
    fcvv[jb][jal[0]][k]=fcvv[jb][jal[0]][k]- stmp0;
    fcvv[kb][0][k]  =   fcvv[kb][0][k]     + stmp0;
    fcvv[kb][kal[0]][k]=fcvv[kb][kal[0]][k]- stmp0;
    fcv0[jb][jal[1]][k]=fcv0[jb][jal[1]][k]- stmp1;
    fcv0[jb][jal[0]][k]=fcv0[jb][jal[0]][k]+ stmp1;
    fcv0[kb][0][k]    = fcv0[kb][0][k]     + stmp1;
    fcv0[kb][kal[0]][k]=fcv0[kb][kal[0]][k]- stmp1;
    fcv1[jb][jal[1]][k]=fcv1[jb][jal[1]][k]- stmp2;
    fcv1[jb][jal[0]][k]=fcv1[jb][jal[0]][k]+ stmp2;
    fcv1[kb][0][k]    = fcv1[kb][0][k]     + stmp2;
    fcv1[kb][kal[0]][k]=fcv1[kb][kal[0]][k]- stmp2;
    fcv2[jb][jal[1]][k]=fcv2[jb][jal[1]][k]- stmp3;
    fcv2[jb][jal[0]][k]=fcv2[jb][jal[0]][k]+ stmp3;
    fcv2[kb][0][k]    = fcv2[kb][0][k]     + stmp3;
    fcv2[kb][kal[0]][k]=fcv2[kb][kal[0]][k]- stmp3;
  }
}

void bonflip(bx,by,bz,ari,rni,eni,en,ar,vol0,iface,iedge,ibon_all,nbon_all,iedb,iedf,ifaed,icnt,ibdmami,idum)
     double  bx[],  by[],  bz[];
     double ari[],rni[NUMBER][3],eni[],en[],ar[],vol0;
     int iface[NF][3],iedge[NE][2];
     int ibon_all[NUMBER][MAXBOND],nbon_all[];
     int iedb[NE][2],iedf[NE][2],ifaed[NF][3];
     int  icnt[],ibdmami[]; 
     long *idum;    
{
  int i,k,ib,jb,it0,it1,i0,imc,if0,if1,iftk,ke,k0,k1,ke0,ke1,kb,kf,kk;
  double r0,eold,enew,esa,esa_ben,esa0,arsa,ar_new,esa_a,ea_new;
  double rs01[4],rsij[4],rsi0[4],rsi1[4],rsj0[4],rsj1[4];
  double ari_new[4],rni_new[4][3],eni_new[4];
  double vasa[2],vo_new,esa_vo,evo_new;
  int iface_new0[3],iface_new1[3];

  for(imc=0;imc<NE;imc++){

		i0=ran2(idum)*NE;
		it0=iedb[i0][0];
		it1=iedb[i0][1];


    /*   avoid overlap of bond and 3 faces   */
    ke0=0;
    for(k=0;k<nbon_all[it0];k++){
      k1=ibon_all[it0][k];
      if(it1==k1){
	ke0=4;
	break;
      }
      for(kk=0;kk<nbon_all[it1];kk++){
	if(k1==ibon_all[it1][kk]){
	  ke0++;
	  break;
	}
      }
    }
    if( ke0>2 )
      break;

    rs01[1]=bx[it0]-bx[it1];
    rs01[2]=by[it0]-by[it1];
    rs01[3]=bz[it0]-bz[it1];
    rs01[0]= rs01[1]*rs01[1]+rs01[2]*rs01[2]+rs01[3]*rs01[3];
    if(rs01[0]<BL_MAXIN2){
      if(rs01[0]>RCUTOU2 ){
	r0=sqrt(rs01[0]);
	enew= exp(1./(RCUTO-r0)+LN80)/(BL_MAX-r0);
      }else{
	enew= 0.;
      }
      
      ib=iedge[i0][0];
      jb=iedge[i0][1];
      if0=iedf[i0][0];
      if1=iedf[i0][1];
      
      
      rsij[1]=bx[ib]-bx[jb];
      rsij[2]=by[ib]-by[jb];
      rsij[3]=bz[ib]-bz[jb];
      rsij[0]= rsij[1]*rsij[1]+rsij[2]*rsij[2]+rsij[3]*rsij[3];
      if(rsij[0]>RCUTOU2 ){
	if(rsij[0]<BL_MAXIN2 )
	  r0=sqrt(rsij[0]);
	else
	  r0=BL_MAXIN;
	eold= exp(1./(RCUTO-r0)+LN80)/(BL_MAX-r0);
	esa= enew -eold;
      }else{
	eold= 0.0;
	esa= enew;
      }
      
      rsi0[1]=bx[ib]-bx[it0];
      rsi0[2]=by[ib]-by[it0];
      rsi0[3]=bz[ib]-bz[it0];
      rsi0[0]= rsi0[1]*rsi0[1]+rsi0[2]*rsi0[2]+rsi0[3]*rsi0[3];
      rsi1[1]=bx[ib]-bx[it1];
      rsi1[2]=by[ib]-by[it1];
      rsi1[3]=bz[ib]-bz[it1];
      rsi1[0]= rsi1[1]*rsi1[1]+rsi1[2]*rsi1[2]+rsi1[3]*rsi1[3];
      rsj0[1]=bx[jb]-bx[it0];
      rsj0[2]=by[jb]-by[it0];
      rsj0[3]=bz[jb]-bz[it0];
      rsj0[0]= rsj0[1]*rsj0[1]+rsj0[2]*rsj0[2]+rsj0[3]*rsj0[3];
      rsj1[1]=bx[jb]-bx[it1];
      rsj1[2]=by[jb]-by[it1];
      rsj1[3]=bz[jb]-bz[it1];
      rsj1[0]= rsj1[1]*rsj1[1]+rsj1[2]*rsj1[2]+rsj1[3]*rsj1[3];
      
      ari_new[0]=ari[ib];
      ari_new[1]=ari[jb];
      ari_new[2]=ari[it0];
      ari_new[3]=ari[it1];
      for(k=0;k<3;k++){
	rni_new[0][k]=rni[ib][k];
	rni_new[1][k]=rni[jb][k];
	rni_new[2][k]=rni[it0][k];
	rni_new[3][k]=rni[it1][k];
      }
      
      for(k=0;k<3;k++){
	if(it0==iface[if0][k]){
	  iftk=k;
	  break;
	}
      }

      iface_new0[0]=it0;
      iface_new1[0]=it0;
      for(k=0;k<3;k++){
	if(ib==iface[if0][k]){
	  if( ((k-iftk)==1)||((k-iftk)==-2) ){
	    iface_new0[1]=ib;
	    iface_new0[2]=it1;
	    iface_new1[1]=it1;
	    iface_new1[2]=jb;
	  }else{
	    iface_new0[1]=it1;
	    iface_new0[2]=ib;
	    iface_new1[1]=jb;
	    iface_new1[2]=it1;
	  }
	  break;
	}
      }
      
      vasa[0]=0.;
      vasa[1]=0.;
      calrn_old(2,ib,jb,it0,rsij,rsi0,rsj0,ari_new,rni_new,vasa,bx,by,bz,iface[if0]);
      calrn_old(3,ib,jb,it1,rsij,rsi1,rsj1,ari_new,rni_new,vasa,bx,by,bz,iface[if1]);
      calrn_new(0,ib,it0,it1,rsi0,rsi1,rs01,ari_new,rni_new,vasa,bx,by,bz,iface_new0);
      calrn_new(1,jb,it0,it1,rsj0,rsj1,rs01,ari_new,rni_new,vasa,bx,by,bz,iface_new1);
      vasa[0]=vasa[0]*0.5;
      vasa[1]=vasa[1]/6.;
      
      esa_ben= -(eni[ib]+eni[jb]+eni[it0]+eni[it1]);
      for(k=0;k<4;k++){
	eni_new[k]= (rni_new[k][0]*rni_new[k][0]+rni_new[k][1]*rni_new[k][1]+rni_new[k][2]*rni_new[k][2])/ari_new[k];
	esa_ben=esa_ben+ eni_new[k];
      }
      esa_ben= esa_ben*KAPPA;
      arsa= -(ari[ib]+ari[jb]+ari[it0]+ari[it1]);
      for(k=0;k<4;k++)
	arsa= arsa+ ari_new[k];
      arsa= arsa*0.125;
      ar_new=ar[0]+ vasa[0];
      ea_new= 0.5*KARE*(ar_new-ARE0)*(ar_new-ARE0);
      esa_a=ea_new-en[4];
      
      vo_new= ar[1]+vasa[1];
      evo_new= 0.5*KVOL*(vo_new-vol0)*(vo_new-vol0);
      esa_vo=evo_new-en[5];
      esa0=esa + esa_ben +esa_a+esa_vo;
      if( (esa0<=0.)||( exp(BETA*esa0)>(1.-ran2(idum)) ) ){
	
	ari[ib]= ari_new[0];
	ari[jb]= ari_new[1];
	ari[it0]=ari_new[2];
	ari[it1]=ari_new[3];
	eni[ib]= eni_new[0];
	eni[jb]= eni_new[1];
	eni[it0]=eni_new[2];
	eni[it1]=eni_new[3];
	en[3]=en[3]+ esa_ben;
	en[1]=en[1]+ esa;
	ar[0]=ar_new;
	ar[1]=vo_new;
	en[4]=ea_new;
	en[5]=evo_new;
	for(k=0;k<3;k++){
	  rni[ib][k]= rni_new[0][k];
	  rni[jb][k]= rni_new[1][k];
	  rni[it0][k]=rni_new[2][k];
	  rni[it1][k]=rni_new[3][k];
	}
	
	iedge[i0][0]=it0;
	iedge[i0][1]=it1;
	iedb[i0][0]=ib;
	iedb[i0][1]=jb;
	for(k=0;k<3;k++){
	  iface[if0][k]=iface_new0[k];
	  iface[if1][k]=iface_new1[k];
	}
	for(k=0;k<3;k++){
	  if(ifaed[if0][k] != i0){
	    ke= ifaed[if0][k];
	    for(kk=0;kk<2;kk++){
	      if(iedge[ke][kk]==ib){
		kb=jb;
		kf=if0;
		break;
	      }else if(iedge[ke][kk]==jb){
		ke0=ke;
		k0=k;
		kb=ib;
		kf=if1;
		break;
	      }
	    }
	    for(kk=0;kk<2;kk++){
	      if(iedb[ke][kk]==kb){
		iedb[ke][kk]=it1;
		iedf[ke][kk]=kf;
		break;
	      }
	    }
	  }
	  if(ifaed[if1][k] != i0){
	    ke= ifaed[if1][k];
	    for(kk=0;kk<2;kk++){
	      if(iedge[ke][kk]==ib){
		ke1=ke;
		k1=k;
		kb=jb;
		kf=if0;
		break;
	      }else if(iedge[ke][kk]==jb){
		kb=ib;
		kf=if1;
		break;
	      }
	    }
	    for(kk=0;kk<2;kk++){
	      if(iedb[ke][kk]==kb){
		iedb[ke][kk]=it0;
		iedf[ke][kk]=kf;
		break;
	      }
	    }
	  }
	}
	ifaed[if0][k0]=ke1;
	ifaed[if1][k1]=ke0;
	icnt[0]++;

	for(k=0;k<nbon_all[ib];k++){
	  if(jb==ibon_all[ib][k]){
	    k1=k;
	    break;
	  }
	}
	for(k=k1;k<nbon_all[ib]-1;k++)
	  ibon_all[ib][k]= ibon_all[ib][k+1];

	for(k=0;k<nbon_all[jb];k++){
	  if(ib==ibon_all[jb][k]){
	    k1=k;
	    break;
	  }
	}
	for(k=k1;k<nbon_all[jb]-1;k++)
	  ibon_all[jb][k]= ibon_all[jb][k+1];

	ibon_all[it0][nbon_all[it0]]=it1;
	ibon_all[it1][nbon_all[it1]]=it0;
	nbon_all[ib]--;
	nbon_all[jb]--;
	nbon_all[it0]++;
	nbon_all[it1]++;
      }
    }
  }

  for(i=0;i<NUMBER+1;i++)
    nbon_all[i]=0;
  for(i=0;i<NE;i++){
    ib=iedge[i][0];
    jb=iedge[i][1];
    ibon_all[ib][nbon_all[ib]]=jb;
    nbon_all[ib]++;
    ibon_all[jb][nbon_all[jb]]=ib;
    nbon_all[jb]++;
  }

  for(i=0;i<NUMBER;i++){
    if(nbon_all[i]>ibdmami[0]) ibdmami[0]=nbon_all[i];
    if(nbon_all[i]<ibdmami[1]) ibdmami[1]=nbon_all[i];
  }

}

void calrn_new(n0,i0,it0,it1,rsi0,rsi1,rs01,ari_new,rni_new,vasa,bx,by,bz,iface0)
     int n0,i0,it0,it1;
     double rsi0[],rsi1[],rs01[],ari_new[],rni_new[4][3];
     double vasa[],bx[],by[],bz[];
     int iface0[];
{
  int k;
  double sgi0[4],sgi1[4],sg01[4];
  double rpdi,rpd0,rpd1,sg0;
  double rij[4],rik[4],gn[4],ygn0;

  rij[1]=bx[iface0[0]]-bx[iface0[1]];
  rij[2]=by[iface0[0]]-by[iface0[1]];
  rij[3]=bz[iface0[0]]-bz[iface0[1]];
  rik[1]=bx[iface0[0]]-bx[iface0[2]];
  rik[2]=by[iface0[0]]-by[iface0[2]];
  rik[3]=bz[iface0[0]]-bz[iface0[2]];
  gn[1]= rij[2]*rik[3] -rij[3]*rik[2];
  gn[2]= rij[3]*rik[1] -rij[1]*rik[3];
  gn[3]= rij[1]*rik[2] -rij[2]*rik[1];
  gn[0]=sqrt(gn[1]*gn[1]+gn[2]*gn[2]+gn[3]*gn[3]);
  ygn0= 1./gn[0];
  vasa[0]=vasa[0] + gn[0];
  vasa[1]=vasa[1] + (bx[i0]*gn[1]+by[i0]*gn[2]+bz[i0]*gn[3]);

  rpdi= rsi0[1]*rsi1[1]+rsi0[2]*rsi1[2]+rsi0[3]*rsi1[3];
  rpd0= -(rsi0[1]*rs01[1]+rsi0[2]*rs01[2]+rsi0[3]*rs01[3]);
  rpd1= rsi1[1]*rs01[1]+rsi1[2]*rs01[2]+rsi1[3]*rs01[3];
  
  sg0= rpdi*ygn0;
  for(k=0;k<4;k++)
    sg01[k]= sg0*rs01[k];
  sg0= rpd0*ygn0;
  for(k=0;k<4;k++)
    sgi1[k]= sg0*rsi1[k];
  sg0= rpd1*ygn0;
  for(k=0;k<4;k++)
    sgi0[k]= sg0*rsi0[k];

  ari_new[n0]=ari_new[n0]+ (sgi0[0]+sgi1[0]);
  for(k=0;k<3;k++)
    rni_new[n0][k]=rni_new[n0][k] +(sgi0[k+1]+sgi1[k+1]);
  ari_new[2]=ari_new[2]+ (sgi0[0]+sg01[0]);
  for(k=0;k<3;k++)
    rni_new[2][k]=rni_new[2][k] +(-sgi0[k+1]+sg01[k+1]);
  ari_new[3]=ari_new[3]+ (sgi1[0]+sg01[0]);
  for(k=0;k<3;k++)
    rni_new[3][k]=rni_new[3][k] -(sgi1[k+1]+sg01[k+1]);
}

void calrn_old(n0,i0,it0,it1,rsi0,rsi1,rs01,ari_new,rni_new,vasa,bx,by,bz,iface0)
     int n0,i0,it0,it1;
     double rsi0[],rsi1[],rs01[],ari_new[],rni_new[4][3];
     double vasa[],bx[],by[],bz[];
     int iface0[];
{
  int k;
  double sgi0[4],sgi1[4],sg01[4];
  double rpdi,rpd0,rpd1,sg0;
  double rij[4],rik[4],gn[4],ygn0;

  rij[1]=bx[iface0[0]]-bx[iface0[1]];
  rij[2]=by[iface0[0]]-by[iface0[1]];
  rij[3]=bz[iface0[0]]-bz[iface0[1]];
  rik[1]=bx[iface0[0]]-bx[iface0[2]];
  rik[2]=by[iface0[0]]-by[iface0[2]];
  rik[3]=bz[iface0[0]]-bz[iface0[2]];
  gn[1]= rij[2]*rik[3] -rij[3]*rik[2];
  gn[2]= rij[3]*rik[1] -rij[1]*rik[3];
  gn[3]= rij[1]*rik[2] -rij[2]*rik[1];
  gn[0]=sqrt(gn[1]*gn[1]+gn[2]*gn[2]+gn[3]*gn[3]);
  ygn0= 1./gn[0];
  vasa[0]=vasa[0] - gn[0];
  vasa[1]=vasa[1] - (bx[i0]*gn[1]+by[i0]*gn[2]+bz[i0]*gn[3]);

  rpdi= rsi0[1]*rsi1[1]+rsi0[2]*rsi1[2]+rsi0[3]*rsi1[3];
  rpd0= -(rsi0[1]*rs01[1]+rsi0[2]*rs01[2]+rsi0[3]*rs01[3]);
  rpd1= rsi1[1]*rs01[1]+rsi1[2]*rs01[2]+rsi1[3]*rs01[3];
  
  sg0= rpdi*ygn0;
  for(k=0;k<4;k++)
    sg01[k]= sg0*rs01[k];
  sg0= rpd0*ygn0;
  for(k=0;k<4;k++)
    sgi1[k]= sg0*rsi1[k];
  sg0= rpd1*ygn0;
  for(k=0;k<4;k++)
    sgi0[k]= sg0*rsi0[k];

  ari_new[0]=ari_new[0]- (sgi0[0]+sgi1[0]);
  for(k=0;k<3;k++)
    rni_new[0][k]=rni_new[0][k] -(sgi0[k+1]+sgi1[k+1]);
  ari_new[1]=ari_new[1]- (sgi0[0]+sg01[0]);
  for(k=0;k<3;k++)
    rni_new[1][k]=rni_new[1][k] +(sgi0[k+1]-sg01[k+1]);
  ari_new[n0]=ari_new[n0]- (sgi1[0]+sg01[0]);
  for(k=0;k<3;k++)
    rni_new[n0][k]=rni_new[n0][k] +(sgi1[k+1]+sg01[k+1]);
}

void cal_vbook(bx,by,bz,ibook,ibooki,ibdmami)
     double  bx[], by[], bz[];
     int ibook[],ibooki[], ibdmami[]; 
{
  int i,j,ibookn;
  double xl,yl,zl,r2;
  double bxi,byi,bzi;

  ibookn=0;
  for(i=0;i<NUMBER;i++){
    ibooki[i]= ibookn;
    bxi=bx[i];
    byi=by[i];
    bzi=bz[i];
    for(j=i+1;j<NUMBER;j++){
      xl=bxi-bx[j];
      yl=byi-by[j];
      zl=bzi-bz[j];
      r2= xl*xl+yl*yl+zl*zl;
      if(r2<RBOOK2){
	ibook[ibookn]=j;
	ibookn++;
      }
    }
  }
  ibooki[NUMBER]=ibookn;
  if(ibookn>ibdmami[2])
    ibdmami[2]=ibookn;
}

void read_init(bx,by,bz,bxv,byv,bzv,bxq,byq,bzq,iface,iedge,ibon_all,nbon_all,iedb,iedf,ifaed,mp,index_pull)
     double  bx[], by[], bz[];
     double  bxv[], byv[], bzv[];
     double  bxq[], byq[], bzq[];
     int iface[NF][3],iedge[NE][2];
		 int index_pull[];
     int ibon_all[NUMBER][MAXBOND],nbon_all[];
     int iedb[NE][2],iedf[NE][2],ifaed[NF][3];
     FILE *mp;
{
	/* read starting configuration */
  int i0,i,ib,jb,k,ned,iboni[MAXBOND],nboni0;
  char dummy[256];  
  FILE *finp;

  finp=fopen("start_redVol08.dat","r");

	/* read positions of the vertices */
  for(i=0;i<NUMBER;i++){
    fscanf(finp,"%lg %lg %lg",bx+i,by+i,bz+i);
	}
  fscanf(finp,"\n"); 

	/* read edges */
  ned=0;
  do{
    fgets(dummy,sizeof(dummy),finp);
    nboni0= sscanf(dummy,"%d %d %d %d %d %d %d %d %d %d %d %d %d",&i0,&iboni[0],&iboni[1],&iboni[2],&iboni[3],&iboni[4],&iboni[5],&iboni[6],&iboni[7],&iboni[8],&iboni[9],&iboni[10],&iboni[11]) -1;
    for(k=0;k<nboni0;k++){
      iedge[ned][0]=i0;
      iedge[ned][1]=iboni[k];
      ned++;
    }
  }while( i0<NUMBER-1);

	/* read faces */
  for(i=0;i<NF;i++){
    fscanf(finp,"%d %d %d", iface[i],iface[i]+1,iface[i]+2);
	}

	/* read velocities of the vertices */
  for(i=0;i<NUMBER;i++){
    fscanf(finp,"%lg %lg %lg", bxv+i,byv+i,bzv+i);
	}
  fclose(finp);

	/* arrange neighbors of the vertices */
  for(i=0;i<NUMBER+1;i++)
    nbon_all[i]=0;
  for(i=0;i<NE;i++){
    ib=iedge[i][0];
    jb=iedge[i][1];
    ibon_all[ib][nbon_all[ib]]=jb;
    nbon_all[ib]++;
    ibon_all[jb][nbon_all[jb]]=ib;
    nbon_all[jb]++;
  }

	/* set positions and velocities to the center-of-mass system */
  set_gg(bx, by,bz);
  set_gg(bxv, byv,bzv);

	/* (re)arrange edges and vertices */
  sort_bd(bx, by, bz,bxv, byv,bzv,bxq,byq,bzq,iface,iedge,ibon_all,nbon_all,iedb,iedf,ifaed,mp,index_pull);
  check_bond(ibon_all,nbon_all,mp);

}

void set_gg(bxv, byv,bzv)
     double  bxv[], byv[], bzv[];
{
	/* shift to center-of-mass system */
  int i;
  double gv[3];

  for(i=0;i<3;i++)
    gv[i]=0.;
  for(i=0;i<NUMBER;i++){
    gv[0]=gv[0]+bxv[i];
    gv[1]=gv[1]+byv[i];
    gv[2]=gv[2]+bzv[i];
  }

  for(i=0;i<3;i++)
    gv[i]=gv[i]/NUMBER;
  for(i=0;i<NUMBER;i++){
    bxv[i]= bxv[i]-gv[0];
    byv[i]= byv[i]-gv[1];
    bzv[i]= bzv[i]-gv[2];
  }

}

double calc_cm(bx, by,bz)
     double  bx[], by[], bz[];
{
	/* calculate the center-of-mass */
  int i;
  double g[3];

  for(i=0;i<3;i++)
    g[i]=0.;
  for(i=0;i<NUMBER;i++){
    g[0]=g[0]+bx[i];
    g[1]=g[1]+by[i];
    g[2]=g[2]+bz[i];
  }

  for(i=0;i<3;i++)
    g[i]=g[i]/NUMBER;

	return g[0], g[1], g[2];

}


void write_cood(bx,by,bz,bxv,byv,bzv,iface,ibon_all,nbon_all,pf,runnumber)
     double  bx[], by[], bz[];
     double  bxv[], byv[], bzv[];
	double pf;
     int iface[NF][3],ibon_all[NUMBER][MAXBOND],nbon_all[];
	int runnumber;
{
  int i,j,k;
  FILE *finp;

	char buf_cood[0x100];
	snprintf(buf_cood, sizeof(buf_cood), "Data/cood_pf%.0f_%d.dat", pf,runnumber);
  finp=fopen(buf_cood,"w");

  for(i=0;i<NUMBER;i++)
    fprintf(finp,"%14.9lg %14.9lg %14.9lg\n",bx[i],by[i],bz[i]);

  for(i=0;i<NUMBER-1;i++){
    if(nbon_all[i]>0){
      for(j=0;j<nbon_all[i];j++){
	if(ibon_all[i][j]>i){
	  fprintf(finp,"%d",i);
	  for(k=0;k<nbon_all[i];k++){
	    if(ibon_all[i][k]>i)
	      fprintf(finp,"  %d",ibon_all[i][k]);
	  }
	  fprintf(finp,"\n");
	  break;
	}
      }
    }
  }
  fprintf(finp,"%d\n",NUMBER-1);

  for(i=0;i<NF;i++)
    fprintf(finp,"%d %d %d\n", iface[i][0],iface[i][1],iface[i][2]);
  for(i=0;i<NUMBER;i++)
    fprintf(finp,"%14.9lg %14.9lg %14.9lg\n",bxv[i],byv[i],bzv[i]);
  fclose(finp);
}

// slightly changed to work for my drawing program
/*
void record_draw(bx,by,bz,ibon_all,nbon_all,istep,dp)
     double  bx[],  by[],  bz[];
     int  ibon_all[NUMBER][MAXBOND],nbon_all[],istep;
     FILE *dp;
{
  int i,j,k;

  fprintf(dp,"%lf\n", istep*DT);

  for(i=0;i<NUMBER;i++)
    fprintf(dp,"%lg %lg %lg\n",bx[i],by[i],bz[i]);
  for(i=0;i<NUMBER-1;i++){
    if(nbon_all[i]>0){
      for(j=0;j<nbon_all[i];j++){
	if(ibon_all[i][j]>i){
	  fprintf(dp,"%d",i);
	  for(k=0;k<nbon_all[i];k++){
	    if(ibon_all[i][k]>i)
	      fprintf(dp,"  %d",ibon_all[i][k]);
	  }
	  fprintf(dp,"\n");
	  break;
	}
      }
    }
  }
  fprintf(dp,"%d\n\n",NUMBER-1);
}
*/

void record_draw(bx,by,bz,iface,istep,dp)
     double  bx[],  by[],  bz[];
     int  iface[NF][3], istep;
     FILE *dp;
{
  int i;

  fprintf(dp,"%lf\n", istep*DT);

  for(i=0;i<NUMBER;i++)
    fprintf(dp,"%lg %lg %lg\n",bx[i],by[i],bz[i]);
  for(i=0;i<NF;i++)
    fprintf(dp,"%d %d %d\n", iface[i][0],iface[i][1],iface[i][2]);

}



void record_drawVMD(bx,by,bz,istep,dpVMD)
     double  bx[],  by[],  bz[];
     int  istep;
     FILE *dpVMD;
{
  /*create output file in .xyz format to display evolution with VMD*/
  int i;

  fprintf(dpVMD,"%d\n", NUMBER);
  fprintf(dpVMD,"%lf\n", istep*DT);

  for(i=0;i<NUMBER;i++){
    fprintf(dpVMD,"H %lg %lg %lg\n",bx[i],by[i],bz[i]);
  }

}


void energy(bx,by,bz,bxv,byv,bzv,rni,gni,iface,iedge,en,ar,rbonl)
     double bx[],by[],bz[],bxv[],byv[],bzv[],rni[NUMBER][3],gni[NUMBER][3];
     int iface[NF][3],iedge[NE][2];
     double en[],ar[],rbonl[];
{
  int i,ib,jb,kb;
  double xlj,ylj,zlj;
  double xlk,ylk,zlk;
  double gn[3];
  double r0,r2,a0,h0;
  double rij[4],ar0,vo0,sar;

  rbonl[0]= 0.;
  rbonl[1]= 0.;
  en[6]=0.;
  vo0=0.;
  ar0=0.;
  for(i=0;i<NUMBER;i++){
    gni[i][0]=0.;
    gni[i][1]=0.;
    gni[i][2]=0.;
  }
  for(i=0;i<NE;i++){
    ib=iedge[i][0];
    jb=iedge[i][1];
    rij[1]= bx[ib]-bx[jb];
    rij[2]= by[ib]-by[jb];
    rij[3]= bz[ib]-bz[jb];
    r0= sqrt(rij[1]*rij[1]+rij[2]*rij[2]+rij[3]*rij[3]);
    rbonl[0]=rbonl[0]+r0;
    rbonl[1]=rbonl[1]+r0*r0;
    if(rbonl[2]<r0)
      rbonl[2]=r0;
    if(rbonl[3]>r0)
      rbonl[3]=r0;
  }
  for(i=0;i<NF;i++){
    ib=iface[i][0];
    jb=iface[i][1];
    kb=iface[i][2];

    xlj=bx[ib]-bx[jb];
    ylj=by[ib]-by[jb];
    zlj=bz[ib]-bz[jb];

    xlk=bx[ib]-bx[kb];
    ylk=by[ib]-by[kb];
    zlk=bz[ib]-bz[kb];

    gn[0]= ylj*zlk - zlj*ylk;
    gn[1]= zlj*xlk - xlj*zlk;
    gn[2]= xlj*ylk - ylj*xlk;

    gni[ib][0]=gni[ib][0] + gn[0];
    gni[ib][1]=gni[ib][1] + gn[1];
    gni[ib][2]=gni[ib][2] + gn[2];
    gni[jb][0]=gni[jb][0] + gn[0];
    gni[jb][1]=gni[jb][1] + gn[1];
    gni[jb][2]=gni[jb][2] + gn[2];
    gni[kb][0]=gni[kb][0] + gn[0];
    gni[kb][1]=gni[kb][1] + gn[1];
    gni[kb][2]=gni[kb][2] + gn[2];

    a0= sqrt(gn[0]*gn[0]+gn[1]*gn[1]+gn[2]*gn[2]);
    h0= gn[0]*bx[ib]+gn[1]*by[ib]+gn[2]*bz[ib];
    vo0=vo0+h0;
		if(h0 < -100){
		printf("test vol %lg face %d xyz %d %d %d\n", h0, i, ib, jb, kb);}
    ar0=ar0+a0;
  }
  rbonl[0]=rbonl[0]/NE;
  rbonl[1]=rbonl[1]/NE;
  vo0=vo0/6.;
  ar0=ar0*0.5;

  ar[2]=0.;
  for(i=0;i<NUMBER;i++){
    r0= rni[i][0]*rni[i][0]+rni[i][1]*rni[i][1]+rni[i][2]*rni[i][2];
    sar = rni[i][0]*gni[i][0]+rni[i][1]*gni[i][1]+rni[i][2]*gni[i][2];
    sar= sar/fabs(sar);
    ar[2]=ar[2]+ sar*sqrt(r0);
  }

  for(i=0;i<NUMBER;i++){
    r2=bxv[i]*bxv[i]+byv[i]*byv[i]+bzv[i]*bzv[i];
    en[6]=en[6] + r2;
  }

  en[6]=en[6]*0.5*BMASSB;

  en[0]=en[1]+en[2]+en[3]+en[4]+en[5];

  if( fabs(ar[1]-vo0)>0.0001)
    printf("vol %lg     %lg %lg\n",ar[1]-vo0,ar[1],vo0);
  if( fabs(ar[0]-ar0)>0.0001)
    printf("area %lg     %lg %lg\n",ar[0]-ar0,ar[0],ar0);

  ar[3]=ar[2]*0.5/sqrt(ar[0]);
  ar[4]= sqrt(ar[0]*0.25/PI);
  ar[4]= ar[1]*3./(4.*PI*ar[4]*ar[4]*ar[4]);
}

void cal_rg(bx,by,bz,rg,gg,ar)
     double  bx[],  by[],  bz[];
     double rg[],gg[],ar[];
{
  int i;
  double r0,Txx[6],rg4,xl,yl,zl,deta;

  for(i=0;i<3;i++)
    gg[i]=0.;
  for(i=0;i<6;i++)
    Txx[i]=0.;
  for(i=0;i<NUMBER;i++){
    gg[0]=gg[0]+bx[i];
    gg[1]=gg[1]+by[i];
    gg[2]=gg[2]+bz[i];
  }
  for(i=0;i<3;i++)
    gg[i]=gg[i]/NUMBER;

  for(i=0;i<NUMBER;i++){
    xl=bx[i]-gg[0];
    yl=by[i]-gg[1];
    zl=bz[i]-gg[2];

    Txx[0]=Txx[0] + xl*xl;
    Txx[1]=Txx[1] + yl*yl;
    Txx[2]=Txx[2] + zl*zl;
    Txx[3]=Txx[3] + xl*yl;
    Txx[4]=Txx[4] + yl*zl;
    Txx[5]=Txx[5] + zl*xl;
  }
  for(i=0;i<6;i++)
    Txx[i]=Txx[i]/NUMBER;
  rg[1]= Txx[0]+Txx[1]+Txx[2];
  rg[0]= sqrt(rg[1]);
  r0=Txx[0]*Txx[1]-Txx[3]*Txx[3]+Txx[1]*Txx[2]-Txx[4]*Txx[4]+Txx[0]*Txx[2]-Txx[5]*Txx[5];
  rg4=rg[1]*rg[1];
  rg[2]=(rg4-3.*r0)/rg4;

  deta= 2.*Txx[3]*Txx[4]*Txx[5]- Txx[0]*Txx[4]*Txx[4]- Txx[1]*Txx[5]*Txx[5]- Txx[2]*Txx[3]*Txx[3]  + Txx[0]*Txx[1]*Txx[2];
  rg[3]= 9.*deta/(rg[1]*r0);

  r0= sqrt(ar[0]*0.25/PI);
  rg[4]= rg[0]/r0;
}

void check_V(bx,by,bz,gg,ar,iface,iedge,mp)
     double bx[],by[],bz[],gg[],ar[];
     int iface[NF][3],iedge[NE][2];
     FILE *mp;
{
  int i,ib,jb,kb;
  double volc[3],gn[3];
  double xlj,ylj,zlj;
  double xlk,ylk,zlk;
  
  for(i=0;i<3;i++)
    volc[i]=0.;
  for(i=0;i<NF;i++){
    ib=iface[i][0];
    jb=iface[i][1];
    kb=iface[i][2];

    xlj=bx[ib]-bx[jb];
    ylj=by[ib]-by[jb];
    zlj=bz[ib]-bz[jb];

    xlk=bx[ib]-bx[kb];
    ylk=by[ib]-by[kb];
    zlk=bz[ib]-bz[kb];

    gn[0]= ylj*zlk - zlj*ylk;
    gn[1]= zlj*xlk - xlj*zlk;
    gn[2]= xlj*ylk - ylj*xlk;

    volc[0]=volc[0]+ gn[0]*(bx[ib]+bx[jb]+bx[kb]);
    volc[1]=volc[1]+ gn[1]*(by[ib]+by[jb]+by[kb]);
    volc[2]=volc[2]+ gn[2]*(bz[ib]+bz[jb]+bz[kb]);
  }
  for(i=0;i<3;i++)
    volc[i]=volc[i]/6.;
  for(i=0;i<3;i++){
    if( ((volc[i]-ar[1])>0.01)||((volc[i]-ar[1])<-0.01) )
      fprintf(mp,"error volume %lg, %lg %lg %lg\n",ar[1],volc[0],volc[1],volc[2]);
  }

}

void mean(en,en_m,en_m2,rg,rgm,ar,arm,rbonl)
     double en[],en_m[],en_m2[],rg[],rgm[],ar[],arm[],rbonl[];
{
  int i;
    
  for(i=0;i<9;i++){
    en_m[i]= en_m[i]+en[i];
    en_m2[i]= en_m2[i]+en[i]*en[i];
  }
  for(i=0;i<5;i++){
    rgm[2*i]=rgm[2*i] + rg[i];
    rgm[2*i+1]=rgm[2*i+1] + rg[i]*rg[i];
  }
  rgm[10]=rgm[10]+ rbonl[0];
  rgm[11]=rgm[11]+ rbonl[1];
  for(i=0;i<5;i++){
    arm[2*i]=arm[2*i] + ar[i];
    arm[2*i+1]=arm[2*i+1] + ar[i]*ar[i];
  }
}
void distmov(bx,by,bz,bxq,byq,bzq,bqmax)
     double  bx[],  by[],  bz[];
     double  bxq[], byq[], bzq[],bqmax[];
{
  int i;
  double r0,xl,yl,zl;

  for(i=0;i<NUMBER;i++){
    xl=bx[i]-bxq[i];
    yl=by[i]-byq[i];
    zl=bz[i]-bzq[i];
    r0= xl*xl+yl*yl+zl*zl;
    if( (r0>bqmax[0]) ) bqmax[0]=r0;
    bxq[i]=bx[i];
    byq[i]=by[i];
    bzq[i]=bz[i];
  }
}

void check_bond(ibon_all,nbon_all,mp)
     int ibon_all[NUMBER][MAXBOND],nbon_all[];
     FILE *mp;
{
	/* check that bond sorting was performed correctly */
  int i,jb,k,kk,kkk,k1,ke0,ke1,k0;
  ke0=0;
  ke1=0;
  for(i=0;i<NUMBER;i++){
    for(k=0;k<nbon_all[i];k++){
      jb=ibon_all[i][k];
      k0=0;
      for(kk=k+1;kk<nbon_all[i];kk++){
	k1=ibon_all[i][kk];
	if(jb==k1 )
	  ke0++;
      }
      for(kk=0;kk<nbon_all[i];kk++){
	k1=ibon_all[i][kk];
	for(kkk=0;kkk<nbon_all[jb];kkk++){
	  if(k1==ibon_all[jb][kkk]){
	    k0++;
	  }
	}
      }
      if( k0 !=2 )
	ke1++;
    }
  }
  if(ke0+ke1>0){
    fprintf(mp,"error ck_bond %d %d\n",ke0,ke1);
    exit(0);
  }
}

void sort_bd(bx, by, bz,bxv, byv, bzv,bxq,byq,bzq,iface,iedge,ibon_all,nbon_all,iedb,iedf,ifaed,mp,index_pull)
     double  bx[], by[], bz[];
     double  bxv[], byv[], bzv[];
     double  bxq[], byq[], bzq[];
     int iface[NF][3],iedge[NE][2];
     int ibon_all[NUMBER][MAXBOND],nbon_all[];
     int iedb[NE][2],iedf[NE][2],ifaed[NF][3];
		 int index_pull[];
     FILE *mp;
{
	/* arrange edges for each vertex */
  int i,j,ib,jb,jc,k,kk,kkk,kb,k0,jbf[3];
  int ib_ck[NUMBER],ib_tmp[NF],ntmp,iface_tmp[NF][3];
  int ibon_tmp[NUMBER][MAXBOND],nbon_tmp[NUMBER];
  double b_tmp[NUMBER];
	int found1 = 1;
	int found2 = 1;

  ntmp=1;
  ib_ck[0]=0;
  ib_tmp[0]=0;
  for(i=1;i<NUMBER;i++){
    ib_ck[i]=-1;
    ib_tmp[i]=NUMBER;
  }
  for(i=0;i<NUMBER;i++){
    ib=ib_tmp[i];
		if(ib == index_pull[0] && found1 == 1){  // the index of the two chosen beads changes!
			index_pull[0] = i;
			found1 = 0;
			//printf("found 1 i ib %d %d \n",i, ib);
			
		}
		if(ib == index_pull[1] && found2 == 1){
			index_pull[1] = i;
			found2 = 0;
			//printf("found 2 i ib %d %d \n",i, ib);
		}
    for(j=0;j<nbon_all[ib];j++){
      if(ib_ck[ibon_all[ib][j]]==-1){
	ib_tmp[ntmp]=ibon_all[ib][j];
	ib_ck[ibon_all[ib][j]]=ntmp;
	ntmp++;
      }
    }
  }

  for(i=0;i<NUMBER;i++){
    if(ib_ck[i]==-1){
      fprintf(mp,"error sort ntmp %d  %d\n",ntmp,i);
      ib_tmp[ntmp]=i;
      ib_ck[i]=ntmp;
      ntmp++;
    }
  }
  if(ntmp !=NUMBER)
    fprintf(mp,"error sort num ntmp %d\n",ntmp);

  for(i=0;i<NUMBER;i++){
    ib=ib_tmp[i];
    nbon_tmp[i]=0;
    for(j=0;j<nbon_all[ib];j++){
      jc=ib_ck[ibon_all[ib][j]];
      if(jc>i){
	k0=0;
	for(k=nbon_tmp[i]-1;k>-1;k--){
	  if(jc>ibon_tmp[i][k]){
	    k0=k+1;
	    break;
	  }
	}
	for(k=nbon_tmp[i];k>k0;k--)
	  ibon_tmp[i][k]=ibon_tmp[i][k-1];
	ibon_tmp[i][k0]=jc;
	nbon_tmp[i]++;
      }
    }
  }
  
  ntmp=0;
  for(i=0;i<NUMBER;i++)
    nbon_all[i]=0;
  for(i=0;i<NUMBER;i++){
    for(j=0;j<nbon_tmp[i];j++){
      iedge[ntmp][0]=i;
      iedge[ntmp][1]=ibon_tmp[i][j];
      ntmp++;
      ibon_all[i][nbon_all[i]]=ibon_tmp[i][j];
      nbon_all[i]++;
      ibon_all[ibon_tmp[i][j]][nbon_all[ibon_tmp[i][j]] ]=i;
      nbon_all[ibon_tmp[i][j]]++;
    }
  }
  for(i=0;i<NUMBER;i++)
    b_tmp[i]=bx[i];
  for(i=0;i<NUMBER;i++)
    bx[i]=b_tmp[ib_tmp[i]];
  for(i=0;i<NUMBER;i++)
    b_tmp[i]=by[i];
  for(i=0;i<NUMBER;i++)
    by[i]=b_tmp[ib_tmp[i]];
  for(i=0;i<NUMBER;i++){
    b_tmp[i]=bz[i];
	}
  for(i=0;i<NUMBER;i++){
    bz[i]=b_tmp[ib_tmp[i]];
	}
  for(i=0;i<NUMBER;i++)
    b_tmp[i]=bxv[i];
  for(i=0;i<NUMBER;i++)
    bxv[i]=b_tmp[ib_tmp[i]];
  for(i=0;i<NUMBER;i++)
    b_tmp[i]=byv[i];
  for(i=0;i<NUMBER;i++)
    byv[i]=b_tmp[ib_tmp[i]];
  for(i=0;i<NUMBER;i++)
    b_tmp[i]=bzv[i];
  for(i=0;i<NUMBER;i++)
    bzv[i]=b_tmp[ib_tmp[i]];
  for(i=0;i<NUMBER;i++){
    bxq[i]=bx[i];
    byq[i]=by[i];
    bzq[i]=bz[i];
  }

  /*****   sort  face    ******/

  for(i=0;i<NF;i++){
    jbf[0]=ib_ck[iface[i][0]];
    jbf[1]=ib_ck[iface[i][1]];
    jbf[2]=ib_ck[iface[i][2]];
    k0=0;
    for(k=1;k<3;k++){
      if(jbf[k0]>jbf[k])
	k0=k;
    }
    iface_tmp[i][0]= jbf[k0];
    iface_tmp[i][1]= jbf[k0+1-((k0+1)/3)*3  ];
    iface_tmp[i][2]= jbf[k0+2-((k0+2)/3)*3  ];
  }

  ntmp=0;
  for(i=0;i<NUMBER;i++){
    for(j=0;j<NF;j++){
      if(iface_tmp[j][0]==i){
	ib_tmp[ntmp]=j;
	ntmp++;
      }
    }
  }
  if(ntmp !=NF)
    fprintf(mp,"error sort nf ntmp %d\n",ntmp);

  for(i=0;i<NF;i++){
    iface[i][0]= iface_tmp[ib_tmp[i]][0];
    iface[i][1]= iface_tmp[ib_tmp[i]][1];
    iface[i][2]= iface_tmp[ib_tmp[i]][2];
  }
  
  /*****   set others    ******/
  
  for(i=0;i<NF;i++)
    ib_tmp[i]=0;
  for(i=0;i<NE;i++){
    ib=iedge[i][0];
    jb=iedge[i][1];
    ntmp=0;
    for(j=0;j<NF;j++){
      for(k=0;k<3;k++){
	if(ib==iface[j][k]){
	  for(kk=0;kk<3;kk++){
	    if(jb==iface[j][kk]){
	      for(kkk=0;kkk<3;kkk++){
		if((k !=kkk)&&(kk !=kkk) ){
		  kb=iface[j][kkk];
		  if( (ntmp==1)&&(iedb[i][0]>kb ) ){
		    iedf[i][1]=iedf[i][0];
		    iedb[i][1]=iedb[i][0];
		    iedf[i][0]=j;
		    iedb[i][0]=kb;
		  }else{
		    iedf[i][ntmp]=j;
		    iedb[i][ntmp]=kb;
		  }
		  ntmp++;

		  ifaed[j][ib_tmp[j]]=i;
		  ib_tmp[j]++;
		}
	      }
	    }
	  }
	}
      }
    }
    if(ntmp !=2){
      fprintf(mp,"error %d ned %d\n",i,ntmp);
      exit(0);
    }
  }

}

void init_para(en_m,en_m2,rgm,arm)
     double en_m[],en_m2[],rgm[],arm[];
{
  int i;
    
  for(i=0;i<9;i++){
    en_m[i]=0.0;
    en_m2[i]=0.0;
  }
  for(i=0;i<12;i++)
    rgm[i]=0.0;
  for(i=0;i<10;i++)
    arm[i]=0.0;

}


double calc_max_dist(bx, by, bz)
	double  bx[], by[], bz[];
{
	double maxdist = 0.0;
	double dist2;
	int i, j;
	for(i=0; i<NUMBER; i++){
		for(j=0; j<=i; j++){
			dist2 = (bx[i] - bx[j])*(bx[i] - bx[j]) + (by[i] - by[j])*(by[i] - by[j]) + (bz[i] - bz[j])*(bz[i] - bz[j]);
			if(dist2 > maxdist){
				maxdist = dist2;
			}
		}
	}
	maxdist = sqrt(maxdist);

	return maxdist;
}


void find_points(index_pull, bx, by, bz, pull,pf)
	int index_pull[2];
	double  bx[], by[], bz[];
	double pull[2][3], pf;
{
	/* finde the two points which are furtest appart, i.e. the ones opposing each other most */

	double maxdist = 0.0;
	double dist2, norm1, norm2, norm;
	int i, j;
	for(i=0; i<NUMBER; i++){
		for(j=0; j<=i; j++){
			dist2 = (bx[i] - bx[j])*(bx[i] - bx[j]) + (by[i] - by[j])*(by[i] - by[j]) + (bz[i] - bz[j])*(bz[i] - bz[j]);
			if(dist2 > maxdist){
				maxdist = dist2;
				index_pull[0] = i;
				index_pull[1] = j;
			}
		}
	}
	maxdist = sqrt(maxdist);

	/* apply force in opposing directions */
	norm = 1.0/sqrt(bx[index_pull[0]]*bx[index_pull[0]] + by[index_pull[0]]*by[index_pull[0]] + bz[index_pull[0]]*bz[index_pull[0]] );
	pull[0][0] = pf*norm*bx[index_pull[0]];
	pull[1][0] = -pf*norm*bx[index_pull[0]];
	pull[0][1] = pf*norm*by[index_pull[0]];
	pull[1][1] = -pf*norm*by[index_pull[0]];
	pull[0][2] = pf*norm*bz[index_pull[0]];
	pull[1][2] = -pf*norm*bz[index_pull[0]];

}


void ext_const_force(extforce, index_pull, pull, istep)
	int index_pull[2];
	int istep;
	double extforce[NUMBER][3], pull[2][3];
{
	/* apply constant external force in two oppsing directions
	   - possible to apply force only for a fixed period of time */
	int i, j;
	for(i=0; i<NUMBER; i++){
		for(j=0; j<3; j++){
			extforce[i][j] = 0.0;
		}
	 }
	if(istep < 2000000){
		extforce[index_pull[0]][0] = pull[0][0];
		extforce[index_pull[0]][1] = pull[0][1];
		extforce[index_pull[0]][2] = pull[0][2];

		extforce[index_pull[1]][0] = pull[1][0];
		extforce[index_pull[1]][1] = pull[1][1];
		extforce[index_pull[1]][2] = pull[1][2];
	}

}

void ext_const_force_single(extforce, index_pull, pull, istep)
	int index_pull[2];
	int istep;
	double extforce[NUMBER][3], pull[2][3];
{
	/* apply constant external force in two oppsing directions
	   - possible to apply force only for a fixed period of time */
	int i, j;
	for(i=0; i<NUMBER; i++){
		for(j=0; j<3; j++){
			extforce[i][j] = 0.0;
		}
	 }
	if(istep < 2000000){
		extforce[index_pull[0]][0] = pull[0][0];
		extforce[index_pull[0]][1] = pull[0][1];
		extforce[index_pull[0]][2] = pull[0][2];

		extforce[index_pull[1]][0] = 0.0;
		extforce[index_pull[1]][1] = 0.0;
		extforce[index_pull[1]][2] = 0.0;
	}

}

void wnoise(idum,wn_f,dif)
     long *idum;
     double wn_f[],dif[];
{
	/* generating white noise */
  double ran2();
  double fac,rsq,v1,v2;
  int i;

  for(i=0;i<(( (NUMBER)*3+1 )/2);i++){
    do {
      v1=2.0*ran2(idum)-1.0;
      v2=2.0*ran2(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    wn_f[2*i]  = v2*fac;
    wn_f[2*i+1]= v1*fac;
  }

  for(i=0;i<(NUMBER*3);i++)
    wn_f[i]  = wn_f[i]*dif[0];
}

double gasdev(idum)
     long *idum;
{
  double ran2();
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;
    
  if  (iset == 0) {
    do {
      v1=2.0*ran2(idum)-1.0;
      v2=2.0*ran2(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(idum)
     long *idum;
{
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  double temp;

  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/* (C) Copr. 1986-92 Numerical Recipes Software 91o)+!|Y". */
