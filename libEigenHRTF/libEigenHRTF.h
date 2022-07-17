// Dear programmer:
// When I wrote this code, only god and
// I knew how it worked.
// Now, only god knows it!
// 
// Therefore, if you are trying to optimize
// this routine and it fails (most surely),
// please increase this counter as a
// warning for the next person:
// 
// total hours wasted here = 13.4


//AUTHOR: Camilo Arevalo - jarevalo29@javerianacali.edu.co
//DATE: 2021-06-24
//LICENSE
//Copyright (c) 2021 Juan Camilo Arvalo Arboleda
//All rights reserved.
//This code was created ONLY for noncommercial research or personal use purpose. 
//Any use, modification, or distribution of this software outside of the previously mentioned 
//purpose or not granted in writing by the author is NOT PERMITTED.
//The use of this software is subject to the following conditions:
//1. The above copyright notice and this permission notice shall be included in all copies or 
//substantial portions of the Software.
//2. You must give appropriate credit to the author of this software.
//3. Derivation of this code cannot be re-licensed under another license without the 
//authorization in writing by the author.
//The software is provided without warranty of any kind. In no event shall 
//the authors or copyright holders be liable for any claim, damages, or other 
//liability, whether in an action of contract, tort or otherwise, arising from, 
//out of or in connection with the software or the use or other dealings 
//in the software.

#include <complex.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#if defined(_WIN64)
#include <windows.h>
#define OPEN_FLAG "rb"
#else
#define OPEN_FLAG "r"
#endif

#ifndef M_PI
#define M_PI (float)3.14159265359
#endif
#define N_DISTANCE 8

//Structure that will hold the offline barycentric interpolation
struct interpolation{
  int16_t idx[1212680][4]; //index of the points in the mesh
  float ws[1212680][4]; //weights of the interpolation
};

//structure that holds the eigen decomposition
struct Eigen{ 
    int sr; //sample rate of the database
    int n_bins; //size of the hrtf
    int16_t n_s[6344]; //number of coeficients for small pinna
    int16_t n_l[6344]; //number of coeficients for large pinna
    int16_t s_delays[6344]; //delays in samples for the small pinna
    int16_t l_delays[6344]; //delays in samples for the large pinna
    float complex *s_mean; //HRTF mean for the small pinna
    float complex **s_V; //Eigen vectors for the small pinna
    float** s_coe; //Coeficientsfor the small pinna
    float complex *l_mean; //HRTF mean for the large pinna
    float complex **l_V; //Eigen vectors for the large pinna
    float** l_coe; //Coeficientsfor the large pinna
};

void read_interpolation(struct interpolation *interpol,char* path ){
    FILE *fp;
	char name[250];
    sprintf(name, "%sinterpolation.dat", path);
    fp = fopen(name, OPEN_FLAG);	
    for (int i = 0; i < 1212680; i++)
	{
		fread(interpol->idx[i],sizeof(int16_t)*4,1,fp);
	}
	for (int i = 0; i < 1212680; i++)
	{
		fread(interpol->ws[i],sizeof(float)*4,1,fp);
	}
    fclose(fp);
}

//sr: sample rate
void read_eigen(struct Eigen *eig,int sr,char* path ){
    char name[512];
    sprintf(name, "%sEigen@%d.dat",path,sr);
    FILE *fp;    
    fp = fopen(name, OPEN_FLAG);
    int n_bins=257;
    if (eig->sr>48000)
    {
        n_bins=513;
    }
    fread(&eig->sr,sizeof(int),1,fp);    
    fread(&eig->n_bins,sizeof(int),1,fp);    
    fread(eig->n_s,sizeof(int16_t)*6344,1,fp);
    fread(eig->n_l,sizeof(int16_t)*6344,1,fp);
    fread(eig->s_delays,sizeof(int16_t)*6344,1,fp);    
    fread(eig->l_delays,sizeof(int16_t)*6344,1,fp);    
    eig->s_mean=(float complex*)malloc(n_bins*sizeof(float complex));
    fread(eig->s_mean,n_bins*sizeof(float complex),1,fp);    
    int max_s=0;
    for(int i=0;i<6344; i++) if(eig->n_s[i]>max_s) max_s=eig->n_s[i];        
    eig->s_V=(float complex **)malloc(n_bins*sizeof(float complex *));        
    for (int i = 0; i < n_bins; i++)
    {
        eig->s_V[i]=(float complex *)malloc(max_s*sizeof(float complex));
        fread(eig->s_V[i],max_s*sizeof(float complex),1,fp);        
    }         
    eig->s_coe=(float**)malloc(6344*sizeof(float*));    
    for (int i = 0; i < 6344; i++)
    {        
        eig->s_coe[i]=(float*)malloc(eig->n_s[i]*sizeof(float));
        fread(eig->s_coe[i],eig->n_s[i]*sizeof(float),1,fp);
    }    
    /////
    eig->l_mean=(float complex*)malloc(n_bins*sizeof(float complex));
    fread(eig->l_mean,n_bins*sizeof(float complex),1,fp);      
    int max_l=0;
    for(int i=0;i<6344; i++) if(eig->n_l[i]>max_l) max_l=eig->n_l[i];
    eig->l_V=(float complex**)malloc(n_bins*sizeof(float complex*));
    for (int i = 0; i < n_bins; i++)
    {
        eig->l_V[i]=(float complex*)malloc(max_l*sizeof(float complex));
        fread(eig->l_V[i],max_l*sizeof(float complex),1,fp);    
    }
    eig->l_coe=(float**)malloc(6344*sizeof(float*));    
    for (int i = 0; i < 6344; i++)
    {        
        eig->l_coe[i]=(float*)malloc(eig->n_l[i]*sizeof(float));        
        fread(eig->l_coe[i],eig->n_l[i]*sizeof(float),1,fp);
    }
    fclose(fp);
}

float deg2rad(float deg) {
    return deg * M_PI / 180.0;
}

//return the index of the most near point in the offline interpolation database 
//for the given coodinates d:distances e:elevation and a:azimuth
int get_index(float d,float e,float a){
	int tmp_e=(int)roundf(e/2)+20;
	int n_e=0;
	int n_d=0;	
	int i=0;
    int distances[]={20,30,40,50,75,100,130,160};
	for (i = 0; i < 7; ++i)
	{
		if (distances[i]<=d && distances[i+1]>=d)
		{			
			n_d=i*170800+(int)roundf((d-distances[i])/((distances[i+1]-distances[i])/10))*17080;
		}
	}
	for (i = 0; i < tmp_e; ++i)
	{
		n_e+=(int)roundf(cosf(deg2rad(i*2-40))*360);
	}	
	return n_d+n_e+(int)roundf(a/(360/roundf(cosf(deg2rad(i*2-40))*360)));
}

//return the index in the original database for the given coordinates
//d:distances e:elevation and a:azimuth.
//if the coordinate asked is not in the original database, 
//it will return -1
int check_coor(float d, float e, float a){
    int distances[]={20,30,40,50,75,100,130,160};
	int n_d=-1;		
	for (int i = 0; i < N_DISTANCE; ++i)
	{
		if (distances[i]==(int)d)
		{
			n_d=i*793;
			break;
		}
	}
	if (n_d<0) return -1;
	if (fmodf(e,10)!=0) return -1;
	if (-40<=e && e<=50)
	{
		if (fmodf(a,5)!=0) return -1;
		return n_d+(((int)e/10)+4)*72+((int)a/5);	
	}
	if (e==60)
	{
		if (fmodf(a,10)!=0) return -1;
		return n_d+720+((int)a/10);	
	}
	if (e==70)
	{
		if (fmodf(a,15)!=0) return -1;
		return n_d+756+((int)a/15);
	}
	if (e==80)
	{
		if (fmodf(a,30)!=0) return -1;
		return n_d+780+((int)a/30);	
	}
	if(e==90) return n_d+792;	
	return -1;
}


/*
Reconstrution of hrtf is done as following: 

lets assume we want to use the small pinna.
                                   | dot product
								   v 
hrtf=s_mean[1..n_bins]+s_coe[idx][1..n_s[idx]]*s_V[1..n_s][1..n_bins]

Interpolation is done as following:

hrtf=s_mean[1..n_bins]+
(sum i from 0 to 4 
	s_coe[inter.idx[idx][i]][1..n_s[inter.idx[idx][i]]]*inter.ws[idx][i])
	*  s_V[1..max(n_s[])][1..n_bins]
	^
	| dot product

*/

//reconstruct from the original database eig the HRTF in the
//idx index for the pinna (0: small, 1: large) and write the
//reconstruction in the out pointer.
void reconstruction(struct Eigen *eig, int idx,int pinna,float complex* out){
	float complex hrtf[eig->n_bins];
	if (pinna)
	{
		memcpy(hrtf,eig->l_mean,sizeof(float complex)*eig->n_bins);
		for (int i = 0; i < eig->n_bins; ++i)
		{
			for (int j = 0; j < eig->n_l[idx]; ++j)			
			{
				hrtf[i]+=eig->l_coe[idx][j]*eig->l_V[i][j];
			}			
		}
	}else{
		memcpy(hrtf,eig->s_mean,sizeof(float complex)*eig->n_bins);
		for (int i = 0; i < eig->n_bins; ++i)
		{
			for (int j = 0; j < eig->n_s[idx]; ++j)			
			{
				hrtf[i]+=eig->s_coe[idx][j]*eig->s_V[i][j];
			}			
		}
	}
	memcpy(out,hrtf,sizeof(float complex)*eig->n_bins);
}

//interpolate and reconstruct from the offline interpolation database inter 
//the HRTF in the idx index (interpolation database) for the pinna (0: small, 1: large) and write the
//interpolated reconstruction in the out pointer.
void interpolation(struct Eigen *eig,struct interpolation *inter,int idx,int pinna,float complex* out){
	int i=0;
	int max=0;
	float complex hrtf[eig->n_bins];	
	if (pinna)
	{
		memcpy(hrtf,eig->l_mean,sizeof(float complex)*eig->n_bins);			
		for (i = 0; i < 4; ++i)
		{
			if (eig->n_l[inter->idx[idx][i]]>max) max=eig->n_l[inter->idx[idx][i]];						
		}		
		for (int j = 0; j < max; ++j)
		{
			float coe=0;
			for (int k = 0; k < 4; ++k)
			{
				if (eig->n_l[inter->idx[idx][k]]>j) coe+=eig->l_coe[inter->idx[idx][k]][j]*inter->ws[idx][k];								
			}			
			for (i = 0; i < eig->n_bins; ++i)
			{				
				hrtf[i]+=eig->l_V[i][j]*coe;
			}			
		}
	}else{
		memcpy(hrtf,eig->s_mean,sizeof(float complex)*eig->n_bins);
		for (i = 0; i < 4; ++i)
		{
			if (eig->n_s[inter->idx[idx][i]]>max) max=eig->n_s[inter->idx[idx][i]];
		}				
		for (int j = 0; j < max; ++j)
		{
			float coe=0;
			for (int k = 0; k < 4; ++k)
			{
				if (eig->n_s[inter->idx[idx][k]]>j) coe+=eig->s_coe[inter->idx[idx][k]][j]*inter->ws[idx][k];
			}
			for (i = 0; i < eig->n_bins; ++i)
			{				
				hrtf[i]+=eig->s_V[i][j]*coe;
			}			
		}
	}	
	memcpy(out,hrtf,sizeof(float complex)*eig->n_bins);	
}

//interpolate the delay from the offline interpolation database inter 
//in the idx index (interpolation database) for the pinna (0: small, 1: large) and return the
//interpolated delay in samples.
int inter_delay(struct Eigen *eig,struct interpolation *inter,int idx,int pinna){
	int i=0;	
	float delay=0;
	if (pinna)
	{
		for (i = 0; i < 4; ++i)
		{
			delay+=(float)eig->l_delays[inter->idx[idx][i]]*inter->ws[idx][i];
		}						
	}else{
		for (i = 0; i < 4; ++i)
		{
			delay+=(float)eig->s_delays[inter->idx[idx][i]]*inter->ws[idx][i];
		}
	}
	return (int)roundf(delay);
}

//get the delay in samples from the original database eig in the 
//index idx for the pinna (0: small, 1: large)
int get_delay(struct Eigen *eig,int idx,int pinna){
	if (pinna)
	{
		return eig->l_delays[idx];
	}else{
		return eig->s_delays[idx];
	}
}

/**
* For the brave souls who get this far: You are the chosen ones,
* the valiant knights of programming who toil away, without rest,
* fixing our most awful code. To you, true saviors, kings of men,
* I say this: never gonna give you up, never gonna let you down,
* never gonna run around and desert you. Never gonna make you cry,
* never gonna say goodbye. Never gonna tell a lie and hurt you.
*/

//get the HRTF by the reconstriction (with or without interpolation using the inter offline interpolated database)
//of the eig database for the given coordinates d:distances e:elevation and a:azimuth
//and for the pinna (0: small, 1: large).
//HRTF reconstruceted for the left pinna will be writen in the filter_l and the rigth pinna in filter_r.
//the delays will be computed as int[2] (0: left pinna, 1: rigth pinna),
//and index of HRTF used in the reconstrucction will be in the idxs (0:original database,1: interpolated database).
void get_filters(struct Eigen *eig, struct interpolation *inter, float d,float e,float a,int pinna, float complex* filter_l,float complex* filter_r,int* delays,int* idxs){
	if (e==90) d=0;
	if (a>=360) a=fmod(a,360);
	if (a<0) a=fabsf(a);
	if (e>90) e=90;
	if (e<-40) e=-40;
	if (d<20) d=20;
	if (d>160) d=160;
	if (!pinna)
	{
		a=fmodf(360-a,360);
	}
	int idx_l=check_coor(d,e,a);
	if (idxs[0]==idx_l)
	{
		return;
	}	
	if (idx_l<0)
	{		
		idx_l=get_index(d,e,a);
		if (idx_l==idxs[1])
		{
			return;
		}
		int idx_r=get_index(d,e,fmodf(360-a,360));
		idxs[1]=idx_l;
		idxs[0]=-2;
		interpolation(eig,inter,idx_l,pinna,filter_l);
		interpolation(eig,inter,idx_r,pinna,filter_r);		
		delays[0]=inter_delay(eig,inter,idx_l,pinna);		
		delays[1]=inter_delay(eig,inter,idx_r,pinna);		
	}
	else{			
		idxs[0]=idx_l;
		idxs[1]=-2;	
		int idx_r=check_coor(d,e,fmodf(360-a,360));
		reconstruction(eig,idx_l,pinna,filter_l);
		reconstruction(eig,idx_r,pinna,filter_r);
		delays[0]=get_delay(eig,idx_l,pinna);
		delays[1]=get_delay(eig,idx_r,pinna);
	}	
}