// pycpptraj note: this file was modified from ClusCo-0.2.1/src/gdt.cpp file
// OPENMP flags were turned off

#include <iostream>
#include <fstream>
#include <math.h>
#include <float.h>
using namespace std;
const double inv3 = (unsigned int)1/3.0;//0.3333333333333333333333333333333333333333333333;
const double root3 = 1.732050807568877193176604123436845839023590087890625;//1.732050807568877294;

inline double dist(float *atoms1, float *atoms2, int i) {
  double xs = atoms2[3*i]-atoms1[3*i];
  double ys = atoms2[3*i+1]-atoms1[3*i+1];
  double zs = atoms2[3*i+2]-atoms1[3*i+2];
  return sqrt(xs*xs + ys*ys + zs*zs);
//  return sqrt(pow(xs,2)+pow(ys,2)+pow(zs,2));
  
}


void rotate(float *atoms1, float *atoms2, bool *indexes, int protlen) {
  	 // calculate rotation matrix:
	// 1) covariance matrix
	unsigned int fraglen = 0;
		float covmat0,covmat1,covmat2,covmat3,covmat4,covmat5,covmat6,covmat7,covmat8; 
                covmat0=covmat1=covmat2=covmat3=covmat4=covmat5=covmat6=covmat7=covmat8=0;
                         
				 fraglen=0;
				 
				 for (int a=0;a<protlen;a++) {
		     if(indexes[a]==true) {
		      fraglen++;
		      
		      float a1x = atoms1[3*a];
		      float a1y = atoms1[3*a+1];
		      float a1z = atoms1[3*a+2];
		      float a2x = atoms2[3*a];
		      float a2y = atoms2[3*a+1];
		      float a2z = atoms2[3*a+2];
		      
		      covmat0 += (a2x)* (a1x);
		      covmat1 += (a2y)* (a1x);
		      covmat2 += (a2z)* (a1x);
		      
		      covmat3 += (a2x)* (a1y);
		      covmat4 += (a2y)* (a1y);
		      covmat5 += (a2z)* (a1y);
		      
		      covmat6 += (a2x)* (a1z);
		      covmat7 += (a2y)* (a1z);
		      covmat8 += (a2z)* (a1z);
		     }
		  }
		  
		  double div = 1.0/fraglen;
		  covmat0 *= div;
		  covmat1 *= div;
		  covmat2 *= div;
		  covmat3 *= div;
		  covmat4 *= div;
		  covmat5 *= div;
		  covmat6 *= div;
		  covmat7 *= div;
		  covmat8 *= div;
		// 2) square of covariance matrix, R2   
		  double r0 = covmat0*covmat0 +   covmat3*covmat3 + covmat6*covmat6;
		  double r1 = covmat0*covmat1 +    covmat3*covmat4 + covmat6*covmat7;
		  double r2 = covmat0*covmat2 +    covmat3*covmat5 + covmat6*covmat8;
		//  double r3 =r1;// covmat1*covmat0 +    covmat4*covmat3 + covmat7*covmat6;
		  double r4 = covmat1*covmat1 +    covmat4*covmat4 + covmat7*covmat7;
		  double r5 = covmat1*covmat2 +    covmat4*covmat5 + covmat7*covmat8;
		  double r8 = covmat2*covmat2 +    covmat5*covmat5 + covmat8*covmat8;
	 		//printf("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",r0,r1,r2,r4,r5,r8);

	// 3) eigenvalues
		    double root0,root1,root2;
		
		double r55 = r5*r5;
		double r22 = r2*r2;
		double r11 = r1*r1;
		    double c0 = r0*(r4*r8-r55) + 2.0*r1*r2*r5  - r4*r22 - r8*r11;
		    double c1 = r0*r4 - r11 + r0*r8 - r22 + r4*r8 - r55;
		    double c2 = r0+r4+r8;
		    double c2Div3 = c2*inv3;
		    double aDiv3 = (c1-c2*c2Div3)*inv3;
		    
		    if(aDiv3>0.0) 
		      aDiv3 = 0.0;
		    double mbDiv2 = 0.5*(c0+c2Div3*(2.0*c2Div3*c2Div3 - c1));
		    double q = mbDiv2*mbDiv2 + aDiv3*aDiv3*aDiv3;
		    if (q>0.0)
		      q=0.0;

		    double magnitude = sqrt(-aDiv3);
		    double angle = atan2(sqrt(-1*q),mbDiv2)*inv3;
		    double sn,cs;//=sin(angle); double cs=cos(angle);
		    sincos(angle,&sn,&cs);
		    double magnitudecs = magnitude*cs;
		    double magnituderoot = magnitude*root3*sn;

		    root0 = c2Div3 + 2.0*magnitudecs ;
		    root1 = c2Div3 - magnitudecs - magnituderoot;
		    root2 = c2Div3 - magnitudecs +magnituderoot;
	   
  if (root0<1e-5)
  root0 = 0;
  if (root1<1e-5)
  root1 = 0;
	if (root2<1e-5)
	 root2 = 0;
		    
		    
		    double minr,maxr,midr;
		    minr = min(root0,min(root1,root2));
		    maxr = max(root0,max(root1,root2));
		   midr = (maxr==root0)?(minr==root2 ? root1 : root2) : (minr==root2 ? root1 : root0);
		    root0=maxr; root1=midr; root2=minr;
		    
	// 4) eigenvectors
		       

		double ev0,ev1,ev2,ev3,ev4,ev5,ev6,ev7,ev8;
		double bf=r1*r5 - r2*r4;
		//first ev
		double clbf = r2*root0 + bf;
		double li = root0-r8;
		double iclbf = 1/clbf;
		ev0 = iclbf*(root0*(li-r4) + r4*r8 - r55 );
		ev1 = iclbf*(r1*li+r2*r5);
		ev2 = 1.0;
		 double len = 1.0/(double)sqrt(ev0*ev0 + ev1*ev1 + 1);
		ev0 *= len; ev1 *= len;ev2 = len;
	
		//second ev
		clbf = r2*root1 + bf;
		iclbf = 1.0/clbf;
		li = root1-r8;
		ev3 = (root1*(li-r4) + r4*r8 - r55 )*iclbf;
		ev4 = (r1*li+r2*r5)/clbf;
		ev5 = 1.0;
		len=1.0/sqrt(ev3*ev3 + ev4*ev4 + 1);
	ev3 *= len; ev4 *= len;ev5 = len;
	
			//third ev

		ev6 = r1*r5 - r2*(r4-root2);
		ev7 = r1*r2 - (r0-root2)*r5;
		ev8 = (r0-root2)*(r4-root2) - r1*r1;
		
		//clbf = r2*root2 + bf;
		//iclbf=(unsigned int)1/clbf;
		//li = root2-r8;
		//ev6 = (root2*(li-r4) + r4*r8 - r55 )*iclbf;
		//ev7 = (r1*li+r2*r5)*iclbf;
		//ev8 = 1.0;
	
	len = 1.0/sqrt(ev6*ev6 + ev7*ev7 + ev8*ev8);
		ev6 *= len;ev7 *= len;ev8 *= len;
	//	printf("%7.3f\n",ev2);

		if ((ev0*(ev4*ev8 - ev5*ev7) - ev1*(ev3*ev8-ev5*ev6) + ev2*(ev3*ev7-ev4*ev6))<0) {
			ev6*=-1.0;
			ev7*=-1.0;
			ev8*=-1.0;
		}
		
			// 5) rotation matrix
		// B = cov_matrix * eigenvectors
   		float B0,B1,B2,B3,B4,B5,B6,B7,B8;

		B0 =  covmat0*ev0 + covmat1*ev1 + covmat2*ev2;
		B1 =  covmat3*ev0 + covmat4*ev1 + covmat5*ev2;
		B2 =  covmat6*ev0 + covmat7*ev1 + covmat8*ev2;
		B3 =  covmat0*ev3 + covmat1*ev4 + covmat2*ev5;
		B4 =  covmat3*ev3 + covmat4*ev4 + covmat5*ev5;
		B5 =  covmat6*ev3 + covmat7*ev4 + covmat8*ev5;
		
	//	len = (unsigned int)1/
len = 1.0/sqrt(B0*B0 + B1*B1 + B2*B2);

		B0 *= len; 		B1 *= len; 		B2 *= len; // normalize
	//	len = (unsigned int)1/
len = 1.0/sqrt(B3*B3 + B4*B4 + B5*B5);

		B3 *= len;		B4 *= len;		B5 *= len;
		// third row of B
		B6 =  B1*B5 - B4*B2; 
		B7 = - B0*B5 + B3*B2;
		B8 =  B0*B4 - B3*B1;
	
// tttttttttttttttttttttttt
		// U = B^T*ev
		double U0,U1,U2,U3,U4,U5,U6,U7,U8;
		
		U0 = B0*ev0 + B3*ev3 + B6*ev6;
		U1 = B0*ev1 + B3*ev4 + B6*ev7;
		U2 = B0*ev2 + B3*ev5 + B6*ev8;
	//	len = (unsigned int)1/
len = 1.0/sqrt(U0*U0 + U1*U1 + U2*U2);
		U0*=len;U1*=len;U2*=len;
		
		U3 = B1*ev0 + B4*ev3 + B7*ev6;
		U4 = B1*ev1 + B4*ev4 + B7*ev7;
		U5 = B1*ev2 + B4*ev5 + B7*ev8;
	//	len = (unsigned int)1/
len = 1.0/sqrt(U3*U3 + U4*U4 + U5*U5);
		U3*=len;U4*=len;U5*=len;
		
		U6 =   U1*U5 - U4*U2; 
		U7 = - U0*U5 + U3*U2;
		U8 =   U0*U4 - U3*U1;
		//len = (unsigned int)1/
len = 1.0/sqrt(U6*U6 + U7*U7 + U8*U8);
		U6*=len;U7*=len;U8*=len;

	// 6) apply rotation
		float xx,xy,xz;

//printf("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n",U0, U1, U2, U3, U4, U5, U6, U7, U8);


		for (int pi = 0; pi < protlen; pi++) {
		   xx = atoms2[3*pi]*U0 + atoms2[3*pi+1]*U1 + atoms2[3*pi+2]*U2;
		   xy = atoms2[3*pi]*U3 + atoms2[3*pi+1]*U4 + atoms2[3*pi+2]*U5;
		   xz = atoms2[3*pi]*U6 + atoms2[3*pi+1]*U7 + atoms2[3*pi+2]*U8;
		   //cout << atoms2[3*pi] << " ";
		   atoms2[3*pi] = xx; atoms2[3*pi+1] = xy; atoms2[3*pi+2] = xz;
		}
}

unsigned short * gdtCPUOneReference(float *reference, float *arr,  int conformers,int protlen,int score) {
    /*
     *  score:
     * 1 = GDT
     * 2 = TM score
     * 3 = MaxSub
    */
	//int triangleLen = conformers*(conformers-1)/2;
	unsigned short* result = NULL;
	result = new unsigned short[conformers];
  // tmscore cutoff
double d0;
if(score==2){
		     d0=1.24*pow((protlen-15),inv3)-1.8;
		          if(protlen<=15 || d0<0.5)
			    	    d0=0.5;
      }
const        int e[] ={protlen,(int)(protlen/4),(int)(protlen/2),4};
//#pragma omp parallel for schedule(dynamic)
   for (int idx=0;idx<conformers;idx++) {

    float* atoms1 = new float[protlen*3]; 
	float* atoms2 = new float[protlen*3];
	
    int s_i = idx*protlen*3; 
//	int s_j =  3*protlen*j;
	int GDT1 = 0,GDT2 = 0,GDT4 = 0,GDT8 = 0;
    double TMScore=-1.0;
	double MaxSub = -1.0;
	
	bool* indexes = new bool[protlen];
                 double* a1_cm = new double[3];
                  double* a2_cm = new double[3];
 int i,si,s,ii,iii,j2,iter,dupa,fragmentlen,ci,cii;	
  for (i = 0; i < protlen*3; i++) { // copy two models into local array
	  atoms1[i] = reference[i];
	  atoms2[i] = arr[s_i+i];
  }  
	
	double cutoff=3.5;
	
	//cutoff += 0.5; //poprawic

	  for (si=1;si<=4;si++) { 

	     s = e[si-1];//protlen/pow(2,(si-1));//tmscore fragments

	    if (s<4)
		continue;	

	    for (ii = 0; ii < protlen-s; ii++) { // make fragments on whole chain
	    	for (iii = 0; iii < protlen; iii++) 
		  indexes[iii] = false;
		for (j2 = 0; j2 < s; j2++)  {
		  indexes[ii+j2]=true;
		  //printf("%d",ii+j2);
	  }

		 iter = 0;

		int alilen=protlen;
		int alilen_old=protlen-1;
		for (dupa=1;dupa<=4;dupa++) { //while(iter<4) {
		  if(alilen==alilen_old)
		    break;
		  else
 		    alilen_old=alilen;
		  
		  fragmentlen=0;
		  iter++;
		  a1_cm[0]=a1_cm[1]=a1_cm[2]=a2_cm[0]=a2_cm[1]=a2_cm[2]=0.0;
		  
		  for (ci = 0; ci < protlen; ci++) {
				 if(indexes[ci]==true) {
				//	 printf("%d",ci);
						fragmentlen++;
						a1_cm[0] += atoms1[3*ci];
						a1_cm[1] += atoms1[3*ci+1];
						a1_cm[2] += atoms1[3*ci+2];
						
						a2_cm[0] += atoms2[3*ci];
						a2_cm[1] += atoms2[3*ci+1];
						a2_cm[2] += atoms2[3*ci+2];
				 }
		 }
		 
		 if(fragmentlen<3) // tego chyba nie 
		   break;
		 
		 float div = 1/(double)fragmentlen;
		 for(int cmi2=0;cmi2<3;cmi2++ ) { 
		 a1_cm[cmi2] *=div;
		 a2_cm[cmi2] *=div;
		 }
		 
		 for (cii=0;cii<protlen;cii++) { // move to the center of mass
		   atoms1[3*cii]   -= a1_cm[0];
		   atoms1[3*cii+1] -= a1_cm[1];
		   atoms1[3*cii+2] -= a1_cm[2];
		      
		   atoms2[3*cii]   -= a2_cm[0];
		   atoms2[3*cii+1] -= a2_cm[1];
		   atoms2[3*cii+2] -= a2_cm[2];   
		
	
			
		 }
		 //delete [] a1_cm; delete [] a2_cm;
		 rotate(atoms1, atoms2, indexes, protlen);
		 
 		int tGDT1=0,tGDT2=0,tGDT4=0,tGDT8=0;
		double tmpTMScore=0.0;
		double tmpMaxSub = 0.0;
		double du;
		double distance;
		for (int pi = 0; pi < protlen; pi++) {
		   // calc distance after superimposition
		   distance = dist(atoms2,atoms1,pi);
if(score==2){
		   du = distance/d0;	   
		   tmpTMScore += 1/(1+du*du);
}
		  
		   if (distance<cutoff) { // jesli jest mniej niz cztery to rozszerz cutoff (tak jak w tmscore)
		     indexes[pi] = true; 
		     alilen++;
		  }
		   else
		     indexes[pi] = false;
double du_ms;
switch(score) {
case 1:
                   if(distance<=8) {
                      tGDT8++;//=tGDT8+1;
                      if(distance<=4) {
                           tGDT4++;//=tGDT4+1;
                           if(distance<=2) {
                               tGDT2++;//=tGDT2+1;
                               if(distance<=1) {
                                          tGDT1++;//=tGDT1+1;
                               }
                           }
                      }
                   }
break;

case 3:	
//if(score==3) {
 du_ms = distance/3.5;
tmpMaxSub += 1/(1.0+du_ms*du_ms);
break;
	}	
		   
}
		if(score==2 && tmpTMScore>TMScore)
		  TMScore = tmpTMScore;
		if(score==3 && tmpMaxSub>MaxSub)
		  MaxSub = tmpMaxSub;
		// update best GDT
		if (score==1) {
		if (tGDT1>GDT1)
		   GDT1 = tGDT1;
		if (tGDT2>GDT2)
		   GDT2 = tGDT2;
	        if (tGDT4>GDT4)
		   GDT4 = tGDT4;
		if (tGDT8>GDT8)
		   GDT8 = tGDT8;
		}

		} //end of while loop
	  }
	    
	  }
	  
	if(score==1) {
        float gdt = (float)(GDT1+GDT2+GDT4+GDT8)/(4.0*protlen);
        result[idx] = (unsigned short)(gdt*1000);
    }
    else if (score==2) {
        result[idx] = (unsigned short)(TMScore*1000.0/protlen);
    }
    else if (score==3) {
        result[idx] = (unsigned short)(MaxSub*1000.0/protlen);
    }
    
    
    
    delete [] atoms1; delete [] atoms2; delete [] indexes;

  }

return result;	
	}

unsigned short * gdtCPU(float *arr, int conformers,int protlen,int score) {

double d0;
if(score==2){
		     d0=1.24*pow((protlen-15),inv3)-1.8;
		          if(protlen<=15 || d0<0.5)
			    	    d0=0.5;
      }
const        int e[] ={4,(int)(protlen/4),(int)(protlen/2),protlen};
int triangleLen = conformers*(conformers-1)/2;
unsigned short* result = NULL;
result = new unsigned short[triangleLen];
		    
  
//#pragma omp parallel for schedule(dynamic)
   for (int idx=0;idx<conformers-1;idx++) {
      int iter000=0;
      
        float* atoms1 = new float[protlen*3]; 
        float* atoms2 = new float[protlen*3];
        bool* indexes = new bool[protlen];
        double* a1_cm = new double[3];
        double* a2_cm = new double[3];
        int i,si,s,ii,iii,j2,iter,fragmentlen,ci,cii;	
        
    for (int j = idx+1;j<conformers;j++) {

	////////////////////////////////////
    int s_i = idx*protlen*3; 
	int s_j =  3*protlen*j;

    int GDT1 = 0,GDT2 = 0,GDT4 = 0,GDT8 = 0;
    double TMScore=-1.0;
	double MaxSub = -1.0;
	

  for (i = 0; i < protlen*3; i++) { // copy two models into local array
  	  atoms1[i] = arr[s_j+i];
	  atoms2[i] = arr[s_i+i];
  }  
	
	double cutoff=3.5;
	
	//cutoff += 0.5; //poprawic

	  for (si=1;si<=4;si++) { 

	     s = e[si-1];//protlen/pow(2,(si-1));//tmscore fragments
	    if (s<4)
		continue;	

	    for (ii = 0; ii < protlen-s; ii++) { // make fragments on whole chain
	    	for (iii = 0; iii < protlen; iii++) 
		  indexes[iii] = false;
		for (j2 = 0; j2 < s; j2++)  
		  indexes[ii+j2]=true;

		 iter = 0;

		int alilen=protlen;
		int alilen_old=protlen-1;
		while(iter<4) {
		  if(alilen==alilen_old)
		    break;
		  else
 		    alilen_old=alilen;
		  
		  fragmentlen=0;
		  iter++;
		  a1_cm[0]=a1_cm[1]=a1_cm[2]=a2_cm[0]=a2_cm[1]=a2_cm[2]=0.0;
		  
		  for (ci = 0; ci < protlen; ci++) {
				 if(indexes[ci]==true) {
						fragmentlen++;
						a1_cm[0] += atoms1[3*ci];
						a1_cm[1] += atoms1[3*ci+1];
						a1_cm[2] += atoms1[3*ci+2];
						
						a2_cm[0] += atoms2[3*ci];
						a2_cm[1] += atoms2[3*ci+1];
						a2_cm[2] += atoms2[3*ci+2];
				 }
		 }
		 
		 if(fragmentlen<3) // tego chyba nie 
		   break;
		 
		 double div = (unsigned int)1/(double)fragmentlen;
		 for(int cmi2=0;cmi2<3;cmi2++ ) { 
		 a1_cm[cmi2] *=div;
		 a2_cm[cmi2] *=div;
		 }
		 
		 for (cii=0;cii<protlen;cii++) { // move to the center of mass
		   atoms1[3*cii]   -= a1_cm[0];
		   atoms1[3*cii+1] -= a1_cm[1];
		   atoms1[3*cii+2] -= a1_cm[2];
		   
		   atoms2[3*cii] -= a2_cm[0];
		   atoms2[3*cii+1] -= a2_cm[1];
		   atoms2[3*cii+2] -= a2_cm[2];  

		 }
		 //delete [] a1_cm; delete [] a2_cm;
		 rotate(atoms1, atoms2, indexes, protlen);
		 
 		int tGDT1=0,tGDT2=0,tGDT4=0,tGDT8=0;
		double tmpTMScore=0.0;
		double tmpMaxSub = 0.0;
		double du;
		double distance;
		for (int pi = 0; pi < protlen; pi++) {
		   // calc distance after superimposition
		   distance = dist(atoms2,atoms1,pi);
if(score==2){
		   du = distance/d0;	   
		   tmpTMScore += 1/(1+du*du);
}
		  
		   if (distance<cutoff) { // jesli jest mniej niz cztery to rozszerz cutoff (tak jak w tmscore)
		     indexes[pi] = true; 
		     alilen++;
		  }
		   else
		     indexes[pi] = false;
double du_ms;
switch(score) {
case 1:
                   if(distance<=8) {
                      tGDT8++;//=tGDT8+1;
                      if(distance<=4) {
                           tGDT4++;//=tGDT4+1;
                           if(distance<=2) {
                               tGDT2++;//=tGDT2+1;
                               if(distance<=1) {
                                          tGDT1++;//=tGDT1+1;
                               }
                           }
                      }
                   }
break;

case 3:	
//if(score==3) {
 du_ms = distance/3.5;
tmpMaxSub += 1/(1+du_ms*du_ms);
break;
	}	
		   
}
		if(score==2 && tmpTMScore>TMScore)
		  TMScore = tmpTMScore;
		if(score==3 && tmpMaxSub>MaxSub)
		  MaxSub = tmpMaxSub;
		// update best GDT
		if (score==1) {
		if (tGDT1>GDT1)
		   GDT1 = tGDT1;
		if (tGDT2>GDT2)
		   GDT2 = tGDT2;
	        if (tGDT4>GDT4)
		   GDT4 = tGDT4;
		if (tGDT8>GDT8)
		   GDT8 = tGDT8;
		}

		} //end of while loop
	  }
	    
	  }
    
	////////////////////////////////////

    int index1d = (int)(iter000+ idx*(conformers-1) - idx*(idx-1)/2);
    if(score==1) {
        float gdt = (float)(GDT1+GDT2+GDT4+GDT8)/(4.0*protlen);
        result[index1d] = (unsigned short)(gdt*1000);
    }
    else if (score==2) {
        result[index1d] = (unsigned short)(TMScore*1000/protlen);
    }
    else if (score==3) {
        result[index1d] = (unsigned short)(MaxSub*1000/protlen);
    }
    
    iter000++;
    

    } // end of pair comparison
  }
  
  return result;
}

unsigned short * gdtCPUOneReferenceExt(float *reference, float *arr,  int conformers,int protlen,int score) {
    /*
     *  score:
     * 1 = GDT
     * 2 = TM score
     * 3 = MaxSub
    */
	//int triangleLen = conformers*(conformers-1)/2;
	unsigned short* result = NULL;
	result = new unsigned short[conformers];
  // tmscore cutoff
double d0;
if(score==2){
		     d0=1.24*pow((protlen-15),inv3)-1.8;
		          if(protlen<=15 || d0<0.5)
			    	    d0=0.5;
      }
const        int e[] ={4,(int)(protlen/4),(int)(protlen/2),protlen};
//#pragma omp parallel for schedule(dynamic)
   for (int idx=0;idx<conformers;idx++) {


    float* atoms1 = new float[protlen*3]; 
	float* atoms2 = new float[protlen*3];
	
    int s_i = idx*protlen*3; 
//	int s_j =  3*protlen*j;
	int GDT1 = 0,GDT2 = 0,GDT4 = 0,GDT8 = 0;
    double TMScore=-1.0;
	double MaxSub = -1.0;
	
	bool* indexes = new bool[protlen];
                 double* a1_cm = new double[3];
                  double* a2_cm = new double[3];
 int i,si,s,ii,iii,j2,iter,dupa,fragmentlen,ci,cii;	
 
	
    for(unsigned int cti=1;cti<=20;cti++) {
        
          for (i = 0; i < protlen*3; i++) { // copy two models into local array
	  atoms1[i] = reference[i];
	  atoms2[i] = arr[s_i+i];
  } 
  
	double cutoff=0.5*cti;
	
	//cutoff += 0.5; //poprawic

	  for (si=1;si<=4;si++) { 

	     s = e[si-1];//protlen/pow(2,(si-1));//tmscore fragments
	

	    for (ii = 0; ii < protlen-s; ii++) { // make fragments on whole chain
	    	for (iii = 0; iii < protlen; iii++) 
		  indexes[iii] = false;
		for (j2 = 0; j2 < s; j2++)  
		  indexes[ii+j2]=true;

		 iter = 0;

		int alilen=protlen;
		int alilen_old=protlen-1;
		for (dupa=1;dupa<=4;dupa++) { //while(iter<4) {
		  if(alilen==alilen_old)
		    break;
		  else
 		    alilen_old=alilen;
		  
		  fragmentlen=0;
		  iter++;
		  a1_cm[0]=a1_cm[1]=a1_cm[2]=a2_cm[0]=a2_cm[1]=a2_cm[2]=0.0;
		  
		  for (ci = 0; ci < protlen; ci++) {
				 if(indexes[ci]==true) {
						fragmentlen++;
						a1_cm[0] += atoms1[3*ci];
						a1_cm[1] += atoms1[3*ci+1];
						a1_cm[2] += atoms1[3*ci+2];
						
						a2_cm[0] += atoms2[3*ci];
						a2_cm[1] += atoms2[3*ci+1];
						a2_cm[2] += atoms2[3*ci+2];
				 }
		 }
		 
		 if(fragmentlen<3) // tego chyba nie 
		   break;
		 
		 double div = (unsigned int)1/(double)fragmentlen;
		 for(int cmi2=0;cmi2<3;cmi2++ ) { 
		 a1_cm[cmi2] *=div;
		 a2_cm[cmi2] *=div;
		 }
		 
		 for (cii=0;cii<protlen;cii++) { // move to the center of mass
		   atoms1[3*cii]   -= a1_cm[0];
		   atoms1[3*cii+1] -= a1_cm[1];
		   atoms1[3*cii+2] -= a1_cm[2];
		   
		   atoms2[3*cii] -= a2_cm[0];
		   atoms2[3*cii+1] -= a2_cm[1];
		   atoms2[3*cii+2] -= a2_cm[2];  

		 }
		 //delete [] a1_cm; delete [] a2_cm;
		 rotate(atoms1, atoms2, indexes, protlen);
		 
 		int tGDT1=0,tGDT2=0,tGDT4=0,tGDT8=0;
		double tmpTMScore=0.0;
		double tmpMaxSub = 0.0;
		double du;
		double distance;
		for (int pi = 0; pi < protlen; pi++) {
		   // calc distance after superimposition
		   distance = dist(atoms2,atoms1,pi);
if(score==2){
		   du = distance/d0;	   
		   tmpTMScore += 1/(1+du*du);
}
		  
		   if (distance<cutoff) { // jesli jest mniej niz cztery to rozszerz cutoff (tak jak w tmscore)
		     indexes[pi] = true; 
		     alilen++;
		  }
		   else
		     indexes[pi] = false;
double du_ms;
switch(score) {
case 1:
                   if(distance<=8) {
                      tGDT8++;//=tGDT8+1;
                      if(distance<=4) {
                           tGDT4++;//=tGDT4+1;
                           if(distance<=2) {
                               tGDT2++;//=tGDT2+1;
                               if(distance<=1) {
                                          tGDT1++;//=tGDT1+1;
                               }
                           }
                      }
                   }
break;

case 3:	
//if(score==3) {
 du_ms = distance/3.5;
tmpMaxSub += 1/(1+du_ms*du_ms);
break;
	}	
		   
}
		if(score==2 && tmpTMScore>TMScore)
		  TMScore = tmpTMScore;
		if(score==3 && tmpMaxSub>MaxSub)
		  MaxSub = tmpMaxSub;
		// update best GDT
		if (score==1) {
		if (tGDT1>GDT1)
		   GDT1 = tGDT1;
		if (tGDT2>GDT2)
		   GDT2 = tGDT2;
	        if (tGDT4>GDT4)
		   GDT4 = tGDT4;
		if (tGDT8>GDT8)
		   GDT8 = tGDT8;
		}

		} //end of while loop
	  }
	    
	  }
	  
  } //end cti
	if(score==1) {
        float gdt = (float)(GDT1+GDT2+GDT4+GDT8)/(4.0*protlen);
        result[idx] = (unsigned short)(gdt*1000);
    }
    else if (score==2) {
        result[idx] = (unsigned short)(TMScore*1000/protlen);
    }
    else if (score==3) {
        result[idx] = (unsigned short)(MaxSub*1000/protlen);
    }
    
    
    
    delete [] atoms1; delete [] atoms2; delete [] indexes;

  }

return result;	
	}
