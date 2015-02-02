#include <cmath>
#include "Matrix_3x3.h"
#include "CpptrajStdio.h"
#include "Constants.h" // PI, RADDEG

// COPY CONSTRUCTOR
Matrix_3x3::Matrix_3x3(const Matrix_3x3& rhs) {
  M_[0] = rhs.M_[0];
  M_[1] = rhs.M_[1];
  M_[2] = rhs.M_[2];
  M_[3] = rhs.M_[3];
  M_[4] = rhs.M_[4];
  M_[5] = rhs.M_[5];
  M_[6] = rhs.M_[6];
  M_[7] = rhs.M_[7];
  M_[8] = rhs.M_[8];
}

/// CONSTRUCTOR - Takes array of 9, row-major
Matrix_3x3::Matrix_3x3(const double *Min) {
  M_[0] = Min[0]; 
  M_[1] = Min[1]; 
  M_[2] = Min[2]; 
  M_[3] = Min[3]; 
  M_[4] = Min[4]; 
  M_[5] = Min[5]; 
  M_[6] = Min[6]; 
  M_[7] = Min[7]; 
  M_[8] = Min[8]; 
}

/// CONSTRUCTOR - Set all elements to xIn 
Matrix_3x3::Matrix_3x3(double xIn) {
  M_[0] = xIn; 
  M_[1] = xIn; 
  M_[2] = xIn; 
  M_[3] = xIn; 
  M_[4] = xIn; 
  M_[5] = xIn; 
  M_[6] = xIn; 
  M_[7] = xIn; 
  M_[8] = xIn; 
}

/// CONSTRUCTOR - Set diagonal
Matrix_3x3::Matrix_3x3(double d1, double d2, double d3) {
  M_[0] = d1;
  M_[1] = 0;
  M_[2] = 0;
  M_[3] = 0;
  M_[4] = d2;
  M_[5] = 0;
  M_[6] = 0;
  M_[7] = 0;
  M_[8] = d3;
}

// Assignment
Matrix_3x3& Matrix_3x3::operator=(const Matrix_3x3& rhs) {
  if (this==&rhs) return *this;
  M_[0] = rhs.M_[0];
  M_[1] = rhs.M_[1];
  M_[2] = rhs.M_[2];
  M_[3] = rhs.M_[3];
  M_[4] = rhs.M_[4];
  M_[5] = rhs.M_[5];
  M_[6] = rhs.M_[6];
  M_[7] = rhs.M_[7];
  M_[8] = rhs.M_[8];
  return *this;
}

void Matrix_3x3::Zero() {
  M_[0] = 0;
  M_[1] = 0;
  M_[2] = 0;
  M_[3] = 0;
  M_[4] = 0;
  M_[5] = 0;
  M_[6] = 0;
  M_[7] = 0;
  M_[8] = 0;
}

// Matrix_3x3::Print()
void Matrix_3x3::Print(const char* Title) const 
{
  mprintf("    %s\n",Title);
  mprintf("     %8.4f %8.4f %8.4f\n", M_[0], M_[1], M_[2]);
  mprintf("     %8.4f %8.4f %8.4f\n", M_[3], M_[4], M_[5]);
  mprintf("     %8.4f %8.4f %8.4f\n", M_[6], M_[7], M_[8]);
}

// -----------------------------------------------------------------------------
/// Max number of iterations to execute Jacobi algorithm
const int Matrix_3x3::MAX_ITERATIONS = 50;

#define ROTATE(ARR,MAJ1,MIN1,MAJ2,MIN2) { \
  dg = ARR[MAJ1 + MIN1]; \
  dh = ARR[MAJ2 + MIN2]; \
  ARR[MAJ1 + MIN1] = dg - ds*(dh+dg*tau); \
  ARR[MAJ2 + MIN2] = dh + ds*(dg-dh*tau); }

// Matrix_3x3::Diagonalize()
/** Diagonalize the matrix using Jacobi method. Eigenvectors are stored in 
  * columns. 
  * \param vecD Output eigenvalues.
  */
int Matrix_3x3::Diagonalize( Vec3& vecD ) 
{
  // Store this matrix
  double mat[9];
  mat[0] = M_[0];
  mat[1] = M_[1];
  mat[2] = M_[2];
  mat[3] = M_[3];
  mat[4] = M_[4];
  mat[5] = M_[5];
  mat[6] = M_[6];
  mat[7] = M_[7];
  mat[8] = M_[8];
  // Create identity matrix
  M_[0] = 1;
  M_[1] = 0;
  M_[2] = 0;
  M_[3] = 0;
  M_[4] = 1;
  M_[5] = 0;
  M_[6] = 0;
  M_[7] = 0;
  M_[8] = 1;
  // Set vectors B and D equal to diagonal of mat. vector Z is 0.
  double vecB[3], vecZ[3];
  vecB[0] = vecD[0] = mat[0]; 
  vecB[1] = vecD[1] = mat[4]; 
  vecB[2] = vecD[2] = mat[8];
  vecZ[0] = 0; 
  vecZ[1] = 0; 
  vecZ[2] = 0;
  // MAIN LOOP
  double tresh = 0;
  //int nrot = 0;
  for (int i = 0; i < MAX_ITERATIONS; ++i) {
    // sm = SUM of UPPER RIGHT TRIANGLE
    double sm = fabs(mat[1]) + fabs(mat[2]) + fabs(mat[5]);
    if (sm == 0) return 0;
    
    if (i < 3)
      tresh = 0.2 * sm / 9;
    else
      tresh = 0;
    // BEGIN INNER LOOP OVER UPPER RIGHT TRIANGLE
    double dt;
    //int p3 = 0;
    int ip, p3;
    for ( ip = p3 = 0; ip < 2; ++ip, p3+=3) {
      for ( int iq = ip + 1; iq < 3; ++iq ) {
        int midx = p3 + iq;
        double dg = 100.0 * fabs(mat[midx]);
        if ( i > 3 && fabs(vecD[ip]) + dg == fabs(vecD[ip]) && 
                      fabs(vecD[iq]) + dg == fabs(vecD[iq]) )
        {
          mat[midx] = 0;
        } else if ( fabs(mat[midx]) > tresh) {
          double dh = vecD[iq] - vecD[ip];
          if (fabs(dh) + dg == fabs(dh))
            dt = mat[p3 + iq] / dh;
          else {
            double theta = 0.5 * dh / mat[midx];
            dt = 1.0 / (fabs(theta)+(double)sqrt(1.0+theta*theta));
            if (theta < 0.0)
              dt = -dt;
          }
          double dc = 1.0 / (double)sqrt(1+dt*dt);
          double ds = dt * dc;
          double tau = ds / (1.0+dc);
          dh = dt * mat[midx];
          vecZ[ip] -= dh;
          vecZ[iq] += dh;
          vecD[ip] -= dh;
          vecD[iq] += dh;
          mat[midx] = 0;
          int j, j3;
          for (j=j3=0; j<=ip-1; j++,j3+=3)
            ROTATE(mat,j3,ip,j3,iq)
          for (int j=ip+1; j<=iq-1; j++)
            ROTATE(mat,p3,j,j*3,iq)
          for (int j=iq+1; j<3; j++)
            ROTATE(mat,p3,j,iq*3,j)

          for (j3=0; j3<9; j3+=3)
            ROTATE(M_,j3,ip,j3,iq)

          //++nrot;
        }
      }
    } // END INNER LOOP OVER UPPER RIGHT TRIANGLE
    vecB[0] += vecZ[0];
    vecD[0] = vecB[0];
    vecZ[0] = 0;
    vecB[1] += vecZ[1];
    vecD[1] = vecB[1];
    vecZ[1] = 0;
    vecB[2] += vecZ[2];
    vecD[2] = vecB[2];
    vecZ[2] = 0;
  }
  mprintf("Too many iterations in routine!\n");
  return 1;
}

// Matrix_3x3::Diagonalize_Sort()
/** Diagonalize the matrix and sort eigenvalues/eigenvectors in 
  * descending order. Eigenvectors will be stored in rows,
  * (V0x, V0y, V0z, V1x, ... V2z).
  * \param EvalOut Output eigenvalues.
  */
int Matrix_3x3::Diagonalize_Sort(Vec3& EvalOut) 
{
  Vec3 Eval;
  if ( Diagonalize( Eval ) ) 
  {
    mprintf("Convergence failed.\n");
    return 1; 
  }
  //printMatrix_3x3("Eigenvector Matrix", Evec);

  if (Eval[0] > Eval[1] && Eval[0] > Eval[2]) { // 0 is max
    if (Eval[1] > Eval[2]) {
      i1_ = 0; i2_ = 1; i3_ = 2;
    } else {
      i1_ = 0; i2_ = 2; i3_ = 1;
    }
  } else if (Eval[1] > Eval[0] && Eval[1] > Eval[2]) { // 1 is max
    if (Eval[0] > Eval[2]) {
      i1_ = 1; i2_ = 0; i3_ = 2;
    } else {
      i1_ = 1; i2_ = 2; i3_ = 0;
    }
  } else if (Eval[0] > Eval[1]) { // 2 is max
    i1_ = 2; i2_ = 0; i3_ = 1;
  } else {
    i1_ = 2; i2_ = 1; i3_ = 0;
  }
  //mprintf("EIGENVALUE ORDER (0=high, 3=med, 6=low): %i %i %i\n",i1_,i2_,i3_);

  // Swap Eigenvectors - place them in rows
  Matrix_3x3 Evec(*this);
  M_[0] = Evec[i1_  ];
  M_[1] = Evec[i1_+3];
  M_[2] = Evec[i1_+6];

  M_[3] = Evec[i2_  ];
  M_[4] = Evec[i2_+3];
  M_[5] = Evec[i2_+6];

  M_[6] = Evec[i3_  ];
  M_[7] = Evec[i3_+3];
  M_[8] = Evec[i3_+6];

  // Swap eigenvalues
  EvalOut[0] = Eval[i1_];
  EvalOut[1] = Eval[i2_];
  EvalOut[2] = Eval[i3_];

  return 0;
}

/**  The jacobi diagonalization procedure can sometimes result
  *  in eigenvectors which when applied to transform the coordinates
  *  result in a a chiral inversion about the Y axis.  This code catches
  *  this case, reversing the offending eigenvectors.
  *  
  *  NOTE: the idea of rotating the coordinate basis vectors came from 
  *  some code posted to the computational chemistry mailing list 
  *  (chemistry@osc) in a summary of methods to perform principal axis 
  *  alignment...
  *
  * It is expected that the eigenvector matrix has eigenvectors in rows.
  */
int Matrix_3x3::jacobiCheckChirality()
{
  Matrix_3x3 points(*this);
  Matrix_3x3 result;
  //points.Print("POINTS"); 

  // rotate vector three into XZ plane
  result.RotationAroundZ( points[2], points[5] ); // Ev0z, Ev1z
  result *= points;
  //result.Print("POINTS1");

  // rotate vector three into Z axis
  points.RotationAroundY( result[2], result[8] ); 
  points *= result;
  //points.Print("POINTS2");

  // rotate vector one into XZ
  result.RotationAroundZ( points[0], points[3] );
  result *= points;
  //result.Print("POINTS3");

  // rotate vector one into X 
  points.RotationAroundY( result[2], result[0] );
  points *= result;
  //points.Print("POINTS4");

  // has Y changed sign? If so, flip Y eigenvector (row 1) 
  if ( points[4] < 0 ) {
    M_[3] = -M_[3];
    M_[4] = -M_[4];
    M_[5] = -M_[5];
    return 1;
  }
  return 0;
}

// Matrix_3x3::Diagonalize_Sort_Chirality()
int Matrix_3x3::Diagonalize_Sort_Chirality(Vec3& EvalOut, int debug)
{
  if ( Diagonalize_Sort( EvalOut ) )
    return 1;

  // Invert eigenvector signs based on ordering to avoid reflections
  if (i1_ == 0 && i2_ == 2 && i3_ == 1) {
    M_[3] = -M_[3];
    M_[4] = -M_[4];
    M_[5] = -M_[5];
  } else if (i1_ == 2 && i2_ == 0 && i3_ == 1) {
    M_[0] = -M_[0];
    M_[1] = -M_[1];
    M_[2] = -M_[2];
    M_[3] = -M_[3];
    M_[4] = -M_[4];
    M_[5] = -M_[5];
    M_[6] = -M_[6];
    M_[7] = -M_[7];
    M_[8] = -M_[8];
  }

  // Flip Y-vector if necessary 
  if (jacobiCheckChirality( ) && debug>0)
    mprintf("Warning: PRINCIPAL: CHECK CHIRALITY: Y eigenvector sign swapped!\n");
   
  return 0;
}

// -----------------------------------------------------------------------------
/** Columns of matrix become rows and vice-versa. */
void Matrix_3x3::Transpose() {
  double U1 = M_[1];
  double U2 = M_[2];
  double U3 = M_[3];
  double U5 = M_[5];
  double U6 = M_[6];
  double U7 = M_[7];
  M_[1] = U3;
  M_[2] = U6;
  M_[3] = U1;
  M_[5] = U7;
  M_[6] = U2;
  M_[7] = U5;
}

// Matrix_3x3::operator*=()
Matrix_3x3& Matrix_3x3::operator*=(const Matrix_3x3& rhs) {
  double Row[9];
  Row[0] = M_[0];
  Row[1] = M_[1];
  Row[2] = M_[2];
  Row[3] = M_[3];
  Row[4] = M_[4];
  Row[5] = M_[5];
  Row[6] = M_[6];
  Row[7] = M_[7];
  Row[8] = M_[8];
  M_[0] = (Row[0] * rhs.M_[0]) + (Row[1] * rhs.M_[3]) + (Row[2] * rhs.M_[6]);
  M_[1] = (Row[0] * rhs.M_[1]) + (Row[1] * rhs.M_[4]) + (Row[2] * rhs.M_[7]);
  M_[2] = (Row[0] * rhs.M_[2]) + (Row[1] * rhs.M_[5]) + (Row[2] * rhs.M_[8]);
  M_[3] = (Row[3] * rhs.M_[0]) + (Row[4] * rhs.M_[3]) + (Row[5] * rhs.M_[6]);
  M_[4] = (Row[3] * rhs.M_[1]) + (Row[4] * rhs.M_[4]) + (Row[5] * rhs.M_[7]);
  M_[5] = (Row[3] * rhs.M_[2]) + (Row[4] * rhs.M_[5]) + (Row[5] * rhs.M_[8]);
  M_[6] = (Row[6] * rhs.M_[0]) + (Row[7] * rhs.M_[3]) + (Row[8] * rhs.M_[6]);
  M_[7] = (Row[6] * rhs.M_[1]) + (Row[7] * rhs.M_[4]) + (Row[8] * rhs.M_[7]);
  M_[8] = (Row[6] * rhs.M_[2]) + (Row[7] * rhs.M_[5]) + (Row[8] * rhs.M_[8]);
  return *this;
}

Matrix_3x3 Matrix_3x3::operator*(Matrix_3x3 const& rhs) const {
  Matrix_3x3 result;
  result.M_[0] = M_[0]*rhs.M_[0] + M_[1]*rhs.M_[3] + M_[2]*rhs.M_[6];
  result.M_[1] = M_[0]*rhs.M_[1] + M_[1]*rhs.M_[4] + M_[2]*rhs.M_[7];
  result.M_[2] = M_[0]*rhs.M_[2] + M_[1]*rhs.M_[5] + M_[2]*rhs.M_[8];
  result.M_[3] = M_[3]*rhs.M_[0] + M_[4]*rhs.M_[3] + M_[5]*rhs.M_[6];
  result.M_[4] = M_[3]*rhs.M_[1] + M_[4]*rhs.M_[4] + M_[5]*rhs.M_[7];
  result.M_[5] = M_[3]*rhs.M_[2] + M_[4]*rhs.M_[5] + M_[5]*rhs.M_[8];
  result.M_[6] = M_[6]*rhs.M_[0] + M_[7]*rhs.M_[3] + M_[8]*rhs.M_[6];
  result.M_[7] = M_[6]*rhs.M_[1] + M_[7]*rhs.M_[4] + M_[8]*rhs.M_[7];
  result.M_[8] = M_[6]*rhs.M_[2] + M_[7]*rhs.M_[5] + M_[8]*rhs.M_[8];
  return result;
}

Matrix_3x3 Matrix_3x3::TransposeMult(Matrix_3x3 const& rhs) const {
  Matrix_3x3 result;
  result.M_[0] = M_[0]*rhs.M_[0] + M_[1]*rhs.M_[1] + M_[2]*rhs.M_[2];
  result.M_[1] = M_[0]*rhs.M_[3] + M_[1]*rhs.M_[4] + M_[2]*rhs.M_[5];
  result.M_[2] = M_[0]*rhs.M_[6] + M_[1]*rhs.M_[7] + M_[2]*rhs.M_[8];
  result.M_[3] = M_[3]*rhs.M_[0] + M_[4]*rhs.M_[1] + M_[5]*rhs.M_[2];
  result.M_[4] = M_[3]*rhs.M_[3] + M_[4]*rhs.M_[4] + M_[5]*rhs.M_[5];
  result.M_[5] = M_[3]*rhs.M_[6] + M_[4]*rhs.M_[7] + M_[5]*rhs.M_[8];
  result.M_[6] = M_[6]*rhs.M_[0] + M_[7]*rhs.M_[1] + M_[8]*rhs.M_[2];
  result.M_[7] = M_[6]*rhs.M_[3] + M_[7]*rhs.M_[4] + M_[8]*rhs.M_[5];
  result.M_[8] = M_[6]*rhs.M_[6] + M_[7]*rhs.M_[7] + M_[8]*rhs.M_[8];
  return result;
}

// Matrix_3x3::RotationAroundZ()
void Matrix_3x3::RotationAroundZ(double a1, double a2) {
  double r = sqrt( a1*a1 + a2*a2 );
  M_[0] = a1 / r; //  cos t
  M_[1] = a2 / r; // -sin t
  M_[2] = 0;
  M_[3] = -M_[1]; //  sin t
  M_[4] = M_[0];  //  cos t
  M_[5] = 0;
  M_[6] = 0;
  M_[7] = 0;
  M_[8] = 1;
}

// Matrix_3x3::RotationAroundY()
void Matrix_3x3::RotationAroundY(double a1, double a2) {
  double r = sqrt( a1*a1 + a2*a2 );
  M_[0] = a2 / r;  //  cos t
  M_[1] = 0;
  M_[2] = -a1 / r; //  sin t
  M_[3] = 0;
  M_[4] = 1;
  M_[5] = 0;
  M_[6] = -M_[2];  // -sin t
  M_[7] = 0;
  M_[8] = M_[0];   //  cos t
}

/** Given an axis of rotation V and a magnitude (radians), calculate a 
  * rotation matrix.
  */
void Matrix_3x3::CalcRotationMatrix(Vec3 const& V, double theta) {
  // Compute all prefactors
  double ux2  = V[0] * V[0];
  double uxuy = V[0] * V[1];
  double uxuz = V[0] * V[2];
  double uy2  = V[1] * V[1];
  double uyuz = V[1] * V[2];
  double uz2  = V[2] * V[2];
  double c    = cos(theta);
  double s    = sin(theta);
  double c1   = 1 - c;
  double uxs  = V[0] * s;
  double uys  = V[1] * s;
  double uzs  = V[2] * s;
  // Store rotation matrix elements
  M_[0] = c + (ux2 * c1);
  M_[3] = (uxuy * c1) + uzs;
  M_[6] = (uxuz * c1) - uys;

  M_[1] = (uxuy * c1) - uzs;
  M_[4] = c + (uy2 * c1);
  M_[7] = (uyuz * c1) + uxs;

  M_[2] = (uxuz * c1) + uys;
  M_[5] = (uyuz * c1) - uxs;
  M_[8] = c + (uz2 * c1);
}

/** Given rotations around the X, Y, and Z axes (radians), calculate a
  * rotation matrix.
  */
void Matrix_3x3::CalcRotationMatrix(double psiX, double psiY, double psiZ) {
  Vec3 V(psiX, psiY, psiZ);
  double Psi = V.Normalize(); 
  //mprintf("\t\tcalcRotationMatrix(%.2lf,%.2lf,%.2lf) Psi=%lf\n",
  //        psiX*Constants::RADDEG,psiY*Constants::RADDEG,psiZ*Constants::RADDEG,Psi*Constants::RADDEG);
  CalcRotationMatrix(V, Psi);
}

/** Return angle of rotation from rotation matrix according to
  * cos(t)=(trace(R)-1)/2
  * Equation taken from :
  *   3D game engine design: a practical approach to real-time Computer Graphics,
  *   Volume 385, By David H. Eberly, 2001, p. 16.
  */
double Matrix_3x3::RotationAngle() {
  double trace = M_[0] + M_[4] + M_[8];
  trace = (trace - 1) / 2;
  return acos( trace );
}

/** If theta is between 0 and pi extract axis of rotation from rotation matrix
  * according to:
  *   R - Rt = (2 * sin(theta)) * S, where S is:
  *     0 -z  y
  *     z  0 -x
  *    -y  x  0
  */
Vec3 Matrix_3x3::AxisOfRotation(double theta) {
  if (theta > 0 && theta < Constants::PI) {
    double dx = 1 / (2 * sin(theta));
    Vec3 result( (M_[5] - M_[7]) * dx,
                 (M_[6] - M_[2]) * dx,
                 (M_[1] - M_[3]) * dx );
    result.Normalize();
    return result;
  } else {
    mprintf("Error: axis_of_rotation: Could not extract axis of rotation, angle is %lf\n",
            Constants::RADDEG*theta);
  }
  return Vec3(0.0, 0.0, 0.0);
}
