//该文件来自于readme
#include "../src/HEAAN.h"

using namespace std;
using namespace NTL;

int main() {

  // Parameters //
  //首先说明模链的定义。模链为一系列逐步变小的密文多项式系数模的集合{qL, qL-1, ..., q1}其中，ql=delta*q_(l-1)。delta是缩放因子。特别地，q1=delta*q0.q0和delta都是用户定义
  long logN = 15;//密文多项式次数的对数
  long logQ = 353;//密文多项式系数模的对数
  long logp = 30; ///缩放因子delta的对数。< Larger logp will give you more correct result (smaller computation noise)
  long slots = 1024; ///< This should be power of two
  long numThread = 8;
	
  // Construct and Generate Public Keys //
  TimeUtils timeutils;
  Ring ring(logN, logQ);//环R_Q=Z_Q[X]/(X^N+1)的两个系数：Q,N
  SecretKey secretKey(ring);
  Scheme scheme(secretKey, ring);
  scheme.addLeftRotKeys(secretKey); ///< When you need left rotation for the vectorized message
  scheme.addRightRotKeys(secretKey); ///< When you need right rotation for the vectorized message
  
  // Make Random Array of Complex //
  complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(slots);
  complex<double>* mvec2 = EvaluatorUtils::randomComplexArray(slots);
  
  // Encrypt Two Arry of Complex //
  Ciphertext cipher1 = scheme.encrypt(mvec1, slots, logp, logQ);
  Ciphertext cipher2 = scheme.encrypt(mvec2, slots, logp, logQ);
  
  //下面的结果都会是353，证明密文中的q等于环中的ring.Q
  cout<<"cipher1 logq="<<cipher1.logq<<endl;
  cout<<"cipher2 logq="<<cipher2.logq<<endl;
  cout<<"ring.logQ="<<ring.logQ<<endl;

  
  // Addition //
  Ciphertext cipherAdd = scheme.add(cipher1, cipher2);
  
  // Multiplication And Rescale //

    //   测得np=13 nprimes
//   long np = ceil((2 + cipher1.logq + cipher2.logq + ring.logN + 2)/59.0);
//   cout<<"np="<<np<<endl;
//   long aa=np<<logN;
//   long bb=13*pow(2,15);
//   cout<<"bb="<<bb<<endl;
//   cout<<"np<<logN="<<aa<<endl;

  Ciphertext cipherMult = scheme.mult(cipher1, cipher2);
  Ciphertext cipherMultAfterReScale = scheme.reScaleBy(cipherMult, logp);//rescale(模交换)，
  
  // Rotation //
  long idx = 1;
  Ciphertext cipherRot = scheme.leftRotate(cipher1, idx);
  
  // Decrypt //
  complex<double>* dvec1 = scheme.decrypt(secretKey, cipher1);
  complex<double>* dvec2 = scheme.decrypt(secretKey, cipher2);
  
  return 0;

}
  