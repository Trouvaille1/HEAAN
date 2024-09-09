#include <iostream>
#include "../src/HEAAN.h"

using namespace std;
using namespace NTL;

int main() {

  // Parameters //
  //首先说明模链的定义。模链为一系列逐步变小的密文多项式系数模的集合{qL, qL-1, ..., q1}其中，ql=delta*q_(l-1)。delta是缩放因子。特别地，q1=delta*q0.q0和delta都是用户定义
  long logN = 15;//密文多项式次数的对数
  long logQ = 353;//密文多项式系数模的对数
  long logp = 30; ///缩放因子delta(精度)的对数。< Larger logp will give you more correct result (smaller computation noise)
  long slots = 8; ///< This should be power of two 消息向量槽数/密文槽数
  long numThread = 8;
	
  // Construct and Generate Public Keys //
  Ring ring(logN, logQ);//环R_Q=Z_Q[X]/(X^N+1)的两个系数：Q,N
  SecretKey secretKey(ring);
  Scheme scheme(secretKey, ring);
  scheme.addLeftRotKeys(secretKey); ///< When you need left rotation for the vectorized message
  scheme.addRightRotKeys(secretKey); ///< When you need right rotation for the vectorized message
  
  // Make Random Array of Complex //
  complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(slots);
  complex<double>* mvec2= EvaluatorUtils::randomComplexArray(slots);

  cout<<"旋转之前的mvec1:"<<endl;
    for(int i=0;i<slots;i++){
        cout<<mvec1[i]<<" ";
    }
    cout<<endl;

  // Encrypt //
  Ciphertext cipher1 = scheme.encrypt(mvec1, slots, logp, logQ);
  Ciphertext cipher2 = scheme.encrypt(mvec2, slots, logp, logQ);
  
  //下面的结果都会是353，证明密文中的q等于环中的ring.Q
  cout<<"cipher1 logq="<<cipher1.logq<<endl;
  cout<<"ring.logQ="<<ring.logQ<<endl;
  
  // Rotation //
  long idx = 2;
  Ciphertext cipherRot = scheme.leftRotate(cipher1, idx);

  
  // Decrypt //
  complex<double>* dvec1 = scheme.decrypt(secretKey, cipherRot);

    cout<<"旋转之后的dvec1:"<<endl;
    for(int i=0;i<slots;i++){
        cout<<dvec1[i]<<" ";
    }
    cout<<endl;

    cout<<"*****************使用leftRotate对密文多项式进行旋转*********************"<<endl;
    cout<<"使用leftRotateFast"<<endl;
    cout<<"旋转之前的mvec2:"<<endl;
    for(int i=0;i<slots;i++){
        cout<<mvec2[i]<<" ";
    }
    cout<<endl;
    Ciphertext cipherRot2 = scheme.leftRotateFast(cipher2, 16);//fast旋转必须是二的幂
    complex<double>* dvec2 = scheme.decrypt(secretKey, cipherRot2);
    cout<<"使用Fast旋转之后的dvec2:"<<endl;
    for(int i=0;i<slots;i++){
        cout<<dvec2[i]<<" ";
    }
    cout<<endl;

    
}
  