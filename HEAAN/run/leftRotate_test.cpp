#include <iostream>
#include "../src/HEAAN.h"

using namespace std;
using namespace NTL;

int main() {

  // Parameters //
  //首先说明模链的定义。模链为一系列逐步变小的密文多项式系数模的集合{qL, qL-1, ..., q1}其中，ql=delta*q_(l-1)。delta是缩放因子。特别地，q1=delta*q0.q0和delta都是用户定义
  long logN = 6;//密文多项式次数的对数 这里试出来至少为6才能运行
  long logQ = 353;//密文多项式系数模的对数
	
  // Construct and Generate Public Keys //
  Ring ring(logN, logQ);//环R_Q=Z_Q[X]/(X^N+1)的两个系数：Q,N
  SecretKey secretKey(ring);
  Scheme scheme(secretKey, ring);
  scheme.addLeftRotKeys(secretKey); ///< When you need left rotation for the vectorized message
  scheme.addRightRotKeys(secretKey); ///< When you need right rotation for the vectorized message
  

  
  // Rotation //
//   Ciphertext cipherRot = scheme.leftRotate(cipher1, idx);
    //初始化多项式数组
    ZZ* a = new ZZ[ring.N];
    for(int i=0;i<ring.N;i++){
        a[i]=i+1;
    }
    cout<<"旋转之前的多项式数组:" <<endl;
    for(int i=0;i<ring.N;i++){
        cout<<a[i]<<" ";
    }
    cout<<endl;

    ZZ* aa = new ZZ[ring.N];
	ring.leftRotate(aa, a, 1);//r=0的时候，pow=1,不会改变数组。r=1的时候，旋转因子为5^1 mod M=5 
    cout<<"旋转之后的多项式数组:" <<endl;
    for(int i=0;i<ring.N;i++){
        cout<<aa[i]<<" ";
    }
    cout<<endl;

    ZZ* conj = new ZZ[ring.N];
	ring.conjugate(conj, a);//r=0的时候，pow=1,不会改变数组。r=1的时候，旋转因子为5^1 mod M=5 
    cout<<"共轭之后的多项式数组:" <<endl;
    for(int i=0;i<ring.N;i++){
        cout<<conj[i]<<" ";
    }
    cout<<endl;
}
  