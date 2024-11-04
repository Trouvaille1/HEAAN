#include <iostream>
#include <fstream>
#include "../src/HEAAN.h"

using namespace std;
using namespace NTL;

int main() {

    // Parameters //
    long logN = 10;//密文多项式次数的对数 这里试出来至少为6才能运行
    long logQ = 353;//密文多项式系数模的对数
	
    // Construct and Generate Public Keys //
    TimeUtils timeutils;
    Ring ring(logN, logQ);//环R_Q=Z_Q[X]/(X^N+1)的两个系数：Q,N

    cout<<"下面输出ring的所有参数"<<endl;
    cout<<"ring.logN="<<ring.logN<<endl;
    cout<<"ring.logQ="<<ring.logQ<<endl;
    cout<<"ring.N="<<ring.N<<endl;
    cout<<"ring.M="<<ring.M<<endl;
    cout<<"ring.Nh="<<ring.Nh<<endl;
    cout<<"ring.logQQ="<<ring.logQQ<<endl;
    cout<<"ring.Q="<<ring.Q<<endl;
    cout<<"ring.QQ="<<ring.QQ<<endl;
    cout<<"ring.qpows为："<<endl;

    //不需要用到所有的qpows
    // for(int i=0;i<ring.logQQ+1;i++){
    //     cout<<ring.qpows[i]<<" ";
    // }
    // cout<<endl;

    long logp = 30; ///缩放因子delta的对数。< Larger logp will give you more correct result (smaller computation noise)
    long slots = 1024; ///< This should be power of two
    long numThread = 8;

    SecretKey secretKey(ring);
    Scheme scheme(secretKey, ring);

    // Make Random Array of Complex //
    complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(slots);
    complex<double>* mvec2 = EvaluatorUtils::randomComplexArray(slots);
  
    // Encrypt Two Arry of Complex //
    Ciphertext cipher1 = scheme.encrypt(mvec1, slots, logp, logQ);
    Ciphertext cipher2 = scheme.encrypt(mvec2, slots, logp, logQ);

    Ciphertext cipherMult = scheme.mult(cipher1, cipher2); //先不考虑完整的乘法，只测试内部NTT的结果是否正确

    ZZ q = ring.qpows[cipher1.logq];//cipher1.logq=Q
    cout<<"q="<<q<<endl;
	ZZ qQ = ring.qpows[cipher1.logq + ring.logQ];//Q^2
    cout<<"qQ="<<qQ<<endl;
	long np = ceil((2 + cipher1.logq + cipher2.logq + ring.logN + 2)/59.0);//中国剩余定理CRT基的余数个数
    cout<<"main函数中np="<<np<<endl;

    //拿到rxi和i的数据
    uint64_t* rx = new uint64_t[np << logN]();//一维数组rx的大小为N*np(13*2^15=425984)。这里的意思是将x表示为中国剩余定理CRT格式，最多有N的系数，每个系数用np个小素数表示
    uint64_t* rxi;
    long i;
	NTL_EXEC_RANGE(np, first, last);//并行循环
	for (i = first; i < last; ++i) {
		rxi = rx + (i << logN);//由于使用了并行循环，这是内部block的位置偏移。rxi指向first到last的小block中的第i行，这一行表示第i个CRT的基，长度为N，表示N个系数都需要和这个基做模运算
		uint64_t pi = ring.multiplier.pVec[i];//first到last的小block中，当前需要处理的第i个CRT的基
		uint64_t pri = ring.multiplier.prVec[i];
		long pTwoki = ring.multiplier.pTwok[i];
		_ntl_general_rem_one_struct* red_ss =ring.multiplier.red_ss_array[i];//预计算的关于模数的信息
		//遍历x系数数组。x中每一个系数都会breakdown为CRT形式
		for (long n = 0; n < ring.N; ++n) {
			rxi[n] = _ntl_general_rem_one_struct_apply(cipher1.ax[n].rep, pi, red_ss);//用该函数快速计算 系数 mod pi 的值，然后存入第i行的第n个数。猜测这个函数为NTL的原语函数，所以要用x[n].rep这种底层的数据结构
		}
		// NTT(rxi, i);
        //拿到rxi和i的数据就跳出循环
        break;
	}
	NTL_EXEC_RANGE_END;////并行循环

    cout<<"i为："<<i<<endl;

    cout<<"下面输出ring.multiplier的所有参数"<<endl;
    
    //ring.multiplier内的nprimes为25  nprimes=ceil(( 2 + logN + 4 * logQ)/ 59.0)
    cout<<"pVec为："<<endl;
    for(int i=0;i<25;i++){
        cout<<ring.multiplier.pVec[i]<<",";
    }
    cout<<endl;

    cout<<"prVec为："<<endl;
    for(int i=0;i<25;i++){
        cout<<ring.multiplier.prVec[i]<<",";
    }
    cout<<endl;

    cout<<"pTwok为："<<endl;
    for(int i=0;i<25;i++){
        cout<<ring.multiplier.pTwok[i]<<",";
    }
    cout<<endl;

    cout<<"pInvVec为："<<endl;
    for(int i=0;i<25;i++){
        cout<<ring.multiplier.pInvVec[i]<<",";
    }
    cout<<endl;

    // 打开文件输出流
    std::ofstream outFile("scaledRootPows0.txt");

    if (!outFile) {
        std::cerr << "无法打开文件进行写入" << std::endl;
        return -1;
    }

    // 写入文件头部，声明数组
    outFile << "const UINT64_T scaledRootPows[N] = {";  // 现在只有一个一维数组

    // 仅遍历第0行的数组并格式化输出
    for (int j = 0; j < ring.N; j++) {
        outFile << ring.multiplier.scaledRootPows[0][j];  // 输出第0行的数组元素
        if (j < ring.N - 1) {
            outFile << ", ";  // 在每个元素后加逗号
        }
    }

    // 写入数组尾部
    outFile << "};" << std::endl;

    // 关闭文件输出流
    outFile.close();
    std::cout << "scaledRootPows 已导出到 scaledRootPows0.txt 文件" << std::endl;

    // cout<<"scaledRootInvPows为："<<endl;
    // for(int i=0;i<25;i++){
    //     for(int j=0;j<ring.N;j++){
    //         cout<<ring.multiplier.scaledRootInvPows[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;

    

    // cout<<"ring.rotGroup为："<<endl;
    // for(int i=0;i<ring.Nh;i++){
    //     cout<<ring.rotGroup[i]<<endl;
    // }
    // cout<<endl;

    // cout<<"ring.ksiPows为："<<endl;
    // for(int i=0;i<ring.M+1;i++){
    //     cout<<ring.ksiPows[i]<<endl;
    // }
    // cout<<endl;

    // cout<<"ring.taylorCoeffsMap为："<<endl;
    // for(auto it=ring.taylorCoeffsMap.begin();it!=ring.taylorCoeffsMap.end();it++){
    //     cout<<it->first<<endl;
    //     for(int i=0;i<ring.N;i++){
    //         cout<<it->second[i]<<" ";
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;



}