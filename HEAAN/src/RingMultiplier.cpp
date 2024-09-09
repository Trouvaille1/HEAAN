/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "RingMultiplier.h"

#include <NTL/BasicThreadPool.h>
#include <NTL/lip.h>
#include <NTL/sp_arith.h>
#include <NTL/tools.h>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <vector>

#include "Primes.h"

RingMultiplier::RingMultiplier(long logN, long logQ) : logN(logN) {
	N = 1 << logN;
	long M = N << 1;

	long bound = 2 + logN + 4 * logQ;//这是对乘法器内部所需模数大小的估计。它结合了环的大小和密文的模数来计算出乘法运算中所需的足够大范围的模数（模数位数的总和）。

	//59是每个质数的位数。NTL库中常用位数为59的质数来处理大整数运算，尤其是乘法。为了支持高效的大数运算，乘法器会使用模数分裂技术，将模数分成多个质数部分来处理。
	//ceil()计算得到需要多少个这样的59位质数。通过将bound除以59.0，得到质数的数量，并使用 ceil 函数向上取整，确保足够的质数。
	long nprimes = ceil(bound / 59.0);//CRT中质数数量。

	pVec = new uint64_t[nprimes];
	prVec = new uint64_t[nprimes];
	pTwok = new long[nprimes];
	pInvVec = new uint64_t[nprimes];
	scaledRootPows = new uint64_t*[nprimes];
	scaledRootInvPows = new uint64_t*[nprimes];
	scaledNInv = new uint64_t[nprimes];

	red_ss_array = new _ntl_general_rem_one_struct*[nprimes];

	for (long i = 0; i < nprimes; ++i) {
		pVec[i] = pPrimesVec[i];
		red_ss_array[i] = _ntl_general_rem_one_struct_build(pVec[i]);//这个函数通过预处理模数 pVec[i] 来加速后续的模数还原操作，返回一个结构体
		pInvVec[i] = inv(pVec[i]);
		pTwok[i] = (2 * ((long) log2(pVec[i]) + 1));
		prVec[i] = (static_cast<unsigned __int128>(1) << pTwok[i]) / pVec[i];
		uint64_t root = findMthRootOfUnity(M, pVec[i]);
		uint64_t rootinv = invMod(root, pVec[i]);
		uint64_t NInv = invMod(N, pVec[i]);
		mulMod(scaledNInv[i], NInv, (1ULL << 32), pVec[i]);
		mulMod(scaledNInv[i], scaledNInv[i], (1ULL << 32), pVec[i]);
		scaledRootPows[i] = new uint64_t[N]();
		scaledRootInvPows[i] = new uint64_t[N]();
		uint64_t power = 1;
		uint64_t powerInv = 1;
		for (long j = 0; j < N; ++j) {
			uint32_t jprime = bitReverse(static_cast<uint32_t>(j)) >> (32 - logN);
			uint64_t rootpow = power;
			mulMod(scaledRootPows[i][jprime], rootpow,(1ULL << 32), pVec[i]);
			mulMod(scaledRootPows[i][jprime], scaledRootPows[i][jprime], (1ULL << 32), pVec[i]);
			uint64_t rootpowInv = powerInv;
			mulMod(scaledRootInvPows[i][jprime], rootpowInv, (1ULL << 32), pVec[i]);
			mulMod(scaledRootInvPows[i][jprime], scaledRootInvPows[i][jprime], (1ULL << 32), pVec[i]);
			mulMod(power, power, root, pVec[i]);
			mulMod(powerInv, powerInv, rootinv, pVec[i]);
		}
	}

	coeffpinv_array = new mulmod_precon_t*[nprimes];
	pProd = new ZZ[nprimes];
	pProdh = new ZZ[nprimes];
	pHat = new ZZ*[nprimes];
	pHatInvModp = new uint64_t*[nprimes];
	for (long i = 0; i < nprimes; ++i) {
		pProd[i] = (i == 0) ? to_ZZ((long) pVec[i]) : pProd[i - 1] * (long) pVec[i];
		pProdh[i] = pProd[i] / 2;
		pHat[i] = new ZZ[i + 1];
		pHatInvModp[i] = new uint64_t[i + 1];
		coeffpinv_array[i] = new mulmod_precon_t[i + 1];
		for (long j = 0; j < i + 1; ++j) {
			pHat[i][j] = ZZ(1);
			for (long k = 0; k < j; ++k) {
				pHat[i][j] *= (long) pVec[k];
			}
			for (long k = j + 1; k < i + 1; ++k) {
				pHat[i][j] *= (long) pVec[k];
			}
			pHatInvModp[i][j] = to_long(pHat[i][j] % (long) pVec[j]);
			pHatInvModp[i][j] = invMod(pHatInvModp[i][j], pVec[j]);
			coeffpinv_array[i][j] = PrepMulModPrecon(pHatInvModp[i][j], pVec[j]);
		}
	}
}

//Cooley-Tukey Radix-2 NTT.
//见论文（这篇论文arxiv版本才有该算法的详细描述，和这个函数的实现一模一样，正式出版的反而没有）P3：https://arxiv.org/pdf/2103.16400
//原算法论文：https://link.springer.com/chapter/10.1007/978-3-319-48965-0_8
void RingMultiplier::NTT(uint64_t* a, long index) {
	long t = N;
	long logt1 = logN + 1;
	uint64_t p = pVec[index];//参数a就是rxi，为rx“矩阵”的第i行。p就是这一行代表的CRT的哪个基
	uint64_t pInv = pInvVec[index];//预计算的该行基的逆元
	//m是当前层每个小蝴蝶结构的长度的一半。m的变化序列：1,2,4,8,...
	for (long m = 1; m < N; m <<= 1) {
		//进入每一层。
		t >>= 1;//t既是当前层每个小蝴蝶结构的两个输入之间的坐标的gap，又是当前层蝴蝶结构的个数。t的变化序列：N/2,N/4,...
		logt1 -= 1;//t1指的是当前层每个小蝴蝶结构两个待处理的数的坐标gap的两倍。即t1=2*t，logt1=logt+1。t1变化序列：N,N/2,N/4,...
		//遍历当前层每个小蝴蝶结构的前半部分即可（即遍历0~m-1），因为蝴蝶运算每次处理两个点
		for (long i = 0; i < m; i++) {
			long j1 = i << logt1;//j1是当前层第一个蝴蝶结构的起始坐标。
			long j2 = j1 + t - 1;//j2是当前层最后一个蝴蝶结构的起始坐标。
			uint64_t W = scaledRootPows[index][m + i];
			//j1到j2的长度为t，t是当前层蝴蝶结构的个数。循环中对每个蝴蝶结构进行运算
			for (long j = j1; j <= j2; j++) {
				//进入当前层当前蝴蝶结构。
				//DIT NTT的蝴蝶运算部分。每次处理a[j]和a[j+t]两个点。公式为：a[j]=(a[j]+a[j+t]*W) mod p，a[j+t]=(a[j]-a[j+t]*W) mod p
				//使用了 Barrett Reduction。见：https://en.wikipedia.org/wiki/Barrett_reduction
				//https://maskray.me/blog/2016-10-03-discrete-fourier-transform
				uint64_t T = a[j + t];
				unsigned __int128
				U = static_cast<unsigned __int128>(T) * W;//U=T*W=a[j + t]*W(U为128位，防止溢出)
				uint64_t U0 = static_cast<uint64_t>(U);//从128位到64位的截断。U0为U的低64位
				uint64_t U1 = U >> 64;//U缩小2^64倍.U1为U的高64位
				uint64_t Q = U0 * pInv;//Q约等于U*p^-1=T*W*p^-1=a[j + t]*W*p^-1
				unsigned __int128
				Hx = static_cast<unsigned __int128>(Q) * p;//Hx约等于Q*p=T*W=a[j + t]*W=U
				uint64_t H = Hx >> 64;//Hx缩小2^64倍。H是Hx的高位部分
				uint64_t V = U1 < H ? U1 + p - H : U1 - H;
				a[j + t] = a[j] < V ? a[j] + p - V : a[j] - V;
				a[j] += V;
				if (a[j] > p)
					a[j] -= p;
			}
		}
	}
}

void RingMultiplier::INTT(uint64_t* a, long index) {
	uint64_t p = pVec[index];
	uint64_t pInv = pInvVec[index];
	long t = 1;
	for (long m = N; m > 1; m >>= 1) {
		long j1 = 0;
		long h = m >> 1;
		for (long i = 0; i < h; i++) {
			long j2 = j1 + t - 1;
			uint64_t W = scaledRootInvPows[index][h + i];
			for (long j = j1; j <= j2; j++) {
				uint64_t U = a[j] + a[j + t];
				if (U > p)
					U -= p;
				uint64_t T =
						a[j] < a[j + t] ? a[j] + p - a[j + t] : a[j] - a[j + t];
				unsigned __int128
				UU = static_cast<unsigned __int128>(T) * W;
				uint64_t U0 = static_cast<uint64_t>(UU);
				uint64_t U1 = UU >> 64;
				uint64_t Q = U0 * pInv;
				unsigned __int128
				Hx = static_cast<unsigned __int128>(Q) * p;
				uint64_t H = Hx >> 64;
				a[j] = U;
				a[j + t] = (U1 < H) ? U1 + p - H : U1 - H;
			}
			j1 += (t << 1);
		}
		t <<= 1;
	}

	uint64_t NScale = scaledNInv[index];
	for (long i = 0; i < N; i++) {
		uint64_t T = a[i];
		unsigned __int128
		U = static_cast<unsigned __int128>(T) * NScale;
		uint64_t U0 = static_cast<uint64_t>(U);
		uint64_t U1 = U >> 64;
		uint64_t Q = U0 * pInv;
		unsigned __int128
		Hx = static_cast<unsigned __int128>(Q) * p;
		uint64_t H = Hx >> 64;
		a[i] = (U1 < H) ? U1 + p - H : U1 - H;
	}
}

//----------------------------------------------------------------------------------
//   FFT
//----------------------------------------------------------------------------------

//将一个多项式的系数数组x，转换为CRT的NTT格式的多项式系数矩阵rx
uint64_t* RingMultiplier::toNTT(ZZ* x, long np) {
	uint64_t* rx = new uint64_t[np << logN]();//一维数组rx的大小为N*np(13*2^15=425984)。这里的意思是将x表示为中国剩余定理CRT格式，最多有N的系数，每个系数用np个小素数表示
	NTL_EXEC_RANGE(np, first, last);//并行循环
	for (long i = first; i < last; ++i) {
		uint64_t* rxi = rx + (i << logN);//由于使用了并行循环，这是内部block的位置偏移。rxi指向irst到last的小block中的第i行，这一行表示第i个CRT的基，长度为N，表示N个系数都需要和这个基做模运算
		uint64_t pi = pVec[i];//first到last的小block中，当前需要处理的第i个CRT的基
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];//预计算的关于模数的信息
		//遍历x系数数组。x中每一个系数都会breakdown为CRT形式
		for (long n = 0; n < N; ++n) {
			rxi[n] = _ntl_general_rem_one_struct_apply(x[n].rep, pi, red_ss);//用该函数快速计算 系数 mod pi 的值，然后存入第i行的第n个数。猜测这个函数为NTL的原语函数，所以要用x[n].rep这种底层的数据结构
		}
		NTT(rxi, i);//第i行做NTT（相当于原来系数在CRT下的一个分量的余数结果做NTT）。注意，rxi中的i指的就是第i行（因为rx是一维数组，只能这么表示），所以index就是1
	}
	NTL_EXEC_RANGE_END;////并行循环
	return rx;
}

void RingMultiplier::addNTTAndEqual(uint64_t* ra, uint64_t* rb, long np) {
	for (long i = 0; i < np; ++i) {
		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logN);
		uint64_t pi = pVec[i];
		for (long n = 0; n < N; ++n) {
			rai[n] += rbi[n];
			if(rai[n] > pi) rai[n] -= pi;
		}
	}
}

//使用CRT将CRT形式下的多项式恢复成普通形式
void RingMultiplier::reconstruct(ZZ* x, uint64_t* rx, long np, ZZ& mod) {
	ZZ* pHatnp = pHat[np - 1];
	uint64_t* pHatInvModpnp = pHatInvModp[np - 1];
	mulmod_precon_t* coeffpinv_arraynp = coeffpinv_array[np - 1];
	ZZ& pProdnp = pProd[np - 1];
	ZZ& pProdhnp = pProdh[np - 1];

	NTL_EXEC_RANGE(N, first, last);
	for (long n = first; n < last; ++n) {
		ZZ& acc = x[n];
		QuickAccumBegin(acc, pProdnp.size());
		for (long i = 0; i < np; i++) {
			long p = pVec[i];
			long tt = pHatInvModpnp[i];
			mulmod_precon_t ttpinv = coeffpinv_arraynp[i];
			long s = MulModPrecon(rx[n + (i << logN)], tt, p, ttpinv);
			QuickAccumMulAdd(acc, pHatnp[i], s);
		}
		QuickAccumEnd(acc);
		QuickRem(x[n], pProdnp);
		if (x[n] > pProdhnp) x[n] -= pProdnp;
		QuickRem(x[n], mod);
	}
	NTL_EXEC_RANGE_END;

}

void RingMultiplier::mult(ZZ* x, ZZ* a, ZZ* b, long np, ZZ& mod) {
	uint64_t* ra = new uint64_t[np << logN]();
	uint64_t* rb = new uint64_t[np << logN]();
	uint64_t* rx = new uint64_t[np << logN]();

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logN);
		uint64_t* rxi = rx + (i << logN);
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];
		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
			rbi[n] = _ntl_general_rem_one_struct_apply(b[n].rep, pi, red_ss);
		}
		NTT(rai, i);
		NTT(rbi, i);
		for (long n = 0; n < N; ++n) {
			mulModBarrett(rxi[n], rai[n], rbi[n], pi, pri, pTwoki);
		}
		INTT(rxi, i);
	}
	NTL_EXEC_RANGE_END;

	reconstruct(x, rx, np, mod);

	delete[] ra;
	delete[] rb;
	delete[] rx;
}

//使用NTT、INTT计算多项式乘法，x=a*b
void RingMultiplier::multNTT(ZZ* x, ZZ* a, uint64_t* rb, long np, ZZ& mod) {
	//变量名带r的表示CRT矩阵形式
	uint64_t* ra = new uint64_t[np << logN]();
	uint64_t* rx = new uint64_t[np << logN]();

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		//第i个CRT的基
		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logN);
		uint64_t* rxi = rx + (i << logN);
		
		uint64_t pi = pVec[i];//first到last的小block中，当前需要处理的第i个CRT的基
		uint64_t pri = prVec[i];//当前需要处理的第i个CRT的基的倒数
		long pTwoki = pTwok[i];//预计算的2的指数
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];//预计算的关于模数的信息
		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
		}
		NTT(rai, i);
		for (long n = 0; n < N; ++n) {
			mulModBarrett(rxi[n], rai[n], rbi[n], pi, pri, pTwoki);
		}
		INTT(rxi, i);
	}
	NTL_EXEC_RANGE_END;

	reconstruct(x, rx, np, mod);

	delete[] ra;
	delete[] rx;
}

//两个NTT形式的多项式相乘，最后返回的结果为非NTT形式的多项式。x=a*b
//mult double NTT
void RingMultiplier::multDNTT(ZZ* x, uint64_t* ra, uint64_t* rb, long np, ZZ& mod) {
	uint64_t* rx = new uint64_t[np << logN]();//NTT形式的结果多项式

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logN);
		uint64_t* rxi = rx + (i << logN);
		uint64_t pi = pVec[i];
		uint64_t pri = prVec[i];
		long pTwoki = pTwok[i];
		_ntl_general_rem_one_struct* red_ss = red_ss_array[i];
		for (long n = 0; n < N; ++n) {
			mulModBarrett(rxi[n], rai[n], rbi[n], pi, pri, pTwoki);//由于是NTT形式的多项式，做element-wise的乘法即可
		}
		INTT(rxi, i);
	}
	NTL_EXEC_RANGE_END;

	reconstruct(x, rx, np, mod);

	delete[] rx;
}

void RingMultiplier::multAndEqual(ZZ* a, ZZ* b, long np, ZZ& mod) {
	uint64_t* ra = new uint64_t[np << logN]();
	uint64_t* rb = new uint64_t[np << logN]();

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logN);
		uint64_t pi = pVec[i];
		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss_array[i]);
			rbi[n] = _ntl_general_rem_one_struct_apply(b[n].rep, pi, red_ss_array[i]);
		}
		NTT(rai, i);
		NTT(rbi, i);
		for (long n = 0; n < N; ++n) {
			mulModBarrett(rai[n], rai[n], rbi[n], pi, prVec[i], pTwok[i]);
		}
		INTT(rai, i);
	}
	NTL_EXEC_RANGE_END;

	ZZ* pHatnp = pHat[np - 1];
	uint64_t* pHatInvModpnp = pHatInvModp[np - 1];

	reconstruct(a, ra, np, mod);

	delete[] ra;
	delete[] rb;
}

void RingMultiplier::multNTTAndEqual(ZZ* a, uint64_t* rb, long np, ZZ& mod) {
	uint64_t* ra = new uint64_t[np << logN]();

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t* rai = ra + (i << logN);
		uint64_t* rbi = rb + (i << logN);
		uint64_t pi = pVec[i];
		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss_array[i]);
		}
		NTT(rai, i);
		for (long n = 0; n < N; ++n) {
			mulModBarrett(rai[n], rai[n], rbi[n], pi, prVec[i], pTwok[i]);
		}
		INTT(rai, i);
	}
	NTL_EXEC_RANGE_END;

	ZZ* pHatnp = pHat[np - 1];
	uint64_t* pHatInvModpnp = pHatInvModp[np - 1];

	reconstruct(a, ra, np, mod);

	delete[] ra;
}


void RingMultiplier::square(ZZ* x, ZZ* a, long np, ZZ& mod) {
	uint64_t* ra = new uint64_t[np << logN]();
	uint64_t* rx = new uint64_t[np << logN]();

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t* rai = ra + (i << logN);
		uint64_t* rxi = rx + (i << logN);
		uint64_t pi = pVec[i];
		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss_array[i]);
		}
		NTT(rai, i);
		for (long n = 0; n < N; ++n) {
			mulModBarrett(rxi[n], rai[n], rai[n], pi, prVec[i], pTwok[i]);
		}
		INTT(rxi, i);
	}
	NTL_EXEC_RANGE_END;

	ZZ* pHatnp = pHat[np - 1];
	uint64_t* pHatInvModpnp = pHatInvModp[np - 1];

	reconstruct(x, rx, np, mod);

	delete[] ra;
	delete[] rx;
}

void RingMultiplier::squareNTT(ZZ* x, uint64_t* ra, long np, ZZ& mod) {
	uint64_t* rx = new uint64_t[np << logN]();

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t* rai = ra + (i << logN);
		uint64_t* rxi = rx + (i << logN);
		uint64_t pi = pVec[i];
		for (long n = 0; n < N; ++n) {
			mulModBarrett(rxi[n], rai[n], rai[n], pi, prVec[i], pTwok[i]);
		}
		INTT(rxi, i);
	}
	NTL_EXEC_RANGE_END;

	reconstruct(x, rx, np, mod);

	delete[] rx;
}

void RingMultiplier::squareAndEqual(ZZ* a, long np, ZZ& mod) {
	uint64_t* ra = new uint64_t[np << logN]();

	NTL_EXEC_RANGE(np, first, last);
	for (long i = first; i < last; ++i) {
		uint64_t* rai = ra + (i << logN);
		uint64_t pi = pVec[i];
		for (long n = 0; n < N; ++n) {
			rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss_array[i]);
		}
		NTT(rai, i);
		for (long n = 0; n < N; ++n) {
			mulModBarrett(rai[n], rai[n], rai[n], pi, prVec[i], pTwok[i]);
		}
		INTT(rai, i);
	}
	NTL_EXEC_RANGE_END;

	reconstruct(a, ra, np, mod);

	delete[] ra;
}

void RingMultiplier::mulMod(uint64_t &r, uint64_t a, uint64_t b, uint64_t m) {
	unsigned __int128 mul = static_cast<unsigned __int128>(a) * b;
	mul %= static_cast<unsigned __int128>(m);
	r = static_cast<uint64_t>(mul);
}

//r = a * b mod p
void RingMultiplier::mulModBarrett(uint64_t& r, uint64_t a, uint64_t b, uint64_t p, uint64_t pr, long twok) {
	unsigned __int128 mul = static_cast<unsigned __int128>(a) * b;
	uint64_t atop, abot;
	abot = static_cast<uint64_t>(mul);
	atop = static_cast<uint64_t>(mul >> 64);
	unsigned __int128 tmp = static_cast<unsigned __int128>(abot) * pr;
	tmp >>= 64;
	tmp += static_cast<unsigned __int128>(atop) * pr;
	tmp >>= twok - 64;
	tmp *= p;
	tmp = mul - tmp;
	r = static_cast<uint64_t>(tmp);
	if(r >= p) {
		r -= p;
	}
}

uint64_t RingMultiplier::invMod(uint64_t x, uint64_t m) {
	return powMod(x, m - 2, m);
}

uint64_t RingMultiplier::powMod(uint64_t x, uint64_t y, uint64_t modulus) {
	uint64_t res = 1;
	while (y > 0) {
		if (y & 1) {
			mulMod(res, res, x, modulus);
		}
		y = y >> 1;
		mulMod(x, x, x, modulus);
	}
	return res;
}

uint64_t RingMultiplier::inv(uint64_t x) {
	return pow(x, static_cast<uint64_t>(-1));
}

uint64_t RingMultiplier::pow(uint64_t x, uint64_t y) {
	uint64_t res = 1;
	while (y > 0) {
		if (y & 1) {
			res *= x;
		}
		y = y >> 1;
		x *= x;
	}
	return res;
}

uint32_t RingMultiplier::bitReverse(uint32_t x) {
	x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));
	x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));
	x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));
	x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));
	return ((x >> 16) | (x << 16));
}

void RingMultiplier::findPrimeFactors(vector<uint64_t> &s, uint64_t number) {
	while (number % 2 == 0) {
		s.push_back(2);
		number /= 2;
	}
	for (uint64_t i = 3; i < sqrt(number); i++) {
		while (number % i == 0) {
			s.push_back(i);
			number /= i;
		}
	}
	if (number > 2) {
		s.push_back(number);
	}
}

uint64_t RingMultiplier::findPrimitiveRoot(uint64_t modulus) {
	vector<uint64_t> s;
	uint64_t phi = modulus - 1;
	findPrimeFactors(s, phi);
	for (uint64_t r = 2; r <= phi; r++) {
		bool flag = false;
		for (auto it = s.begin(); it != s.end(); it++) {
			if (powMod(r, phi / (*it), modulus) == 1) {
				flag = true;
				break;
			}
		}
		if (flag == false) {
			return r;
		}
	}
	return -1;
}

// Algorithm to find m-th primitive root in Z_mod
uint64_t RingMultiplier::findMthRootOfUnity(uint64_t M, uint64_t mod) {
    uint64_t res;
    res = findPrimitiveRoot(mod);
    if((mod - 1) % M == 0) {
        uint64_t factor = (mod - 1) / M;
        res = powMod(res, factor, mod);
        return res;
    }
    else {
        return -1;
    }
}


