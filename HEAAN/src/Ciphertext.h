/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/
#ifndef HEAAN_CIPHERTEXT_H_
#define HEAAN_CIPHERTEXT_H_

#include <NTL/ZZ.h>

#include <fstream>

using namespace std;
using namespace NTL;

/**
 * Ciphertext is an RLWE instance (ax, bx = mx + ex - ax * sx) in the ring Z_q[X] / (X^N + 1)
 */

//ct和pt都要存储p、q、N、n（消息槽数和密文槽数）。pt存储mx（多项式的大整数系数数组），ct存储ax和bx（两个多项式的大整数系数数组）。
class Ciphertext {
public:

	ZZ* ax; ///< a(x) - part of RLWE instance
	ZZ* bx; ///< b(x) - part of RLWE instance

	//在同态加密（特别是CKKS方案）中，明文通常不是整数，而是浮点数或小数。这些数在加密前需要通过量化来转化为整数，因为加密方案本质上只能处理整数运算。
	long logp; ///< number of message quantized bits   缩放因子delta的对数（精度），即小数部分需要保留的精度（小数位数）
	long logq; ///< number of modulus bits    密文模的对数（可以看作是SEAL中的scale。在SEAL中，密文的scale必须相同才能做乘、加运算）

	long N; ///< degree of RLWE
	long n; ///< number of slots

	Ciphertext(ZZ* ax = NULL, ZZ* bx = NULL, long logp = 0, long logq = 0, long N = 0, long n = 0);

	Ciphertext(const Ciphertext& o);

	Ciphertext& operator=(const Ciphertext &o);

	~Ciphertext();

};

#endif
