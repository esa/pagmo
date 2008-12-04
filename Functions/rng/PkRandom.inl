
//Mersenne Twister pseudorandom number generator of 32 and 64-bit.
//2005, Diego Park <diegopark@gmail.com>
//
//A C-program for MT19937, with initialization improved 2002/1/26.
//Coded by Takuji Nishimura and Makoto Matsumoto.
//
//Before using, initialize the state by using init_genrand(seed)
//or init_by_array(init_key, key_length).
//
//Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
//All rights reserved.
//
//A C-program for MT19937-64 (2004/9/29 version).
//Coded by Takuji Nishimura and Makoto Matsumoto.
//
//This is a 64-bit version of Mersenne Twister pseudorandom number
//generator.
//
//Before using, initialize the state by using init_genrand64(seed)
//or init_by_array64(init_key, key_length).
//
//Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
//All rights reserved.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions
//are met:
//
//1. Redistributions of source code must retain the above copyright
//notice, this list of conditions and the following disclaimer.
//
//2. Redistributions in binary form must reproduce the above copyright
//notice, this list of conditions and the following disclaimer in the
//documentation and/or other materials provided with the distribution.
//
//3. The names of its contributors may not be used to endorse or promote
//products derived from this software without specific prior written
//permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//References:
//T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
//ACM Transactions on Modeling and
//Computer Simulation 10. (2000) 348--357.
//M. Matsumoto and T. Nishimura,
//``Mersenne Twister: a 623-dimensionally equidistributed
//uniform pseudorandom number generator''
//ACM Transactions on Modeling and
//Computer Simulation 8. (Jan. 1998) 3--30.
//
//Any feedback is very welcome.
//http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
//email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)

template <class Integer>
Integer Random<Integer>::next() {
	if (index >= MaxState) {
		seed();
	}

	return hashFunction(state[index++]);
}

template <class Integer>
void Random<Integer>::seed() {
	index = 0;

	size_t i;
	for (i = 0; i < MaxState - Period; i++) {
		state[i] = twist(state[i + Period], state[i], state[i + 1]);
	}
	for (; i < MaxState - 1; i++) {
		state[i] = twist(state[i + (Period - MaxState)], state[i], state[i + 1]);
	}
	state[MaxState - 1] = twist(state[Period - 1], state[MaxState - 1], state[0]);
}
