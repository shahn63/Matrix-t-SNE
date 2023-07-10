#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix jointprob(NumericMatrix mat1, NumericMatrix mat2, int n1, int n2) {

  // 결과를 저장할 행렬
	NumericMatrix result(n1, n1), out(n1,n1), tresult(n1,n1), num1(n1,n1);
	NumericMatrix mat3(n1,2);

  // 벡터의 원소를 제곱하고 더하기
	for (int jj = 0; jj < n2; jj++){
		NumericVector sk2(n1);
		mat3(_,0)=mat1(_,jj);
		mat3(_,1)=mat2(_,jj);
		for (int i = 0; i < n1; i++){
			double value1 = mat1(i, jj);
			double value2 = mat2(i, jj);
			sk2(i) = pow(value1, 2) + pow(value2, 2);
			for (int j = 0; j < n1; j++) {
				double sum = 0.0;
				for (int k = 0; k < 2; ++k) {
					sum += mat3(i, k) * mat3(j, k);  // A의 i번째 열과 j번째 열의 내적
				}
				out(i, j) = -2 * sum;
			}
		}


  // 전치행렬 계산
		for (int i = 0; i < n1; i++) {
			result(_,i)=sk2;
			tresult(i,_)=sk2;
		}

		num1 += result + out + tresult;
	}
	return num1;
}

